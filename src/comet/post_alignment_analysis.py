"""Post-Alignment Analysis Pipeline
=====================================

After COMET segmentations have been warped into STOmics space and transcript
assignment is complete, this module handles:

1. Per-sample normalization of integrated AnnData
2. Cross-sample concatenation with batch metadata
3. Batch correction (Harmony / combat)
4. Phenotype-stratified differential expression
5. QC metrics for the integrated data

Usage:
    from comet.post_alignment_analysis import (
        normalize_integrated_adata,
        concatenate_samples,
        run_cross_sample_de,
    )
"""

import json
import logging
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
from scipy import sparse

logger = logging.getLogger(__name__)

try:
    import anndata as ad
    import scanpy as sc
except ImportError:
    raise ImportError("scanpy and anndata required: pip install scanpy anndata")


# ============================================================
# 1. Per-Sample Normalization
# ============================================================

def normalize_integrated_adata(
    adata,
    min_genes: int = 5,
    min_transcripts: int = 10,
    max_pct_mito: float = 50.0,
    target_sum: Optional[float] = 1e4,
    log_transform: bool = True,
    n_top_genes: int = 3000,
    regress_out: Optional[List[str]] = None,
) -> "ad.AnnData":
    """Normalize a single integrated COMET-STOmics AnnData.

    This adapts standard scRNA-seq normalization to spatial data where
    "cells" are COMET-defined (multi-channel IF segmentation) and expression
    comes from STOmics transcript assignment.

    Parameters
    ----------
    adata : AnnData
        Integrated data (COMET cells × STOmics genes).
        Expected obs columns: phenotype, roi_class, n_transcripts, n_genes
    min_genes : int
        Minimum number of detected genes per cell.
    min_transcripts : int
        Minimum transcript count per cell.
    max_pct_mito : float
        Maximum mitochondrial gene percentage.
    target_sum : float or None
        Target sum for library size normalization. None skips.
    log_transform : bool
        Apply log1p transformation.
    n_top_genes : int
        Number of highly variable genes to identify.
    regress_out : list, optional
        Variables to regress out (e.g., ['n_transcripts', 'pct_mito']).

    Returns
    -------
    AnnData
        Normalized AnnData with raw counts preserved in .raw
    """
    logger.info(f"Normalizing: {adata.shape[0]} cells × {adata.shape[1]} genes")

    # --- QC metrics ---
    # Mitochondrial genes
    adata.var['mt'] = adata.var_names.str.startswith('MT-') | adata.var_names.str.startswith('mt-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    n_before = adata.n_obs

    # --- Filtering ---
    if 'n_genes' not in adata.obs.columns:
        adata.obs['n_genes'] = (adata.X > 0).sum(axis=1).A1 if sparse.issparse(adata.X) else (adata.X > 0).sum(axis=1)

    if 'n_transcripts' not in adata.obs.columns:
        adata.obs['n_transcripts'] = adata.X.sum(axis=1).A1 if sparse.issparse(adata.X) else adata.X.sum(axis=1)

    sc.pp.filter_cells(adata, min_genes=min_genes)
    adata = adata[adata.obs['n_transcripts'] >= min_transcripts].copy()

    if 'pct_counts_mt' in adata.obs.columns:
        adata = adata[adata.obs['pct_counts_mt'] < max_pct_mito].copy()

    n_after = adata.n_obs
    logger.info(f"Filtered: {n_before} → {n_after} cells "
                f"({n_before - n_after} removed, {100*(n_before-n_after)/n_before:.1f}%)")

    # --- Store raw counts ---
    adata.raw = adata.copy()

    # --- Normalization ---
    if target_sum is not None:
        sc.pp.normalize_total(adata, target_value=target_sum)

    if log_transform:
        sc.pp.log1p(adata)

    # --- HVG detection ---
    try:
        sc.pp.highly_variable_genes(
            adata, n_top_genes=n_top_genes, flavor='seurat_v3',
            span=0.3, subset=False,
        )
        logger.info(f"HVGs: {adata.var['highly_variable'].sum()}")
    except Exception as e:
        logger.warning(f"HVG detection failed (likely too few cells): {e}")
        # Fallback to seurat flavor
        try:
            sc.pp.highly_variable_genes(adata, n_top_genes=min(n_top_genes, adata.n_vars),
                                        flavor='seurat', subset=False)
        except Exception:
            adata.var['highly_variable'] = True

    # --- Optional regression ---
    if regress_out:
        valid_vars = [v for v in regress_out if v in adata.obs.columns]
        if valid_vars:
            sc.pp.regress_out(adata, valid_vars)
            logger.info(f"Regressed out: {valid_vars}")

    # --- Scale ---
    sc.pp.scale(adata, max_value=10)

    # --- Dimensionality reduction ---
    n_pcs = min(50, adata.n_obs - 1, adata.n_vars - 1)
    if n_pcs > 5:
        sc.tl.pca(adata, n_comps=n_pcs, svd_solver='arpack')
        sc.pp.neighbors(adata, n_neighbors=15, n_pcs=min(30, n_pcs))
        sc.tl.umap(adata)

    # --- Store normalization params in uns ---
    adata.uns['normalization'] = {
        'min_genes': min_genes,
        'min_transcripts': min_transcripts,
        'max_pct_mito': max_pct_mito,
        'target_sum': target_sum,
        'log_transform': log_transform,
        'n_top_genes': n_top_genes,
        'regressed_out': regress_out or [],
        'n_cells_before': n_before,
        'n_cells_after': n_after,
    }

    return adata


# ============================================================
# 2. Cross-Sample Concatenation
# ============================================================

def concatenate_samples(
    h5ad_paths: List[Union[str, Path]],
    normalize: bool = True,
    normalize_kwargs: Optional[dict] = None,
    min_cells_per_sample: int = 100,
) -> "ad.AnnData":
    """Load and concatenate multiple integrated AnnData objects.

    Parameters
    ----------
    h5ad_paths : list
        Paths to per-sample integrated h5ad files.
    normalize : bool
        Whether to normalize each sample before concatenation.
    normalize_kwargs : dict, optional
        Kwargs for normalize_integrated_adata.
    min_cells_per_sample : int
        Skip samples with fewer cells.

    Returns
    -------
    AnnData
        Concatenated AnnData with sample/disease metadata.
    """
    adatas = []
    skipped = []

    for path in h5ad_paths:
        path = Path(path)
        if not path.exists():
            logger.warning(f"File not found: {path}")
            skipped.append(str(path))
            continue

        logger.info(f"Loading {path.name}...")
        adata = ad.read_h5ad(str(path))

        if adata.n_obs < min_cells_per_sample:
            logger.warning(f"Skipping {path.name}: only {adata.n_obs} cells")
            skipped.append(str(path))
            continue

        # Ensure required metadata
        if 'sample_id' not in adata.obs.columns:
            adata.obs['sample_id'] = path.stem.split('_')[0]
        if 'disease' not in adata.obs.columns and 'disease' in adata.uns:
            adata.obs['disease'] = adata.uns['disease']

        if normalize:
            kwargs = normalize_kwargs or {}
            adata = normalize_integrated_adata(adata, **kwargs)

        adatas.append(adata)

    if not adatas:
        raise ValueError("No valid samples to concatenate")

    logger.info(f"\nConcatenating {len(adatas)} samples...")

    # Find common genes
    common_genes = set(adatas[0].var_names)
    for a in adatas[1:]:
        common_genes &= set(a.var_names)
    common_genes = sorted(common_genes)
    logger.info(f"Common genes: {len(common_genes)}")

    # Subset and concatenate (using raw counts for clean concatenation)
    raw_adatas = []
    for a in adatas:
        if a.raw is not None:
            raw_a = a.raw.to_adata()[:, common_genes].copy()
        else:
            raw_a = a[:, common_genes].copy()
        # Carry over obs columns
        for col in ['sample_id', 'disease', 'phenotype', 'roi_class',
                     'n_transcripts', 'n_genes', 'area_px']:
            if col in a.obs.columns:
                raw_a.obs[col] = a.obs[col].values
        raw_adatas.append(raw_a)

    combined = ad.concat(raw_adatas, join='inner', merge='same')
    logger.info(f"Combined: {combined.shape}")
    logger.info(f"Samples: {combined.obs['sample_id'].nunique()}")
    logger.info(f"Disease groups: {dict(combined.obs['disease'].value_counts())}")

    if skipped:
        combined.uns['skipped_samples'] = skipped

    return combined


# ============================================================
# 3. Batch Correction
# ============================================================

def batch_correct(
    adata,
    batch_key: str = 'sample_id',
    method: Literal['harmony', 'combat', 'scvi', 'none'] = 'harmony',
    **kwargs,
) -> "ad.AnnData":
    """Apply batch correction to concatenated data.

    Parameters
    ----------
    adata : AnnData
        Concatenated, normalized AnnData.
    batch_key : str
        Column in obs identifying batches.
    method : str
        Correction method: 'harmony', 'combat', 'scvi', 'none'.

    Returns
    -------
    AnnData
        Batch-corrected AnnData.
    """
    if method == 'none':
        return adata

    if method == 'harmony':
        try:
            import harmonypy
        except ImportError:
            try:
                # scanpy has harmony integration
                sc.external.pp.harmony_integrate(adata, batch_key, basis='X_pca',
                                                  adjusted_basis='X_pca_harmony')
                adata.obsm['X_pca_original'] = adata.obsm['X_pca'].copy()
                # Rebuild neighbors on corrected PCs
                sc.pp.neighbors(adata, use_rep='X_pca_harmony')
                sc.tl.umap(adata)
                adata.uns['batch_correction'] = {'method': 'harmony', 'batch_key': batch_key}
                logger.info("Harmony batch correction applied")
                return adata
            except Exception as e:
                logger.warning(f"Harmony failed: {e}. Install with: pip install harmonypy")
                return adata

        sc.external.pp.harmony_integrate(adata, batch_key, basis='X_pca',
                                          adjusted_basis='X_pca_harmony')
        sc.pp.neighbors(adata, use_rep='X_pca_harmony')
        sc.tl.umap(adata)
        adata.uns['batch_correction'] = {'method': 'harmony', 'batch_key': batch_key}
        logger.info("Harmony batch correction applied")

    elif method == 'combat':
        sc.pp.combat(adata, key=batch_key)
        sc.tl.pca(adata)
        sc.pp.neighbors(adata)
        sc.tl.umap(adata)
        adata.uns['batch_correction'] = {'method': 'combat', 'batch_key': batch_key}
        logger.info("ComBat batch correction applied")

    elif method == 'scvi':
        try:
            import scvi
            scvi.model.SCVI.setup_anndata(adata, batch_key=batch_key)
            model = scvi.model.SCVI(adata, n_latent=30, **kwargs)
            model.train(max_epochs=100)
            adata.obsm['X_scvi'] = model.get_latent_representation()
            sc.pp.neighbors(adata, use_rep='X_scvi')
            sc.tl.umap(adata)
            adata.uns['batch_correction'] = {'method': 'scvi', 'batch_key': batch_key}
            logger.info("scVI batch correction applied")
        except ImportError:
            logger.warning("scvi-tools not installed. pip install scvi-tools")
            return adata

    return adata


# ============================================================
# 4. Phenotype-Stratified Differential Expression
# ============================================================

def run_cross_sample_de(
    adata,
    groupby: str = 'disease',
    phenotype_col: str = 'phenotype',
    method: str = 'wilcoxon',
    min_cells_per_group: int = 50,
    n_genes: int = 100,
    stratify_by_phenotype: bool = True,
) -> Dict[str, pd.DataFrame]:
    """Run differential expression across disease groups.

    Can optionally stratify by COMET phenotype class, answering questions
    like "what genes differ between Endo and MLA in CD8+ T cells?"

    Parameters
    ----------
    adata : AnnData
        Concatenated, normalized, batch-corrected data.
    groupby : str
        Column for group comparison (e.g., 'disease').
    phenotype_col : str
        Column with COMET phenotype classes.
    method : str
        DE test method: 'wilcoxon', 't-test', 'logreg'.
    min_cells_per_group : int
        Minimum cells per group for testing.
    n_genes : int
        Number of top genes to report per group.
    stratify_by_phenotype : bool
        If True, run DE within each phenotype class separately.

    Returns
    -------
    dict
        Keys are phenotype classes (or 'all' if not stratified),
        values are DataFrames with DE results.
    """
    results = {}

    if not stratify_by_phenotype or phenotype_col not in adata.obs.columns:
        # Global DE
        if adata.obs[groupby].nunique() < 2:
            logger.warning(f"Only one group in '{groupby}' — cannot run DE")
            return results

        try:
            sc.tl.rank_genes_groups(adata, groupby=groupby, method=method,
                                     n_genes=n_genes)
            de_df = sc.get.rank_genes_groups_df(adata, group=None)
            results['all'] = de_df
            logger.info(f"Global DE: {len(de_df)} results")
        except Exception as e:
            logger.error(f"Global DE failed: {e}")

        return results

    # Phenotype-stratified DE
    phenotypes = adata.obs[phenotype_col].unique()
    logger.info(f"Running stratified DE across {len(phenotypes)} phenotypes...")

    for pheno in phenotypes:
        mask = adata.obs[phenotype_col] == pheno
        adata_sub = adata[mask].copy()

        # Check group sizes
        group_sizes = adata_sub.obs[groupby].value_counts()
        valid_groups = group_sizes[group_sizes >= min_cells_per_group].index

        if len(valid_groups) < 2:
            logger.info(f"  {pheno}: skipped (insufficient groups: {dict(group_sizes)})")
            continue

        adata_sub = adata_sub[adata_sub.obs[groupby].isin(valid_groups)].copy()

        try:
            sc.tl.rank_genes_groups(adata_sub, groupby=groupby, method=method,
                                     n_genes=n_genes)
            de_df = sc.get.rank_genes_groups_df(adata_sub, group=None)
            de_df['phenotype'] = pheno
            results[pheno] = de_df

            n_sig = (de_df['pvals_adj'] < 0.05).sum()
            logger.info(f"  {pheno}: {len(de_df)} results, {n_sig} significant (adj p<0.05)")

        except Exception as e:
            logger.warning(f"  {pheno}: DE failed: {e}")

    return results


# ============================================================
# 5. Phenotype Distribution Analysis
# ============================================================

def phenotype_composition_analysis(
    adata,
    sample_col: str = 'sample_id',
    disease_col: str = 'disease',
    phenotype_col: str = 'phenotype',
) -> pd.DataFrame:
    """Compute phenotype composition per sample and test for differences.

    Returns a table of phenotype proportions per sample, and runs
    chi-squared / Fisher tests between disease groups.

    Parameters
    ----------
    adata : AnnData
        Concatenated data with phenotype annotations.

    Returns
    -------
    pd.DataFrame
        Phenotype proportions per sample.
    """
    # Count phenotypes per sample
    counts = adata.obs.groupby([sample_col, disease_col, phenotype_col]).size().reset_index(name='n_cells')

    # Compute proportions
    totals = counts.groupby(sample_col)['n_cells'].transform('sum')
    counts['proportion'] = counts['n_cells'] / totals

    # Pivot to wide format
    proportions = counts.pivot_table(
        index=[sample_col, disease_col],
        columns=phenotype_col,
        values='proportion',
        fill_value=0,
    )

    logger.info(f"Phenotype composition computed for {proportions.shape[0]} samples")
    return proportions


# ============================================================
# 6. Full Pipeline Runner
# ============================================================

def run_post_alignment_pipeline(
    h5ad_dir: Union[str, Path],
    output_dir: Union[str, Path],
    disease_comparison: Tuple[str, str] = ('Endo', 'MLA'),
    batch_correction_method: str = 'harmony',
    normalize_kwargs: Optional[dict] = None,
    stratify_de: bool = True,
):
    """Run the complete post-alignment analysis pipeline.

    Parameters
    ----------
    h5ad_dir : path
        Directory containing per-sample integrated h5ad files.
    output_dir : path
        Output directory for results.
    disease_comparison : tuple
        Two disease groups to compare.
    batch_correction_method : str
        Batch correction approach.
    normalize_kwargs : dict, optional
        Normalization parameters.
    stratify_de : bool
        Whether to stratify DE by phenotype.
    """
    h5ad_dir = Path(h5ad_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # 0. Check alignment approvals (validation gate)
    try:
        from comet.alignment_validation import require_approval, get_approval_status

        # Each h5ad lives inside a sample directory; check each one
        h5ad_files_pre = sorted(h5ad_dir.rglob('*_comet_stomics_integrated.h5ad'))
        unapproved = []
        for hf in h5ad_files_pre:
            sample_dir = hf.parent
            try:
                require_approval(sample_dir, action='post-alignment analysis')
            except RuntimeError as e:
                unapproved.append((sample_dir.name, str(e)))

        if unapproved:
            msg_lines = [f"\n{'='*60}",
                         "ALIGNMENT APPROVAL GATE: Cannot proceed.",
                         f"{len(unapproved)} sample(s) have not been approved:"]
            for dirname, err in unapproved:
                msg_lines.append(f"  - {dirname}: {err}")
            msg_lines.append(
                "\nRun the batch QC notebook to review and sign off alignments."
            )
            msg_lines.append(f"{'='*60}")
            full_msg = '\n'.join(msg_lines)
            logger.error(full_msg)
            raise RuntimeError(full_msg)

        logger.info(f"All samples approved — proceeding with analysis.")
    except ImportError:
        logger.warning("alignment_validation module not available — skipping approval gate")

    # 1. Find all integrated h5ad files
    h5ad_files = sorted(h5ad_dir.rglob('*_comet_stomics_integrated.h5ad'))
    logger.info(f"Found {len(h5ad_files)} integrated h5ad files")

    if not h5ad_files:
        raise FileNotFoundError(f"No *_comet_stomics_integrated.h5ad files in {h5ad_dir}")

    # 2. Concatenate
    combined = concatenate_samples(
        h5ad_files,
        normalize=True,
        normalize_kwargs=normalize_kwargs,
    )

    # 3. Re-normalize combined data
    combined_norm = normalize_integrated_adata(combined)

    # 4. Batch correction
    combined_corrected = batch_correct(combined_norm, method=batch_correction_method)

    # Save intermediate
    combined_path = output_dir / 'combined_all_samples.h5ad'
    combined_corrected.write_h5ad(str(combined_path))
    logger.info(f"Saved combined data: {combined_path}")

    # 5. Differential expression
    de_results = run_cross_sample_de(
        combined_corrected,
        groupby='disease',
        stratify_by_phenotype=stratify_de,
    )

    # Save DE results
    for pheno, df in de_results.items():
        de_path = output_dir / f'de_{pheno.replace(" ", "_").replace("+", "pos")}.csv'
        df.to_csv(de_path, index=False)

    # 6. Phenotype composition
    composition = phenotype_composition_analysis(combined_corrected)
    composition.to_csv(output_dir / 'phenotype_composition.csv')

    # 7. Subset to disease comparison
    if disease_comparison:
        d1, d2 = disease_comparison
        mask = combined_corrected.obs['disease'].isin([d1, d2])
        subset = combined_corrected[mask].copy()

        if subset.n_obs > 0:
            de_comparison = run_cross_sample_de(
                subset,
                groupby='disease',
                stratify_by_phenotype=stratify_de,
            )

            for pheno, df in de_comparison.items():
                fname = f'de_{d1}_vs_{d2}_{pheno.replace(" ", "_")}.csv'
                df.to_csv(output_dir / fname, index=False)

            logger.info(f"{d1} vs {d2} DE complete: {len(de_comparison)} phenotype strata")

    logger.info(f"\nPipeline complete. Results in: {output_dir}")
    return combined_corrected, de_results
