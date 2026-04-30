# Consolidated Action Plan: COMET-STOmics-MSI Pipeline Completion

**Date:** 2026-04-29
**Status:** 27/41 samples processed, key gaps identified

---

## Status Reconciliation: What the Tracker Says vs. What's in This Repo

The status tracker (2026-04-24) references several files that exist only on the workstation's `feat/enhanced-phenotype-roi-microbiome` branch and have NOT been pushed to this repo:

| Referenced File | In This Repo? | Notes |
|---|---|---|
| `comet/registration.py` | NO | Multi-method registration (QuPath+VALIS, QuPath+ECC) |
| `config/qupath_transforms.json` | NO | Manual affine seeds for 29 samples |
| `script08_alignment_qc_snapshots.py` | NO | 4-panel zoomed QC visualization |
| `script09_*` | NO | Unknown |
| `script10_cellbin_comparison_zoomed.py` | NO | Tight-zoom 3-panel comparison |

**Action required:** Commit and push (or merge) the `feat/enhanced-phenotype-roi-microbiome` branch so all code is in one place. This is blocking collaboration and reproducibility.

---

## Gap Analysis: What Remains

### Tier 1: Blocking Progress (Do First)

#### A. Process Remaining 14 Samples
- **SO32** (Endo), **SO33** (MLA) — directories exist, registration incomplete
- **SO57** — directory exists but empty
- **SO8, SO11, SO40** (Ovarian) — have cellbin, need QuPath transforms + alignment
- **SO36, SO39, SO67** (Ovarian) — missing ALL raw data (no FASTQ/mask/cellbin)
- **SO64** (USC) — not started, has cellbin
- Requires: QuPath manual affine seeds → registration → transcript assignment

#### B. MSI → H&E Registration (NEW)
- **Current approach:** MSI composite → DAPI (script07_registration.ipynb)
- **Better approach:** MSI composite → H&E (morphologically most similar)
- **New code:** `msi/utils/registration.py` — MSItoHERegistrar class with:
  - Metabolite subset selection (tissue_contrast / variance / spatial / manual)
  - Composite image creation (PCA / SNR / MIP / weighted_sum)
  - VALIS registration with fallback cascade
  - Batch processing support
  - QC visualization

#### C. Merge Feature Branch
- 12 modified + 47 untracked files on workstation
- `registration.py`, QuPath transforms, scripts 08-10 only exist there
- Critical for anyone else to reproduce or extend the work

### Tier 2: Analysis Pipeline (Do After Alignment)

#### D. Post-Alignment Normalization (NEW)
- **New code:** `comet/post_alignment_analysis.py` with:
  - Per-sample normalization (library size, log, HVG, PCA, UMAP)
  - Spatial-aware QC (min genes, min transcripts, mito filtering)
  - Cross-sample concatenation with common gene intersection
  - Batch correction (Harmony / ComBat / scVI)
  - Phenotype-stratified differential expression
  - Phenotype composition analysis across disease groups

#### E. Cross-Sample Analysis
- Concatenate all 27+ integrated AnnData objects
- Endo vs MLA comparison (primary analysis goal)
- Phenotype distribution analysis across disease groups
- Protein-gene correlation validation at scale

#### F. Cellbin Comparison Scale-Up
- Currently only SO4 and SO34 have cellbin comparison results
- Need to run for all 25 samples with cellbin masks
- Expression concordance analysis (per-cell Pearson r)

### Tier 3: Integration (Final)

#### G. MSI-COMET-STOmics Triple Integration
- Chain: MSI → H&E → (existing H&E → COMET DAPI → STOmics DAPI)
- Assign MSI features to COMET-defined cells
- Multi-modal AnnData: genes + proteins + metabolites per cell

#### H. QC and Documentation
- QC snapshots for remaining 6 processed samples
- Publication-quality figures
- Methods documentation

---

## New Code Created in This Session

### 1. `msi/utils/registration.py` — MSI-to-H&E Registration

Key classes and functions:

- **`MSItoHERegistrar`** — Main registration class
  - `select_channels()` — Auto-select metabolites best for registration
  - `create_composite()` — Build composite image from selected channels
  - `register()` — Run VALIS with fallback cascade (PCA → SNR → MIP)
  - `transform_coordinates()` — Warp MSI coordinates to H&E space
  - `transform_adata()` — Apply transform to AnnData, store in `obsm['spatial_he']`
  - `plot_qc()` — 4-panel QC visualization

- **`select_metabolite_subset()`** — Channel selection strategies:
  - `tissue_contrast`: Channels with highest tissue-vs-background difference
  - `variance`: Highest spatial variance
  - `spatial`: Highest Moran's I
  - `manual`: User-specified channel list

- **`batch_register_msi_to_he()`** — Process multiple samples

**Usage example:**
```python
from msi.utils.registration import MSItoHERegistrar

reg = MSItoHERegistrar(
    msi_adata=adata_glycans,
    he_image_path='T:/Sammy Data/HE/SO24_HE.ndpi',
    output_dir='T:/Sammy Data/projects/out/msi_registration/SO24',
    channel_selection_method='tissue_contrast',
    composite_method='pca',
    n_channels=15,
)

# Register with automatic fallback
result = reg.register()

# Apply to all MSI modalities (same coordinate space)
adata_glycans = reg.transform_adata(adata_glycans)
adata_peptides = reg.transform_adata(adata_peptides)
adata_metabolites = reg.transform_adata(adata_metabolites)

# QC
reg.plot_qc()
```

### 2. `comet/post_alignment_analysis.py` — Post-Alignment Pipeline

Key functions:

- **`normalize_integrated_adata()`** — Per-sample normalization
  - QC metrics (mito%, gene count, transcript count)
  - Library size normalization → log1p → HVG → PCA → UMAP
  - Configurable filtering thresholds

- **`concatenate_samples()`** — Cross-sample integration
  - Loads multiple h5ad files
  - Finds common gene intersection
  - Preserves phenotype/disease metadata

- **`batch_correct()`** — Batch effect removal
  - Harmony (default, fast, works well for spatial)
  - ComBat (parametric, good for known batch structure)
  - scVI (deep learning, best for complex batch effects)

- **`run_cross_sample_de()`** — Differential expression
  - Global DE: all cells, disease group comparison
  - **Phenotype-stratified DE**: "What genes differ between Endo and MLA specifically in CD8+ T cells?" — runs DE within each of the 24 COMET phenotype classes

- **`phenotype_composition_analysis()`** — Cell type proportions per sample

- **`run_post_alignment_pipeline()`** — End-to-end runner

**Usage example:**
```python
from comet.post_alignment_analysis import run_post_alignment_pipeline

combined, de_results = run_post_alignment_pipeline(
    h5ad_dir='T:/Sammy Data/projects/out/comet_stomics_alignment',
    output_dir='T:/Sammy Data/projects/out/cross_sample_analysis',
    disease_comparison=('Endo', 'MLA'),
    batch_correction_method='harmony',
    stratify_de=True,
)
```

---

## VALIS Robustness: Design Decisions

The existing pipeline already has strong robustness features (QuPath affine seed → VALIS refinement, ECC fallback). The new MSI registration adds:

1. **Fallback cascade**: If PCA composite fails VALIS, automatically tries SNR, then MIP
2. **Channel selection**: Poor channels (uniform, noisy) excluded before composite creation
3. **Grayscale H&E inversion**: H&E images are inverted so tissue is bright (matching DAPI/MSI contrast)
4. **Error-based method selection**: The composite method with lowest VALIS error wins
5. **Affine seed support**: Manual transforms can be provided for difficult cases
6. **QC gating**: Registration results with error > threshold trigger fallback

For the COMET-STOmics pipeline, the QuPath manual affine + VALIS refinement approach is already the gold standard. The main improvement needed is: **ensure all 41 QuPath transforms are saved in `config/qupath_transforms.json`** and committed.

---

## Execution Order

```
Week 1:
  [x] Merge feature branch (push registration.py, scripts 08-10, transforms)
  [ ] Run remaining 14 samples through COMET-STOmics alignment
      (need QuPath transforms for SO8, SO11, SO32, SO33, SO40, SO57, SO64)
  [ ] MSI → H&E registration for 6 ovarian samples (SO24, SO32-34, SO44-45)

Week 2:
  [ ] Post-alignment normalization on all 27+ samples
  [ ] Cross-sample concatenation + batch correction
  [ ] Endo vs MLA differential expression (global + phenotype-stratified)
  [ ] Cellbin comparison scale-up (remaining 23 samples)

Week 3:
  [ ] MSI-COMET-STOmics triple integration (6 samples with MSI)
  [ ] Phenotype composition analysis
  [ ] Publication figures and methods documentation
```

---

## Files in This Repo (Updated 2026-04-30)

| File | Purpose | Status |
|------|---------|--------|
| `comet/alignment_utils.py` | Core COMET-STOmics alignment (now creates validation manifest on completion) | Production |
| `comet/alignment_validation.py` | **NEW** Validation & signoff system — file fingerprinting, QC visualization, approval gates | Ready |
| `comet/run_batch_alignment.py` | CLI batch runner for alignment | Ready |
| `comet/post_alignment_analysis.py` | Normalization + DE + batch correction (now gates on alignment approval) | Ready |
| `comet/script06_alignment_qc_signoff.ipynb` | **NEW** Interactive QC review & signoff notebook | Ready |
| `msi/utils/registration.py` | MSI→H&E registration | Ready |
| `msi/utils/normalization.py` | MSI-specific normalization | Production |
| `msi/utils/extraction.py` | OME-TIFF intensity extraction | Production |
| `msi/utils/spatial.py` | Moran's I, spatial weights | Production |
| `msi/utils/clustering.py` | PCA, UMAP, Leiden | Production |
| `msi/config/default_params.yaml` | Pipeline configuration | Production |
| `stomics/saw_data_audit.py` | STOmics data inventory audit | Complete |
| `CONSOLIDATED_ACTION_PLAN.md` | This document | Current |

### Validation System Architecture

The alignment validation system prevents the "wrong warp file" problem through three interlocking mechanisms:

1. **Automatic manifest creation** — `run_alignment_pipeline()` now creates a `validation_manifest.json` at the end of every alignment run, recording SHA-256 hashes of all input/output files.

2. **Approval gate** — `run_post_alignment_pipeline()` checks every sample's manifest before proceeding. If any sample is pending, rejected, or has files that changed since approval, the pipeline stops with a clear error message.

3. **Interactive review** — `script06_alignment_qc_signoff.ipynb` provides a structured walkthrough: status dashboard → QC image generation → per-sample review with global overview + ROI zoom panels → sign-off/reject with reason → final status confirmation.
