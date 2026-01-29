"""I/O utilities for MSI data.

This module provides functions for reading and writing MSI data in various
formats including AnnData (.h5ad) and parquet/CSV.
"""

import logging
from pathlib import Path
from typing import Optional, Union, List, Dict

import numpy as np
import pandas as pd

try:
    import anndata as ad
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False

logger = logging.getLogger(__name__)


def load_parquet_channels(
    input_dir: Union[str, Path],
    sample_id: Optional[str] = None,
) -> Dict[str, pd.DataFrame]:
    """Load all parquet channel files from a directory.

    Parameters
    ----------
    input_dir : str or Path
        Directory containing parquet files.
    sample_id : str, optional
        Sample ID subdirectory to load from.

    Returns
    -------
    channels : dict
        Dictionary mapping channel names to DataFrames.
    """
    input_dir = Path(input_dir)

    if sample_id:
        input_dir = input_dir / sample_id

    parquet_files = list(input_dir.glob("*.parquet"))
    logger.info(f"Found {len(parquet_files)} parquet files in {input_dir}")

    channels = {}
    for filepath in parquet_files:
        channel_name = filepath.stem
        df = pd.read_parquet(filepath)
        channels[channel_name] = df
        logger.debug(f"Loaded {channel_name}: {len(df)} rows")

    return channels


def merge_channels_to_matrix(
    channels: Dict[str, pd.DataFrame],
    intensity_col: Optional[str] = None,
    join: str = 'outer',
    fill_value: float = 0.0,
) -> pd.DataFrame:
    """Merge multiple channel DataFrames into a single matrix.

    Parameters
    ----------
    channels : dict
        Dictionary mapping channel names to DataFrames with 'x', 'y', and intensity columns.
    intensity_col : str, optional
        Name of intensity column. If None, uses the third column.
    join : str, default='outer'
        Join type: 'outer' (union of pixels), 'inner' (intersection).
    fill_value : float, default=0.0
        Value for missing intensities.

    Returns
    -------
    merged : pd.DataFrame
        DataFrame with columns ['x', 'y'] + channel names.
    """
    if not channels:
        raise ValueError("No channels provided")

    # Process each channel
    dfs = []
    channel_names = []

    for name, df in channels.items():
        df = df.copy()

        # Determine intensity column
        if intensity_col and intensity_col in df.columns:
            int_col = intensity_col
        else:
            # Use the column that's not x or y
            other_cols = [c for c in df.columns if c not in ('x', 'y')]
            if not other_cols:
                raise ValueError(f"No intensity column found in channel {name}")
            int_col = other_cols[0]

        # Rename intensity column to channel name
        df = df[['x', 'y', int_col]].rename(columns={int_col: name})
        dfs.append(df)
        channel_names.append(name)

    # Merge all dataframes
    merged = dfs[0]
    for df in dfs[1:]:
        merged = pd.merge(merged, df, on=['x', 'y'], how=join)

    # Fill missing values
    merged = merged.fillna(fill_value)

    logger.info(f"Merged {len(channel_names)} channels: {merged.shape[0]} pixels x {len(channel_names)} channels")

    return merged


def create_anndata_from_merged(
    merged_df: pd.DataFrame,
    sample_id: Optional[str] = None,
    modality: Optional[str] = None,
) -> 'ad.AnnData':
    """Create an AnnData object from merged channel data.

    Parameters
    ----------
    merged_df : pd.DataFrame
        DataFrame with 'x', 'y' and channel columns.
    sample_id : str, optional
        Sample identifier.
    modality : str, optional
        Data modality (e.g., 'glycans', 'metabolites', 'peptides').

    Returns
    -------
    adata : AnnData
        Annotated data matrix.
    """
    if not HAS_ANNDATA:
        raise ImportError("anndata is required. Install with: pip install anndata")

    # Separate coordinates from intensities
    coord_cols = ['x', 'y']
    channel_cols = [c for c in merged_df.columns if c not in coord_cols]

    # Create intensity matrix
    X = merged_df[channel_cols].values.astype(np.float32)

    # Create observation metadata
    obs = pd.DataFrame({
        'x': merged_df['x'].values,
        'y': merged_df['y'].values,
    })

    if sample_id:
        obs['sample_id'] = sample_id

    # Create variable metadata
    var = pd.DataFrame(index=channel_cols)
    var['channel_name'] = channel_cols

    if modality:
        var['modality'] = modality

    # Create AnnData
    adata = ad.AnnData(
        X=X,
        obs=obs,
        var=var,
    )

    # Add spatial coordinates
    adata.obsm['spatial'] = merged_df[['x', 'y']].values.astype(np.float32)

    # Add metadata
    adata.uns['msi_metadata'] = {
        'sample_id': sample_id,
        'modality': modality,
        'n_pixels': X.shape[0],
        'n_channels': X.shape[1],
    }

    logger.info(f"Created AnnData: {adata.n_obs} pixels x {adata.n_vars} channels")

    return adata


def save_anndata(
    adata: 'ad.AnnData',
    filepath: Union[str, Path],
    compression: str = 'gzip',
) -> None:
    """Save an AnnData object to h5ad format.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    filepath : str or Path
        Output file path (.h5ad).
    compression : str, default='gzip'
        Compression method.
    """
    filepath = Path(filepath)
    filepath.parent.mkdir(parents=True, exist_ok=True)

    adata.write_h5ad(filepath, compression=compression)
    logger.info(f"Saved AnnData to {filepath}")


def load_anndata(filepath: Union[str, Path]) -> 'ad.AnnData':
    """Load an AnnData object from h5ad format.

    Parameters
    ----------
    filepath : str or Path
        Input file path (.h5ad).

    Returns
    -------
    adata : AnnData
        Loaded annotated data matrix.
    """
    if not HAS_ANNDATA:
        raise ImportError("anndata is required")

    adata = ad.read_h5ad(filepath)
    logger.info(f"Loaded AnnData from {filepath}: {adata.n_obs} x {adata.n_vars}")

    return adata


def anndata_to_parquet(
    adata: 'ad.AnnData',
    output_dir: Union[str, Path],
    sample_id: str,
    include_coordinates: bool = True,
    include_clusters: bool = True,
) -> List[str]:
    """Export AnnData to parquet format.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    output_dir : str or Path
        Output directory.
    sample_id : str
        Sample identifier.
    include_coordinates : bool, default=True
        Include x, y coordinates.
    include_clusters : bool, default=True
        Include cluster assignments if available.

    Returns
    -------
    saved_files : list
        List of saved file paths.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved_files = []

    # Export intensity matrix
    df = pd.DataFrame(
        adata.X,
        columns=adata.var_names,
    )

    if include_coordinates and 'spatial' in adata.obsm:
        df['x'] = adata.obsm['spatial'][:, 0]
        df['y'] = adata.obsm['spatial'][:, 1]

    if include_clusters and 'leiden' in adata.obs.columns:
        df['cluster'] = adata.obs['leiden'].values

    intensities_path = output_dir / f'{sample_id}_intensities.parquet'
    df.to_parquet(intensities_path, index=False)
    saved_files.append(str(intensities_path))

    # Export observation metadata
    obs_path = output_dir / f'{sample_id}_obs.parquet'
    adata.obs.to_parquet(obs_path)
    saved_files.append(str(obs_path))

    # Export variable metadata
    var_path = output_dir / f'{sample_id}_var.parquet'
    adata.var.to_parquet(var_path)
    saved_files.append(str(var_path))

    logger.info(f"Exported {len(saved_files)} parquet files to {output_dir}")

    return saved_files


def export_for_visualization(
    adata: 'ad.AnnData',
    output_path: Union[str, Path],
    channels: Optional[List[str]] = None,
    include_umap: bool = True,
) -> None:
    """Export data for visualization tools.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with clustering results.
    output_path : str or Path
        Output CSV file path.
    channels : list, optional
        Specific channels to include. If None, includes top 10 by variance.
    include_umap : bool, default=True
        Include UMAP coordinates.
    """
    df = pd.DataFrame()

    # Add spatial coordinates
    if 'spatial' in adata.obsm:
        df['x'] = adata.obsm['spatial'][:, 0]
        df['y'] = adata.obsm['spatial'][:, 1]

    # Add cluster assignments
    if 'leiden' in adata.obs.columns:
        df['cluster'] = adata.obs['leiden'].values

    # Add UMAP coordinates
    if include_umap and 'X_umap' in adata.obsm:
        df['umap_1'] = adata.obsm['X_umap'][:, 0]
        df['umap_2'] = adata.obsm['X_umap'][:, 1]

    # Add channel intensities
    if channels is None:
        # Select top channels by variance
        variances = np.var(adata.X, axis=0)
        top_idx = np.argsort(variances)[-10:][::-1]
        channels = [adata.var_names[i] for i in top_idx]

    for channel in channels:
        if channel in adata.var_names:
            idx = list(adata.var_names).index(channel)
            df[channel] = adata.X[:, idx]

    # Save
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, index=False)

    logger.info(f"Exported visualization data to {output_path}")


def list_samples(input_dir: Union[str, Path], pattern: str = "*.h5ad") -> List[str]:
    """List available sample files in a directory.

    Parameters
    ----------
    input_dir : str or Path
        Directory to search.
    pattern : str, default="*.h5ad"
        Glob pattern for files.

    Returns
    -------
    samples : list
        List of sample IDs (filenames without extension).
    """
    input_dir = Path(input_dir)
    files = list(input_dir.glob(pattern))

    samples = [f.stem for f in files]
    logger.info(f"Found {len(samples)} samples in {input_dir}")

    return samples
