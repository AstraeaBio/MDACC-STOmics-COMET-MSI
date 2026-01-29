"""Spatial statistics utilities for MSI data.

This module provides functions for spatial autocorrelation analysis,
including Moran's I calculation and spatial weights construction.
"""

import logging
from typing import Optional, Tuple, Union, Literal

import numpy as np
import pandas as pd
from scipy import spatial

try:
    from libpysal import weights as pysal_weights
    from esda.moran import Moran
    HAS_SPATIAL_LIBS = True
except ImportError:
    HAS_SPATIAL_LIBS = False

try:
    import anndata as ad
except ImportError:
    ad = None

logger = logging.getLogger(__name__)


def create_knn_weights(
    coordinates: np.ndarray,
    k: int = 6,
    ids: Optional[np.ndarray] = None,
) -> 'pysal_weights.KNN':
    """Create K-nearest neighbors spatial weights matrix.

    Parameters
    ----------
    coordinates : np.ndarray
        Array of shape (n, 2) with x, y coordinates.
    k : int, default=6
        Number of nearest neighbors.
    ids : np.ndarray, optional
        Array of observation IDs.

    Returns
    -------
    w : KNN
        PySAL KNN weights object.

    Raises
    ------
    ImportError
        If libpysal is not installed.
    """
    if not HAS_SPATIAL_LIBS:
        raise ImportError("libpysal is required for spatial weights. "
                         "Install with: pip install libpysal")

    if ids is not None:
        w = pysal_weights.KNN.from_array(coordinates, k=k, ids=ids)
    else:
        w = pysal_weights.KNN.from_array(coordinates, k=k)

    return w


def create_distance_weights(
    coordinates: np.ndarray,
    threshold: float,
    binary: bool = True,
    ids: Optional[np.ndarray] = None,
) -> 'pysal_weights.DistanceBand':
    """Create distance-based spatial weights matrix.

    Parameters
    ----------
    coordinates : np.ndarray
        Array of shape (n, 2) with x, y coordinates.
    threshold : float
        Distance threshold for neighbors.
    binary : bool, default=True
        If True, weights are binary (0/1). If False, weights are 1/distance.
    ids : np.ndarray, optional
        Array of observation IDs.

    Returns
    -------
    w : DistanceBand
        PySAL DistanceBand weights object.
    """
    if not HAS_SPATIAL_LIBS:
        raise ImportError("libpysal is required for spatial weights.")

    w = pysal_weights.DistanceBand.from_array(
        coordinates,
        threshold=threshold,
        binary=binary,
        ids=ids,
    )

    return w


def calculate_morans_i(
    values: np.ndarray,
    weights: 'pysal_weights.W',
    permutations: int = 999,
) -> dict:
    """Calculate Moran's I statistic for spatial autocorrelation.

    Parameters
    ----------
    values : np.ndarray
        Array of values for which to calculate spatial autocorrelation.
    weights : W
        PySAL weights object defining spatial relationships.
    permutations : int, default=999
        Number of permutations for significance testing.

    Returns
    -------
    result : dict
        Dictionary containing:
        - 'I': Moran's I statistic
        - 'p_value': p-value from permutation test
        - 'z_score': z-score
        - 'expected_I': expected value under null hypothesis
    """
    if not HAS_SPATIAL_LIBS:
        raise ImportError("esda is required for Moran's I. "
                         "Install with: pip install esda")

    moran = Moran(values, weights, permutations=permutations)

    return {
        'I': moran.I,
        'p_value': moran.p_sim,
        'z_score': moran.z_sim,
        'expected_I': moran.EI,
    }


def morans_i_per_channel(
    data: np.ndarray,
    coordinates: np.ndarray,
    k: int = 6,
    permutations: int = 999,
    channel_names: Optional[list] = None,
    n_jobs: int = 1,
) -> pd.DataFrame:
    """Calculate Moran's I for each channel/feature in the data.

    Parameters
    ----------
    data : np.ndarray
        Data matrix with shape (n_observations, n_channels).
    coordinates : np.ndarray
        Spatial coordinates with shape (n_observations, 2).
    k : int, default=6
        Number of nearest neighbors for weights.
    permutations : int, default=999
        Number of permutations for significance testing.
    channel_names : list, optional
        Names for each channel.
    n_jobs : int, default=1
        Number of parallel jobs. Currently only supports 1.

    Returns
    -------
    results : pd.DataFrame
        DataFrame with Moran's I statistics for each channel.
    """
    n_channels = data.shape[1]

    if channel_names is None:
        channel_names = [f'channel_{i}' for i in range(n_channels)]

    # Create spatial weights
    logger.info(f"Creating KNN weights (k={k})...")
    w = create_knn_weights(coordinates, k=k)

    results = []
    for i, name in enumerate(channel_names):
        logger.info(f"Computing Moran's I for {name} ({i+1}/{n_channels})...")
        values = data[:, i]

        # Skip if all values are the same (no variance)
        if np.std(values) == 0:
            logger.warning(f"Skipping {name}: no variance")
            results.append({
                'channel': name,
                'I': np.nan,
                'p_value': np.nan,
                'z_score': np.nan,
                'expected_I': np.nan,
            })
            continue

        try:
            mi_result = calculate_morans_i(values, w, permutations=permutations)
            mi_result['channel'] = name
            results.append(mi_result)
        except Exception as e:
            logger.error(f"Error computing Moran's I for {name}: {e}")
            results.append({
                'channel': name,
                'I': np.nan,
                'p_value': np.nan,
                'z_score': np.nan,
                'expected_I': np.nan,
            })

    return pd.DataFrame(results)


def filter_spatially_coherent(
    morans_results: pd.DataFrame,
    threshold: float = 0.2,
    p_value_threshold: float = 0.05,
) -> pd.DataFrame:
    """Filter channels based on spatial coherence (Moran's I).

    Parameters
    ----------
    morans_results : pd.DataFrame
        Results from morans_i_per_channel().
    threshold : float, default=0.2
        Minimum Moran's I threshold for inclusion.
    p_value_threshold : float, default=0.05
        Maximum p-value for significance.

    Returns
    -------
    filtered : pd.DataFrame
        Filtered results containing only spatially coherent channels.
    """
    mask = (
        (morans_results['I'] >= threshold) &
        (morans_results['p_value'] <= p_value_threshold)
    )

    filtered = morans_results[mask].copy()
    logger.info(f"Filtered to {len(filtered)} spatially coherent channels "
               f"(threshold={threshold}, p<{p_value_threshold})")

    return filtered


def get_modality_threshold(modality: str) -> float:
    """Get recommended Moran's I threshold for a given modality.

    Parameters
    ----------
    modality : str
        Data modality: 'glycans', 'metabolites', or 'peptides'.

    Returns
    -------
    threshold : float
        Recommended Moran's I threshold.
    """
    thresholds = {
        'glycans': 0.2,
        'metabolites': 0.05,
        'peptides': 0.2,
    }

    if modality.lower() not in thresholds:
        logger.warning(f"Unknown modality '{modality}', using default threshold 0.1")
        return 0.1

    return thresholds[modality.lower()]


def compute_spatial_stats_anndata(
    adata: 'ad.AnnData',
    k: int = 6,
    permutations: int = 999,
    layer: Optional[str] = None,
    key_added: str = 'morans_i',
) -> 'ad.AnnData':
    """Compute Moran's I for all features in an AnnData object.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix with spatial coordinates in .obsm['spatial'].
    k : int, default=6
        Number of nearest neighbors for weights.
    permutations : int, default=999
        Number of permutations for significance testing.
    layer : str, optional
        Layer to use. If None, uses .X
    key_added : str, default='morans_i'
        Key for storing results in .var and .uns.

    Returns
    -------
    adata : AnnData
        Modified AnnData with Moran's I results in .var[key_added].
    """
    if ad is None:
        raise ImportError("anndata is required")

    # Get coordinates
    if 'spatial' not in adata.obsm:
        raise ValueError("AnnData must have spatial coordinates in .obsm['spatial']")

    coordinates = adata.obsm['spatial']

    # Get data
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    if hasattr(X, 'toarray'):
        X = X.toarray()

    # Compute Moran's I
    results = morans_i_per_channel(
        X,
        coordinates,
        k=k,
        permutations=permutations,
        channel_names=list(adata.var_names),
    )

    # Store results
    results = results.set_index('channel')
    for col in ['I', 'p_value', 'z_score', 'expected_I']:
        adata.var[f'{key_added}_{col}'] = results[col]

    adata.uns[key_added] = {
        'k': k,
        'permutations': permutations,
    }

    return adata


def compute_nearest_neighbor_distances(
    coordinates: np.ndarray,
    k: int = 1,
) -> np.ndarray:
    """Compute distances to k-nearest neighbors for each point.

    Parameters
    ----------
    coordinates : np.ndarray
        Array of shape (n, 2) with x, y coordinates.
    k : int, default=1
        Number of nearest neighbors.

    Returns
    -------
    distances : np.ndarray
        Array of shape (n, k) with distances to k nearest neighbors.
    """
    tree = spatial.cKDTree(coordinates)
    distances, _ = tree.query(coordinates, k=k+1)  # +1 because point is its own neighbor

    return distances[:, 1:]  # Exclude self-distance
