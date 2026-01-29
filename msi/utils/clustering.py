"""Clustering utilities for MSI data.

This module provides functions for dimensionality reduction, clustering,
and visualization of MSI data using PCA, UMAP, and Leiden algorithms.
"""

import logging
from typing import Optional, Tuple, Union, Literal, List

import numpy as np
import pandas as pd
from scipy import sparse

try:
    import scanpy as sc
    import anndata as ad
    HAS_SCANPY = True
except ImportError:
    HAS_SCANPY = False

try:
    import umap
    HAS_UMAP = True
except ImportError:
    HAS_UMAP = False

logger = logging.getLogger(__name__)


def run_pca(
    data: np.ndarray,
    n_components: int = 21,
    random_state: int = 42,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Run PCA dimensionality reduction.

    Parameters
    ----------
    data : np.ndarray
        Data matrix with observations in rows.
    n_components : int, default=21
        Number of principal components.
    random_state : int, default=42
        Random seed for reproducibility.

    Returns
    -------
    pca_result : np.ndarray
        PCA-transformed data of shape (n_obs, n_components).
    explained_variance : np.ndarray
        Explained variance ratio for each component.
    components : np.ndarray
        Principal component loadings.
    """
    from sklearn.decomposition import PCA

    # Handle NaN values
    data = np.nan_to_num(data, nan=0.0)

    pca = PCA(n_components=n_components, random_state=random_state)
    pca_result = pca.fit_transform(data)

    return pca_result, pca.explained_variance_ratio_, pca.components_


def run_umap(
    data: np.ndarray,
    n_neighbors: int = 15,
    n_components: int = 2,
    min_dist: float = 0.1,
    metric: str = 'euclidean',
    random_state: int = 42,
) -> np.ndarray:
    """Run UMAP dimensionality reduction.

    Parameters
    ----------
    data : np.ndarray
        Data matrix (often PCA-reduced) with observations in rows.
    n_neighbors : int, default=15
        Number of neighbors for UMAP.
    n_components : int, default=2
        Number of UMAP dimensions.
    min_dist : float, default=0.1
        Minimum distance between points in embedding.
    metric : str, default='euclidean'
        Distance metric.
    random_state : int, default=42
        Random seed for reproducibility.

    Returns
    -------
    umap_result : np.ndarray
        UMAP embedding of shape (n_obs, n_components).
    """
    if not HAS_UMAP:
        raise ImportError("umap-learn is required. Install with: pip install umap-learn")

    reducer = umap.UMAP(
        n_neighbors=n_neighbors,
        n_components=n_components,
        min_dist=min_dist,
        metric=metric,
        random_state=random_state,
    )

    # Handle NaN values
    data = np.nan_to_num(data, nan=0.0)

    return reducer.fit_transform(data)


def build_knn_graph(
    data: np.ndarray,
    n_neighbors: int = 15,
    metric: str = 'euclidean',
) -> Tuple[np.ndarray, np.ndarray]:
    """Build a k-nearest neighbors graph.

    Parameters
    ----------
    data : np.ndarray
        Data matrix with observations in rows.
    n_neighbors : int, default=15
        Number of neighbors.
    metric : str, default='euclidean'
        Distance metric.

    Returns
    -------
    indices : np.ndarray
        Neighbor indices of shape (n_obs, n_neighbors).
    distances : np.ndarray
        Neighbor distances of shape (n_obs, n_neighbors).
    """
    from sklearn.neighbors import NearestNeighbors

    nn = NearestNeighbors(n_neighbors=n_neighbors, metric=metric)
    nn.fit(data)
    distances, indices = nn.kneighbors(data)

    return indices, distances


def run_leiden(
    adjacency: np.ndarray,
    resolution: float = 0.6,
    random_state: int = 42,
) -> np.ndarray:
    """Run Leiden community detection algorithm.

    Parameters
    ----------
    adjacency : np.ndarray
        Adjacency matrix (can be sparse).
    resolution : float, default=0.6
        Resolution parameter for clustering.
    random_state : int, default=42
        Random seed for reproducibility.

    Returns
    -------
    labels : np.ndarray
        Cluster labels for each observation.
    """
    try:
        import leidenalg
        import igraph as ig
    except ImportError:
        raise ImportError("leidenalg and igraph are required. "
                         "Install with: pip install leidenalg igraph")

    # Convert to sparse if not already
    if not sparse.issparse(adjacency):
        adjacency = sparse.csr_matrix(adjacency)

    # Create igraph from adjacency
    sources, targets = adjacency.nonzero()
    weights = adjacency[sources, targets].A1

    g = ig.Graph(directed=False)
    g.add_vertices(adjacency.shape[0])
    g.add_edges(list(zip(sources, targets)))
    g.es['weight'] = weights

    # Run Leiden
    partition = leidenalg.find_partition(
        g,
        leidenalg.RBConfigurationVertexPartition,
        resolution_parameter=resolution,
        seed=random_state,
        weights='weight',
    )

    return np.array(partition.membership)


def cluster_anndata(
    adata: 'ad.AnnData',
    n_pcs: int = 21,
    n_neighbors: int = 15,
    resolution: float = 0.6,
    layer: Optional[str] = None,
    use_highly_variable: bool = True,
    random_state: int = 42,
) -> 'ad.AnnData':
    """Run full clustering pipeline on AnnData object.

    Performs PCA, neighbor graph construction, Leiden clustering, and UMAP.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    n_pcs : int, default=21
        Number of principal components.
    n_neighbors : int, default=15
        Number of neighbors for graph construction.
    resolution : float, default=0.6
        Resolution for Leiden clustering.
    layer : str, optional
        Layer to use. If None, uses .X
    use_highly_variable : bool, default=True
        Whether to subset to highly variable features.
    random_state : int, default=42
        Random seed for reproducibility.

    Returns
    -------
    adata : AnnData
        Modified AnnData with clustering results.
    """
    if not HAS_SCANPY:
        raise ImportError("scanpy is required. Install with: pip install scanpy")

    # Use scanpy's pipeline
    logger.info("Running PCA...")
    sc.tl.pca(adata, n_comps=n_pcs, random_state=random_state)

    logger.info(f"Computing neighbors (n={n_neighbors})...")
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=random_state)

    logger.info(f"Running Leiden clustering (resolution={resolution})...")
    sc.tl.leiden(adata, resolution=resolution, random_state=random_state)

    logger.info("Computing UMAP...")
    sc.tl.umap(adata, random_state=random_state)

    # Store parameters
    adata.uns['clustering_params'] = {
        'n_pcs': n_pcs,
        'n_neighbors': n_neighbors,
        'resolution': resolution,
        'random_state': random_state,
    }

    return adata


def select_highly_variable(
    data: np.ndarray,
    n_top: int = 50,
    method: Literal['variance', 'dispersion', 'cv'] = 'variance',
    feature_names: Optional[list] = None,
) -> Tuple[np.ndarray, list]:
    """Select highly variable features.

    Parameters
    ----------
    data : np.ndarray
        Data matrix with observations in rows and features in columns.
    n_top : int, default=50
        Number of top features to select.
    method : str, default='variance'
        Selection method: 'variance', 'dispersion', or 'cv' (coefficient of variation).
    feature_names : list, optional
        Names of features.

    Returns
    -------
    selected_data : np.ndarray
        Data subset to highly variable features.
    selected_names : list
        Names of selected features.
    """
    if feature_names is None:
        feature_names = [f'feature_{i}' for i in range(data.shape[1])]

    if method == 'variance':
        scores = np.nanvar(data, axis=0)
    elif method == 'dispersion':
        means = np.nanmean(data, axis=0)
        variances = np.nanvar(data, axis=0)
        scores = variances / (means + 1e-10)
    elif method == 'cv':
        means = np.nanmean(data, axis=0)
        stds = np.nanstd(data, axis=0)
        scores = stds / (means + 1e-10)
    else:
        raise ValueError(f"Unknown method: {method}")

    # Select top features
    top_indices = np.argsort(scores)[-n_top:][::-1]
    selected_data = data[:, top_indices]
    selected_names = [feature_names[i] for i in top_indices]

    logger.info(f"Selected {n_top} highly variable features using {method}")

    return selected_data, selected_names


def cluster_metrics(
    labels: np.ndarray,
    data: Optional[np.ndarray] = None,
) -> dict:
    """Compute clustering quality metrics.

    Parameters
    ----------
    labels : np.ndarray
        Cluster labels.
    data : np.ndarray, optional
        Data used for clustering (for silhouette score).

    Returns
    -------
    metrics : dict
        Dictionary of clustering metrics.
    """
    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)

    metrics = {
        'n_clusters': n_clusters,
        'cluster_sizes': {int(l): int(np.sum(labels == l)) for l in unique_labels},
        'min_cluster_size': int(min(np.sum(labels == l) for l in unique_labels)),
        'max_cluster_size': int(max(np.sum(labels == l) for l in unique_labels)),
    }

    # Compute silhouette score if data provided
    if data is not None and n_clusters > 1:
        from sklearn.metrics import silhouette_score
        try:
            metrics['silhouette_score'] = float(silhouette_score(data, labels))
        except Exception as e:
            logger.warning(f"Could not compute silhouette score: {e}")

    return metrics


def save_clustering_results(
    adata: 'ad.AnnData',
    output_dir: str,
    sample_id: str,
    save_umap: bool = True,
    save_clusters: bool = True,
) -> List[str]:
    """Save clustering results to CSV files.

    Parameters
    ----------
    adata : AnnData
        AnnData with clustering results.
    output_dir : str
        Output directory.
    sample_id : str
        Sample identifier for filenames.
    save_umap : bool, default=True
        Whether to save UMAP coordinates.
    save_clusters : bool, default=True
        Whether to save cluster assignments.

    Returns
    -------
    saved_files : list
        List of saved file paths.
    """
    from pathlib import Path

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    saved_files = []

    if save_clusters and 'leiden' in adata.obs.columns:
        clusters_df = pd.DataFrame({
            'cluster': adata.obs['leiden'].values,
        })

        # Add spatial coordinates if available
        if 'spatial' in adata.obsm:
            clusters_df['x'] = adata.obsm['spatial'][:, 0]
            clusters_df['y'] = adata.obsm['spatial'][:, 1]

        clusters_path = output_dir / f'{sample_id}_clusters.csv'
        clusters_df.to_csv(clusters_path, index=False)
        saved_files.append(str(clusters_path))
        logger.info(f"Saved clusters to {clusters_path}")

    if save_umap and 'X_umap' in adata.obsm:
        umap_df = pd.DataFrame({
            'umap_1': adata.obsm['X_umap'][:, 0],
            'umap_2': adata.obsm['X_umap'][:, 1],
        })

        if 'leiden' in adata.obs.columns:
            umap_df['cluster'] = adata.obs['leiden'].values

        umap_path = output_dir / f'{sample_id}_umap.csv'
        umap_df.to_csv(umap_path, index=False)
        saved_files.append(str(umap_path))
        logger.info(f"Saved UMAP to {umap_path}")

    return saved_files
