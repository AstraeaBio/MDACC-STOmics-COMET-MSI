"""Normalization utilities for MSI data.

This module provides normalization methods commonly used in Mass Spectrometry
Imaging analysis, including RobustScaler and z-score normalization.
"""

import logging
from typing import Optional, Tuple, Union, Literal

import numpy as np
import pandas as pd
from scipy import stats

try:
    import anndata as ad
except ImportError:
    ad = None

logger = logging.getLogger(__name__)


def robust_scale(
    data: np.ndarray,
    quantile_range: Tuple[float, float] = (0.05, 0.95),
    with_centering: bool = True,
    axis: int = 0,
) -> np.ndarray:
    """Apply robust scaling using quantiles.

    This method is less sensitive to outliers than standard scaling.
    It centers data using median and scales using the interquantile range.

    Parameters
    ----------
    data : np.ndarray
        Data to scale. Can be 1D or 2D.
    quantile_range : tuple, default=(0.05, 0.95)
        Quantile range used for scaling. Default uses 5th to 95th percentile.
    with_centering : bool, default=True
        Whether to center the data by subtracting the median.
    axis : int, default=0
        Axis along which to compute statistics. 0 for column-wise, 1 for row-wise.

    Returns
    -------
    scaled_data : np.ndarray
        Scaled data with same shape as input.
    """
    data = np.asarray(data, dtype=np.float64)

    # Calculate quantiles
    q_low = np.percentile(data, quantile_range[0] * 100, axis=axis, keepdims=True)
    q_high = np.percentile(data, quantile_range[1] * 100, axis=axis, keepdims=True)
    iqr = q_high - q_low

    # Avoid division by zero
    iqr = np.where(iqr == 0, 1.0, iqr)

    if with_centering:
        median = np.median(data, axis=axis, keepdims=True)
        scaled = (data - median) / iqr
    else:
        scaled = data / iqr

    return scaled


def zscore_normalize(
    data: np.ndarray,
    axis: int = 0,
    ddof: int = 0,
) -> np.ndarray:
    """Apply z-score normalization (standardization).

    Centers data to zero mean and unit variance.

    Parameters
    ----------
    data : np.ndarray
        Data to normalize.
    axis : int, default=0
        Axis along which to compute statistics.
    ddof : int, default=0
        Degrees of freedom for standard deviation calculation.

    Returns
    -------
    normalized_data : np.ndarray
        Z-score normalized data.
    """
    return stats.zscore(data, axis=axis, ddof=ddof, nan_policy='omit')


def log_transform(
    data: np.ndarray,
    pseudocount: float = 1.0,
    base: Optional[float] = None,
) -> np.ndarray:
    """Apply log transformation with optional pseudocount.

    Parameters
    ----------
    data : np.ndarray
        Data to transform.
    pseudocount : float, default=1.0
        Value added before log to handle zeros.
    base : float, optional
        Log base. If None, uses natural log. Common values: 2, 10.

    Returns
    -------
    transformed : np.ndarray
        Log-transformed data.
    """
    data = np.asarray(data, dtype=np.float64)

    if base is None:
        return np.log(data + pseudocount)
    elif base == 2:
        return np.log2(data + pseudocount)
    elif base == 10:
        return np.log10(data + pseudocount)
    else:
        return np.log(data + pseudocount) / np.log(base)


def tic_normalize(
    data: np.ndarray,
    target_sum: Optional[float] = None,
) -> np.ndarray:
    """Total Ion Current (TIC) normalization.

    Normalizes each observation (row) to sum to a constant value.

    Parameters
    ----------
    data : np.ndarray
        Data matrix with observations in rows.
    target_sum : float, optional
        Target sum for each row. If None, uses median of row sums.

    Returns
    -------
    normalized : np.ndarray
        TIC-normalized data.
    """
    data = np.asarray(data, dtype=np.float64)

    row_sums = data.sum(axis=1, keepdims=True)

    if target_sum is None:
        target_sum = np.median(row_sums[row_sums > 0])

    # Avoid division by zero
    row_sums = np.where(row_sums == 0, 1.0, row_sums)

    return data * (target_sum / row_sums)


def normalize_anndata(
    adata: 'ad.AnnData',
    method: Literal['robust', 'zscore', 'log', 'tic', 'none'] = 'robust',
    layer: Optional[str] = None,
    copy: bool = False,
    **kwargs,
) -> Optional['ad.AnnData']:
    """Normalize an AnnData object using specified method.

    Parameters
    ----------
    adata : AnnData
        Annotated data matrix.
    method : str, default='robust'
        Normalization method: 'robust', 'zscore', 'log', 'tic', or 'none'.
    layer : str, optional
        Layer to normalize. If None, uses .X
    copy : bool, default=False
        If True, return a copy instead of modifying in place.
    **kwargs
        Additional arguments passed to the normalization function.

    Returns
    -------
    adata : AnnData or None
        Returns AnnData if copy=True, otherwise modifies in place.
    """
    if ad is None:
        raise ImportError("anndata is required for normalize_anndata")

    if copy:
        adata = adata.copy()

    # Get data matrix
    if layer is not None:
        X = adata.layers[layer]
    else:
        X = adata.X

    # Convert to dense if sparse
    if hasattr(X, 'toarray'):
        X = X.toarray()
    else:
        X = np.asarray(X)

    # Store raw data before normalization
    if adata.raw is None and method != 'none':
        adata.raw = adata.copy()

    # Apply normalization
    if method == 'robust':
        quantile_range = kwargs.get('quantile_range', (0.05, 0.95))
        X_norm = robust_scale(X, quantile_range=quantile_range, axis=0)
    elif method == 'zscore':
        X_norm = zscore_normalize(X, axis=0)
    elif method == 'log':
        pseudocount = kwargs.get('pseudocount', 1.0)
        base = kwargs.get('base', None)
        X_norm = log_transform(X, pseudocount=pseudocount, base=base)
    elif method == 'tic':
        target_sum = kwargs.get('target_sum', None)
        X_norm = tic_normalize(X, target_sum=target_sum)
    elif method == 'none':
        X_norm = X
    else:
        raise ValueError(f"Unknown normalization method: {method}")

    # Store normalized data
    if layer is not None:
        adata.layers[layer] = X_norm
    else:
        adata.X = X_norm

    # Log the normalization
    adata.uns['normalization'] = {
        'method': method,
        'parameters': kwargs,
    }

    logger.info(f"Applied {method} normalization")

    if copy:
        return adata


def compute_qc_metrics(
    data: np.ndarray,
    feature_names: Optional[list] = None,
) -> pd.DataFrame:
    """Compute quality control metrics for each feature (channel).

    Parameters
    ----------
    data : np.ndarray
        Data matrix with observations in rows and features in columns.
    feature_names : list, optional
        Names for each feature.

    Returns
    -------
    metrics : pd.DataFrame
        DataFrame with QC metrics for each feature.
    """
    data = np.asarray(data, dtype=np.float64)
    n_features = data.shape[1]

    if feature_names is None:
        feature_names = [f'feature_{i}' for i in range(n_features)]

    metrics = pd.DataFrame({
        'feature': feature_names,
        'mean': np.nanmean(data, axis=0),
        'std': np.nanstd(data, axis=0),
        'median': np.nanmedian(data, axis=0),
        'min': np.nanmin(data, axis=0),
        'max': np.nanmax(data, axis=0),
        'pct_nonzero': (data > 0).sum(axis=0) / data.shape[0] * 100,
        'pct_5': np.nanpercentile(data, 5, axis=0),
        'pct_95': np.nanpercentile(data, 95, axis=0),
    })

    return metrics


def filter_by_intensity(
    data: np.ndarray,
    min_total: Optional[float] = None,
    max_total: Optional[float] = None,
    min_percentile: Optional[float] = None,
    max_percentile: Optional[float] = None,
) -> np.ndarray:
    """Create a boolean mask for filtering observations by total intensity.

    Parameters
    ----------
    data : np.ndarray
        Data matrix with observations in rows.
    min_total : float, optional
        Minimum total intensity threshold.
    max_total : float, optional
        Maximum total intensity threshold.
    min_percentile : float, optional
        Minimum percentile threshold (e.g., 1 for 1st percentile).
    max_percentile : float, optional
        Maximum percentile threshold (e.g., 99 for 99th percentile).

    Returns
    -------
    mask : np.ndarray
        Boolean array indicating which observations pass the filter.
    """
    total_intensity = data.sum(axis=1)
    mask = np.ones(len(total_intensity), dtype=bool)

    if min_percentile is not None:
        min_val = np.percentile(total_intensity, min_percentile)
        mask &= total_intensity >= min_val

    if max_percentile is not None:
        max_val = np.percentile(total_intensity, max_percentile)
        mask &= total_intensity <= max_val

    if min_total is not None:
        mask &= total_intensity >= min_total

    if max_total is not None:
        mask &= total_intensity <= max_total

    return mask
