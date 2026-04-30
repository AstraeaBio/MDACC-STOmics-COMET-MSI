"""MSI-to-H&E Registration Utilities
====================================

Registers Mass Spectrometry Imaging data to H&E histology images, which are
morphologically the closest match to MSI tissue images. From H&E space, the
transformation chain to COMET/STOmics is well-established via VALIS.

Registration chain:
    MSI → H&E → COMET DAPI → STOmics DAPI

Design choices:
    - MSI→H&E is the primary alignment step (morphological similarity)
    - Metabolite subset selection improves composite image quality for registration
    - Multiple composite strategies with automatic QC-based selection
    - Fallback cascade: PCA → SNR → MIP if registration fails
    - Affine seed support for difficult cases (manual or QuPath-derived)

Usage:
    from msi.utils.registration import MSItoHERegistrar

    reg = MSItoHERegistrar(
        msi_adata=adata_glycans,
        he_image_path='path/to/sample_HE.ndpi',
        output_dir='path/to/output',
    )
    reg.register()
    warped_coords = reg.transform_coordinates(adata_glycans.obsm['spatial'])
"""

import json
import logging
import shutil
import tempfile
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd
import tifffile
from scipy import ndimage
from skimage import exposure, transform

logger = logging.getLogger(__name__)

try:
    from valis import registration as valis_reg
    HAS_VALIS = True
except ImportError:
    HAS_VALIS = False
    logger.warning("VALIS not installed. Install with: pip install valis-wsi")

try:
    import cv2
    HAS_CV2 = True
except ImportError:
    HAS_CV2 = False

try:
    import anndata as ad
except ImportError:
    ad = None


# ============================================================
# Composite Image Creation
# ============================================================

def select_metabolite_subset(
    adata,
    method: Literal['tissue_contrast', 'variance', 'spatial', 'manual'] = 'tissue_contrast',
    n_channels: int = 15,
    manual_channels: Optional[List[str]] = None,
    spatial_threshold: float = 0.1,
) -> List[str]:
    """Select a subset of MSI channels best suited for registration.

    The key insight is that not all MSI channels are useful for alignment —
    some are noisy, some are uniform, and some highlight tissue boundaries
    (which is what we need for morphological registration to H&E).

    Parameters
    ----------
    adata : AnnData
        MSI AnnData with expression matrix and spatial coordinates.
    method : str
        Selection strategy:
        - 'tissue_contrast': Channels with highest tissue-vs-background contrast
        - 'variance': Highest spatial variance channels
        - 'spatial': Channels with highest Moran's I (spatial autocorrelation)
        - 'manual': Use provided channel list
    n_channels : int
        Number of channels to select.
    manual_channels : list, optional
        Channel names for manual selection.
    spatial_threshold : float
        Minimum Moran's I for spatial method.

    Returns
    -------
    list
        Selected channel names.
    """
    if method == 'manual' and manual_channels:
        valid = [c for c in manual_channels if c in adata.var_names]
        if len(valid) < len(manual_channels):
            missing = set(manual_channels) - set(valid)
            logger.warning(f"Missing channels: {missing}")
        return valid

    X = adata.X.toarray() if hasattr(adata.X, 'toarray') else adata.X

    if method == 'tissue_contrast':
        # Channels where tissue regions differ most from background
        # Approximate tissue mask via total intensity threshold
        total_per_pixel = X.sum(axis=1)
        threshold = np.percentile(total_per_pixel[total_per_pixel > 0], 25)
        tissue_mask = total_per_pixel > threshold

        if tissue_mask.sum() < 10 or (~tissue_mask).sum() < 10:
            logger.warning("Cannot distinguish tissue/background, falling back to variance")
            method = 'variance'
        else:
            tissue_mean = X[tissue_mask].mean(axis=0)
            bg_mean = X[~tissue_mask].mean(axis=0)
            tissue_std = X[tissue_mask].std(axis=0) + 1e-10
            contrast = np.abs(tissue_mean - bg_mean) / tissue_std
            top_idx = np.argsort(contrast)[-n_channels:]
            return [adata.var_names[i] for i in top_idx]

    if method == 'variance':
        variance = X.var(axis=0)
        top_idx = np.argsort(variance)[-n_channels:]
        return [adata.var_names[i] for i in top_idx]

    if method == 'spatial':
        # Use pre-computed Moran's I if available
        if 'morans_i' in adata.var.columns:
            mi = adata.var['morans_i'].values
            above_threshold = mi > spatial_threshold
            if above_threshold.sum() >= n_channels:
                top_idx = np.argsort(mi)[-n_channels:]
            else:
                top_idx = np.argsort(mi)[-n_channels:]
            return [adata.var_names[i] for i in top_idx]
        else:
            logger.warning("No morans_i in var; falling back to variance")
            return select_metabolite_subset(adata, method='variance', n_channels=n_channels)

    raise ValueError(f"Unknown method: {method}")


def create_composite_image(
    adata,
    channels: Optional[List[str]] = None,
    method: Literal['pca', 'snr', 'mip', 'weighted_sum'] = 'pca',
    pca_n_components: int = 3,
    snr_top_n: int = 10,
    output_shape: Optional[Tuple[int, int]] = None,
) -> np.ndarray:
    """Create a composite image from MSI data for registration.

    Parameters
    ----------
    adata : AnnData
        MSI data with obsm['spatial'] coordinates.
    channels : list, optional
        Subset of channels to use. If None, uses all.
    method : str
        Composite creation method.
    pca_n_components : int
        Number of PCA components (for 'pca' method).
    output_shape : tuple, optional
        (height, width) of output. If None, inferred from coordinates.

    Returns
    -------
    np.ndarray
        Composite image (uint8), either (H, W) or (H, W, C).
    """
    coords = adata.obsm['spatial']

    if channels:
        valid_channels = [c for c in channels if c in adata.var_names]
        X = adata[:, valid_channels].X
    else:
        X = adata.X

    X = X.toarray() if hasattr(X, 'toarray') else X.copy()

    if output_shape is None:
        x_max = int(coords[:, 0].max()) + 1
        y_max = int(coords[:, 1].max()) + 1
    else:
        y_max, x_max = output_shape

    ix = np.clip(coords[:, 0].astype(int), 0, x_max - 1)
    iy = np.clip(coords[:, 1].astype(int), 0, y_max - 1)

    if method == 'pca':
        from sklearn.decomposition import PCA
        n_comp = min(pca_n_components, X.shape[1])
        pca = PCA(n_components=n_comp)
        pca_result = pca.fit_transform(X)

        composite = np.zeros((y_max, x_max, n_comp), dtype=np.float32)
        composite[iy, ix, :] = pca_result

        for c in range(n_comp):
            ch = composite[:, :, c]
            composite[:, :, c] = exposure.rescale_intensity(ch, out_range=(0, 255))

        return composite.astype(np.uint8)

    elif method == 'snr':
        means = X.mean(axis=0)
        stds = X.std(axis=0) + 1e-10
        snr = means / stds
        top_idx = np.argsort(snr)[-snr_top_n:]
        X_top = X[:, top_idx]
        weights = snr[top_idx]
        weights = weights / weights.sum()
        weighted = X_top @ weights

        composite = np.zeros((y_max, x_max), dtype=np.float32)
        composite[iy, ix] = weighted
        composite = exposure.rescale_intensity(composite, out_range=(0, 255))
        return composite.astype(np.uint8)

    elif method == 'mip':
        max_intensity = X.max(axis=1)
        composite = np.zeros((y_max, x_max), dtype=np.float32)
        composite[iy, ix] = max_intensity
        composite = exposure.rescale_intensity(composite, out_range=(0, 255))
        return composite.astype(np.uint8)

    elif method == 'weighted_sum':
        # Weight by coefficient of variation — channels with clear
        # tissue structure contribute more
        cv = X.std(axis=0) / (X.mean(axis=0) + 1e-10)
        weights = cv / cv.sum()
        weighted = X @ weights

        composite = np.zeros((y_max, x_max), dtype=np.float32)
        composite[iy, ix] = weighted
        composite = exposure.rescale_intensity(composite, out_range=(0, 255))
        return composite.astype(np.uint8)

    else:
        raise ValueError(f"Unknown method: {method}")


# ============================================================
# H&E Image Loading
# ============================================================

def load_he_image(
    he_path: Union[str, Path],
    max_dim: int = 4096,
    grayscale: bool = True,
) -> np.ndarray:
    """Load H&E image, optionally downsample for registration.

    Supports NDPI (via tifffile/openslide), standard TIFF, and PNG/JPG.

    Parameters
    ----------
    he_path : path
        Path to H&E image file.
    max_dim : int
        Maximum dimension (for memory-safe loading of large images).
    grayscale : bool
        Convert to grayscale for registration (recommended).

    Returns
    -------
    np.ndarray
        H&E image (uint8).
    """
    he_path = Path(he_path)

    if he_path.suffix.lower() in ('.ndpi', '.svs', '.mrxs'):
        # Whole-slide image — use openslide if available, else tifffile thumbnail
        try:
            import openslide
            slide = openslide.OpenSlide(str(he_path))
            # Get a reasonable resolution level
            dims = slide.dimensions  # (width, height)
            downsample = max(dims[0], dims[1]) / max_dim
            level = slide.get_best_level_for_downsample(downsample)
            img = np.array(slide.read_region((0, 0), level, slide.level_dimensions[level]))
            img = img[:, :, :3]  # drop alpha
            slide.close()
        except ImportError:
            logger.info("openslide not available, using tifffile for WSI thumbnail")
            with tifffile.TiffFile(str(he_path)) as tif:
                # Read the smallest pyramid level that's ≥ max_dim
                for page in reversed(tif.pages):
                    if max(page.shape[:2]) >= max_dim or page == tif.pages[0]:
                        img = page.asarray()
                        break
                    img = page.asarray()
    else:
        img = tifffile.imread(str(he_path))

    # Handle various channel arrangements
    if img.ndim == 2:
        return img.astype(np.uint8) if img.dtype != np.uint8 else img

    if img.ndim == 3:
        if img.shape[0] in (3, 4) and img.shape[0] < img.shape[2]:
            img = np.moveaxis(img, 0, -1)
        if img.shape[2] == 4:
            img = img[:, :, :3]

    # Downsample if needed
    if max(img.shape[:2]) > max_dim:
        scale = max_dim / max(img.shape[:2])
        img = transform.rescale(img, scale, channel_axis=2 if img.ndim == 3 else None,
                                preserve_range=True, anti_aliasing=True).astype(np.uint8)

    if grayscale and img.ndim == 3:
        # Convert H&E to grayscale — invert so tissue is bright (like DAPI)
        gray = np.mean(img.astype(np.float32), axis=2)
        gray = 255 - gray  # invert: tissue bright, background dark
        return exposure.rescale_intensity(gray, out_range=(0, 255)).astype(np.uint8)

    return img.astype(np.uint8) if img.dtype != np.uint8 else img


# ============================================================
# Registration Engine
# ============================================================

class MSItoHERegistrar:
    """Register MSI data to H&E histology image using VALIS.

    The registration chain is MSI → H&E. From H&E space, you can chain
    to COMET/STOmics via the existing QuPath+VALIS pipeline.

    Parameters
    ----------
    msi_adata : AnnData
        MSI data with obsm['spatial'] and expression matrix.
    he_image_path : path
        Path to H&E image (NDPI, TIFF, PNG, etc.).
    output_dir : path
        Output directory for registration artifacts.
    channel_subset : list, optional
        Specific MSI channels to use for composite. If None, auto-selects.
    channel_selection_method : str
        How to auto-select channels: 'tissue_contrast', 'variance', 'spatial', 'manual'
    n_channels : int
        Number of channels for auto-selection.
    composite_method : str
        Composite creation: 'pca', 'snr', 'mip', 'weighted_sum'
    max_image_dim : int
        Maximum image dimension for registration (memory control).
    affine_seed : np.ndarray, optional
        Initial 3×3 affine matrix (from QuPath or manual alignment).
    """

    def __init__(
        self,
        msi_adata,
        he_image_path: Union[str, Path],
        output_dir: Union[str, Path],
        channel_subset: Optional[List[str]] = None,
        channel_selection_method: str = 'tissue_contrast',
        n_channels: int = 15,
        composite_method: str = 'pca',
        max_image_dim: int = 4096,
        affine_seed: Optional[np.ndarray] = None,
    ):
        self.msi_adata = msi_adata
        self.he_path = Path(he_image_path)
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        self.channel_subset = channel_subset
        self.channel_selection_method = channel_selection_method
        self.n_channels = n_channels
        self.composite_method = composite_method
        self.max_image_dim = max_image_dim
        self.affine_seed = affine_seed

        # State
        self.selected_channels = None
        self.composite_image = None
        self.he_image = None
        self.registrar = None
        self.transform_matrix = None
        self.quality_metrics = None

    def select_channels(self) -> List[str]:
        """Auto-select or validate MSI channels for registration."""
        if self.channel_subset:
            self.selected_channels = [c for c in self.channel_subset
                                      if c in self.msi_adata.var_names]
        else:
            self.selected_channels = select_metabolite_subset(
                self.msi_adata,
                method=self.channel_selection_method,
                n_channels=self.n_channels,
            )

        logger.info(f"Selected {len(self.selected_channels)} channels for composite")
        return self.selected_channels

    def create_composite(self) -> np.ndarray:
        """Create composite MSI image for registration."""
        if self.selected_channels is None:
            self.select_channels()

        self.composite_image = create_composite_image(
            self.msi_adata,
            channels=self.selected_channels,
            method=self.composite_method,
        )

        # Save composite
        composite_path = self.output_dir / 'msi_composite.tif'
        tifffile.imwrite(str(composite_path), self.composite_image)
        logger.info(f"Saved composite ({self.composite_method}): {composite_path}")

        return self.composite_image

    def load_he(self) -> np.ndarray:
        """Load and preprocess H&E image."""
        self.he_image = load_he_image(
            self.he_path,
            max_dim=self.max_image_dim,
            grayscale=True,
        )

        he_path = self.output_dir / 'he_grayscale.tif'
        tifffile.imwrite(str(he_path), self.he_image)
        logger.info(f"Loaded H&E: {self.he_image.shape}")

        return self.he_image

    def register(
        self,
        fallback_methods: Optional[List[str]] = None,
    ) -> dict:
        """Run full registration pipeline with fallback cascade.

        Parameters
        ----------
        fallback_methods : list, optional
            Composite methods to try if primary fails.
            Default: ['pca', 'snr', 'mip']

        Returns
        -------
        dict
            Registration results with quality metrics.
        """
        if not HAS_VALIS:
            raise ImportError("VALIS required. Install with: pip install valis-wsi")

        if fallback_methods is None:
            fallback_methods = ['pca', 'snr', 'mip']

        # Ensure primary method is first
        if self.composite_method in fallback_methods:
            fallback_methods.remove(self.composite_method)
        methods_to_try = [self.composite_method] + fallback_methods

        if self.selected_channels is None:
            self.select_channels()

        if self.he_image is None:
            self.load_he()

        best_result = None
        best_error = float('inf')

        for method in methods_to_try:
            logger.info(f"\nAttempting registration with {method} composite...")

            try:
                self.composite_method = method
                composite = create_composite_image(
                    self.msi_adata,
                    channels=self.selected_channels,
                    method=method,
                )

                result = self._run_valis(composite)

                # Evaluate quality
                error = result.get('mean_error', float('inf'))
                logger.info(f"  {method}: mean_error={error:.2f}")

                if error < best_error:
                    best_error = error
                    best_result = result
                    best_result['composite_method'] = method
                    self.composite_image = composite

                # If error is acceptable, stop trying
                if error < 50:  # pixels
                    logger.info(f"  Acceptable registration with {method}")
                    break

            except Exception as e:
                logger.warning(f"  {method} failed: {e}")
                continue

        if best_result is None:
            raise RuntimeError("All registration methods failed")

        self.registrar = best_result.get('registrar')
        self.quality_metrics = best_result
        self._save_results(best_result)

        return best_result

    def _run_valis(self, composite: np.ndarray) -> dict:
        """Execute VALIS registration between composite and H&E."""
        # VALIS needs files on disk in the same directory
        tmp_dir = self.output_dir / '_valis_input'
        tmp_dir.mkdir(exist_ok=True)

        # Ensure both images are comparable
        # Resize composite to match H&E aspect ratio approximately
        he_h, he_w = self.he_image.shape[:2]

        if composite.ndim == 3:
            # For PCA composites, take first component as grayscale for VALIS
            comp_gray = composite[:, :, 0] if composite.shape[2] >= 1 else composite
        else:
            comp_gray = composite

        # Save images
        msi_path = tmp_dir / 'msi_composite.tif'
        he_path = tmp_dir / 'he_reference.tif'
        tifffile.imwrite(str(msi_path), comp_gray)
        tifffile.imwrite(str(he_path), self.he_image)

        # Run VALIS
        valis_out = self.output_dir / '_valis_output'

        registrar = valis_reg.Valis(
            str(tmp_dir),
            str(valis_out),
            reference_img_f='he_reference.tif',
            align_to_reference=True,
        )

        rigid_registrar, non_rigid_registrar, error_df = registrar.register()

        # Extract error metrics
        result = {
            'registrar': registrar,
            'error_df': error_df,
            'mean_error': error_df['mean'].mean() if len(error_df) > 0 else float('inf'),
        }

        # Extract affine matrix from the MSI slide
        msi_slide = registrar.get_slide('msi_composite.tif')
        if msi_slide is not None:
            try:
                M = msi_slide.M
                result['affine_matrix'] = M.tolist() if hasattr(M, 'tolist') else M
            except Exception:
                pass

        return result

    def transform_coordinates(
        self,
        coords: np.ndarray,
    ) -> np.ndarray:
        """Transform MSI pixel coordinates to H&E space.

        Parameters
        ----------
        coords : np.ndarray
            (N, 2) array of (x, y) coordinates in MSI pixel space.

        Returns
        -------
        np.ndarray
            (N, 2) array of (x, y) coordinates in H&E pixel space.
        """
        if self.registrar is None:
            raise RuntimeError("Must call register() first")

        msi_slide = self.registrar.get_slide('msi_composite.tif')
        if msi_slide is None:
            raise RuntimeError("MSI slide not found in registrar")

        warped = msi_slide.warp_xy(coords)
        return warped

    def transform_adata(
        self,
        adata,
        coord_key: str = 'spatial',
        warped_key: str = 'spatial_he',
    ):
        """Transform AnnData spatial coordinates to H&E space.

        Parameters
        ----------
        adata : AnnData
            MSI data (can be different modality than what was registered,
            as long as coordinates are in the same MSI pixel space).
        coord_key : str
            Key in obsm for input coordinates.
        warped_key : str
            Key in obsm for warped coordinates.

        Returns
        -------
        AnnData
            Updated AnnData with warped coordinates.
        """
        coords = adata.obsm[coord_key]
        warped = self.transform_coordinates(coords)
        adata.obsm[warped_key] = warped
        adata.obs['x_he'] = warped[:, 0]
        adata.obs['y_he'] = warped[:, 1]

        # Store registration metadata
        adata.uns['registration_msi_to_he'] = {
            'method': self.quality_metrics.get('composite_method', 'unknown'),
            'mean_error': float(self.quality_metrics.get('mean_error', -1)),
            'n_channels_used': len(self.selected_channels) if self.selected_channels else 0,
            'he_image': str(self.he_path),
        }

        return adata

    def _save_results(self, result: dict):
        """Save registration results and QC artifacts."""
        # Save quality metrics
        metrics = {k: v for k, v in result.items()
                   if k not in ('registrar', 'error_df')}
        metrics_path = self.output_dir / 'registration_metrics.json'
        with open(metrics_path, 'w') as f:
            json.dump(metrics, f, indent=2, default=str)

        # Save selected channels
        if self.selected_channels:
            channels_path = self.output_dir / 'selected_channels.json'
            with open(channels_path, 'w') as f:
                json.dump(self.selected_channels, f, indent=2)

        logger.info(f"Saved registration results to {self.output_dir}")

    def plot_qc(self, save: bool = True):
        """Generate QC visualization of the registration."""
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(14, 14))

        # Panel 1: H&E reference
        ax = axes[0, 0]
        ax.imshow(self.he_image, cmap='gray')
        ax.set_title('H&E Reference (inverted grayscale)')
        ax.axis('off')

        # Panel 2: MSI composite
        ax = axes[0, 1]
        if self.composite_image.ndim == 3:
            ax.imshow(self.composite_image)
        else:
            ax.imshow(self.composite_image, cmap='gray')
        ax.set_title(f'MSI Composite ({self.composite_method})')
        ax.axis('off')

        # Panel 3: Overlay (if registered)
        ax = axes[1, 0]
        if self.registrar is not None:
            try:
                msi_slide = self.registrar.get_slide('msi_composite.tif')
                if msi_slide is not None and hasattr(msi_slide, 'warped_img'):
                    warped = msi_slide.warped_img
                    if warped is not None:
                        # Blend
                        he_norm = self.he_image.astype(np.float32) / 255
                        w_norm = warped.astype(np.float32) / 255 if warped.ndim == 2 else warped[:, :, 0].astype(np.float32) / 255
                        # Resize to match if needed
                        if he_norm.shape != w_norm.shape:
                            w_norm = transform.resize(w_norm, he_norm.shape, preserve_range=True)
                        blend = np.stack([he_norm, w_norm, np.zeros_like(he_norm)], axis=2)
                        ax.imshow(np.clip(blend, 0, 1))
                        ax.set_title('Overlay (magenta=H&E, green=MSI)')
                    else:
                        ax.text(0.5, 0.5, 'Warped image not available', ha='center', va='center', transform=ax.transAxes)
                else:
                    ax.text(0.5, 0.5, 'MSI slide not found', ha='center', va='center', transform=ax.transAxes)
            except Exception as e:
                ax.text(0.5, 0.5, f'Error: {e}', ha='center', va='center', transform=ax.transAxes)
        else:
            ax.text(0.5, 0.5, 'Not registered yet', ha='center', va='center', transform=ax.transAxes)
        ax.axis('off')

        # Panel 4: Coordinate scatter (original vs warped)
        ax = axes[1, 1]
        coords = self.msi_adata.obsm['spatial']
        ax.scatter(coords[:, 0], coords[:, 1], s=0.5, alpha=0.3, c='blue', label='MSI pixels')
        if self.registrar is not None:
            try:
                warped_coords = self.transform_coordinates(coords)
                ax.scatter(warped_coords[:, 0], warped_coords[:, 1], s=0.5, alpha=0.3, c='red', label='Warped to H&E')
            except Exception:
                pass
        ax.legend()
        ax.set_aspect('equal')
        ax.set_title('Coordinate Transform')

        plt.suptitle(f'MSI→H&E Registration QC', fontsize=14, fontweight='bold')
        plt.tight_layout()

        if save:
            fig.savefig(self.output_dir / 'registration_qc.png', dpi=150, bbox_inches='tight')

        return fig


# ============================================================
# Batch Registration
# ============================================================

def batch_register_msi_to_he(
    sample_configs: List[Dict],
    output_base: Union[str, Path],
    channel_selection_method: str = 'tissue_contrast',
    composite_method: str = 'pca',
    n_channels: int = 15,
    max_image_dim: int = 4096,
) -> pd.DataFrame:
    """Register multiple MSI samples to their corresponding H&E images.

    Parameters
    ----------
    sample_configs : list of dict
        Each dict must have:
        - 'sample_id': str
        - 'msi_adata_path': path to MSI h5ad
        - 'he_image_path': path to H&E image
        Optional:
        - 'channel_subset': list of channel names
        - 'affine_seed': 3x3 matrix
    output_base : path
        Base output directory (per-sample subdirectories created).
    channel_selection_method, composite_method, n_channels, max_image_dim:
        Passed to MSItoHERegistrar.

    Returns
    -------
    pd.DataFrame
        Summary of registration results per sample.
    """
    output_base = Path(output_base)
    results = []

    for i, config in enumerate(sample_configs):
        sid = config['sample_id']
        logger.info(f"\n{'='*60}")
        logger.info(f"Registering {sid} ({i+1}/{len(sample_configs)})")
        logger.info(f"{'='*60}")

        try:
            adata = ad.read_h5ad(config['msi_adata_path'])
            sample_out = output_base / sid
            sample_out.mkdir(parents=True, exist_ok=True)

            reg = MSItoHERegistrar(
                msi_adata=adata,
                he_image_path=config['he_image_path'],
                output_dir=sample_out,
                channel_subset=config.get('channel_subset'),
                channel_selection_method=channel_selection_method,
                n_channels=n_channels,
                composite_method=composite_method,
                max_image_dim=max_image_dim,
                affine_seed=config.get('affine_seed'),
            )

            result = reg.register()

            # Apply transform to AnnData and save
            adata = reg.transform_adata(adata)
            warped_path = sample_out / f'{sid}_msi_he_aligned.h5ad'
            adata.write_h5ad(str(warped_path))

            # Generate QC plot
            reg.plot_qc()

            results.append({
                'sample_id': sid,
                'status': 'success',
                'composite_method': result.get('composite_method', 'unknown'),
                'mean_error': result.get('mean_error', -1),
                'n_channels': len(reg.selected_channels) if reg.selected_channels else 0,
                'output_path': str(warped_path),
            })

        except Exception as e:
            logger.error(f"Failed for {sid}: {e}")
            import traceback
            traceback.print_exc()
            results.append({
                'sample_id': sid,
                'status': f'error: {str(e)[:80]}',
                'composite_method': '',
                'mean_error': -1,
                'n_channels': 0,
                'output_path': '',
            })

    summary = pd.DataFrame(results)
    summary.to_csv(output_base / 'batch_registration_summary.csv', index=False)
    logger.info(f"\nBatch registration complete. Summary saved to {output_base}")

    return summary
