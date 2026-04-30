"""Alignment Validation & Signoff System
=========================================

Provides visual QC and a gating mechanism so that downstream analysis only
proceeds with validated, signed-off alignments. Prevents the silent use of
stale or wrong VALIS warp files.

Design:
    1. Each alignment produces a `validation_manifest.json` that records
       exactly which files were used (registrar, affine, source images,
       warp outputs) plus a SHA-256 fingerprint of each.
    2. QC visualizations are generated automatically: a global overview
       and N local ROI zoom panels.
    3. The manifest starts with `"status": "pending"`. A human reviews
       the QC images and calls `sign_off()` (or the notebook helper
       `interactive_signoff()`) to flip it to `"approved"`.
    4. Downstream functions (normalization, DE, export) check
       `is_approved()` and refuse to proceed on unapproved samples.

Usage:
    # After alignment
    val = AlignmentValidator(sample_dir='path/to/SO34_A03979E2')
    val.generate_qc()          # produces PNG panels
    val.interactive_signoff()   # opens images, asks for approval

    # Before downstream
    if val.is_approved():
        run_normalization(...)
    else:
        print(f"Alignment not signed off for {val.sample_id}")
"""

import hashlib
import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Literal, Optional, Tuple, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

try:
    import tifffile
except ImportError:
    tifffile = None

try:
    import matplotlib
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection
except ImportError:
    plt = None


# ============================================================
# File Fingerprinting
# ============================================================

def file_hash(path: Union[str, Path], quick: bool = True) -> str:
    """SHA-256 hash of a file. If quick=True, hashes first+last 1MB only."""
    path = Path(path)
    if not path.exists():
        return "MISSING"

    h = hashlib.sha256()
    size = path.stat().st_size

    if quick and size > 2 * 1024 * 1024:
        with open(path, 'rb') as f:
            h.update(f.read(1024 * 1024))
            f.seek(max(0, size - 1024 * 1024))
            h.update(f.read(1024 * 1024))
        h.update(str(size).encode())
    else:
        with open(path, 'rb') as f:
            for chunk in iter(lambda: f.read(65536), b''):
                h.update(chunk)

    return h.hexdigest()[:16]


# ============================================================
# Validation Manifest
# ============================================================

MANIFEST_NAME = 'validation_manifest.json'


def create_manifest(
    sample_dir: Union[str, Path],
    sample_id: str,
    chip_id: str,
    files_used: Dict[str, Union[str, Path]],
    registration_method: str = 'unknown',
    extra_metadata: Optional[dict] = None,
) -> dict:
    """Create a validation manifest that locks down which files were used.

    Parameters
    ----------
    sample_dir : path
        Output directory for this sample.
    sample_id : str
        e.g. 'SO34'
    chip_id : str
        e.g. 'A03979E2'
    files_used : dict
        Keys are role names, values are paths. Typical keys:
        'comet_dapi', 'stomics_dapi', 'geojson_source', 'registrar_dir',
        'warped_geojson', 'integrated_h5ad', 'comet_cell_mask', 'affine_json'
    registration_method : str
        e.g. 'qupath_valis', 'qupath_ecc', 'valis_auto'
    extra_metadata : dict, optional
        Any additional info to record.

    Returns
    -------
    dict
        The manifest contents.
    """
    sample_dir = Path(sample_dir)

    manifest = {
        'schema_version': 1,
        'sample_id': sample_id,
        'chip_id': chip_id,
        'registration_method': registration_method,
        'created': datetime.now().isoformat(),
        'status': 'pending',  # pending | approved | rejected
        'approved_by': None,
        'approved_at': None,
        'rejection_reason': None,
        'files': {},
    }

    for role, path in files_used.items():
        p = Path(path)
        manifest['files'][role] = {
            'path': str(p),
            'exists': p.exists(),
            'hash': file_hash(p) if p.exists() else None,
            'size_mb': round(p.stat().st_size / (1024 * 1024), 1) if p.exists() else None,
        }

    if extra_metadata:
        manifest['metadata'] = extra_metadata

    # Save
    manifest_path = sample_dir / MANIFEST_NAME
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)

    logger.info(f"Created manifest: {manifest_path}")
    return manifest


def load_manifest(sample_dir: Union[str, Path]) -> Optional[dict]:
    """Load the validation manifest for a sample directory."""
    path = Path(sample_dir) / MANIFEST_NAME
    if not path.exists():
        return None
    with open(path) as f:
        return json.load(f)


def verify_manifest(manifest: dict) -> Dict[str, bool]:
    """Check that all files in the manifest still match their recorded hashes.

    This catches the case where someone re-ran VALIS and overwrote the warp
    files but the old manifest (possibly approved) is still sitting there.

    Returns
    -------
    dict
        Per-file verification results. True = matches, False = changed/missing.
    """
    results = {}
    for role, info in manifest.get('files', {}).items():
        path = Path(info['path'])
        if not path.exists():
            results[role] = False
            continue

        current_hash = file_hash(path)
        results[role] = (current_hash == info.get('hash'))

    return results


# ============================================================
# QC Visualization
# ============================================================

def _load_geojson(path):
    """Load GeoJSON file."""
    with open(path) as f:
        return json.load(f)


def _extract_centroids(features, max_n=None):
    """Extract centroids from GeoJSON features."""
    centroids = []
    for feat in features[:max_n] if max_n else features:
        coords = np.array(feat['geometry']['coordinates'][0])
        cx = np.mean(coords[:-1, 0]) if len(coords) > 1 else coords[0, 0]
        cy = np.mean(coords[:-1, 1]) if len(coords) > 1 else coords[0, 1]
        centroids.append([cx, cy])
    return np.array(centroids) if centroids else np.empty((0, 2))


def generate_global_qc(
    sample_dir: Union[str, Path],
    stomics_dapi_path: Optional[Union[str, Path]] = None,
    comet_dapi_path: Optional[Union[str, Path]] = None,
    warped_geojson_path: Optional[Union[str, Path]] = None,
    stomics_cellbin_mask_path: Optional[Union[str, Path]] = None,
    comet_cell_mask_path: Optional[Union[str, Path]] = None,
    sample_id: str = '',
    max_cells: int = 10000,
    dpi: int = 150,
) -> Path:
    """Generate a global overview QC panel.

    4-panel layout:
        [1] DAPI blend: magenta=COMET, green=STOmics (registration quality)
        [2] STOmics DAPI + warped COMET cell outlines (cyan)
        [3] COMET cell mask vs STOmics cellbin mask coverage
        [4] Centroid distance histogram (alignment error distribution)

    Parameters
    ----------
    sample_dir : path
        Output directory for this sample.
    stomics_dapi_path, comet_dapi_path : path, optional
        DAPI images. Auto-discovered from sample_dir if not provided.
    warped_geojson_path : path, optional
        Warped COMET segmentations.

    Returns
    -------
    Path to saved QC image.
    """
    if plt is None:
        raise ImportError("matplotlib required for QC visualization")

    sample_dir = Path(sample_dir)

    # Auto-discover files if not provided
    if warped_geojson_path is None:
        candidates = list(sample_dir.glob('*_warped_label.geojson')) + \
                     list(sample_dir.glob('*_warped_segmentations.geojson'))
        warped_geojson_path = candidates[0] if candidates else None

    if stomics_dapi_path is None:
        candidates = list(sample_dir.glob('*stomics*dapi*.tif')) + \
                     list(sample_dir.glob('ssDNA_*_regist.tif'))
        stomics_dapi_path = candidates[0] if candidates else None

    if comet_dapi_path is None:
        candidates = list(sample_dir.glob('*_dapi.tif')) + \
                     list(sample_dir.glob('*comet*dapi*.tif'))
        comet_dapi_path = candidates[0] if candidates else None

    fig, axes = plt.subplots(2, 2, figsize=(18, 18))

    # --- Panel 1: DAPI Blend ---
    ax = axes[0, 0]
    has_blend = False
    if stomics_dapi_path and comet_dapi_path and tifffile:
        try:
            st_dapi = tifffile.imread(str(stomics_dapi_path))
            co_dapi = tifffile.imread(str(comet_dapi_path))

            # Normalize both to 0-1
            st_norm = st_dapi.astype(np.float32)
            st_norm = (st_norm - st_norm.min()) / (st_norm.max() - st_norm.min() + 1e-10)
            co_norm = co_dapi.astype(np.float32)
            co_norm = (co_norm - co_norm.min()) / (co_norm.max() - co_norm.min() + 1e-10)

            # If different sizes, note it (registration should have handled this)
            if st_norm.shape != co_norm.shape:
                from skimage.transform import resize
                co_norm = resize(co_norm, st_norm.shape, preserve_range=True)

            # Blend: magenta=COMET, green=STOmics
            blend = np.stack([co_norm, st_norm, co_norm], axis=2)
            ax.imshow(np.clip(blend, 0, 1))
            ax.set_title('DAPI Blend (magenta=COMET, green=STOmics)', fontsize=11)
            has_blend = True
        except Exception as e:
            ax.text(0.5, 0.5, f'DAPI blend failed:\n{e}',
                    ha='center', va='center', transform=ax.transAxes, fontsize=10)
    else:
        ax.text(0.5, 0.5, 'DAPI images not available',
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
    ax.axis('off')

    # --- Panel 2: STOmics DAPI + COMET outlines ---
    ax = axes[0, 1]
    if stomics_dapi_path and warped_geojson_path and tifffile:
        try:
            if 'st_dapi' not in dir() or st_dapi is None:
                st_dapi = tifffile.imread(str(stomics_dapi_path))
            ax.imshow(st_dapi, cmap='gray', alpha=0.8)

            geojson = _load_geojson(warped_geojson_path)
            features = geojson['features']
            # Subsample for rendering speed
            if len(features) > max_cells:
                idx = np.random.choice(len(features), max_cells, replace=False)
                features_sub = [features[i] for i in sorted(idx)]
            else:
                features_sub = features

            patches = []
            for feat in features_sub:
                coords = feat['geometry']['coordinates']
                if coords and coords[0]:
                    exterior = np.array(coords[0])
                    patches.append(MplPolygon(exterior, closed=True))

            if patches:
                pc = PatchCollection(patches, alpha=0.15, edgecolor='cyan',
                                     facecolor='none', linewidth=0.3)
                ax.add_collection(pc)
            ax.set_title(f'STOmics DAPI + COMET Outlines (n={len(features):,})', fontsize=11)
        except Exception as e:
            ax.text(0.5, 0.5, f'Overlay failed:\n{e}',
                    ha='center', va='center', transform=ax.transAxes, fontsize=10)
    else:
        ax.text(0.5, 0.5, 'Missing DAPI or warped GeoJSON',
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
    ax.axis('off')

    # --- Panel 3: Mask coverage comparison ---
    ax = axes[1, 0]
    if comet_cell_mask_path and stomics_cellbin_mask_path and tifffile:
        try:
            comet_mask = tifffile.imread(str(comet_cell_mask_path)) > 0
            st_mask = tifffile.imread(str(stomics_cellbin_mask_path)) > 0

            if comet_mask.shape != st_mask.shape:
                from skimage.transform import resize
                st_mask = resize(st_mask, comet_mask.shape, preserve_range=True) > 0.5

            # Color: red=COMET only, blue=STOmics only, yellow=overlap
            overlay = np.zeros((*comet_mask.shape, 3), dtype=np.float32)
            overlay[comet_mask & ~st_mask] = [1, 0.3, 0.3]     # red: COMET only
            overlay[~comet_mask & st_mask] = [0.3, 0.3, 1]     # blue: STOmics only
            overlay[comet_mask & st_mask] = [1, 1, 0.3]        # yellow: overlap

            # Compute Jaccard
            intersection = (comet_mask & st_mask).sum()
            union = (comet_mask | st_mask).sum()
            jaccard = intersection / union if union > 0 else 0

            ax.imshow(overlay)
            ax.set_title(f'Mask Coverage (Jaccard={jaccard:.3f})\n'
                         f'Red=COMET only, Blue=STOmics only, Yellow=overlap', fontsize=10)
        except Exception as e:
            ax.text(0.5, 0.5, f'Mask comparison failed:\n{e}',
                    ha='center', va='center', transform=ax.transAxes, fontsize=10)
    else:
        ax.text(0.5, 0.5, 'Cell masks not available\n(generate with rasterize step)',
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
    ax.axis('off')

    # --- Panel 4: Centroid distance histogram ---
    ax = axes[1, 1]
    if warped_geojson_path:
        try:
            if 'geojson' not in dir() or geojson is None:
                geojson = _load_geojson(warped_geojson_path)

            # Try to load integrated h5ad for STOmics cell coords
            h5ad_candidates = list(sample_dir.glob('*_comet_stomics_integrated.h5ad'))
            if h5ad_candidates:
                import anndata as ad
                integrated = ad.read_h5ad(str(h5ad_candidates[0]))
                if 'spatial' in integrated.obsm:
                    # Use a subsample of COMET centroids and find nearest STOmics cell
                    from scipy.spatial import cKDTree
                    comet_cents = _extract_centroids(geojson['features'], max_n=50000)
                    st_coords = integrated.obsm['spatial']
                    tree = cKDTree(st_coords)
                    distances, _ = tree.query(comet_cents, k=1)

                    ax.hist(distances, bins=100, range=(0, 200), color='steelblue',
                            edgecolor='none', alpha=0.8)
                    median_d = np.median(distances)
                    pct_30 = (distances < 30).sum() / len(distances) * 100
                    pct_50 = (distances < 50).sum() / len(distances) * 100

                    # Quality classification
                    if median_d < 15:
                        quality, color = 'EXCELLENT', 'green'
                    elif median_d < 30:
                        quality, color = 'GOOD', 'forestgreen'
                    elif median_d < 50:
                        quality, color = 'ACCEPTABLE', 'orange'
                    else:
                        quality, color = 'POOR', 'red'

                    ax.axvline(median_d, color='red', linestyle='--', linewidth=2,
                               label=f'Median={median_d:.1f}px')
                    ax.axvline(30, color='green', linestyle=':', alpha=0.6, label='30px')
                    ax.axvline(50, color='orange', linestyle=':', alpha=0.6, label='50px')
                    ax.set_title(f'COMET↔STOmics Distance: {quality}\n'
                                 f'Median={median_d:.1f}px, <30px={pct_30:.0f}%, <50px={pct_50:.0f}%',
                                 fontsize=10, color=color)
                    ax.legend(fontsize=9)
                else:
                    ax.text(0.5, 0.5, 'No spatial coords in integrated h5ad',
                            ha='center', va='center', transform=ax.transAxes)
            else:
                ax.text(0.5, 0.5, 'Integrated h5ad not found\n(run alignment pipeline first)',
                        ha='center', va='center', transform=ax.transAxes)
        except Exception as e:
            ax.text(0.5, 0.5, f'Distance computation failed:\n{e}',
                    ha='center', va='center', transform=ax.transAxes, fontsize=10)
    else:
        ax.text(0.5, 0.5, 'Warped GeoJSON not available',
                ha='center', va='center', transform=ax.transAxes, fontsize=12)
    ax.set_xlabel('Distance (pixels)')
    ax.set_ylabel('Count')

    plt.suptitle(f'Alignment QC: {sample_id}', fontsize=16, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    qc_path = sample_dir / f'{sample_id}_alignment_qc_global.png'
    fig.savefig(qc_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved global QC: {qc_path}")

    return qc_path


def generate_roi_qc(
    sample_dir: Union[str, Path],
    stomics_dapi_path: Optional[Union[str, Path]] = None,
    warped_geojson_path: Optional[Union[str, Path]] = None,
    sample_id: str = '',
    n_rois: int = 4,
    roi_size: int = 500,
    dpi: int = 150,
) -> Path:
    """Generate local ROI zoom panels for detailed alignment inspection.

    Auto-selects N regions with high cell density and produces tight-zoom
    panels showing DAPI + overlaid COMET cell boundaries.

    Parameters
    ----------
    sample_dir : path
        Output directory.
    n_rois : int
        Number of ROI panels.
    roi_size : int
        Half-width of each ROI window in pixels.

    Returns
    -------
    Path to saved QC image.
    """
    if plt is None:
        raise ImportError("matplotlib required")

    sample_dir = Path(sample_dir)

    # Auto-discover
    if warped_geojson_path is None:
        candidates = list(sample_dir.glob('*_warped_label.geojson')) + \
                     list(sample_dir.glob('*_warped_segmentations.geojson'))
        warped_geojson_path = candidates[0] if candidates else None

    if stomics_dapi_path is None:
        candidates = list(sample_dir.glob('*stomics*dapi*.tif')) + \
                     list(sample_dir.glob('ssDNA_*_regist.tif'))
        stomics_dapi_path = candidates[0] if candidates else None

    if not warped_geojson_path:
        logger.warning("No warped GeoJSON found — skipping ROI QC")
        return None

    geojson = _load_geojson(warped_geojson_path)
    features = geojson['features']
    centroids = _extract_centroids(features)

    if len(centroids) == 0:
        logger.warning("No cells in GeoJSON — skipping ROI QC")
        return None

    # Select ROI centers: find regions of highest cell density
    # Use a coarse grid to find hot spots
    from scipy.stats import binned_statistic_2d
    x_range = (centroids[:, 0].min(), centroids[:, 0].max())
    y_range = (centroids[:, 1].min(), centroids[:, 1].max())
    n_bins = 20

    counts, x_edges, y_edges, _ = binned_statistic_2d(
        centroids[:, 0], centroids[:, 1],
        values=None, statistic='count',
        bins=n_bins, range=[x_range, y_range],
        expand_binnumbers=True,
    )

    # Get top-N bins by cell count, spread them out spatially
    flat_idx = np.argsort(counts.ravel())[::-1]
    selected_centers = []
    min_spacing = roi_size * 2  # Don't overlap

    for idx in flat_idx:
        if len(selected_centers) >= n_rois:
            break
        ix, iy = np.unravel_index(idx, counts.shape)
        cx = (x_edges[ix] + x_edges[ix + 1]) / 2
        cy = (y_edges[iy] + y_edges[iy + 1]) / 2

        # Check spacing from already-selected centers
        too_close = False
        for sc in selected_centers:
            if np.sqrt((cx - sc[0])**2 + (cy - sc[1])**2) < min_spacing:
                too_close = True
                break

        if not too_close:
            selected_centers.append((cx, cy))

    n_actual = len(selected_centers)
    if n_actual == 0:
        logger.warning("Could not find suitable ROI centers")
        return None

    # Load DAPI for background
    dapi = None
    if stomics_dapi_path and tifffile:
        try:
            dapi = tifffile.imread(str(stomics_dapi_path))
        except Exception as e:
            logger.warning(f"Could not load DAPI: {e}")

    # Generate panels: 2 columns per ROI (DAPI+outlines, outlines only)
    n_cols = 2
    fig, axes = plt.subplots(n_actual, n_cols, figsize=(7 * n_cols, 7 * n_actual))
    if n_actual == 1:
        axes = axes.reshape(1, -1)

    for i, (cx, cy) in enumerate(selected_centers):
        x_lo, x_hi = cx - roi_size, cx + roi_size
        y_lo, y_hi = cy - roi_size, cy + roi_size

        # Filter features in this ROI
        roi_features = []
        for feat in features:
            coords = np.array(feat['geometry']['coordinates'][0])
            fc = np.mean(coords[:-1], axis=0) if len(coords) > 1 else coords[0]
            if x_lo <= fc[0] <= x_hi and y_lo <= fc[1] <= y_hi:
                roi_features.append(feat)

        # Panel A: DAPI + COMET outlines
        ax = axes[i, 0]
        if dapi is not None:
            # Crop DAPI to ROI
            y_lo_i = max(0, int(y_lo))
            y_hi_i = min(dapi.shape[0], int(y_hi))
            x_lo_i = max(0, int(x_lo))
            x_hi_i = min(dapi.shape[1], int(x_hi))
            roi_dapi = dapi[y_lo_i:y_hi_i, x_lo_i:x_hi_i]
            ax.imshow(roi_dapi, cmap='gray', extent=[x_lo_i, x_hi_i, y_hi_i, y_lo_i])

        patches = []
        colors = []
        for feat in roi_features:
            coords = np.array(feat['geometry']['coordinates'][0])
            patches.append(MplPolygon(coords, closed=True))
            # Color by phenotype if available
            pheno = feat.get('properties', {}).get('phenotype', '')
            if 'TUMOR' in pheno.upper():
                colors.append('red')
            elif 'T CELL' in pheno.upper() or 'CD8' in pheno.upper() or 'CD4' in pheno.upper():
                colors.append('lime')
            elif 'MACROPHAGE' in pheno.upper() or 'CD68' in pheno.upper():
                colors.append('yellow')
            else:
                colors.append('cyan')

        if patches:
            pc = PatchCollection(patches, alpha=0.25, edgecolors=colors,
                                 facecolor='none', linewidths=0.8)
            ax.add_collection(pc)

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_hi, y_lo)  # Inverted Y for image coords
        ax.set_title(f'ROI {i+1}: DAPI + COMET outlines ({len(roi_features)} cells)', fontsize=10)
        ax.set_aspect('equal')
        ax.axis('off')

        # Panel B: COMET outlines only (for shape assessment)
        ax = axes[i, 1]
        if patches:
            patches2 = []
            for feat in roi_features:
                coords = np.array(feat['geometry']['coordinates'][0])
                patches2.append(MplPolygon(coords, closed=True))

            pc2 = PatchCollection(patches2, alpha=0.4, edgecolors=colors,
                                  facecolors=[(c, 0.15) for c in colors] if len(colors) > 0 else 'cyan',
                                  linewidths=0.8)
            # Use individual colors
            for j, (patch, col) in enumerate(zip(patches2, colors)):
                patch.set_edgecolor(col)
                patch.set_facecolor((*matplotlib.colors.to_rgb(col), 0.15))
                ax.add_patch(patch)

        ax.set_xlim(x_lo, x_hi)
        ax.set_ylim(y_hi, y_lo)
        ax.set_title(f'ROI {i+1}: Cell boundaries colored by phenotype', fontsize=10)
        ax.set_aspect('equal')
        ax.set_facecolor('black')
        ax.axis('off')

    plt.suptitle(f'ROI Alignment QC: {sample_id}  (window={roi_size*2}px)',
                 fontsize=14, fontweight='bold', y=0.995)
    plt.tight_layout(rect=[0, 0, 1, 0.98])

    qc_path = sample_dir / f'{sample_id}_alignment_qc_rois.png'
    fig.savefig(qc_path, dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved ROI QC: {qc_path}")

    return qc_path


# ============================================================
# The Validator Class
# ============================================================

class AlignmentValidator:
    """Manages alignment validation and signoff for a single sample.

    Typical workflow:
        val = AlignmentValidator('path/to/SO34_A03979E2')
        val.create_manifest(files_used={...})  # record what was used
        val.generate_qc()                       # produce QC images
        val.interactive_signoff()                # human reviews + approves
        # --- later ---
        assert val.is_approved()                # gate downstream steps
    """

    def __init__(self, sample_dir: Union[str, Path]):
        self.sample_dir = Path(sample_dir)
        self.manifest = load_manifest(self.sample_dir)

        # Try to infer sample_id and chip_id from directory name
        name = self.sample_dir.name
        parts = name.split('_', 1)
        self.sample_id = parts[0] if parts else name
        self.chip_id = parts[1] if len(parts) > 1 else ''

    def create_manifest(self, files_used: Dict[str, Union[str, Path]],
                        registration_method: str = 'unknown',
                        extra_metadata: Optional[dict] = None):
        """Create a fresh manifest for this sample."""
        self.manifest = create_manifest(
            self.sample_dir,
            sample_id=self.sample_id,
            chip_id=self.chip_id,
            files_used=files_used,
            registration_method=registration_method,
            extra_metadata=extra_metadata,
        )

    def generate_qc(self, **kwargs) -> Tuple[Optional[Path], Optional[Path]]:
        """Generate both global and ROI QC panels."""
        global_path = generate_global_qc(
            self.sample_dir,
            sample_id=self.sample_id,
            **{k: v for k, v in kwargs.items()
               if k in ('stomics_dapi_path', 'comet_dapi_path', 'warped_geojson_path',
                         'stomics_cellbin_mask_path', 'comet_cell_mask_path',
                         'max_cells', 'dpi')},
        )

        roi_path = generate_roi_qc(
            self.sample_dir,
            sample_id=self.sample_id,
            **{k: v for k, v in kwargs.items()
               if k in ('stomics_dapi_path', 'warped_geojson_path',
                         'n_rois', 'roi_size', 'dpi')},
        )

        return global_path, roi_path

    def verify_files(self) -> Dict[str, bool]:
        """Check that all manifest files still match their recorded hashes."""
        if self.manifest is None:
            logger.warning("No manifest found — run create_manifest() first")
            return {}
        return verify_manifest(self.manifest)

    def is_approved(self, verify_hashes: bool = True) -> bool:
        """Check whether this alignment has been signed off.

        Parameters
        ----------
        verify_hashes : bool
            If True, also verifies that files haven't changed since approval.
            This catches the "wrong warp file" problem — if someone re-ran
            VALIS and the warped GeoJSON changed, the old approval is invalid.
        """
        if self.manifest is None:
            return False

        if self.manifest.get('status') != 'approved':
            return False

        if verify_hashes:
            checks = self.verify_files()
            if not all(checks.values()):
                changed = [k for k, v in checks.items() if not v]
                logger.warning(
                    f"APPROVAL INVALIDATED for {self.sample_id}: "
                    f"files changed since signoff: {changed}. "
                    f"Re-run QC and re-approve."
                )
                return False

        return True

    def sign_off(self, approved: bool = True, reviewer: str = 'unknown',
                 reason: str = ''):
        """Record approval or rejection.

        Parameters
        ----------
        approved : bool
            True to approve, False to reject.
        reviewer : str
            Who is signing off.
        reason : str
            Rejection reason (if rejected).
        """
        if self.manifest is None:
            raise RuntimeError("No manifest — run create_manifest() first")

        # Re-verify file hashes at signoff time
        checks = self.verify_files()
        if not all(checks.values()):
            changed = [k for k, v in checks.items() if not v]
            logger.warning(f"WARNING: Some files changed since manifest creation: {changed}")

        self.manifest['status'] = 'approved' if approved else 'rejected'
        self.manifest['approved_by'] = reviewer
        self.manifest['approved_at'] = datetime.now().isoformat()
        if not approved:
            self.manifest['rejection_reason'] = reason

        # Re-hash all files at approval time (so future verification catches changes)
        for role, info in self.manifest.get('files', {}).items():
            p = Path(info['path'])
            if p.exists():
                info['hash'] = file_hash(p)
                info['exists'] = True

        # Save
        manifest_path = self.sample_dir / MANIFEST_NAME
        with open(manifest_path, 'w') as f:
            json.dump(self.manifest, f, indent=2)

        status_str = 'APPROVED' if approved else f'REJECTED ({reason})'
        logger.info(f"Alignment {status_str} for {self.sample_id} by {reviewer}")

    def interactive_signoff(self, reviewer: str = 'Trevor'):
        """Interactive signoff in a Jupyter notebook.

        Shows the QC images inline and prompts for approval.
        """
        from IPython.display import display, Image, HTML

        # Show QC images
        global_qc = self.sample_dir / f'{self.sample_id}_alignment_qc_global.png'
        roi_qc = self.sample_dir / f'{self.sample_id}_alignment_qc_rois.png'

        if global_qc.exists():
            display(HTML(f'<h3>Global QC: {self.sample_id}</h3>'))
            display(Image(filename=str(global_qc), width=900))
        else:
            print(f"No global QC image found. Run generate_qc() first.")

        if roi_qc.exists():
            display(HTML(f'<h3>ROI QC: {self.sample_id}</h3>'))
            display(Image(filename=str(roi_qc), width=900))

        # File verification
        checks = self.verify_files()
        if checks:
            display(HTML('<h4>File Integrity</h4>'))
            for role, ok in checks.items():
                icon = '✓' if ok else '✗ CHANGED'
                color = 'green' if ok else 'red'
                display(HTML(f'<span style="color:{color}">{icon} {role}</span>'))

        # Prompt
        print('\n' + '='*60)
        response = input(f'Approve alignment for {self.sample_id}? [y/n/reason]: ').strip()

        if response.lower() in ('y', 'yes', 'approve'):
            self.sign_off(approved=True, reviewer=reviewer)
            print(f'APPROVED by {reviewer}')
        else:
            reason = response if response.lower() not in ('n', 'no') else 'Rejected by reviewer'
            self.sign_off(approved=False, reviewer=reviewer, reason=reason)
            print(f'REJECTED: {reason}')


# ============================================================
# Batch Validation Helpers
# ============================================================

def get_approval_status(
    alignment_dir: Union[str, Path],
    sample_dirs: Optional[List[str]] = None,
) -> pd.DataFrame:
    """Get approval status for all samples in an alignment output directory.

    Parameters
    ----------
    alignment_dir : path
        Base directory containing per-sample subdirectories.
    sample_dirs : list, optional
        Specific subdirectory names. If None, scans all.

    Returns
    -------
    pd.DataFrame
        Status summary per sample.
    """
    alignment_dir = Path(alignment_dir)

    if sample_dirs is None:
        sample_dirs = [d.name for d in alignment_dir.iterdir() if d.is_dir()]

    rows = []
    for dirname in sorted(sample_dirs):
        d = alignment_dir / dirname
        manifest = load_manifest(d)

        if manifest is None:
            rows.append({
                'directory': dirname,
                'sample_id': dirname.split('_')[0],
                'status': 'NO_MANIFEST',
                'registration_method': '',
                'approved_by': '',
                'files_ok': '',
            })
            continue

        # Verify files
        checks = verify_manifest(manifest)
        all_ok = all(checks.values()) if checks else False

        rows.append({
            'directory': dirname,
            'sample_id': manifest.get('sample_id', ''),
            'status': manifest.get('status', 'unknown'),
            'registration_method': manifest.get('registration_method', ''),
            'approved_by': manifest.get('approved_by', ''),
            'approved_at': manifest.get('approved_at', ''),
            'files_ok': 'YES' if all_ok else 'CHANGED' if checks else 'UNCHECKED',
            'n_files': len(manifest.get('files', {})),
        })

    return pd.DataFrame(rows)


def require_approval(sample_dir: Union[str, Path], action: str = 'downstream analysis'):
    """Gate function — raises RuntimeError if alignment is not approved.

    Use this at the start of any downstream step to enforce the signoff gate.

    Parameters
    ----------
    sample_dir : path
        Sample output directory.
    action : str
        Description of what's being attempted (for error message).

    Raises
    ------
    RuntimeError
        If alignment is not approved or files have changed since approval.
    """
    val = AlignmentValidator(sample_dir)

    if not val.is_approved():
        status = val.manifest.get('status', 'NO_MANIFEST') if val.manifest else 'NO_MANIFEST'

        if status == 'pending':
            msg = (f"Cannot proceed with {action} for {val.sample_id}: "
                   f"alignment has not been reviewed. "
                   f"Run QC and sign_off() first.")
        elif status == 'rejected':
            reason = val.manifest.get('rejection_reason', 'unknown')
            msg = (f"Cannot proceed with {action} for {val.sample_id}: "
                   f"alignment was REJECTED ({reason}). Fix and re-register.")
        elif status == 'approved':
            msg = (f"Cannot proceed with {action} for {val.sample_id}: "
                   f"alignment was previously approved but FILES HAVE CHANGED. "
                   f"Someone may have re-run VALIS. Re-review and re-approve.")
        else:
            msg = (f"Cannot proceed with {action} for {val.sample_id}: "
                   f"no validation manifest found. Run the alignment pipeline first.")

        raise RuntimeError(msg)
