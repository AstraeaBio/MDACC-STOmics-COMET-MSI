"""
Alignment QC Snapshot Generator
================================

Generates high-resolution zoomed snapshots for visual QC of COMET-STOmics
alignment. Each snapshot shows multiple panels at a selected ROI:

  1. DAPI blend (magenta=COMET, green=STOmics) for registration QC
  2. STOmics DAPI + COMET cell outlines (cyan, QuPath-affine warped)
  3. STOmics DAPI + STOmics cellbin outlines (yellow, from mask)
  4. Transcript density + both segmentations

Uses the QuPath affine transform to warp COMET polygons into STOmics
space on-the-fly (fast matrix multiply, no VALIS needed). For visual QC
at the ROI window scale this is sufficient — the VALIS non-rigid
correction is small relative to the window size.

Usage:
    # Single sample
    python script08_alignment_qc_snapshots.py --sample SO4

    # Custom ROI center (in STOmics pixel coords)
    python script08_alignment_qc_snapshots.py --sample SO4 --roi 12000,8000

    # Batch: all samples with qupath_valis data
    python script08_alignment_qc_snapshots.py --batch

    # Adjust window size and number of auto-ROIs
    python script08_alignment_qc_snapshots.py --sample SO4 --window 2000 --n-rois 6
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import cv2
import h5py
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from matplotlib.collections import PatchCollection
from matplotlib.colors import Normalize
from matplotlib.patches import Polygon as MplPolygon
from scipy.ndimage import gaussian_filter

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from comet.alignment_utils import ALL_SAMPLES, DefaultPaths
from comet.registration import load_qupath_transforms

logger = logging.getLogger(__name__)

WORK_DIR = DefaultPaths.OUTPUT_BASE / 'comet_stomics_alignment'
QC_OUTPUT_DIR = DefaultPaths.OUTPUT_BASE / 'alignment_qc_snapshots'
TRANSFORMS_PATH = DefaultPaths.REPO_DIR / 'config' / 'qupath_transforms.json'
DEFAULT_WINDOW = 1500
DEFAULT_N_ROIS = 4
TRANSCRIPT_BIN_SIZE = 10
DENSITY_SIGMA = 3.0


# ── I/O Helpers (Windows network drive compatible) ───────────────────────────

def _memmap_tiff(path):
    """Memory-map an uncompressed single-strip TIFF. Returns memmap array."""
    path = Path(path)
    with tifffile.TiffFile(str(path)) as tif:
        page = tif.pages[0]
        if page.compression != 1 or len(page.dataoffsets) != 1:
            # Compressed/multi-strip: use tifffile (works for small strips)
            return tifffile.imread(str(path))
        return np.memmap(str(path), dtype=page.dtype, mode='r',
                         offset=page.dataoffsets[0], shape=page.shape)


def _load_tiff_roi(path, y0, y1, x0, x1):
    """Load a rectangular ROI from a TIFF. Uses memmap for uncompressed,
    tifffile for compressed (reads full then crops)."""
    path = Path(path)
    with tifffile.TiffFile(str(path)) as tif:
        page = tif.pages[0]
        h, w = page.shape[:2]

    y0, y1 = max(0, y0), min(h, y1)
    x0, x1 = max(0, x0), min(w, x1)

    try:
        mm = _memmap_tiff(path)
        crop = np.array(mm[y0:y1, x0:x1])
        if isinstance(mm, np.memmap):
            del mm
        return crop
    except (OSError, ValueError):
        # Fallback: read full image (works for compressed multi-strip TIFFs)
        img = tifffile.imread(str(path))
        return img[y0:y1, x0:x1].copy()


def _read_json_chunked(path):
    """Read a large JSON file in chunks (avoids Windows read limit)."""
    chunks = []
    with open(path, 'r', encoding='utf-8') as f:
        while True:
            chunk = f.read(64 * 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
    return json.loads(''.join(chunks))


# ── Data Loading ─────────────────────────────────────────────────────────────

def load_stomics_dapi(sample_dir, chip_id):
    """Load STOmics DAPI as memmap (lazy access, fast)."""
    for candidate in [
        sample_dir / 'qupath_valis' / 'input_images' / f'ssDNA_{chip_id}_regist.tif',
        sample_dir / 'input_images' / f'ssDNA_{chip_id}_regist.tif',
        DefaultPaths.stomics_dapi_path(chip_id),
    ]:
        if candidate.exists():
            logger.info(f"STOmics DAPI: {candidate.name}")
            return _memmap_tiff(candidate)
    logger.error(f"STOmics DAPI not found for {chip_id}")
    return None


def load_comet_prewarped_dapi(sample_dir, sample_id):
    """Load the QuPath pre-warped COMET DAPI (already in STOmics space)."""
    path = sample_dir / 'qupath_valis' / 'input_images' / f'{sample_id}_prewarped.tif'
    if not path.exists():
        return None
    logger.info(f"COMET pre-warped DAPI: {path.name}")
    return _memmap_tiff(path)


def load_original_geojson(sample_id, aligned=False):
    """Load original COMET GeoJSON (in COMET pixel space)."""
    path = DefaultPaths.geojson_path(sample_id, aligned=aligned)
    if not path.exists():
        logger.warning(f"Original GeoJSON not found: {path}")
        return None
    logger.info(f"Loading original GeoJSON: {path.name}")
    return _read_json_chunked(path)


def load_qupath_affine(sample_id, modality='comet'):
    """Load the forward QuPath affine matrix (2x3) for a sample."""
    try:
        sample_cfg = load_qupath_transforms(str(TRANSFORMS_PATH), sample_id)
        matrix = sample_cfg['transforms'][modality]['matrix']
        logger.info(f"QuPath affine loaded for {sample_id}")
        return np.array(matrix, dtype=np.float64)
    except (KeyError, FileNotFoundError) as e:
        logger.warning(f"No QuPath transform for {sample_id}: {e}")
        return None


def load_transcript_positions(chip_id):
    """Load raw DNB transcript positions from bin1 tissue GEF."""
    gef_path = DefaultPaths.stomics_tissue_gef_path(chip_id)
    if not gef_path.exists():
        logger.warning(f"Tissue GEF not found: {gef_path}")
        return None, None, None

    logger.info(f"Loading bin1 transcripts: {gef_path.name}")
    with h5py.File(str(gef_path), 'r') as f:
        bin1 = f['geneExp']['bin1']
        expression = bin1['expression']
        n_records = expression.shape[0]
        n_genes = len(bin1['gene'])
        logger.info(f"  {n_genes} genes, {n_records:,} DNB records")

        chunk_size = 5_000_000
        all_x = np.empty(n_records, dtype=np.int32)
        all_y = np.empty(n_records, dtype=np.int32)
        all_counts = np.empty(n_records, dtype=np.int32)
        for i in range(0, n_records, chunk_size):
            end = min(i + chunk_size, n_records)
            chunk = expression[i:end]
            all_x[i:end] = chunk['x']
            all_y[i:end] = chunk['y']
            all_counts[i:end] = chunk['count']

    logger.info(f"  Loaded {n_records:,} transcript positions")
    return all_x, all_y, all_counts


# ── Affine Warp of GeoJSON ───────────────────────────────────────────────────

def warp_polygons_in_roi(geojson_data, matrix_2x3, roi_x, roi_y, window):
    """Apply affine to GeoJSON polygons and return those within the ROI.

    Only warps polygons whose centroids (after transform) fall in the ROI,
    skipping the rest. This is much faster than warping the full GeoJSON.
    """
    M = np.array(matrix_2x3, dtype=np.float64)
    x0, x1 = roi_x - window, roi_x + window
    y0, y1 = roi_y - window, roi_y + window

    # Expand ROI slightly for centroid check (cells at edge)
    margin = window * 0.1
    polys = []
    for feat in geojson_data.get('features', []):
        coords = np.array(feat['geometry']['coordinates'][0], dtype=np.float64)
        # Quick centroid check in original space isn't helpful since
        # we don't know where it maps to. Warp centroid first.
        cx, cy = np.mean(coords[:, 0]), np.mean(coords[:, 1])
        warped_c = M @ np.array([cx, cy, 1.0])
        if not (x0 - margin <= warped_c[0] <= x1 + margin and
                y0 - margin <= warped_c[1] <= y1 + margin):
            continue
        # Warp all vertices
        ones = np.ones((len(coords), 1))
        aug = np.hstack([coords, ones])
        warped = (M @ aug.T).T
        polys.append(warped)

    return polys


# ── Mask Contour Extraction ─────────────────────────────────────────────────

def extract_contours_from_mask(mask_path, roi_x, roi_y, window):
    """Extract cell contours from a labeled mask TIFF within an ROI."""
    if mask_path is None or not Path(mask_path).exists():
        return [], 0

    crop = _load_tiff_roi(mask_path, roi_y - window, roi_y + window,
                          roi_x - window, roi_x + window)

    with tifffile.TiffFile(str(mask_path)) as tif:
        h, w = tif.pages[0].shape[:2]
    actual_y0 = max(0, roi_y - window)
    actual_x0 = max(0, roi_x - window)

    binary = (crop > 0).astype(np.uint8)
    n_cells = len(np.unique(crop)) - 1

    contours, _ = cv2.findContours(binary, cv2.RETR_EXTERNAL,
                                    cv2.CHAIN_APPROX_SIMPLE)
    abs_contours = []
    for cnt in contours:
        if len(cnt) < 3:
            continue
        pts = cnt.squeeze()
        if pts.ndim < 2:
            continue
        abs_pts = pts.astype(np.float64)
        abs_pts[:, 0] += actual_x0
        abs_pts[:, 1] += actual_y0
        abs_contours.append(abs_pts)

    return abs_contours, n_cells


# ── Transcript Density ───────────────────────────────────────────────────────

def build_transcript_density(tx, ty, t_counts, roi_x, roi_y, window,
                             bin_size=TRANSCRIPT_BIN_SIZE):
    """Build 2D transcript density heatmap for an ROI window."""
    x0, x1 = roi_x - window, roi_x + window
    y0, y1 = roi_y - window, roi_y + window

    mask = (tx >= x0) & (tx < x1) & (ty >= y0) & (ty < y1)
    lx, ly, lc = tx[mask], ty[mask], t_counts[mask]

    n_bx = max(1, int((x1 - x0) / bin_size))
    n_by = max(1, int((y1 - y0) / bin_size))

    if len(lx) == 0:
        return np.zeros((n_by, n_bx)), (x0, x1, y0, y1)

    density, _, _ = np.histogram2d(
        ly.astype(np.float64), lx.astype(np.float64),
        bins=[n_by, n_bx], range=[[y0, y1], [x0, x1]],
        weights=lc.astype(np.float64),
    )
    density = gaussian_filter(density, sigma=DENSITY_SIGMA)
    return density, (x0, x1, y0, y1)


# ── ROI Selection ────────────────────────────────────────────────────────────

def select_rois(dapi_shape, tx, ty, t_counts, n_rois=DEFAULT_N_ROIS,
                window=DEFAULT_WINDOW):
    """Auto-select ROI centers from transcript density."""
    h, w = dapi_shape[:2]
    rois = []
    coarse_bin = 200

    if tx is not None and len(tx) > 0:
        n_bx = max(1, w // coarse_bin)
        n_by = max(1, h // coarse_bin)
        density, _, _ = np.histogram2d(
            ty.astype(np.float64), tx.astype(np.float64),
            bins=[n_by, n_bx], range=[[0, h], [0, w]],
            weights=t_counts.astype(np.float64),
        )

        # Pick top-N density peaks with minimum separation
        flat_order = np.argsort(density.ravel())[::-1]
        for idx in flat_order:
            if len(rois) >= n_rois:
                break
            iy, ix = np.unravel_index(idx, density.shape)
            rx = int(np.clip((ix + 0.5) * coarse_bin, window, w - window))
            ry = int(np.clip((iy + 0.5) * coarse_bin, window, h - window))
            if all(abs(rx - ex) > window or abs(ry - ey) > window
                   for _, ex, ey in rois):
                label = 'peak' if len(rois) == 0 else f'region_{len(rois)}'
                rois.append((label, rx, ry))

    # Fill remaining with random positions
    if len(rois) < n_rois:
        rng = np.random.default_rng(42)
        for i in range(n_rois - len(rois)):
            rx = int(rng.uniform(window, max(window + 1, w - window)))
            ry = int(rng.uniform(window, max(window + 1, h - window)))
            rois.append((f'random_{len(rois)}', rx, ry))

    return rois


# ── Rendering ────────────────────────────────────────────────────────────────

def _normalize_crop(img, y0, y1, x0, x1):
    """Extract and normalize an image ROI to float [0,1]."""
    crop = np.array(img[y0:y1, x0:x1], dtype=np.float32)
    if crop.ndim == 3:
        crop = np.mean(crop, axis=-1)
    vals = crop[crop > 0]
    if len(vals) > 0:
        p1, p99 = np.percentile(vals, [1, 99.5])
        crop = np.clip((crop - p1) / max(p99 - p1, 1), 0, 1)
    return crop


def render_snapshot(stomics_dapi, comet_prewarped_dapi,
                    comet_mask_path, cellbin_mask_path,
                    tx, ty, t_counts,
                    roi_name, roi_x, roi_y, window,
                    sample_id, chip_id, output_dir):
    """Render a 4-panel QC snapshot for one ROI."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    h, w = stomics_dapi.shape[:2]
    y0 = max(0, roi_y - window)
    y1 = min(h, roi_y + window)
    x0 = max(0, roi_x - window)
    x1 = min(w, roi_x + window)
    extent = [x0, x1, y0, y1]

    # Normalize DAPI crops
    st_crop = _normalize_crop(stomics_dapi, y0, y1, x0, x1)

    # DAPI blend (magenta=COMET, green=STOmics)
    has_blend = comet_prewarped_dapi is not None
    if has_blend:
        ch, cw = comet_prewarped_dapi.shape[:2]
        cy0 = max(0, min(y0, ch))
        cy1 = max(0, min(y1, ch))
        cx0 = max(0, min(x0, cw))
        cx1 = max(0, min(x1, cw))
        if cy1 > cy0 and cx1 > cx0:
            comet_crop = _normalize_crop(comet_prewarped_dapi, cy0, cy1, cx0, cx1)
            # Ensure same size
            target_h, target_w = st_crop.shape[:2]
            if comet_crop.shape != st_crop.shape:
                comet_crop = cv2.resize(comet_crop, (target_w, target_h))
        else:
            has_blend = False

    # COMET cell contours from rasterized mask (proper shapes, already in STOmics space)
    comet_contours, n_comet = extract_contours_from_mask(
        comet_mask_path, roi_x, roi_y, window)

    # STOmics cellbin contours from rasterized mask
    cellbin_contours, n_cellbin = extract_contours_from_mask(
        cellbin_mask_path, roi_x, roi_y, window)

    # Transcript density
    has_tx = tx is not None and len(tx) > 0
    density = None
    n_umi, n_dnb = 0, 0
    if has_tx:
        density, _ = build_transcript_density(tx, ty, t_counts,
                                              roi_x, roi_y, window)
        t_mask = (tx >= x0) & (tx < x1) & (ty >= y0) & (ty < y1)
        n_umi = int(t_counts[t_mask].sum())
        n_dnb = int(t_mask.sum())

    # ── Figure ────────────────────────────────────────────────────────
    fig, axes = plt.subplots(2, 2, figsize=(20, 20), facecolor='black')
    fig.subplots_adjust(hspace=0.08, wspace=0.08)
    for ax in axes.flat:
        ax.set_facecolor('black')
        ax.tick_params(colors='white', labelsize=7)
        for spine in ax.spines.values():
            spine.set_color('white')

    def _draw_polys(ax, polys, color, lw=0.6, alpha=0.8):
        if not polys:
            return
        patches = [MplPolygon(p, closed=True) for p in polys]
        ax.add_collection(PatchCollection(
            patches, facecolor='none', edgecolor=color,
            linewidth=lw, alpha=alpha))

    def _draw_contours(ax, contours, color, lw=0.6, alpha=0.8):
        if not contours:
            return
        patches = [MplPolygon(c, closed=True) for c in contours]
        ax.add_collection(PatchCollection(
            patches, facecolor='none', edgecolor=color,
            linewidth=lw, alpha=alpha))

    # Panel A: DAPI blend
    ax = axes[0, 0]
    if has_blend:
        blend = np.zeros((*st_crop.shape, 3))
        blend[..., 0] = comet_crop   # R = COMET (magenta)
        blend[..., 1] = st_crop      # G = STOmics
        blend[..., 2] = comet_crop   # B = COMET (magenta)
        ax.imshow(blend, extent=extent, origin='upper', aspect='equal')
        ax.set_title('DAPI blend (magenta=COMET, green=STOmics)',
                     color='white', fontsize=11, fontweight='bold', pad=8)
    else:
        ax.imshow(st_crop, cmap='gray', extent=extent, origin='upper',
                  aspect='equal', vmin=0, vmax=1)
        ax.set_title('STOmics DAPI (no COMET pre-warp available)',
                     color='white', fontsize=11, fontweight='bold', pad=8)
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)

    # Panel B: DAPI + COMET cells (from rasterized mask in STOmics space)
    ax = axes[0, 1]
    ax.imshow(st_crop, cmap='gray', extent=extent, origin='upper',
              aspect='equal', vmin=0, vmax=1)
    _draw_contours(ax, comet_contours, 'cyan')
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_title(f'COMET cells (n={n_comet:,}, from mask)',
                 color='cyan', fontsize=12, fontweight='bold', pad=8)

    # Panel C: DAPI + STOmics cellbin
    ax = axes[1, 0]
    ax.imshow(st_crop, cmap='gray', extent=extent, origin='upper',
              aspect='equal', vmin=0, vmax=1)
    _draw_contours(ax, cellbin_contours, 'yellow')
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_title(f'STOmics cellbin (n={n_cellbin:,})',
                 color='yellow', fontsize=12, fontweight='bold', pad=8)

    # Panel D: Combined with transcript density
    ax = axes[1, 1]
    ax.imshow(st_crop, cmap='gray', extent=extent, origin='upper',
              aspect='equal', vmin=0, vmax=1, alpha=0.7)
    if density is not None and density.max() > 0:
        ax.imshow(density, cmap='hot', extent=[x0, x1, y1, y0],
                  origin='upper', aspect='equal', alpha=0.35,
                  norm=Normalize(vmin=0, vmax=np.percentile(density, 99)))
    _draw_contours(ax, comet_contours, 'cyan', lw=0.5, alpha=0.7)
    _draw_contours(ax, cellbin_contours, 'yellow', lw=0.5, alpha=0.7)
    ax.set_xlim(x0, x1); ax.set_ylim(y1, y0)
    ax.set_title(f'Combined ({n_dnb:,} DNBs, {n_umi:,} UMIs)',
                 color='white', fontsize=11, fontweight='bold', pad=8)

    fig.suptitle(
        f'{sample_id} / {chip_id}  |  ROI: {roi_name}  '
        f'(center={roi_x},{roi_y}  window={window*2}px)',
        color='white', fontsize=14, fontweight='bold', y=0.995)

    fname = f'{sample_id}_{roi_name}_qc.png'
    out_path = output_dir / fname
    fig.savefig(str(out_path), dpi=200, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close(fig)
    logger.info(f"  Saved: {fname} ({n_comet} COMET, {n_cellbin} cellbin, "
                f"{n_umi:,} UMIs)")
    return out_path


# ── Per-Sample Driver ────────────────────────────────────────────────────────

def process_sample(sample_id, chip_id=None, work_dir=WORK_DIR,
                   output_dir=QC_OUTPUT_DIR, n_rois=DEFAULT_N_ROIS,
                   window=DEFAULT_WINDOW, manual_rois=None):
    """Generate QC snapshots for one sample."""
    if chip_id is None:
        if sample_id not in ALL_SAMPLES:
            logger.error(f"Unknown sample: {sample_id}")
            return []
        chip_id = ALL_SAMPLES[sample_id]['chip_id']
        aligned = ALL_SAMPLES[sample_id].get('aligned', False)
    else:
        aligned = False

    sample_dir = Path(work_dir) / f'{sample_id}_{chip_id}'
    if not sample_dir.exists():
        logger.error(f"Sample directory not found: {sample_dir}")
        return []

    sample_out = Path(output_dir) / f'{sample_id}_{chip_id}'
    sample_out.mkdir(parents=True, exist_ok=True)

    logger.info(f"{'='*60}")
    logger.info(f"Processing {sample_id} ({chip_id})")
    logger.info(f"{'='*60}")

    # 1. Load STOmics DAPI (memmap, instant)
    stomics_dapi = load_stomics_dapi(sample_dir, chip_id)
    if stomics_dapi is None:
        return []

    # 2. Load COMET pre-warped DAPI for blend (memmap, instant)
    comet_dapi = load_comet_prewarped_dapi(sample_dir, sample_id)

    # 3. Find mask paths (contours extracted per-ROI, no full load)
    comet_mask_path = sample_dir / f'{sample_id}_comet_cell_mask.tif'
    if not comet_mask_path.exists():
        comet_mask_path = None
        logger.info("No COMET cell mask found")
    else:
        logger.info(f"COMET mask: {comet_mask_path.name}")

    cellbin_mask_path = sample_dir / f'{sample_id}_stomics_cellbin_mask.tif'
    if not cellbin_mask_path.exists():
        cellbin_mask_path = None
        logger.info("No cellbin mask found")
    else:
        logger.info(f"Cellbin mask: {cellbin_mask_path.name}")

    # 5. Load transcripts
    logger.info("Loading transcript positions...")
    tx, ty, t_counts = load_transcript_positions(chip_id)

    # 6. Select ROIs
    if manual_rois:
        rois = [(f'manual_{i}', x, y) for i, (x, y) in enumerate(manual_rois)]
    else:
        rois = select_rois(stomics_dapi.shape, tx, ty, t_counts,
                           n_rois=n_rois, window=window)

    logger.info(f"Selected {len(rois)} ROIs:")
    for name, rx, ry in rois:
        logger.info(f"  {name}: ({rx}, {ry})")

    # 7. Render
    outputs = []
    for roi_name, roi_x, roi_y in rois:
        out = render_snapshot(
            stomics_dapi=stomics_dapi,
            comet_prewarped_dapi=comet_dapi,
            comet_mask_path=comet_mask_path,
            cellbin_mask_path=cellbin_mask_path,
            tx=tx, ty=ty, t_counts=t_counts,
            roi_name=roi_name, roi_x=roi_x, roi_y=roi_y,
            window=window,
            sample_id=sample_id, chip_id=chip_id,
            output_dir=sample_out,
        )
        outputs.append(out)

    # 8. Overview
    _render_overview(stomics_dapi, rois, window,
                     sample_id, chip_id, sample_out)

    logger.info(f"Done: {sample_id} — {len(outputs)} snapshots -> {sample_out}")
    return outputs


def _render_overview(dapi, rois, window, sample_id, chip_id, output_dir):
    """Low-res overview showing ROI locations on tissue."""
    h, w = dapi.shape[:2]
    ds = max(1, max(h, w) // 2000)
    dapi_small = np.array(dapi[::ds, ::ds], dtype=np.float32)
    vals = dapi_small[dapi_small > 0]
    if len(vals) > 0:
        p1, p99 = np.percentile(vals, [1, 99])
        dapi_small = np.clip((dapi_small - p1) / max(p99 - p1, 1), 0, 1)

    fig, ax = plt.subplots(1, 1, figsize=(12, 12), facecolor='black')
    ax.set_facecolor('black')
    ax.imshow(dapi_small, cmap='gray', extent=[0, w, h, 0],
              aspect='equal', vmin=0, vmax=1)

    colors = plt.cm.Set1(np.linspace(0, 1, max(len(rois), 1)))
    for i, (name, rx, ry) in enumerate(rois):
        rect = plt.Rectangle((rx - window, ry - window), 2 * window, 2 * window,
                              linewidth=2, edgecolor=colors[i], facecolor='none')
        ax.add_patch(rect)
        ax.text(rx - window + 30, ry - window + 80, name,
                color=colors[i], fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', fc='black', alpha=0.7))

    ax.set_title(f'{sample_id} / {chip_id} — ROI Overview',
                 color='white', fontsize=14, fontweight='bold')
    ax.tick_params(colors='white', labelsize=8)
    for spine in ax.spines.values():
        spine.set_color('white')

    out_path = output_dir / f'{sample_id}_roi_overview.png'
    fig.savefig(str(out_path), dpi=150, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close(fig)
    logger.info(f"  Saved overview: {out_path.name}")


# ── Batch + CLI ──────────────────────────────────────────────────────────────

def run_batch(work_dir=WORK_DIR, output_dir=QC_OUTPUT_DIR,
              n_rois=DEFAULT_N_ROIS, window=DEFAULT_WINDOW):
    """Run QC for all samples with qupath_valis data."""
    results = {}
    for sd in sorted(Path(work_dir).iterdir()):
        if not sd.is_dir():
            continue
        parts = sd.name.split('_', 1)
        if len(parts) != 2:
            continue
        sample_id, chip_id = parts
        if not (sd / 'qupath_valis').exists():
            continue
        try:
            outs = process_sample(sample_id, chip_id, work_dir=work_dir,
                                  output_dir=output_dir, n_rois=n_rois,
                                  window=window)
            results[sample_id] = len(outs)
        except Exception as e:
            logger.error(f"Failed {sample_id}: {e}", exc_info=True)
            results[sample_id] = 0

    ok = sum(1 for v in results.values() if v > 0)
    logger.info(f"\nBatch: {ok}/{len(results)} samples succeeded")
    for sid, n in sorted(results.items()):
        logger.info(f"  {sid}: {n} snapshots" if n > 0 else f"  {sid}: FAILED")
    return results


def main():
    parser = argparse.ArgumentParser(
        description='Alignment QC snapshots for COMET-STOmics pipeline')
    parser.add_argument('--sample', type=str)
    parser.add_argument('--batch', action='store_true')
    parser.add_argument('--roi', type=str,
                        help='"x,y" or "x1,y1;x2,y2" in STOmics pixels')
    parser.add_argument('--window', type=int, default=DEFAULT_WINDOW)
    parser.add_argument('--n-rois', type=int, default=DEFAULT_N_ROIS)
    parser.add_argument('--work-dir', type=str, default=str(WORK_DIR))
    parser.add_argument('--output-dir', type=str, default=str(QC_OUTPUT_DIR))
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%H:%M:%S')

    if not args.batch and not args.sample:
        parser.error('Specify --sample SOXX or --batch')

    manual_rois = None
    if args.roi:
        manual_rois = []
        for pair in args.roi.split(';'):
            x, y = pair.strip().split(',')
            manual_rois.append((int(x.strip()), int(y.strip())))

    if args.batch:
        run_batch(args.work_dir, args.output_dir, args.n_rois, args.window)
    else:
        process_sample(args.sample, work_dir=args.work_dir,
                       output_dir=args.output_dir, n_rois=args.n_rois,
                       window=args.window, manual_rois=manual_rois)


if __name__ == '__main__':
    main()
