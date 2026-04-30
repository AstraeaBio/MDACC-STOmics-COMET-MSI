"""
Zoomed Cell-Level Comparison: COMET Labels vs STOmics CellBin
==============================================================

Generates tight-zoom (500px window) panels showing individual cell
boundaries from both segmentations overlaid on STOmics DAPI. Uses
polygon outlines directly from GeoJSON (not rasterized masks) for
crisp cell boundary visualization.

Panels per ROI:
  1. COMET label cells (cyan) on DAPI — with phenotype color coding
  2. STOmics cellbin cells (yellow) on DAPI
  3. Side-by-side overlay — both segmentations on same DAPI crop

Auto-selects ROIs where BOTH segmentations have high cell density.

Usage:
    python script10_cellbin_comparison_zoomed.py --sample SO4
    python script10_cellbin_comparison_zoomed.py --sample SO4 --window 300
    python script10_cellbin_comparison_zoomed.py --sample SO4 --roi 7000,5500
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import cv2
import matplotlib.pyplot as plt
import numpy as np
import tifffile
from matplotlib.collections import PatchCollection
from matplotlib.patches import Polygon as MplPolygon

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from comet.alignment_utils import ALL_SAMPLES, DefaultPaths

logger = logging.getLogger(__name__)

WORK_DIR = DefaultPaths.OUTPUT_BASE / 'comet_stomics_alignment'
OUTPUT_DIR = DefaultPaths.OUTPUT_BASE / 'cellbin_comparison_zoomed'
DEFAULT_WINDOW = 500
DEFAULT_N_ROIS = 4

# Phenotype color map (subset of common phenotypes)
PHENOTYPE_COLORS = {
    'TUMOR CELLS CD56+': '#e6194b',
    'TUMOR CELLS PanCK+': '#e6194b',
    'CD8+ T CELLS': '#3cb44b',
    'CD4+ T CELLS': '#4363d8',
    'CD68+ MACROPHAGES': '#f58231',
    'CD163+ MACROPHAGES': '#911eb4',
    'CD20+ B CELLS': '#42d4f4',
    'CD45RO+ MEMORY': '#f032e6',
    'FOXP3+ TREGS': '#bfef45',
    'CD11C+ DCs': '#fabebe',
    'VESSELS CD31+': '#469990',
    'FIBROBLASTS aSMA+': '#dcbeff',
    'STROMA': '#aaffc3',
}
DEFAULT_CELL_COLOR = '#00ffff'  # cyan fallback


def load_tiff_roi(path, y0, y1, x0, x1):
    """Load a rectangular ROI from a TIFF."""
    path = Path(path)
    with tifffile.TiffFile(str(path)) as tif:
        h, w = tif.pages[0].shape[:2]
    y0, y1 = max(0, y0), min(h, y1)
    x0, x1 = max(0, x0), min(w, x1)
    img = tifffile.imread(str(path))
    return img[y0:y1, x0:x1].copy()


def normalize_dapi(crop):
    """Normalize DAPI crop to float [0,1] with contrast stretch."""
    crop = crop.astype(np.float32)
    if crop.ndim == 3:
        crop = np.mean(crop, axis=-1)
    vals = crop[crop > 0]
    if len(vals) > 0:
        p1, p99 = np.percentile(vals, [1, 99.5])
        crop = np.clip((crop - p1) / max(p99 - p1, 1), 0, 1)
    return crop


def read_json_chunked(path):
    """Read large JSON in chunks (Windows network drive compatible)."""
    chunks = []
    with open(path, 'r', encoding='utf-8') as f:
        while True:
            chunk = f.read(64 * 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
    return json.loads(''.join(chunks))


def extract_polygons_in_roi(geojson_data, x0, y0, x1, y1, margin=50):
    """Extract polygon coordinates and metadata for features within ROI."""
    polys = []
    for feat in geojson_data.get('features', []):
        geom = feat.get('geometry', {})
        props = feat.get('properties', {})
        coords_list = geom.get('coordinates', [])
        if not coords_list:
            continue

        # Get exterior ring
        exterior = np.array(coords_list[0], dtype=np.float64)
        if len(exterior) < 3:
            continue

        # Centroid check
        cx, cy = np.mean(exterior[:, 0]), np.mean(exterior[:, 1])
        if not (x0 - margin <= cx <= x1 + margin and
                y0 - margin <= cy <= y1 + margin):
            continue

        # Get phenotype from classification
        phenotype = None
        cls = props.get('classification', {})
        if isinstance(cls, dict):
            phenotype = cls.get('name')

        polys.append({
            'coords': exterior,
            'phenotype': phenotype,
            'label': props.get('object_index', props.get('label', '')),
        })

    return polys


def extract_cellbin_polygons_in_roi(geojson_data, x0, y0, x1, y1, margin=50):
    """Extract STOmics cellbin polygons within ROI."""
    polys = []
    for feat in geojson_data.get('features', []):
        geom = feat.get('geometry', {})
        coords_list = geom.get('coordinates', [])
        if not coords_list:
            continue

        exterior = np.array(coords_list[0], dtype=np.float64)
        if len(exterior) < 3:
            continue

        cx, cy = np.mean(exterior[:, 0]), np.mean(exterior[:, 1])
        if not (x0 - margin <= cx <= x1 + margin and
                y0 - margin <= cy <= y1 + margin):
            continue

        polys.append(exterior)

    return polys


def find_dense_rois(comet_geojson, cellbin_geojson, dapi_shape,
                    window=DEFAULT_WINDOW, n_rois=DEFAULT_N_ROIS):
    """Find ROI centers where both segmentations have high cell density."""
    h, w = dapi_shape[:2]
    step = window * 3  # coarse grid

    # Build centroid arrays for both
    comet_centroids = []
    for feat in comet_geojson.get('features', []):
        coords = feat.get('geometry', {}).get('coordinates', [[]])
        ext = coords[0] if coords else []
        if len(ext) >= 3:
            arr = np.array(ext, dtype=np.float64)
            comet_centroids.append([np.mean(arr[:, 0]), np.mean(arr[:, 1])])
    comet_centroids = np.array(comet_centroids) if comet_centroids else np.empty((0, 2))

    cellbin_centroids = []
    for feat in cellbin_geojson.get('features', []):
        coords = feat.get('geometry', {}).get('coordinates', [[]])
        ext = coords[0] if coords else []
        if len(ext) >= 3:
            arr = np.array(ext, dtype=np.float64)
            cellbin_centroids.append([np.mean(arr[:, 0]), np.mean(arr[:, 1])])
    cellbin_centroids = np.array(cellbin_centroids) if cellbin_centroids else np.empty((0, 2))

    logger.info(f"COMET centroids: {len(comet_centroids):,}")
    logger.info(f"Cellbin centroids: {len(cellbin_centroids):,}")

    candidates = []
    for cy in range(window + 100, h - window - 100, step):
        for cx in range(window + 100, w - window - 100, step):
            x0, x1_ = cx - window, cx + window
            y0, y1_ = cy - window, cy + window

            if len(comet_centroids) > 0:
                cm_mask = ((comet_centroids[:, 0] >= x0) &
                           (comet_centroids[:, 0] < x1_) &
                           (comet_centroids[:, 1] >= y0) &
                           (comet_centroids[:, 1] < y1_))
                n_comet = cm_mask.sum()
            else:
                n_comet = 0

            if len(cellbin_centroids) > 0:
                cb_mask = ((cellbin_centroids[:, 0] >= x0) &
                           (cellbin_centroids[:, 0] < x1_) &
                           (cellbin_centroids[:, 1] >= y0) &
                           (cellbin_centroids[:, 1] < y1_))
                n_cellbin = cb_mask.sum()
            else:
                n_cellbin = 0

            if n_comet > 20 and n_cellbin > 20:
                candidates.append((cx, cy, n_comet, n_cellbin,
                                   min(n_comet, n_cellbin)))

    # Sort by min(comet, cellbin) — want both to be dense
    candidates.sort(key=lambda r: r[4], reverse=True)

    # Pick top N with minimum separation
    rois = []
    min_sep = window * 2.5
    for cx, cy, nc, ns, _ in candidates:
        if len(rois) >= n_rois:
            break
        if all(abs(cx - rx) > min_sep or abs(cy - ry) > min_sep
               for _, rx, ry, _, _ in rois):
            rois.append((f'roi_{len(rois)}', cx, cy, nc, ns))

    return rois


def render_zoomed_comparison(stomics_dapi_path, comet_polys, cellbin_polys,
                             roi_name, roi_x, roi_y, window,
                             sample_id, chip_id, output_dir):
    """Render 3-panel zoomed comparison for one ROI."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    x0, x1 = roi_x - window, roi_x + window
    y0, y1 = roi_y - window, roi_y + window
    extent = [x0, x1, y0, y1]

    # Load and normalize DAPI crop
    dapi_crop = load_tiff_roi(stomics_dapi_path, y0, y1, x0, x1)
    dapi_norm = normalize_dapi(dapi_crop)

    def draw_comet_polys(ax, polys, use_phenotype_colors=True):
        if not polys:
            return
        patches = []
        colors = []
        for p in polys:
            pts = p['coords']
            patches.append(MplPolygon(pts, closed=True))
            if use_phenotype_colors and p['phenotype']:
                colors.append(PHENOTYPE_COLORS.get(p['phenotype'],
                                                    DEFAULT_CELL_COLOR))
            else:
                colors.append(DEFAULT_CELL_COLOR)
        ax.add_collection(PatchCollection(
            patches, facecolor='none', edgecolors=colors,
            linewidths=0.8, alpha=0.9))

    def draw_cellbin_polys(ax, polys, color='#ffff00', lw=0.8):
        if not polys:
            return
        patches = [MplPolygon(p, closed=True) for p in polys]
        ax.add_collection(PatchCollection(
            patches, facecolor='none', edgecolor=color,
            linewidth=lw, alpha=0.9))

    # ── 3-panel figure ──
    fig, axes = plt.subplots(1, 3, figsize=(30, 10), facecolor='black')
    fig.subplots_adjust(wspace=0.05)

    for ax in axes:
        ax.set_facecolor('black')
        ax.tick_params(colors='white', labelsize=7)
        for spine in ax.spines.values():
            spine.set_color('white')
        ax.set_xlim(x0, x1)
        ax.set_ylim(y1, y0)

    # Panel 1: COMET label cells (phenotype-colored) on DAPI
    ax = axes[0]
    ax.imshow(dapi_norm, cmap='gray', extent=extent, origin='upper',
              aspect='equal', vmin=0, vmax=1)
    draw_comet_polys(ax, comet_polys, use_phenotype_colors=True)
    ax.set_title(f'COMET label cells (n={len(comet_polys):,})',
                 color='cyan', fontsize=13, fontweight='bold', pad=8)

    # Panel 2: STOmics cellbin (yellow) on DAPI
    ax = axes[1]
    ax.imshow(dapi_norm, cmap='gray', extent=extent, origin='upper',
              aspect='equal', vmin=0, vmax=1)
    draw_cellbin_polys(ax, cellbin_polys)
    ax.set_title(f'STOmics cellbin (n={len(cellbin_polys):,})',
                 color='yellow', fontsize=13, fontweight='bold', pad=8)

    # Panel 3: Both overlaid
    ax = axes[2]
    ax.imshow(dapi_norm, cmap='gray', extent=extent, origin='upper',
              aspect='equal', vmin=0, vmax=1)
    draw_comet_polys(ax, comet_polys, use_phenotype_colors=False)
    draw_cellbin_polys(ax, cellbin_polys, color='#ffff00', lw=0.6)
    ax.set_title(f'Overlay (cyan=COMET, yellow=cellbin)',
                 color='white', fontsize=13, fontweight='bold', pad=8)

    fig.suptitle(
        f'{sample_id} / {chip_id}  |  {roi_name}  '
        f'(center={roi_x},{roi_y}  window={window*2}px)',
        color='white', fontsize=15, fontweight='bold', y=1.02)

    # Legend for phenotypes present
    present_phenos = set(p['phenotype'] for p in comet_polys
                         if p['phenotype'])
    if present_phenos:
        from matplotlib.lines import Line2D
        legend_entries = []
        for ph in sorted(present_phenos):
            c = PHENOTYPE_COLORS.get(ph, DEFAULT_CELL_COLOR)
            legend_entries.append(
                Line2D([0], [0], color=c, lw=2, label=ph))
        if len(legend_entries) <= 12:
            axes[0].legend(handles=legend_entries, loc='lower left',
                           fontsize=6, framealpha=0.7,
                           facecolor='black', labelcolor='white',
                           edgecolor='gray')

    fname = f'{sample_id}_{roi_name}_zoomed.png'
    out_path = output_dir / fname
    fig.savefig(str(out_path), dpi=250, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close(fig)
    logger.info(f"  Saved: {fname} (COMET={len(comet_polys)}, "
                f"cellbin={len(cellbin_polys)})")
    return out_path


def process_sample(sample_id, chip_id=None, work_dir=WORK_DIR,
                   output_dir=OUTPUT_DIR, n_rois=DEFAULT_N_ROIS,
                   window=DEFAULT_WINDOW, manual_rois=None):
    """Generate zoomed comparison for one sample."""
    if chip_id is None:
        if sample_id not in ALL_SAMPLES:
            logger.error(f"Unknown sample: {sample_id}")
            return []
        chip_id = ALL_SAMPLES[sample_id]['chip_id']

    sample_dir = Path(work_dir) / f'{sample_id}_{chip_id}'
    if not sample_dir.exists():
        logger.error(f"Sample directory not found: {sample_dir}")
        return []

    sample_out = Path(output_dir) / f'{sample_id}_{chip_id}'

    logger.info(f"{'='*60}")
    logger.info(f"Zoomed comparison: {sample_id} ({chip_id})")
    logger.info(f"{'='*60}")

    # 1. Find DAPI
    stomics_dapi_path = None
    for candidate in [
        sample_dir / 'qupath_valis' / 'input_images' / f'ssDNA_{chip_id}_regist.tif',
        sample_dir / 'input_images' / f'ssDNA_{chip_id}_regist.tif',
        DefaultPaths.stomics_dapi_path(chip_id),
    ]:
        if candidate.exists():
            stomics_dapi_path = candidate
            break
    if stomics_dapi_path is None:
        logger.error("STOmics DAPI not found")
        return []
    logger.info(f"DAPI: {stomics_dapi_path}")

    with tifffile.TiffFile(str(stomics_dapi_path)) as tif:
        dapi_shape = tif.pages[0].shape

    # 2. Load COMET warped label GeoJSON (cells with phenotypes)
    comet_geojson_path = sample_dir / f'{sample_id}_warped_label.geojson'
    if not comet_geojson_path.exists():
        # Fallback to warped_segmentations (masks-level, no phenotype)
        comet_geojson_path = sample_dir / f'{sample_id}_warped_segmentations.geojson'
        if not comet_geojson_path.exists():
            comet_geojson_path = sample_dir / f'{sample_id}_affine_warped.geojson'

    if not comet_geojson_path.exists():
        logger.error(f"No COMET GeoJSON found in {sample_dir}")
        return []
    logger.info(f"COMET GeoJSON: {comet_geojson_path.name}")
    comet_geojson = read_json_chunked(str(comet_geojson_path))
    logger.info(f"  {len(comet_geojson.get('features', []))} features")

    # 3. Load STOmics cellbin GeoJSON
    cellbin_geojson_path = sample_dir / f'{sample_id}_stomics_cellbin_borders.geojson'
    if not cellbin_geojson_path.exists():
        logger.error(f"STOmics cellbin GeoJSON not found: {cellbin_geojson_path}")
        logger.info("Run export_cellbin_borders_to_geojson() first")
        return []
    logger.info(f"Cellbin GeoJSON: {cellbin_geojson_path.name}")
    cellbin_geojson = read_json_chunked(str(cellbin_geojson_path))
    logger.info(f"  {len(cellbin_geojson.get('features', []))} features")

    # 4. Select ROIs
    if manual_rois:
        rois = [(f'manual_{i}', x, y, 0, 0)
                for i, (x, y) in enumerate(manual_rois)]
        logger.info(f"Using {len(rois)} manual ROI(s)")
    else:
        logger.info("Auto-selecting ROIs with dual coverage...")
        rois = find_dense_rois(comet_geojson, cellbin_geojson,
                               dapi_shape, window=window, n_rois=n_rois)
        if not rois:
            logger.error("No ROIs found with dual coverage")
            return []

    for name, rx, ry, nc, ns in rois:
        logger.info(f"  {name}: ({rx},{ry}) COMET~{nc} cellbin~{ns}")

    # 5. Render each ROI
    outputs = []
    for roi_name, roi_x, roi_y, _, _ in rois:
        x0, x1 = roi_x - window, roi_x + window
        y0, y1 = roi_y - window, roi_y + window

        comet_polys = extract_polygons_in_roi(
            comet_geojson, x0, y0, x1, y1, margin=50)
        cellbin_polys = extract_cellbin_polygons_in_roi(
            cellbin_geojson, x0, y0, x1, y1, margin=50)

        logger.info(f"  {roi_name}: COMET={len(comet_polys)}, "
                     f"cellbin={len(cellbin_polys)}")

        out = render_zoomed_comparison(
            stomics_dapi_path, comet_polys, cellbin_polys,
            roi_name, roi_x, roi_y, window,
            sample_id, chip_id, sample_out)
        outputs.append(out)

    # 6. Summary figure (all ROIs on DAPI overview)
    render_roi_overview(stomics_dapi_path, dapi_shape, rois, window,
                        sample_id, chip_id, sample_out)

    logger.info(f"Done: {len(outputs)} zoomed comparisons -> {sample_out}")
    return outputs


def render_roi_overview(dapi_path, dapi_shape, rois, window,
                        sample_id, chip_id, output_dir):
    """Low-res overview showing ROI locations."""
    h, w = dapi_shape[:2]
    ds = max(1, max(h, w) // 2000)

    dapi = tifffile.imread(str(dapi_path))
    dapi_small = dapi[::ds, ::ds].astype(np.float32)
    if dapi_small.ndim == 3:
        dapi_small = np.mean(dapi_small, axis=-1)
    vals = dapi_small[dapi_small > 0]
    if len(vals) > 0:
        p1, p99 = np.percentile(vals, [1, 99])
        dapi_small = np.clip((dapi_small - p1) / max(p99 - p1, 1), 0, 1)

    fig, ax = plt.subplots(1, 1, figsize=(12, 12), facecolor='black')
    ax.set_facecolor('black')
    ax.imshow(dapi_small, cmap='gray', extent=[0, w, h, 0],
              aspect='equal', vmin=0, vmax=1)

    colors = plt.cm.Set1(np.linspace(0, 1, max(len(rois), 1)))
    for i, (name, rx, ry, nc, ns) in enumerate(rois):
        rect = plt.Rectangle((rx - window, ry - window),
                              2 * window, 2 * window,
                              linewidth=2, edgecolor=colors[i],
                              facecolor='none')
        ax.add_patch(rect)
        ax.text(rx - window + 10, ry - window + 40, name,
                color=colors[i], fontsize=9, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.2', fc='black', alpha=0.7))

    ax.set_title(f'{sample_id} / {chip_id} — Zoomed ROI Overview',
                 color='white', fontsize=14, fontweight='bold')
    ax.tick_params(colors='white', labelsize=8)
    for spine in ax.spines.values():
        spine.set_color('white')

    out_path = Path(output_dir) / f'{sample_id}_zoomed_overview.png'
    fig.savefig(str(out_path), dpi=150, bbox_inches='tight',
                facecolor='black', edgecolor='none')
    plt.close(fig)
    logger.info(f"  Saved overview: {out_path.name}")


def main():
    parser = argparse.ArgumentParser(
        description='Zoomed COMET label vs STOmics cellbin comparison')
    parser.add_argument('--sample', type=str, required=True)
    parser.add_argument('--chip', type=str, default=None)
    parser.add_argument('--roi', type=str,
                        help='"x,y" or "x1,y1;x2,y2" in STOmics pixels')
    parser.add_argument('--window', type=int, default=DEFAULT_WINDOW)
    parser.add_argument('--n-rois', type=int, default=DEFAULT_N_ROIS)
    parser.add_argument('--work-dir', type=str, default=str(WORK_DIR))
    parser.add_argument('--output-dir', type=str, default=str(OUTPUT_DIR))
    parser.add_argument('-v', '--verbose', action='store_true')
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s %(levelname)-8s %(message)s', datefmt='%H:%M:%S')

    manual_rois = None
    if args.roi:
        manual_rois = []
        for pair in args.roi.split(';'):
            x, y = pair.strip().split(',')
            manual_rois.append((int(x.strip()), int(y.strip())))

    process_sample(args.sample, chip_id=args.chip,
                   work_dir=args.work_dir, output_dir=args.output_dir,
                   n_rois=args.n_rois, window=args.window,
                   manual_rois=manual_rois)


if __name__ == '__main__':
    main()
