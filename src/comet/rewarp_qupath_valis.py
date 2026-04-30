"""
Re-warp COMET GeoJSON using QuPath affine + VALIS non-rigid registration.
==========================================================================

The original pipeline may have used plain VALIS registration (without the
QuPath affine pre-warp), producing poorly aligned warped GeoJSON files.
This script re-generates the warped GeoJSON by:

  1. Loading the original COMET GeoJSON (in COMET pixel space)
  2. Applying the QuPath affine transform (coarse alignment)
  3. Loading the saved VALIS registrar from qupath_valis/
  4. Applying VALIS non-rigid warp to the affine-warped GeoJSON
  5. Saving the corrected warped GeoJSON

The qupath_valis registrar pickle must already exist (from a prior pipeline
run with method='qupath_valis'). This script does NOT re-run registration.

Usage:
    # Must run in maitra_spatial conda env (has VALIS installed)
    conda run -n maitra_spatial python script09_rewarp_qupath_valis.py --sample SO6
    conda run -n maitra_spatial python script09_rewarp_qupath_valis.py --batch
"""

import argparse
import json
import logging
import pickle
import sys
from copy import deepcopy
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from comet.alignment_utils import (
    ALL_SAMPLES,
    DefaultPaths,
    warp_geojson_with_valis as _warp_geojson_with_valis,
)
from comet.registration import (
    register_qupath_affine,
    apply_affine_to_geojson,
)

logger = logging.getLogger(__name__)

WORK_DIR = DefaultPaths.OUTPUT_BASE / 'comet_stomics_alignment'
TRANSFORMS_PATH = DefaultPaths.REPO_DIR / 'config' / 'qupath_transforms.json'


# ── VALIS Registrar Loading ─────────────────────────────────────────────────

def load_qupath_valis_registrar(sample_dir):
    """Load the saved VALIS registrar from qupath_valis/ subdirectory."""
    pickle_path = (sample_dir / 'qupath_valis' / 'valis_output' /
                   'input_images' / 'data' / 'input_images_registrar.pickle')
    if not pickle_path.exists():
        logger.error(f"Registrar pickle not found: {pickle_path}")
        return None

    logger.info(f"Loading VALIS registrar: {pickle_path} "
                f"({pickle_path.stat().st_size / 1e6:.0f} MB)")

    # Copy to local temp to avoid network drive read limits
    import tempfile, shutil
    with tempfile.NamedTemporaryFile(suffix='.pickle', delete=False) as tmp:
        tmp_path = tmp.name
    try:
        logger.info(f"  Copying to local temp: {tmp_path}")
        shutil.copy2(str(pickle_path), tmp_path)
        with open(tmp_path, 'rb') as f:
            reg = pickle.load(f)
    finally:
        Path(tmp_path).unlink(missing_ok=True)

    slides = list(reg.slide_dict.values())
    logger.info(f"Registrar slides: {[s.name for s in slides]}")
    return reg


# ── Load GeoJSON with chunked reading ────────────────────────────────────────

def load_geojson_chunked(path):
    """Load a large GeoJSON file using chunked reads (Windows compat)."""
    path = Path(path)
    chunks = []
    with open(path, 'r', encoding='utf-8') as f:
        while True:
            chunk = f.read(64 * 1024 * 1024)
            if not chunk:
                break
            chunks.append(chunk)
    text = ''.join(chunks)
    return json.loads(text)


# ── Per-Sample Driver ────────────────────────────────────────────────────────

def rewarp_sample(sample_id, chip_id=None, work_dir=WORK_DIR,
                  transforms_path=TRANSFORMS_PATH, dry_run=False):
    """Re-warp GeoJSON for one sample using QuPath affine + VALIS.

    Args:
        sample_id: e.g., 'SO6'
        chip_id: If None, looked up from ALL_SAMPLES
        work_dir: Pipeline output root
        transforms_path: Path to qupath_transforms.json
        dry_run: If True, only check prerequisites without warping

    Returns:
        Path to the new warped GeoJSON, or None on failure
    """
    if chip_id is None:
        if sample_id not in ALL_SAMPLES:
            logger.error(f"Unknown sample: {sample_id}")
            return None
        chip_id = ALL_SAMPLES[sample_id]['chip_id']

    sample_dir = Path(work_dir) / f'{sample_id}_{chip_id}'
    aligned = ALL_SAMPLES.get(sample_id, {}).get('aligned', False)

    logger.info(f"{'='*60}")
    logger.info(f"Re-warping {sample_id} ({chip_id})")
    logger.info(f"{'='*60}")

    # Check prerequisites
    geojson_path = DefaultPaths.geojson_path(sample_id, aligned=aligned)
    if not geojson_path.exists():
        logger.error(f"Original GeoJSON not found: {geojson_path}")
        return None

    try:
        qupath_matrix = register_qupath_affine(sample_id, transforms_path)
    except (KeyError, FileNotFoundError) as e:
        logger.error(f"QuPath transform error: {e}")
        return None
    logger.info(f"QuPath affine matrix:\n{qupath_matrix}")

    registrar = load_qupath_valis_registrar(sample_dir)
    if registrar is None:
        return None

    if dry_run:
        logger.info(f"Dry run: prerequisites OK for {sample_id}")
        return geojson_path

    # Step 1: Load original GeoJSON
    logger.info(f"Loading original GeoJSON: {geojson_path.name}")
    geojson_data = load_geojson_chunked(geojson_path)
    n_features = len(geojson_data.get('features', []))
    logger.info(f"  {n_features} features loaded")

    # Step 2: Apply QuPath affine
    logger.info("Applying QuPath affine transform...")
    affine_path = sample_dir / f'{sample_id}_affine_warped.geojson'
    apply_affine_to_geojson(geojson_data, qupath_matrix,
                            output_path=affine_path)
    del geojson_data

    # Step 3: Apply VALIS non-rigid warp
    logger.info("Applying VALIS non-rigid warp...")
    warped_path = sample_dir / f'{sample_id}_warped_segmentations.geojson'

    # Back up old warped GeoJSON if it exists
    if warped_path.exists():
        backup = sample_dir / f'{sample_id}_warped_segmentations.OLD.geojson'
        if not backup.exists():
            import shutil
            shutil.move(str(warped_path), str(backup))
            logger.info(f"Backed up old warped GeoJSON -> {backup.name}")
        else:
            logger.info("Old backup already exists, overwriting warped GeoJSON")

    warped_geojson = _warp_geojson_with_valis(registrar, affine_path,
                                               warped_path)

    logger.info(f"Done: {sample_id} -> {warped_path}")
    return warped_path


# ── Batch Driver ─────────────────────────────────────────────────────────────

def run_batch(work_dir=WORK_DIR, transforms_path=TRANSFORMS_PATH,
              dry_run=False):
    """Re-warp all samples that have qupath_valis registrar data."""
    results = {}

    for sample_dir in sorted(Path(work_dir).iterdir()):
        if not sample_dir.is_dir():
            continue
        parts = sample_dir.name.split('_', 1)
        if len(parts) != 2:
            continue
        sample_id, chip_id = parts

        # Check if qupath_valis registrar exists
        pickle_path = (sample_dir / 'qupath_valis' / 'valis_output' /
                       'input_images' / 'data' / 'input_images_registrar.pickle')
        if not pickle_path.exists():
            logger.info(f"Skipping {sample_id}: no qupath_valis registrar")
            continue

        try:
            out = rewarp_sample(sample_id, chip_id,
                                work_dir=work_dir,
                                transforms_path=transforms_path,
                                dry_run=dry_run)
            results[sample_id] = 'OK' if out else 'FAILED'
        except Exception as e:
            logger.error(f"Failed {sample_id}: {e}", exc_info=True)
            results[sample_id] = f'ERROR: {e}'

    # Summary
    logger.info(f"\n{'='*60}")
    logger.info("Batch re-warp summary:")
    ok = sum(1 for v in results.values() if v == 'OK')
    logger.info(f"  {ok}/{len(results)} succeeded")
    for sid, status in sorted(results.items()):
        logger.info(f"  {sid}: {status}")

    return results


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Re-warp COMET GeoJSON using QuPath + VALIS registration',
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--sample', type=str,
                        help='Sample ID (e.g., SO6)')
    parser.add_argument('--batch', action='store_true',
                        help='Process all samples with qupath_valis data')
    parser.add_argument('--dry-run', action='store_true',
                        help='Check prerequisites only')
    parser.add_argument('--work-dir', type=str, default=str(WORK_DIR))
    parser.add_argument('--transforms', type=str, default=str(TRANSFORMS_PATH))
    parser.add_argument('-v', '--verbose', action='store_true')

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%H:%M:%S',
    )

    if not args.batch and not args.sample:
        parser.error('Specify --sample SOXX or --batch')

    if args.batch:
        run_batch(work_dir=args.work_dir,
                  transforms_path=args.transforms,
                  dry_run=args.dry_run)
    else:
        rewarp_sample(args.sample,
                      work_dir=args.work_dir,
                      transforms_path=args.transforms,
                      dry_run=args.dry_run)


if __name__ == '__main__':
    main()
