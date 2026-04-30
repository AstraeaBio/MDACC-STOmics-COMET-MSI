#!/usr/bin/env python3
"""
run_batch_alignment.py - Standalone batch COMET-STOmics alignment
=================================================================

Runs the COMET-STOmics alignment pipeline across all samples with
available data. Designed for autonomous overnight execution.

Can be run from command line in a screen/tmux session:

    tmux new -s alignment
    python run_batch_alignment.py --sample-set all_with_data 2>&1 | tee alignment.log

Or for specific samples:
    python run_batch_alignment.py --samples SO24 SO25 SO30

Or just Endo/MLA:
    python run_batch_alignment.py --sample-set endo_mla

Options:
    --sample-set:   endo_mla | endo_only | mla_only | all_with_data | all_samples
    --samples:      Specific sample IDs to process
    --output-dir:   Override output directory
    --force:        Re-process samples that already have results
    --method:       Aggregation method: 'nearest' (fast) or 'polygon' (precise)
    --max-distance: Max distance for cell matching (default: 50 pixels)
    --dry-run:      Show what would be processed without running
"""

import argparse
import json
import logging
import sys
import traceback
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# Add repo root to path
REPO_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO_DIR))

from comet.alignment_utils import (
    DefaultPaths, ENDO_SAMPLES, MLA_SAMPLES, ALL_ENDO_MLA,
    ALL_SAMPLES, OVARIAN_WITH_STOMICS, SAMPLES_WITH_CELLBIN,
    check_sample_data, run_alignment_pipeline,
)

# ============================================================
# Logging
# ============================================================

def setup_logging(output_dir: Path, verbose: bool = False):
    """Configure logging to both console and file."""
    log_file = output_dir / f'batch_alignment_{datetime.now():%Y%m%d_%H%M%S}.log'

    handlers = [
        logging.StreamHandler(sys.stdout),
        logging.FileHandler(str(log_file)),
    ]

    logging.basicConfig(
        level=logging.DEBUG if verbose else logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%H:%M:%S',
        handlers=handlers,
    )
    return log_file


# ============================================================
# Sample Selection
# ============================================================

SAMPLE_SETS = {
    'endo_mla': ALL_ENDO_MLA,
    'endo_only': ENDO_SAMPLES,
    'mla_only': MLA_SAMPLES,
    'all_with_data': {**ALL_ENDO_MLA, **OVARIAN_WITH_STOMICS},
    'all_samples': ALL_SAMPLES,
}


def get_target_samples(sample_set: str = None, sample_ids: list = None) -> dict:
    """Resolve which samples to target."""
    if sample_ids:
        return {s: ALL_SAMPLES[s] for s in sample_ids if s in ALL_SAMPLES}

    if sample_set and sample_set in SAMPLE_SETS:
        return SAMPLE_SETS[sample_set]

    return ALL_ENDO_MLA


def check_all_data(target_samples: dict) -> pd.DataFrame:
    """Check data availability for all target samples."""
    rows = []
    for sid, info in target_samples.items():
        status = check_sample_data(sid, info['chip_id'], info.get('aligned', False))
        status['disease'] = info.get('disease', 'Unknown')
        status['can_register'] = (
            status['geojson'] and status['comet_bs'] and status['stomics_dapi']
        )
        status['can_integrate'] = (
            status['can_register'] and
            (status['stomics_cellbin'] or status['stomics_h5ad'])
        )
        rows.append(status)

    return pd.DataFrame(rows)


def build_run_list(status_df: pd.DataFrame, output_dir: Path,
                   skip_existing: bool = True) -> list:
    """Build list of samples ready to process."""
    samples_to_run = []

    for _, row in status_df.iterrows():
        sid = row['sample_id']
        chip = row['chip_id']

        if not row['can_register']:
            logging.info(f'SKIP {sid}: missing required data '
                         f'(geojson={row["geojson"]}, '
                         f'comet_bs={row["comet_bs"]}, '
                         f'stomics_dapi={row["stomics_dapi"]})')
            continue

        if skip_existing:
            metrics_file = output_dir / f'{sid}_{chip}' / f'{sid}_{chip}_alignment_metrics.json'
            if metrics_file.exists():
                logging.info(f'SKIP {sid}: already processed (use --force to re-run)')
                continue

        samples_to_run.append({
            'sample_id': sid,
            'chip_id': chip,
            'disease': row['disease'],
            'has_stomics': row['stomics_cellbin'] or row['stomics_h5ad'],
        })

    return samples_to_run


# ============================================================
# Batch Runner
# ============================================================

def run_batch(samples_to_run: list, output_dir: Path,
              aggregation_method: str = 'nearest',
              max_distance: int = 50,
              save_registered_slides: bool = False) -> pd.DataFrame:
    """Run the alignment pipeline on all samples in the list."""
    batch_results = []
    start_time = datetime.now()

    total = len(samples_to_run)
    for i, sample in enumerate(samples_to_run):
        sid = sample['sample_id']
        chip = sample['chip_id']

        logging.info(f'\n{"#" * 70}')
        logging.info(f'Processing {i+1}/{total}: {sid} ({sample["disease"]}, {chip})')
        logging.info(f'{"#" * 70}')

        sample_start = datetime.now()

        try:
            result = run_alignment_pipeline(
                sample_id=sid,
                chip_id=chip,
                work_dir=output_dir,
                aggregation_method=aggregation_method,
                max_distance=max_distance,
                save_registered_slides=save_registered_slides,
            )

            elapsed = (datetime.now() - sample_start).total_seconds() / 60

            summary = {
                'sample_id': sid,
                'chip_id': chip,
                'disease': sample['disease'],
                'status': result.get('status', 'unknown'),
                'n_comet_cells': result.get('n_comet_cells', 0),
                'elapsed_min': round(elapsed, 1),
            }

            if 'alignment_metrics' in result:
                m = result['alignment_metrics']
                summary.update({
                    'quality': m['quality'],
                    'median_distance': round(m['median_distance'], 1),
                    'pct_within_50px': round(m['pct_within_50px'], 1),
                    'n_stomics_cells': m.get('n_stomics_cells', 0),
                })

            if 'integrated_adata' in result:
                summary['n_integrated_genes'] = result['integrated_adata'].n_vars

            batch_results.append(summary)
            logging.info(f'  -> {result.get("status", "unknown")} ({elapsed:.1f} min)')

        except Exception as e:
            elapsed = (datetime.now() - sample_start).total_seconds() / 60
            logging.error(f'ERROR processing {sid}: {e}')
            traceback.print_exc()
            batch_results.append({
                'sample_id': sid,
                'chip_id': chip,
                'disease': sample['disease'],
                'status': f'error: {str(e)[:80]}',
                'elapsed_min': round(elapsed, 1),
            })

    total_elapsed = (datetime.now() - start_time).total_seconds() / 60
    logging.info(f'\nBatch complete in {total_elapsed:.1f} min')

    return pd.DataFrame(batch_results)


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser(
        description='Batch COMET-STOmics alignment pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_batch_alignment.py --sample-set endo_mla
  python run_batch_alignment.py --sample-set all_with_data --force
  python run_batch_alignment.py --samples SO24 SO25 SO30 SO34
  python run_batch_alignment.py --dry-run --sample-set all_samples
        """
    )

    parser.add_argument('--sample-set', choices=list(SAMPLE_SETS.keys()),
                        default='endo_mla',
                        help='Predefined sample set (default: endo_mla)')
    parser.add_argument('--samples', nargs='+',
                        help='Specific sample IDs (overrides --sample-set)')
    parser.add_argument('--output-dir',
                        default='T:/Sammy Data/projects/out/comet_stomics_alignment',
                        help='Output directory')
    parser.add_argument('--method', choices=['nearest', 'polygon'],
                        default='nearest',
                        help='Aggregation method (default: nearest)')
    parser.add_argument('--max-distance', type=int, default=50,
                        help='Max distance for cell matching in pixels (default: 50)')
    parser.add_argument('--force', action='store_true',
                        help='Re-process samples that already have results')
    parser.add_argument('--save-slides', action='store_true',
                        help='Save warped slide images (large files)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Show what would be processed without running')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Enable debug logging')

    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    log_file = setup_logging(output_dir, args.verbose)

    logging.info('=' * 70)
    logging.info('COMET-STOmics Batch Alignment Pipeline')
    logging.info(f'Started: {datetime.now():%Y-%m-%d %H:%M:%S}')
    logging.info(f'Output: {output_dir}')
    logging.info(f'Log: {log_file}')
    logging.info('=' * 70)

    # Resolve target samples
    target = get_target_samples(args.sample_set, args.samples)
    logging.info(f'\nTarget: {len(target)} samples '
                 f'({"custom" if args.samples else args.sample_set})')

    # Check data availability
    logging.info('\nChecking data availability...')
    status_df = check_all_data(target)

    can_integrate = status_df['can_integrate'].sum()
    can_register = status_df['can_register'].sum()
    cannot = (~status_df['can_register']).sum()

    logging.info(f'  Can run full pipeline: {can_integrate}')
    logging.info(f'  Can register only:     {can_register - can_integrate}')
    logging.info(f'  Cannot process:        {cannot}')

    # Build run list
    samples_to_run = build_run_list(
        status_df, output_dir, skip_existing=not args.force
    )

    logging.info(f'\nSamples to process: {len(samples_to_run)}')

    if not samples_to_run:
        logging.info('Nothing to process. Use --force to re-run existing results.')
        return

    for s in samples_to_run:
        stomics_str = 'full pipeline' if s['has_stomics'] else 'WARP ONLY'
        logging.info(f'  {s["sample_id"]:6s} ({s["disease"]:8s}, {s["chip_id"]}) - {stomics_str}')

    if args.dry_run:
        logging.info('\n** DRY RUN - no processing performed **')
        return

    # Run batch
    logging.info(f'\nStarting batch processing ({args.method} method, '
                 f'max_distance={args.max_distance}px)...')

    results_df = run_batch(
        samples_to_run, output_dir,
        aggregation_method=args.method,
        max_distance=args.max_distance,
        save_registered_slides=args.save_slides,
    )

    # Save results summary
    summary_path = output_dir / f'batch_summary_{datetime.now():%Y%m%d_%H%M%S}.csv'
    results_df.to_csv(summary_path, index=False)

    # Print final summary
    logging.info('\n' + '=' * 70)
    logging.info('BATCH RESULTS')
    logging.info('=' * 70)

    n_complete = (results_df['status'] == 'complete').sum()
    n_failed = (results_df['status'].str.startswith('error')).sum() if 'status' in results_df.columns else 0
    n_other = len(results_df) - n_complete - n_failed

    logging.info(f'  Complete: {n_complete}')
    logging.info(f'  Failed:   {n_failed}')
    logging.info(f'  Other:    {n_other}')

    if 'elapsed_min' in results_df.columns:
        total_min = results_df['elapsed_min'].sum()
        logging.info(f'  Total processing time: {total_min:.0f} min ({total_min/60:.1f} hours)')

    if 'median_distance' in results_df.columns:
        completed = results_df[results_df['status'] == 'complete']
        if len(completed) > 0:
            logging.info(f'\nAlignment Quality:')
            for _, row in completed.iterrows():
                quality = row.get('quality', '?')
                median = row.get('median_distance', '?')
                pct50 = row.get('pct_within_50px', '?')
                logging.info(f'  {row["sample_id"]:6s} ({row["disease"]:8s}): '
                             f'{quality} (median={median}px, within_50px={pct50}%)')

    if n_failed > 0:
        logging.info(f'\nFailed samples:')
        failed = results_df[results_df['status'].str.startswith('error')]
        for _, row in failed.iterrows():
            logging.info(f'  {row["sample_id"]}: {row["status"]}')

    logging.info(f'\nResults saved: {summary_path}')
    logging.info(f'Log saved: {log_file}')
    logging.info('=' * 70)


if __name__ == '__main__':
    main()
