#!/usr/bin/env python3
"""
saw_data_audit.py - Audit STOmics chips for SAW pipeline readiness
=================================================================

Scans T: drive locations to determine which chips have:
- Existing cellbin GEF outputs (already processed)
- Raw FASTQ files (needed for SAW)
- Mask files (needed for SAW)
- Image files / IPR records (needed for CellBin segmentation)
- GEM expression files (needed for standalone CellBin)

Run on machine with access to T: drive.

Usage:
    python saw_data_audit.py
    python saw_data_audit.py --output saw_audit_results.csv
"""

import argparse
import csv
import os
import sys
from datetime import datetime
from pathlib import Path


# ============================================================
# Complete sample-to-chip mapping from alignment_utils.py
# ============================================================
ALL_CHIPS = {
    # Ovarian
    'SO1':  'D03453A6',
    'SO4':  'C03027C4',
    'SO5':  'C03027F5',
    'SO6':  'C03036D6',
    'SO8':  'C03030F4',
    'SO11': 'C03137E3',
    'SO14': 'C03137D4',
    'SO15': 'C03137E6',
    'SO20': 'D04164A1',
    'SO22': 'A03979G5',
    'SO26': 'C03449E3',
    'SO36': 'C04139G6',
    'SO39': 'C04140E2',
    'SO40': 'C04143D3',
    'SO67': 'B04106E6',
    # Endo
    'SO24': 'D03451C5',
    'SO30': 'D03453C2',
    'SO32': 'D04160C4',
    'SO44': 'B04101A3',
    'SO46': 'C04139D3',
    'SO58': 'C04139G2',
    'SO59': 'C04143G2',
    # MLA
    'SO25': 'C03450G3',
    'SO33': 'D04165G2',
    'SO34': 'A03979E2',
    'SO45': 'A04100A6',
    'SO47': 'C04138E4',
    'SO48': 'B04106C3',
    'SO50': 'C04138A5',
    'SO51': 'C04140G4',
    # CAH
    'SO52': 'C04142G4',
    'SO53': 'C04140A3',
    'SO54': 'C04143F3',
    'SO55': 'C04141G2',
    'SO56': 'C04143C6',
    # USC
    'SO64': 'B04101E6',
}

DISEASE_MAP = {
    'SO1': 'Ovarian', 'SO4': 'Ovarian', 'SO5': 'Ovarian', 'SO6': 'Ovarian',
    'SO8': 'Ovarian', 'SO11': 'Ovarian', 'SO14': 'Ovarian', 'SO15': 'Ovarian',
    'SO20': 'Ovarian', 'SO22': 'Ovarian', 'SO26': 'Ovarian', 'SO36': 'Ovarian',
    'SO39': 'Ovarian', 'SO40': 'Ovarian', 'SO67': 'Ovarian',
    'SO24': 'Endo', 'SO30': 'Endo', 'SO32': 'Endo', 'SO44': 'Endo',
    'SO46': 'Endo', 'SO58': 'Endo', 'SO59': 'Endo',
    'SO25': 'MLA', 'SO33': 'MLA', 'SO34': 'MLA', 'SO45': 'MLA',
    'SO47': 'MLA', 'SO48': 'MLA', 'SO50': 'MLA', 'SO51': 'MLA',
    'SO52': 'CAH', 'SO53': 'CAH', 'SO54': 'CAH+Endo', 'SO55': 'CAH',
    'SO56': 'CAH',
    'SO64': 'USC',
}

# Known STOmics data locations
STOMICS_LOCATIONS = [
    Path('T:/Sammy Data'),
    Path('T:/Data/215_Sammy_Data/6_Files_Mar3'),
    Path('T:/Data/215_Sammy_Data/Sammy Data'),
]

# Additional places to search for FASTQ / mask files
EXTRA_SEARCH_DIRS = [
    Path('T:/Sammy Data/FASTQ'),
    Path('T:/Data/215_Sammy_Data/FASTQ'),
    Path('T:/Data/FASTQ'),
    Path('T:/Data/215_Sammy_Data/masks'),
]


def audit_chip(chip_id: str) -> dict:
    """Check what data exists for a given chip ID."""
    result = {
        'chip_id': chip_id,
        'has_cellbin_gef': False,
        'has_adjusted_cellbin': False,
        'has_regist_tif': False,
        'has_gem_gz': False,
        'has_fastq': False,
        'has_mask': False,
        'has_ipr': False,
        'has_image_tar': False,
        'has_tile_images': False,
        'cellbin_location': '',
        'fastq_r1': '',
        'fastq_r2': '',
        'mask_location': '',
        'gem_location': '',
        'ipr_location': '',
        'image_tar_location': '',
        'chip_dirs_found': [],
    }

    for base in STOMICS_LOCATIONS:
        chip_dir = base / chip_id
        if not chip_dir.exists():
            continue

        result['chip_dirs_found'].append(str(chip_dir))

        # Check cellbin outputs in 03.ssDNA_analysis
        ssdna_dir = chip_dir / '03.ssDNA_analysis'
        if ssdna_dir.exists():
            try:
                for f in ssdna_dir.iterdir():
                    name = f.name.lower()
                    if name.endswith('.cellbin.gef') and 'adjusted' not in name:
                        result['has_cellbin_gef'] = True
                        result['cellbin_location'] = str(ssdna_dir)
                    if 'adjusted.cellbin.gef' in name:
                        result['has_adjusted_cellbin'] = True
                    if 'regist.tif' in name:
                        result['has_regist_tif'] = True
            except PermissionError:
                pass

        # Check for GEM expression files
        try:
            for f in chip_dir.rglob('*.gem.gz'):
                result['has_gem_gz'] = True
                result['gem_location'] = str(f)
                break
            if not result['has_gem_gz']:
                for f in chip_dir.rglob('*.gem'):
                    result['has_gem_gz'] = True
                    result['gem_location'] = str(f)
                    break
        except PermissionError:
            pass

        # Check for FASTQ files
        try:
            fq_files = list(chip_dir.rglob('*.fq.gz')) + list(chip_dir.rglob('*.fastq.gz'))
            if fq_files:
                result['has_fastq'] = True
                r1_files = [f for f in fq_files if '_1.fq' in f.name or '_R1' in f.name or 'read_1' in f.name]
                r2_files = [f for f in fq_files if '_2.fq' in f.name or '_R2' in f.name or 'read_2' in f.name]
                if r1_files:
                    result['fastq_r1'] = ','.join(str(f) for f in r1_files)
                if r2_files:
                    result['fastq_r2'] = ','.join(str(f) for f in r2_files)
                if not r1_files and not r2_files:
                    result['fastq_r1'] = str(fq_files[0].parent)  # Just record the directory
        except PermissionError:
            pass

        # Check for mask files
        try:
            for f in chip_dir.rglob('*.h5'):
                fname = f.name.lower()
                if 'mask' in fname or 'barcodetopos' in fname or chip_id.lower() in fname:
                    result['has_mask'] = True
                    result['mask_location'] = str(f)
                    break
            if not result['has_mask']:
                for f in chip_dir.rglob('*.bin'):
                    if 'mask' in f.name.lower():
                        result['has_mask'] = True
                        result['mask_location'] = str(f)
                        break
        except PermissionError:
            pass

        # Check for image files
        try:
            for f in chip_dir.rglob('*.ipr'):
                result['has_ipr'] = True
                result['ipr_location'] = str(f)
                break
            for f in chip_dir.rglob('*image*.tar.gz'):
                result['has_image_tar'] = True
                result['image_tar_location'] = str(f)
                break
            # Check for tile images directory
            for d in chip_dir.rglob('*tiles*'):
                if d.is_dir():
                    result['has_tile_images'] = True
                    break
        except PermissionError:
            pass

    # Search additional FASTQ / mask locations
    if not result['has_fastq']:
        for search_dir in EXTRA_SEARCH_DIRS:
            if search_dir.exists():
                try:
                    matches = list(search_dir.glob(f'*{chip_id}*'))
                    if matches:
                        result['has_fastq'] = True
                        result['fastq_r1'] = str(matches[0])
                        break
                except PermissionError:
                    pass

    # Convert chip_dirs_found list to string for CSV
    result['chip_dirs_found'] = '; '.join(result['chip_dirs_found'])

    return result


def classify_readiness(r: dict) -> str:
    """Classify a chip's readiness for processing."""
    if r['has_cellbin_gef']:
        return 'DONE'
    if r['has_fastq'] and r['has_mask']:
        if r['has_ipr'] or r['has_image_tar']:
            return 'READY_FULL'      # Can run full SAW with cellbin
        elif r['has_gem_gz'] and r['has_tile_images']:
            return 'READY_CELLBIN'   # Can run standalone CellBin
        else:
            return 'READY_NO_IMAGE'  # Can run SAW without cellbin, or need images
    if r['has_gem_gz']:
        if r['has_tile_images']:
            return 'CELLBIN_ONLY'    # Can run standalone CellBin
        else:
            return 'NEED_IMAGES'     # Have expression but need images
    return 'NEED_DATA'               # Missing fundamental inputs


def main():
    parser = argparse.ArgumentParser(description='Audit STOmics chips for SAW pipeline readiness')
    parser.add_argument('--output', '-o', default='saw_data_audit.csv',
                        help='Output CSV file path')
    parser.add_argument('--verbose', '-v', action='store_true',
                        help='Print detailed per-chip info')
    args = parser.parse_args()

    print("=" * 80)
    print(f"SAW Pipeline Data Audit Report - {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    print("=" * 80)

    # Check if locations are accessible
    accessible = []
    for loc in STOMICS_LOCATIONS:
        if loc.exists():
            accessible.append(str(loc))
    print(f"\nAccessible data locations: {len(accessible)}/{len(STOMICS_LOCATIONS)}")
    for a in accessible:
        print(f"  OK: {a}")
    for loc in STOMICS_LOCATIONS:
        if str(loc) not in accessible:
            print(f"  MISSING: {loc}")

    if not accessible:
        print("\nERROR: No data locations accessible. Are you on a machine with T: drive?")
        sys.exit(1)

    # Audit all chips
    results = []
    for sample_id in sorted(ALL_CHIPS.keys(), key=lambda x: int(x[2:])):
        chip_id = ALL_CHIPS[sample_id]
        r = audit_chip(chip_id)
        r['sample_id'] = sample_id
        r['disease'] = DISEASE_MAP.get(sample_id, 'Unknown')
        r['readiness'] = classify_readiness(r)
        results.append(r)

    # Summary statistics
    total = len(results)
    by_status = {}
    for r in results:
        by_status.setdefault(r['readiness'], []).append(r)

    print(f"\n{'─' * 60}")
    print("SUMMARY")
    print(f"{'─' * 60}")
    print(f"Total chips: {total}")
    for status in ['DONE', 'READY_FULL', 'READY_CELLBIN', 'READY_NO_IMAGE',
                    'CELLBIN_ONLY', 'NEED_IMAGES', 'NEED_DATA']:
        chips = by_status.get(status, [])
        if chips:
            print(f"  {status:20s}: {len(chips):3d}  ({', '.join(r['sample_id'] for r in chips)})")

    # Detailed listing
    print(f"\n{'─' * 60}")
    print("DETAILED STATUS")
    print(f"{'─' * 60}")
    print(f"{'Sample':8s} {'Chip':12s} {'Disease':10s} {'Status':18s} {'CellBin':8s} {'FASTQ':6s} {'Mask':5s} {'IPR':4s} {'GEM':4s}")
    print(f"{'─'*8} {'─'*12} {'─'*10} {'─'*18} {'─'*8} {'─'*6} {'─'*5} {'─'*4} {'─'*4}")

    for r in results:
        print(f"{r['sample_id']:8s} {r['chip_id']:12s} {r['disease']:10s} {r['readiness']:18s} "
              f"{'YES' if r['has_cellbin_gef'] else 'no':8s} "
              f"{'YES' if r['has_fastq'] else 'no':6s} "
              f"{'YES' if r['has_mask'] else 'no':5s} "
              f"{'YES' if r['has_ipr'] else 'no':4s} "
              f"{'YES' if r['has_gem_gz'] else 'no':4s}")

    # Action items
    needs_processing = [r for r in results if r['readiness'] != 'DONE']
    if needs_processing:
        print(f"\n{'─' * 60}")
        print(f"ACTION ITEMS: {len(needs_processing)} chips need processing")
        print(f"{'─' * 60}")

        ready = [r for r in needs_processing if r['readiness'].startswith('READY')]
        need_data = [r for r in needs_processing if not r['readiness'].startswith('READY')]

        if ready:
            print(f"\nReady to process ({len(ready)} chips):")
            for r in ready:
                print(f"  {r['sample_id']} ({r['chip_id']}) - {r['readiness']}")

        if need_data:
            print(f"\nNeed additional data ({len(need_data)} chips):")
            for r in need_data:
                missing = []
                if not r['has_fastq']: missing.append('FASTQ')
                if not r['has_mask']: missing.append('MASK')
                if not r['has_ipr'] and not r['has_image_tar']: missing.append('IMAGES')
                print(f"  {r['sample_id']} ({r['chip_id']}) - Need: {', '.join(missing)}")

    # Save CSV
    output_path = Path(args.output)
    fieldnames = ['sample_id', 'chip_id', 'disease', 'readiness',
                  'has_cellbin_gef', 'has_adjusted_cellbin', 'has_regist_tif',
                  'has_gem_gz', 'has_fastq', 'has_mask', 'has_ipr',
                  'has_image_tar', 'has_tile_images',
                  'cellbin_location', 'fastq_r1', 'fastq_r2',
                  'mask_location', 'gem_location', 'ipr_location',
                  'image_tar_location', 'chip_dirs_found']

    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\nDetailed audit saved to: {output_path.resolve()}")


if __name__ == '__main__':
    main()
