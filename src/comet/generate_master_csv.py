#!/usr/bin/env python3
"""Generate master lookup CSV cross-referencing COMET, MSI, and STOmics data."""

import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
from comet.alignment_utils import ALL_SAMPLES, SAMPLES_WITH_CELLBIN, DefaultPaths


def check_file(p):
    """Return filename if path exists, else empty string."""
    return p.name if p.exists() else ''


def find_he(sample_id):
    """Search all H&E subdirectories for a sample."""
    he_base = DefaultPaths.BASE / 'HE'
    # Try multiple naming patterns in multiple locations
    locations = [
        ('HE', he_base),
        ('HE/second set', he_base / 'second set'),
        ('HE/third set', he_base / 'third set'),
    ]
    so_num = sample_id[2:] if sample_id.startswith('SO') else sample_id
    padded = so_num.zfill(3)
    patterns = [
        f'{sample_id}.ndpi',                                    # SO24.ndpi
        f'S0{so_num}.ndpi',                                     # S024.ndpi
        f'S{padded}.ndpi',                                      # S024.ndpi
        f'S00{so_num}.ndpi' if len(so_num) == 1 else None,     # S001.ndpi
    ]
    for loc_name, loc_path in locations:
        for pattern in patterns:
            if pattern and (loc_path / pattern).exists():
                return pattern, loc_name
    return '', ''


def find_msi(sample_id):
    """Check MSI data for all three modalities.

    MSI files are in T:/Sammy Data/Third set results/{modality} OME-TIFF/
    with naming like 'S024 ovary manual glycans.ome.tif' (zero-padded)
    or 'SO8 ovary manual glycans.ome.tif' (SO prefix).
    """
    msi_base = DefaultPaths.BASE / 'Third set results'
    dirs = {
        'glycan': msi_base / 'glycan OME-TIFF',
        'metabolite': msi_base / 'metabolite OME-TIFF',
        'peptide': msi_base / 'peptide OME-TIFF',
    }
    mod_suffix = {'glycan': 'glycans', 'metabolite': 'mets', 'peptide': 'peptides'}

    results = {}
    for mod, mod_dir in dirs.items():
        if not mod_dir.exists():
            results[mod] = ''
            continue

        so_num = sample_id[2:] if sample_id.startswith('SO') else sample_id
        suffix = mod_suffix[mod]

        # Build candidate filenames in order of likelihood
        # Zero-padded 3-digit: S024, S001
        padded = so_num.zfill(3)
        candidates = [
            f'{sample_id} ovary manual {suffix}.ome.tif',       # SO24 ovary...
            f'{sample_id}ovary manual {suffix}.ome.tif',        # SO51ovary... (no space)
            f'S0{so_num} ovary manual {suffix}.ome.tif',        # S024 ovary...
            f'S{padded} ovary manual {suffix}.ome.tif',         # S024 ovary...
            f'S00{so_num} ovary manual {suffix}.ome.tif',       # S001 ovary... (for single digit)
        ]

        found = ''
        for p in candidates:
            if (mod_dir / p).exists():
                found = p
                break

        # Fallback: glob for any file containing the number
        if not found:
            for f in mod_dir.glob(f'*.ome.tif'):
                # Match on the numeric part
                fname_lower = f.name.lower()
                if f'o{so_num} ' in fname_lower or f'o{so_num}o' in fname_lower or f'{padded} ' in fname_lower:
                    found = f.name
                    break

        results[mod] = found
    return results


def find_stomics_location(chip_id):
    """Determine which location holds the STOmics data."""
    primary = DefaultPaths.BASE / chip_id / '03.ssDNA_analysis'
    if primary.exists():
        return f'T:/Sammy Data/{chip_id}/'

    for root in DefaultPaths.STOMICS_ALT_ROOTS:
        alt = root / chip_id / '03.ssDNA_analysis'
        if alt.exists():
            return str(root / chip_id).replace('\\', '/')
    return ''


def main():
    output_path = Path('T:/0_Organizational/Git/MDACC-STOmics-COMET-MSI/DATA_INVENTORY_MASTER.csv')

    # Sort samples by numeric ID
    def sort_key(item):
        num = ''.join(c for c in item[0] if c.isdigit())
        return int(num) if num else 0

    rows = []
    for sample_id, info in sorted(ALL_SAMPLES.items(), key=sort_key):
        chip_id = info['chip_id']
        disease = info['disease']
        aligned = info['aligned']

        # COMET BS
        suffix = '(Aligned)' if aligned else ''
        # Check for single and dual BS images
        bs_single = DefaultPaths.BS_DIR / f'{sample_id}_BS.ome.tiff'
        bs_1 = DefaultPaths.BS_DIR / f'{sample_id}(1)_BS.ome.tiff'
        bs_2 = DefaultPaths.BS_DIR / f'{sample_id}(2)_BS.ome.tiff'
        if bs_single.exists():
            comet_bs = bs_single.name
        elif bs_1.exists():
            parts = []
            if bs_1.exists(): parts.append(bs_1.name)
            if bs_2.exists(): parts.append(bs_2.name)
            comet_bs = '; '.join(parts)
        else:
            comet_bs = ''

        # Visiopharm MLD
        mld = DefaultPaths.mld_path(sample_id, aligned)
        mld_file = check_file(mld)

        # GeoJSON
        geojson = DefaultPaths.geojson_path(sample_id, aligned)
        geojson_file = check_file(geojson)

        # H&E
        he_file, he_loc = find_he(sample_id)

        # MSI
        msi = find_msi(sample_id)

        # STOmics - always check all locations (don't rely on SAMPLES_WITH_CELLBIN)
        stomics_loc = find_stomics_location(chip_id)

        # STOmics file checks
        stomics_dapi = ''
        stomics_cellbin = ''
        stomics_raw_cellbin = ''
        if stomics_loc:
            dapi_path = DefaultPaths.stomics_dapi_path(chip_id)
            cellbin_path = DefaultPaths.stomics_cellbin_path(chip_id)
            raw_cellbin_path = DefaultPaths.stomics_raw_cellbin_path(chip_id)
            stomics_dapi = check_file(dapi_path)
            stomics_cellbin = check_file(cellbin_path)
            stomics_raw_cellbin = check_file(raw_cellbin_path)

        rows.append({
            'sample_id': sample_id,
            'chip_id': chip_id,
            'disease': disease,
            'aligned_bs': aligned,
            'comet_bs': comet_bs,
            'visiopharm_mld': mld_file,
            'geojson': geojson_file,
            'he_image': he_file,
            'he_location': he_loc,
            'msi_glycan': msi['glycan'],
            'msi_metabolite': msi['metabolite'],
            'msi_peptide': msi['peptide'],
            'stomics_dapi': stomics_dapi,
            'stomics_cellbin': stomics_cellbin,
            'stomics_raw_cellbin': stomics_raw_cellbin,
            'stomics_location': stomics_loc,
        })

    # Write CSV
    fieldnames = list(rows[0].keys())
    with open(output_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} samples to {output_path}")

    # Summary stats
    has = lambda field: sum(1 for r in rows if r[field])
    print(f"\nData availability (out of {len(rows)} samples):")
    print(f"  COMET BS:           {has('comet_bs')}")
    print(f"  Visiopharm MLD:     {has('visiopharm_mld')}")
    print(f"  GeoJSON:            {has('geojson')}")
    print(f"  H&E:                {has('he_image')}")
    print(f"  MSI Glycan:         {has('msi_glycan')}")
    print(f"  MSI Metabolite:     {has('msi_metabolite')}")
    print(f"  MSI Peptide:        {has('msi_peptide')}")
    print(f"  STOmics DAPI:       {has('stomics_dapi')}")
    print(f"  STOmics Cellbin:    {has('stomics_cellbin')}")


if __name__ == '__main__':
    main()
