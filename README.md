# MDACC-STOmics-COMET-MSI

Multi-modal spatial biology analysis pipeline integrating three imaging modalities for endometriosis and MLA tissue characterization at MD Anderson Cancer Center.

## Modalities

- **COMET** — Multi-channel immunofluorescence protein imaging (48-channel, 24 phenotype classes via Visiopharm)
- **STOmics** — Stereo-seq spatial transcriptomics (subcellular resolution gene expression)
- **MSI** — Mass spectrometry imaging (glycans, peptides, metabolites)

## Repository Structure

```
src/                   Library code (pip install -e .)
  comet/               COMET-STOmics alignment, validation, analysis
  msi/                 MSI extraction, normalization, spatial analysis, registration
  stomics/             STOmics data auditing and processing

notebooks/             Pipeline notebooks (ordered by step)
  comet/               00-10: GeoJSON export → alignment → QC → cellbin comparison
  msi/                 01-07: extraction → AnnData → QC → clustering → DE → registration
  stomics/             01-04: preprocessing → QC → clustering → scRNA integration

config/                Pipeline configuration and alignment transforms
docs/                  Action plans, data inventory, reference figures
archive/               Superseded notebook variants (preserved for reference)
```

## Setup

```bash
# Clone
git clone https://github.com/AstraeaBio/MDACC-STOmics-COMET-MSI.git
cd MDACC-STOmics-COMET-MSI

# Install as editable package
pip install -e ".[dev,registration]"

# Set up automatic notebook output stripping
nbstripout --install
```

## Pipeline Overview

### COMET-STOmics Alignment
1. Export Visiopharm MLD segmentations to GeoJSON (`src/comet/mld_to_geojson.py`)
2. Register COMET DAPI → STOmics DAPI via VALIS with QuPath affine seeds
3. Warp cell segmentations to STOmics coordinate space
4. Aggregate gene expression per COMET-defined cell → integrated AnnData
5. QC review and sign-off (validation gates downstream analysis)

### MSI Analysis
1. Extract intensities from OME-TIFF files
2. Build AnnData per modality (glycans, peptides, metabolites)
3. Normalize, QC, spatial analysis (Moran's I)
4. Clustering (PCA → UMAP → Leiden)
5. Differential analysis across disease groups
6. Register MSI → H&E for multi-modal integration

### Cross-Modal Integration
COMET cells serve as the integration scaffold: each cell gets gene expression (STOmics), protein levels (COMET), and metabolite features (MSI) through chained spatial registration.

## Alignment Validation

The pipeline includes a validation and sign-off system that prevents downstream analysis from running on unapproved alignments. After each alignment run, a `validation_manifest.json` is created with SHA-256 file hashes. The manifest must be reviewed and approved (via `notebooks/comet/06_alignment_qc_signoff.ipynb`) before normalization or differential expression will proceed.

## License

See [LICENSE](LICENSE).
