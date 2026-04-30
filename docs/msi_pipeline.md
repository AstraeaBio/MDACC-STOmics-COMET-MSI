# MSI Analysis Pipeline

Comprehensive Mass Spectrometry Imaging (MSI) analysis pipeline supporting:
- All three modalities: glycans, peptides, metabolites
- Clustering/segmentation, differential analysis, and multi-modal integration
- Optional VALIS registration for alignment with DAPI/H&E
- Dual output formats: AnnData (.h5ad) and parquet/CSV

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                    MSI Analysis Pipeline                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                 │
│  OME-TIFF Files                                                 │
│       │                                                         │
│       ▼                                                         │
│  ┌─────────────────────────┐                                    │
│  │ Script 01               │                                    │
│  │ Intensity Extraction    │ ──► Per-channel parquet files      │
│  │ (CLI tool)              │                                    │
│  └─────────────────────────┘                                    │
│       │                                                         │
│       ▼                                                         │
│  ┌─────────────────────────┐                                    │
│  │ Script 02               │                                    │
│  │ Create AnnData          │ ──► .h5ad files                    │
│  └─────────────────────────┘                                    │
│       │                                                         │
│       ▼                                                         │
│  ┌─────────────────────────┐                                    │
│  │ Script 03               │                                    │
│  │ Preprocessing & QC      │ ──► Normalized data, QC reports    │
│  └─────────────────────────┘                                    │
│       │                                                         │
│       ▼                                                         │
│  ┌─────────────────────────┐                                    │
│  │ Script 04               │                                    │
│  │ Spatial Analysis        │ ──► Moran's I, coherent channels   │
│  └─────────────────────────┘                                    │
│       │                                                         │
│       ▼                                                         │
│  ┌─────────────────────────┐                                    │
│  │ Script 05               │                                    │
│  │ Clustering              │ ──► PCA, UMAP, Leiden clusters     │
│  └─────────────────────────┘                                    │
│       │                                                         │
│       ▼                                                         │
│  ┌─────────────────────────┐                                    │
│  │ Script 06               │                                    │
│  │ Differential Analysis   │ ──► Group comparisons, heatmaps    │
│  └─────────────────────────┘                                    │
│       │                                                         │
│       ▼                                                         │
│  ┌─────────────────────────┐                                    │
│  │ Script 07 (Optional)    │                                    │
│  │ Registration            │ ──► Warped coordinates             │
│  └─────────────────────────┘                                    │
│                                                                 │
└─────────────────────────────────────────────────────────────────┘
```

## Scripts

### Script 01: Intensity Extraction (CLI)

Extracts pixel intensities from multi-channel OME-TIFF files.

```bash
python script01_intensity_extraction.py \
    --input-dir /path/to/tiffs \
    --output-dir /path/to/output \
    --modality glycans \
    --pattern "*.ome.tif" \
    --min-intensity 1.0 \
    --workers-images 4 \
    --workers-channels 8
```

**Output:** Per-channel parquet files with columns `(x, y, intensity)`

### Script 02: Create AnnData Objects

Merges per-channel parquet files into AnnData objects.

**Input:** Parquet files from Script 01
**Output:** `.h5ad` files with:
- `.X`: intensity matrix (pixels × channels)
- `.obs`: x, y, sample_id
- `.obsm['spatial']`: coordinate array
- `.var`: channel metadata

### Script 03: Preprocessing & QC

Quality control, normalization, and filtering.

**Features:**
- Normalization: RobustScaler (default), z-score, log, TIC, or none
- QC metrics: total intensity, channel statistics, outliers
- Filtering: low-quality pixels and poor-signal channels

### Script 04: Spatial Analysis

Computes Moran's I spatial autocorrelation per channel.

**Features:**
- KNN spatial weights (k=6 default)
- Significance testing with permutations
- Modality-specific thresholds:
  - Glycans: > 0.2
  - Metabolites: > 0.05
  - Peptides: > 0.2

### Script 05: Clustering

PCA, UMAP, and Leiden clustering.

**Default Parameters:**
- n_pcs: 21
- n_neighbors: 15
- leiden_resolution: 0.6
- umap_min_dist: 0.1

### Script 06: Differential Analysis

Group comparisons and variance-based analysis.

**Features:**
- Per-group variance calculation
- Peak classification (group-specific vs shared)
- Statistical testing with FDR correction
- Volcano plots and clustermaps

### Script 07: Registration (Optional)

VALIS-based image registration to reference images.

**Features:**
- Composite image creation (PCA, SNR, MIP)
- Registration to DAPI/H&E reference
- Coordinate transformation
- Quality assessment

## Output Structure

```
out/out_msi/
├── {modality}_parquet/           # Raw extracted intensities
│   └── {sample_id}/
│       └── channel_{idx}.parquet
├── {modality}_h5ad/              # AnnData objects
│   └── {sample_id}.h5ad
├── qc_reports/                   # QC metrics and plots
│   └── {sample_id}_qc.csv
├── spatial_stats/                # Moran's I results
│   └── {sample_id}_morans_i.csv
├── clustering_results/           # Cluster assignments
│   ├── {sample_id}_clusters.csv
│   └── {sample_id}_umap.csv
└── differential/                 # Group comparisons
    └── {comparison_name}/
        ├── differential_results.csv
        └── clustermap.png
```

## Configuration

Default parameters are in `../config/default_params.yaml`. Key settings:

| Parameter | Default | Description |
|-----------|---------|-------------|
| min_intensity | 1.0 | Threshold for pixel inclusion |
| normalization | "robust" | Normalization method |
| robust_quantile | (0.05, 0.95) | Quantile range for RobustScaler |
| n_pcs | 21 | PCA components |
| n_neighbors | 15 | KNN neighbors for clustering |
| leiden_resolution | 0.6 | Leiden clustering resolution |
| morans_i_knn | 6 | KNN for spatial weights |

## Dependencies

```
numpy
pandas
scipy
scikit-learn
scanpy
anndata
tifffile
pyarrow
umap-learn
leidenalg
esda            # for Moran's I
libpysal        # for spatial weights
matplotlib
seaborn
pyyaml
valis-wsi       # optional, for registration
```

Install with:
```bash
pip install numpy pandas scipy scikit-learn scanpy anndata tifffile pyarrow umap-learn leidenalg esda libpysal matplotlib seaborn pyyaml
pip install valis-wsi  # optional, requires libvips
```

## Usage Example

```python
# 1. Extract intensities (CLI or Python)
from utils import extraction
results = extraction.batch_extract_images(
    input_dir="data/msi/glycans",
    output_dir="out/glycans_parquet",
    pattern="*.ome.tif",
    min_intensity=1.0,
)

# 2. Create AnnData
from utils import io as msi_io
channels = msi_io.load_parquet_channels("out/glycans_parquet", sample_id="SO1")
merged = msi_io.merge_channels_to_matrix(channels)
adata = msi_io.create_anndata_from_merged(merged, sample_id="SO1", modality="glycans")

# 3. Preprocessing
from utils import normalization as msi_norm
msi_norm.normalize_anndata(adata, method='robust')

# 4. Spatial analysis
from utils import spatial as msi_spatial
morans = msi_spatial.morans_i_per_channel(
    adata.X, adata.obsm['spatial'], k=6
)

# 5. Clustering
from utils import clustering as msi_cluster
msi_cluster.cluster_anndata(adata, n_pcs=21, n_neighbors=15, resolution=0.6)
```

## Utility Modules

The `utils/` directory contains reusable functions:

- **extraction.py**: TIFF reading, intensity extraction
- **normalization.py**: RobustScaler, z-score, QC metrics
- **spatial.py**: Moran's I, spatial weights, autocorrelation
- **clustering.py**: PCA, UMAP, Leiden, marker analysis
- **io.py**: AnnData/parquet I/O utilities

## Notes

- The pipeline preserves existing `script01a_*.ipynb` files for reference
- Memory-intensive operations (spatial analysis) support subsampling
- All notebooks have batch processing sections for multiple samples
- Random seeds are configurable for reproducibility
