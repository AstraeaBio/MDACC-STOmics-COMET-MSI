# Sammy Scripts

## Project Overview
This repository contains scientific analysis scripts for multimodal spatial biology data, including:
- Mass Spectrometry Imaging (MSI) analysis
- COMET imaging analysis
- Spatial transcriptomics (STOmics) data processing
- Image registration and spatial integration

## Directory Structure

### `/comet/`
COMET imaging analysis scripts
- `script01_comet_analysis.py` - Main COMET analysis pipeline
- `script01_comet_analysis.ipynb` - Interactive notebook version
- `script01_comet_analysis_mem_optimized.ipynb` - Memory-optimized version
- `script02_image_check_and_conversions.ipynb` - Image preprocessing and format conversions

### `/msi/`
Mass Spectrometry Imaging (MSI) intensity extraction
- `script01a_intensity_extraction_msi_peptides.ipynb` - Peptide intensity extraction
- `script01a_intensity_extraction_msi_metab.ipynb` - Metabolite intensity extraction
- `script01a_intensity_extraction_msi_glycans.ipynb` - Glycan intensity extraction

### `/stomics/`
Spatial transcriptomics (STOmics) analysis pipeline
- `script01_pre-processing.ipynb` - Initial data preprocessing
- `script01_pre-processing_new.ipynb` - Updated preprocessing pipeline
- `script02_quality_control.ipynb` - Quality control metrics and filtering
- `script03_spatial_clustering.ipynb` - Spatial clustering analysis
- `script04_sc_RNA_seq.ipynb` - Single-cell RNA-seq analysis
- `script05_msi_analysis.ipynb` - MSI integration

### `/testing/`
Experimental and testing scripts
- `script01_image_registration.ipynb` - Image registration algorithms
- `script01_image_registration_stomics_comet.ipynb` - STOmics-COMET registration
- `script01_coordinate_plots_test_for_matches.ipynb` - Coordinate matching validation

## Technologies & Tools
- Python (Jupyter notebooks and scripts)
- Likely dependencies: numpy, pandas, scipy, scikit-image, scanpy, matplotlib
- Spatial omics data processing
- Image analysis and registration
- Single-cell analysis

## Analysis Workflow
1. **Data Preprocessing** - Quality control and normalization
2. **Spatial Analysis** - Clustering and spatial statistics
3. **Single-cell Analysis** - Cell type identification and expression patterns
4. **MSI Integration** - Multi-modal data integration
5. **Image Registration** - Aligning different imaging modalities

## Expected Input Folder Structure

### Base Directory Structure
```
~/ext_hd_sammy/
├── data/
│   ├── COMET/
│   │   ├── original/              # Original COMET IF images (.tiff)
│   │   └── visio_output/          # Visio segmentation masks
│   │       └── {analysis_folder}/
│   │           └── {analysis_folder}/  # .tif segmentation masks
│   ├── stomics/
│   │   └── gene_exp/
│   │       └── {Chip_ID}/
│   │           └── 03.ssDNA_analysis/
│   │               └── ssDNA_{Chip_ID}_regist.tif  # STOmics DAPI images
│   └── msi/
│       ├── Glycan OMEtif files/   # Glycan MSI data (.ome.tif)
│       └── he/                    # H&E stained slides (.ndpi)
└── projects/
    └── out/
        ├── out_comet/             # COMET analysis outputs
        ├── out_stomics/           # STOmics analysis outputs
        │   ├── script01_output/   # Preprocessing outputs
        │   ├── script02_output/   # QC outputs
        │   └── script03_output/   # Clustering outputs
        └── matched_images_stomics_COMET/  # Matched images for registration
            └── {Sample_ID}_{Chip_ID}/
```

### Sample Naming Convention
- **Sample IDs**: SO1, SO4, SO5, SO6, etc.
- **Chip IDs**: D03453A6, C03027C4, C03036D6, etc.
- **Whitelist chips**: C04143D3, B04101E6, C04143G2, B04101A3, A04100A6, C04139D3

### COMET Input Files
- **Original images**: `{data}/COMET/original/{Sample_ID}_BS.ome.tiff`
  - Multi-channel TIFF images with protein markers
  - 23+ channels including: DAPI, CY5, TRITC, CD4, CD31, CD8, CD45RO, CD163, FOXp3, TOX1, CD20, CD86, CD68, CD10, CD133, CD11C, GZMB, Periostin, CD66B, COL1A1, CD56, KERATIN8, CD45, aSMA
- **Segmentation masks**: `{data}/COMET/visio_output/{folder}/{Sample_ID}_BS.tif`
  - Multi-channel TIFF files from Visiopharm analysis
  - **Channel 0** (`v_img[0, ...]`) contains the cell segmentation masks
  - Pixel values represent cell IDs (labels)
  - Binary masks processed with `label(v_img[0, ...] > 0)` to create labeled regions

### STOmics Input Files
- **Gene expression data**: Located in `{Chip_ID}/` folders
- **DAPI images**: `{Chip_ID}/03.ssDNA_analysis/ssDNA_{Chip_ID}_regist.tif`
- Expected format: Cellbin output from STOmics pipeline

### MSI Input Files
- **Glycan images**: `{data}/msi/Glycan OMEtif files/{Sample_ID} ovary manual glycans.ome.tif`
- **H&E images**: `{data}/msi/he/{Sample_ID}.ndpi`

## Intermediate Files and Outputs

### COMET Pipeline Outputs

#### script01_comet_analysis.py
**Output location**: `~/ext_hd_sammy/projects/out/out_comet/`

**Files generated**:
- `{Sample_ID}_BS_batch{start}-{end}_protein.parquet` (temporary batch files)
  - Columns: Protein markers + ['centroid_x', 'centroid_y']
  - Contains mean protein expression per cell for batch of cells
- `{Sample_ID}_BS_protein_combined.parquet` (final output)
  - Combined DataFrame with all cells
  - Rows: Individual cells (labeled by segmentation mask)
  - Columns: Mean intensity per protein channel + spatial coordinates
  - Index: Cell label IDs from segmentation

**Processing parameters**:
- Batch size: 150 cells per batch
- Multiprocessing: 4 processes with maxtasksperchild=1

#### script01_comet_analysis.ipynb (Visiopharm Segmentation Analysis)
**Output location**: `~/ext_hd_sammy/projects/out/out_comet_geojson/` and `~/ext_hd_sammy/projects/out/out_comet_csv_polygons/`

**Visiopharm Segmentation Processing**:

The Visiopharm segmentation TIFF files are processed to extract cell boundaries and properties:

1. **Loading segmentation masks**:
   ```python
   v_img = tiff.imread(visio_path)  # Multi-channel TIFF
   single_cell = label(v_img[0, ...] > 0)  # Extract and label channel 0
   props = regionprops(single_cell)  # Extract cell properties
   ```

2. **Cell properties extracted** (using scikit-image `regionprops`):
   - `.coords` - All pixel coordinates belonging to each cell
   - `.centroid` - Cell centroid as (y, x) tuple
   - `.label` - Unique cell ID number from segmentation
   - `.image` - Binary mask of the individual cell
   - `.bbox` - Bounding box [min_row, min_col, max_row, max_col]
   - `.area` - Cell area in pixels

3. **Polygon/boundary extraction**:
   - Uses `find_contours()` from scikit-image to extract cell boundaries
   - Simplifies polygons with `approximate_polygon(tolerance=1.5)` to reduce vertex count
   - Handles complex shapes with exterior rings and interior holes
   - **Important**: X-coordinates are flipped to match coordinate systems: `x_flipped = (image_width - 1) - x`
   - Ensures proper polygon orientation (counter-clockwise for exterior, clockwise for holes)

4. **Output files generated**:

   **GeoJSON format** (`{Sample_ID}_BS_masks.geojson`):
   - FeatureCollection with one Feature per cell
   - Geometry: Polygon with coordinates as [exterior_ring, hole1, hole2, ...]
   - Properties for each cell:
     - `label`: Cell ID from segmentation
     - `area_px`: Cell area in pixels
     - `bbox`: Bounding box coordinates
   - Compatible with GIS tools and image registration software (e.g., VALIS)

   **CSV format** (`{Sample_ID}_BS_vertices.csv`):
   - Columns: `vertex_x`, `vertex_y`, `cell_id`
   - Each row contains arrays of polygon vertices for one cell
   - Useful for custom spatial analysis workflows

5. **Coordinate system transformations**:
   - Original Visiopharm output uses image pixel coordinates
   - X-axis is flipped: `orig_img = np.flip(orig_img, axis=2)` when overlaying
   - Ensures alignment between segmentation masks and original COMET images

6. **Protein intensity extraction per cell**:
   - For each labeled cell: `mask = single_cell == label_id`
   - Extracts mean intensity across all 24 protein channels using the cell mask
   - Cell centroids stored as (centroid_x, centroid_y) for spatial mapping

### STOmics Pipeline Outputs

#### script01_pre-processing (Cellbin output)
**Output location**: `~/ext_hd_sammy/projects/out/out_stomics/script01_output/{folder_name}/`

**Files generated**:
- `{Chip_ID}_{folder_name}.h5ad` - AnnData objects with raw spatial transcriptomics data
  - `.obs`: Cell metadata with x, y coordinates
  - `.var`: Gene information
  - `.X`: Gene expression matrix
  - `.obsm['spatial']`: Spatial coordinates array

#### script02_quality_control.ipynb
**Output location**: `~/ext_hd_sammy/projects/out/out_stomics/script02_output/{folder_name}/`

**Files generated**:
- `{Chip_ID}_{folder_name}.h5ad` - Filtered AnnData objects
  - `.obs['ncounts']`: Total transcript counts per cell
  - `.obs['ngenes']`: Number of genes detected per cell
  - `.obs['sample_id']`: Chip ID
  - `.raw`: Unfiltered data backup
  - Filtered by: min_genes (1st percentile), max_counts (99th percentile), min_counts (1st percentile)

- `metadata_transcripts.csv` - QC statistics
  - Sample ID
  - Median ncounts, median ngenes
  - 99th percentile ncounts/ngenes
  - 1st percentile ncounts/ngenes

- `distribution_ncounts_ngenes.png` - QC visualization plots

#### script03_spatial_clustering.ipynb
**Output location**: `~/ext_hd_sammy/projects/out/out_stomics/script03_output/`

**Files generated**:
- `banksy_results_{Sample_ID}.pkl` - BANKSY dictionary with spatial features
- `banksy_results_{Sample_ID}_results_df.pkl` - Clustering results DataFrame
- `banksy_results_{Sample_ID}_max_num_labels.pkl` - Maximum cluster labels
- `banksy_results_{Sample_ID}_p_names.pkl` - Parameter names
- `banksy_results_{Sample_ID}_adata.h5ad` - AnnData with clustering results
  - `.obs['labels_{method}_pc{dims}_nc{lambda}_r{resolution}']`: Cluster assignments
  - `.obsm['pc{dims}_banksy']`: BANKSY-derived PCA embeddings
  - Highly variable genes subset for analysis

**BANKSY parameters**:
- k_geom: 30 (nearest neighbors)
- max_m: 1 (azimuthal transform order)
- nbr_weight_decay: "scaled_gaussian"
- pca_dims: [21]
- lambda_list: [0.8] (spatial weight)
- resolutions: [0.15] (Leiden clustering)
- num_nn: 25 (for clustering)

### Image Registration Outputs

#### script02_image_check_and_conversions.ipynb
**Output location**: `~/ext_hd_sammy/projects/out/matched_images_stomics_COMET/{Sample_ID}_{Chip_ID}/`

**Files organized**:
- `{Sample_ID}_BS.ome.tiff` - COMET IF image
- `ssDNA_{Chip_ID}_regist.tif` - STOmics DAPI image
- `{Sample_ID}_BS.tif` - COMET segmentation mask
- `{Sample_ID}_glycans.ome.tif` - Single-channel glycan MSI data

#### script01_image_registration_stomics_comet.ipynb
**Output location**: `~/ext_hd_sammy/projects/out/registration_dst_stomics_COMET/{Sample_ID}_{Chip_ID}/`

**Files generated**:
- `registered_slides/` - Warped and registered images
  - `{Sample_ID}_BS.ome.tiff` - Registered COMET image
  - `ssDNA_{Chip_ID}_regist.ome.tiff` - Registered STOmics DAPI
- `processed/` - Intermediate processed images (PNG format)
- Registration metadata and transformation matrices from VALIS

**Registration parameters**:
- align_to_reference: True
- reference_img_f: `ssDNA_{Chip_ID}_regist.tif` (STOmics DAPI as reference)
- Both rigid and non-rigid registration performed

## Notes for AI Assistant
- This is a scientific analysis repository focused on spatial multi-omics
- Scripts are numbered sequentially to indicate workflow order
- Multiple versions of scripts may exist (optimized, new, etc.)
- HTML exports likely represent executed notebook outputs
- Be mindful of memory constraints when working with large imaging datasets
- Consider computational efficiency for memory-intensive operations
- File paths use both Linux (`/home/shamini_pathomics_io/`) and home directory (`~/`) conventions
- Sample-Chip ID matching is critical - maintained in metadata CSV file
- Processing is typically done on whitelisted chips or matched sample sets only
