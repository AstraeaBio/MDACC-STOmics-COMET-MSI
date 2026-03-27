---
title: "feat: COMET vs STOmics Cellbin Segmentation Comparison"
type: feat
status: completed
date: 2026-03-27
---

# feat: COMET vs STOmics Cellbin Segmentation Comparison

## Overview

Build an analysis pipeline to compare COMET cell segmentation masks (from Visiopharm, warped via VALIS) against STOmics DAPI-derived cell segmentation (from cellbin GEF) in the same coordinate space. The goal is to quantify overlap, match cells between the two segmentations, compare expression profiles of matched cells, and produce visualizations to inform which segmentation is more accurate or how to integrate them.

## Problem Frame

We have two independent cell segmentation approaches for the same tissue:
1. **COMET**: Multi-channel immunofluorescence (23+ protein markers) segmented by Visiopharm, warped to STOmics coordinate space via VALIS registration. ~521K cells for SO34.
2. **STOmics**: DAPI-only nuclear staining segmented by the STOmics pipeline, natively in STOmics coordinate space. ~504K cells for SO34 with 32-vertex polygon borders.

It is unknown which segmentation is "better" or whether they capture different cell populations. Comparing them spatially and by expression concordance will inform whether to use one mask, merge them, or use each for different downstream analyses.

## Requirements Trace

- R0. Validate prerequisites: STOmics cellbin GEF exists and is valid (>1 cell); COMET warped GeoJSON exists in STOmics coordinate space; visual alignment inspection possible via QuPath-compatible GeoJSON export
- R1. Load STOmics cellbin borders from GEF and rasterize into a labeled mask
- R2. Rasterize COMET warped GeoJSON into a labeled mask (reuse existing `rasterize_geojson_to_mask`)
- R3. Compute pixel-level overlap metrics (agreement, Jaccard) between the two masks
- R4. Match cells between segmentations via centroid KD-tree (1:1, fragmented, merged)
- R5. Compare expression profiles of matched cell pairs (gene count correlation)
- R6. Compare cell morphology (size distributions, shape)
- R7. Produce visualizations: overlap maps, histograms, scatter plots, concordance, zoomed ROI alignment overlays

## Scope Boundaries

- SO34 (A03979E2) only as pilot — batch extension is not in scope
- No re-registration — uses existing VALIS-warped COMET GeoJSON
- No mask merging or consensus building — this is analysis/comparison only
- No modification to existing pipeline functions (additive only)

## Context & Research

### Relevant Code and Patterns

- `comet/alignment_utils.py:357` — `rasterize_geojson_to_mask()`: rasterizes COMET GeoJSON → uint32 mask, returns `(mask, phenotype_map)`. Reusable for COMET mask.
- `comet/alignment_utils.py:942` — `load_stomics_cellbin_gef()`: loads cellbin into AnnData with spatial coordinates. Does NOT load border polygons.
- `comet/alignment_utils.py:1315` — `map_comet_to_stomics_cells()`: KD-tree 1:1 nearest-neighbor mapping with max_distance threshold. Returns mapping DataFrame.
- `comet/alignment_utils.py:1348` — `compute_alignment_metrics()`: distance-based alignment quality metrics.
- `comet/alignment_utils.py:1382` — `plot_alignment_validation()`: existing COMET vs STOmics overlay plots.
- `comet/script03_comet_stomics_alignment.ipynb` — Pilot alignment notebook pattern to follow.
- `comet/script06_comet_cellbin_generation.ipynb` — Enhanced pipeline notebook pattern.

### STOmics Cellbin GEF Structure (A03979E2)

- `/cellBin/cell`: structured array with `(id:u4, x:i4, y:i4, offset:u4, geneCount:u2, expCount:u2, dnbCount:u2, area:u2, cellTypeID:u2, clusterID:u2)` — 504,086 cells
- `/cellBin/cellBorder`: int16 array shape `(504086, 32, 2)` — 32-vertex polygon per cell as offsets from centroid `(x, y)`
- `/cellBin/gene`: gene names and offset/count for expression lookup
- `/cellBin/cellExp` or `/cellBin/geneExp`: per-cell expression data
- Coordinate range: X:[3064, 22176], Y:[3531, 22407] — within STOmics DAPI (23520x23520)
- Cell area median: 348 pixels, geneCount median: 44 genes

### Key Data Paths (SO34)

- STOmics cellbin: `T:/Sammy Data/A03979E2/03.ssDNA_analysis/A03979E2.adjusted.cellbin.gef`
- STOmics DAPI: `T:/Sammy Data/A03979E2/03.ssDNA_analysis/ssDNA_A03979E2_regist.tif` (23520x23520)
- COMET warped GeoJSON: `{output_dir}/SO34_warped_segmentations.geojson` (from VALIS pipeline)
- COMET masks GeoJSON: `T:/Sammy Data/Vis_Analysis_All_Samples/geojson_output/SO34_BS_masks.geojson`

## Key Technical Decisions

- **Border polygon reconstruction**: STOmics cellBorder stores int16 offsets from centroid. Reconstruct absolute coordinates as `centroid + offset` for each of 32 vertices. This is straightforward and well-documented in the GEF spec.
- **Rasterization approach**: Use cv2.fillPoly (same as existing `rasterize_geojson_to_mask`) on reconstructed STOmics polygons. Output uint32 mask with cell IDs.
- **Cell matching strategy**: Extend beyond simple 1:1 nearest-neighbor. Use mutual KD-tree queries with distance thresholds to identify:
  - **1:1 matches**: reciprocal nearest neighbors within threshold
  - **Fragmented**: one COMET cell maps to multiple STOmics cells (COMET over-segmented or STOmics under-segmented)
  - **Merged**: multiple COMET cells map to one STOmics cell (opposite)
  - **Unmatched**: cells with no partner within threshold
- **Expression comparison**: For matched cell pairs, compare STOmics-native gene counts (from cellbin GEF) vs COMET-aggregated gene counts (from `aggregate_transcripts_by_mask`). Pearson correlation per matched pair.
- **Notebook as primary deliverable**: Interactive notebook with inline visualizations, consistent with project pattern.

## Open Questions

### Resolved During Planning

- **Which STOmics cellbin to use?** Use `adjusted.cellbin.gef` (post-alignment), not raw `.cellbin.gef`.
- **How to handle coordinate systems?** Both are already in STOmics pixel space — COMET via VALIS warp, STOmics natively. No additional transforms needed.
- **What about SO1?** D03453A6 cellbin is broken (1 super-cell). SO34 is the correct pilot.

### Deferred to Implementation

- **Optimal distance threshold for matching**: Start with 30px (cell diameter ~20px at STOmics resolution), tune based on observed distribution.
- **Whether to use IoU vs centroid-distance for matching**: Start with centroid KD-tree (fast), add polygon IoU for matched pairs as a quality metric.
- **Visualization aesthetics**: Determine colormaps and subplot layouts during implementation.

## Implementation Units

- [x] **Unit 0: Prerequisite Validation & QuPath Alignment Inspection**

**Goal:** Add a preflight validation function that confirms both data sources are present and usable before running the comparison. Also export a QuPath-compatible GeoJSON of the STOmics cellbin borders so that both COMET warped GeoJSON and STOmics cellbin GeoJSON can be loaded side-by-side in QuPath on the DAPI image for visual alignment inspection.

**Requirements:** R0

**Dependencies:** None

**Files:**
- Modify: `comet/alignment_utils.py`

**Approach:**
- Add `validate_cellbin_comparison_inputs(sample_id, chip_id, output_dir)` that checks:
  1. **STOmics cellbin GEF exists** — file present at `DefaultPaths.stomics_cellbin_path(chip_id)`
  2. **STOmics cellbin is valid** — open with h5py, read `/cellBin/cell` count, reject if n_cells <= 1 (catches the D03453A6 "super-cell" case)
  3. **STOmics cellbin has borders** — verify `/cellBin/cellBorder` dataset exists
  4. **STOmics DAPI exists** — needed for mask dimensions and visualization
  5. **COMET warped GeoJSON exists** — check `{output_dir}/{sample_id}_warped_segmentations.geojson`; if missing, check if source GeoJSON + registration data exist so the user knows what step to run first
  6. Returns a status dict with pass/fail per check, overall go/no-go, and informative messages
- Add `export_cellbin_borders_to_geojson(gef_path, output_path)` that:
  1. Reads `/cellBin/cell` (centroids, IDs, areas) and `/cellBin/cellBorder` (int16 offsets)
  2. Reconstructs absolute polygon vertices per cell
  3. Exports as a QuPath-compatible GeoJSON FeatureCollection (same format as COMET GeoJSON: Polygon geometry, properties with `label`, `area_px`, optionally `geneCount`)
  4. This enables loading both the COMET warped GeoJSON and the STOmics cellbin GeoJSON as overlays on the STOmics DAPI in QuPath for direct visual comparison of alignment quality
  5. Subsample option (e.g., every Nth cell) for faster QuPath loading on 500K+ cell datasets

**Patterns to follow:**
- `check_sample_data()` at line 249 for validation pattern
- `export_geojson_from_tif()` at line 1628 for GeoJSON export pattern
- COMET warped GeoJSON format for QuPath compatibility (FeatureCollection with Polygon features)

**Test scenarios:**
- SO34 (A03979E2): all checks should pass (valid cellbin, borders present, DAPI exists)
- SO1 (D03453A6): should fail cellbin validity check (1 cell)
- A sample without warped GeoJSON: should fail with informative message about running alignment first
- Exported STOmics GeoJSON loads correctly in QuPath alongside COMET warped GeoJSON

**Verification:**
- Validation function returns clear pass/fail for each prerequisite
- Exported GeoJSON opens in QuPath and polygons overlay correctly on STOmics DAPI
- Visual inspection in QuPath confirms alignment quality before proceeding

- [x] **Unit 1: STOmics Cellbin Border Loading & Rasterization**

**Goal:** Add utility functions to load STOmics cellbin border polygons from GEF and rasterize them into a labeled mask.

**Requirements:** R1

**Dependencies:** Unit 0 (shares border loading logic from `export_cellbin_borders_to_geojson`)

**Files:**
- Modify: `comet/alignment_utils.py`

**Approach:**
- Add `load_stomics_cellbin_borders(gef_path)` function that reads `/cellBin/cell` (centroids, areas, IDs) and `/cellBin/cellBorder` (int16 offsets), reconstructs absolute polygon vertices as `centroid + offset`, returns dict with cell metadata and polygon list.
- Add `rasterize_cellbin_to_mask(cellbin_data, mask_shape, output_path=None)` that takes the output of the border loader and uses cv2.fillPoly to create a uint32 labeled mask, same pattern as `rasterize_geojson_to_mask`.
- Handle edge cases: zero-area offset polygons (all 32 vertices at centroid), coordinates outside mask bounds.

**Patterns to follow:**
- `rasterize_geojson_to_mask()` at line 357 for cv2.fillPoly pattern and int32 workaround
- `_load_cellbin_h5py()` at line 1046 for h5py GEF reading pattern

**Test scenarios:**
- Load A03979E2 cellbin: expect 504,086 cells with 32-vertex polygons
- Reconstructed polygons should have coordinates in [3064, 22176] x [3531, 22407] range
- Rasterized mask should have ~504K unique non-zero labels
- Mask shape matches STOmics DAPI (23520, 23520)

**Verification:**
- `load_stomics_cellbin_borders` returns data for all 504K cells
- Rasterized mask pixel coverage is reasonable (not 0%, not 100%)
- No errors on SO34 data

- [x] **Unit 2: Enhanced Cell Matching (1:1, Fragmented, Merged)**

**Goal:** Build a cell matching function that goes beyond simple nearest-neighbor to classify match types between two segmentations.

**Requirements:** R4

**Dependencies:** Unit 1

**Files:**
- Modify: `comet/alignment_utils.py`

**Approach:**
- Add `compare_segmentation_cells(comet_centroids, stomics_centroids, max_distance=30)` that:
  1. Builds KD-trees for both sets of centroids
  2. Queries each COMET centroid's nearest STOmics neighbor AND vice versa
  3. Classifies matches:
     - **1:1**: mutual nearest neighbors within threshold
     - **Fragmented** (COMET→many STOmics): one COMET cell is nearest for multiple STOmics cells
     - **Merged** (many COMET→one STOmics): multiple COMET cells share same nearest STOmics cell
     - **COMET-only**: no STOmics neighbor within threshold
     - **STOmics-only**: no COMET neighbor within threshold
  4. Returns summary dict + detailed DataFrames for each match type

**Patterns to follow:**
- `map_comet_to_stomics_cells()` at line 1315 for KD-tree usage pattern

**Test scenarios:**
- With SO34 data (~521K COMET, ~504K STOmics), expect majority 1:1 matches
- Summary should report counts for all 5 match categories
- Match distances should have a reasonable distribution (peaked near 0 for good alignment)

**Verification:**
- All cells classified into exactly one category
- Sum of all categories equals total cells from both segmentations
- 1:1 matches have reciprocal nearest-neighbor property

- [x] **Unit 3: Pixel-Level Overlap Metrics**

**Goal:** Compute spatial overlap metrics between the two rasterized masks.

**Requirements:** R3

**Dependencies:** Unit 1

**Files:**
- Modify: `comet/alignment_utils.py`

**Approach:**
- Add `compute_mask_overlap_metrics(comet_mask, stomics_mask)` that computes:
  1. **Pixel agreement**: fraction of pixels where both masks agree (both 0 or both >0)
  2. **Foreground overlap**: Jaccard index of foreground regions (any cell vs background)
  3. **Per-matched-cell IoU**: For 1:1 matched cell pairs (from Unit 2), compute IoU of their individual polygon areas
  4. **Coverage stats**: fraction of tissue covered by each mask independently
  5. **COMET-only pixels**: pixels assigned in COMET but not STOmics
  6. **STOmics-only pixels**: pixels assigned in STOmics but not COMET
- Return metrics dict with all computed values

**Patterns to follow:**
- `compute_alignment_metrics()` at line 1348 for metrics dict pattern

**Test scenarios:**
- Pixel agreement should be between 0 and 1
- Per-cell IoU distribution should be informative (not all 0 or all 1)
- Coverage fractions should be reasonable given tissue area

**Verification:**
- All metrics are finite, in valid ranges
- Per-cell IoU computed for a meaningful number of matched pairs

- [x] **Unit 4: Expression Concordance Analysis**

**Goal:** Compare gene expression profiles between matched COMET and STOmics cells.

**Requirements:** R5

**Dependencies:** Unit 2

**Files:**
- Modify: `comet/alignment_utils.py`

**Approach:**
- Add `compare_matched_cell_expression(comet_adata, stomics_adata, match_df)` that:
  1. For each 1:1 matched pair, extracts expression vectors from both AnnDatas
  2. Computes per-cell metrics: total gene count, total transcript count, number of shared genes
  3. Computes per-gene correlation across all matched cells (for highly expressed genes)
  4. Computes aggregate metrics: median Pearson r, fraction of cells with r > 0.5
  5. Returns summary dict + per-pair DataFrame
- The COMET AnnData comes from `aggregate_transcripts_by_mask` (existing pipeline output or recomputed).
- The STOmics AnnData comes from `load_stomics_cellbin_gef`.
- Need to align on common gene set (intersection of var_names).

**Patterns to follow:**
- `aggregate_transcripts_by_mask()` at line 1865 for how COMET AnnData is structured
- `load_stomics_cellbin_gef()` at line 942 for STOmics AnnData structure

**Test scenarios:**
- Common gene set should have thousands of genes (both are from same tissue)
- Per-cell correlations should vary (some high, some low depending on cell type/matching quality)
- Aggregate correlation should be positive for well-matched cells

**Verification:**
- Expression comparison runs without errors on SO34 matched pairs
- Results include per-cell and per-gene correlation metrics

- [x] **Unit 5: Comparison Notebook with Visualizations**

**Goal:** Create the interactive notebook that runs the full comparison pipeline and produces all visualizations.

**Requirements:** R6, R7 (and integrates R1-R5)

**Dependencies:** Units 0-4

**Files:**
- Create: `comet/script07_cellbin_comparison.ipynb`

**Approach:**
- Follow notebook structure from `script03_comet_stomics_alignment.ipynb` and `script06_comet_cellbin_generation.ipynb`
- Notebook sections:
  1. **Setup**: Imports, paths, sample config (SO34 / A03979E2)
  2. **Preflight validation**: Run `validate_cellbin_comparison_inputs()`, display status table, halt if any critical check fails
  3. **QuPath alignment export**: Run `export_cellbin_borders_to_geojson()` to produce STOmics cellbin GeoJSON; print paths for both COMET warped + STOmics cellbin GeoJSON so user can load in QuPath for visual inspection
  4. **Load COMET data**: Load warped GeoJSON, rasterize to mask, extract centroids
  5. **Load STOmics cellbin**: Load borders, rasterize to mask, extract centroids/metadata
  6. **Visual alignment overlay**: Side-by-side mask contour overlays on DAPI — full tissue view + 3 zoomed ROIs (user-selectable crop coordinates) showing both segmentations in different colors. These serve as in-notebook visual alignment quality metrics.
  7. **Cell matching**: Run `compare_segmentation_cells`, display summary table
  8. **Pixel overlap**: Run `compute_mask_overlap_metrics`, display metrics
  9. **Cell morphology**: Size distribution histograms (both segmentations overlaid), scatter of matched cell areas
  10. **Expression concordance**: Run `compare_matched_cell_expression`, correlation plots
  11. **Spatial maps**: Color-coded maps showing match type (1:1/fragmented/merged/unmatched) across tissue
  12. **Summary**: Key findings table, recommendations

**Visualizations to include:**
- Preflight status table (green/red per check)
- Dual-color mask contour overlay on DAPI (full tissue + 3 zoomed ROIs with COMET=cyan, STOmics=magenta)
- Cell size histogram (COMET blue, STOmics orange, overlaid)
- Matched cell area scatter (COMET area vs STOmics area for 1:1 pairs)
- Match distance histogram
- Match type spatial map (cells colored by match category)
- Per-cell gene count correlation scatter
- Per-gene correlation bar chart (top N genes)
- Summary metrics table

**Patterns to follow:**
- `script06_comet_cellbin_generation.ipynb` for notebook cell organization
- `plot_alignment_validation()` at line 1382 for matplotlib subplot pattern

**Test scenarios:**
- Notebook runs end-to-end on SO34 without errors
- All visualization cells produce plots
- Summary metrics are printed clearly

**Verification:**
- All notebook cells execute successfully
- Visualizations are informative and correctly labeled
- Output directory contains saved figures

## System-Wide Impact

- **Interaction graph:** New functions in `alignment_utils.py` are additive — no existing function signatures change.
- **Error propagation:** New functions should handle missing data gracefully (e.g., cellbin GEF without borders dataset).
- **State lifecycle risks:** None — all outputs are fresh computations, no shared state modification.
- **API surface parity:** N/A — notebook-only analysis.

## Risks & Dependencies

- **VALIS alignment quality**: The comparison assumes the COMET→STOmics VALIS registration is reasonably good. If registration is poor, most cells will appear unmatched. **Mitigated by Unit 0**: QuPath GeoJSON export enables visual inspection of alignment quality before running the comparison, and the zoomed ROI overlays in the notebook provide in-notebook visual diagnostics.
- **Broken cellbin data**: Some STOmics cellbin GEFs are broken (e.g., D03453A6 has only 1 "super-cell"). **Mitigated by Unit 0**: preflight validation rejects cellbin with <=1 cell and verifies border data exists.
- **Missing warped GeoJSON**: Comparison requires COMET cells already warped to STOmics space. **Mitigated by Unit 0**: preflight checks for warped GeoJSON and provides actionable guidance if missing.
- **Memory**: Two uint32 masks (23520x23520) = ~4.4 GB total. Manageable on the analysis machine but worth noting.
- **Expression comparison validity**: COMET AnnData aggregates STOmics transcripts using the COMET mask, while STOmics AnnData uses the STOmics mask. They share the same underlying transcript data, so correlations reflect mask agreement more than biological concordance.

## Sources & References

- Related code: `comet/alignment_utils.py`, `comet/script03_comet_stomics_alignment.ipynb`
- STOmics cellbin format: GEF HDF5 spec with cellBorder offsets
- DEVLOG: `DEVLOG_comet_cellbin_generation.md` — documents STOmics cellbin survey results
