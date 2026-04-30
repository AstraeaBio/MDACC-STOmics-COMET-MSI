# Data Completeness Report
**Survey Date:** 2026-03-05

## Executive Summary
Comprehensive inventory of all spatial multi-omics data across COMET, MSI, and STOmics modalities. A total of 46 unique samples identified with varying levels of data completeness across modalities.

---

## File Inventory Summary

### COMET Imaging Data
- **Location:** `T:/Sammy Data/BS/`
- **Total files:** 47 COMET BS OME-TIFF images
- **Sample coverage:** 41 unique sample IDs (some have duplicates: SO20, SO22, SO25, SO30, SO47, SO59)
- **File size:** ~48GB per image
- **Status:** Complete coverage for all analyzed samples

### Visiopharm Segmentation (MLD Files)
- **Location:** `T:/Sammy Data/Vis_Analysis_All_Samples/`
- **Total files:** 40 MLD binary files
- **Sample coverage:** 40 unique sample IDs
- **Status:** 6 files marked as "Aligned" (SO20, SO22, SO25, SO30, SO47, SO59)
- **Gap:** SO13, SO23, SO55, SO56, SO58, SO59, SO64, SO67 have reduced COMET coverage

### GeoJSON Cell Segmentation Masks
- **Location:** `T:/Sammy Data/Vis_Analysis_All_Samples/geojson_output/`
- **Total files:** 31 GeoJSON files
- **Sample coverage:** 31 unique sample IDs
- **Status:** INCOMPLETE - missing GeoJSON for 9 samples: SO13, SO23, SO55, SO56, SO58, SO59, SO64, SO67
- **Gap Analysis:** 31/40 = 77.5% coverage
  - Missing: SO23, SO55, SO56, SO58, SO59, SO64, SO67 (7 total)
  - Note: SO13 has MLD but no GeoJSON

---

## H&E Histology Images
- **Total files:** 32 NDPI images
- **Locations:**
  - Top-level `T:/Sammy Data/HE/`: 6 files (S024, S032, S033, S034, S044, S045)
  - `T:/Sammy Data/HE/second set/`: 8 files (SO1, SO19, SO22, SO23, SO30, SO4, SO5, SO6)
  - `T:/Sammy Data/HE/third set/`: 18 files (SO8, SO14, SO15, SO36, SO39, SO40, SO48, SO50, SO51, SO52, SO53, SO54, SO55, SO56, SO57, SO58, SO59, SO64)
- **Naming convention:** Mixed (S0XX in top-level, SO format in subdirectories)
- **Coverage:** 32/41 = 78% of unique samples
- **Missing:** SO2, SO9, SO10, SO11, SO12, SO13, SO26, SO46, SO47, SO67

---

## Mass Spectrometry Imaging (MSI)
- **Location:** `T:/Data/215_Sammy_Data/6_Files_Mar3/`
- **Modalities:** Glycans, Metabolites, Peptides
- **Sample coverage:** 6 ovarian samples with all three modalities

### Glycan OME-TIFF
- **Location:** `T:/Data/215_Sammy_Data/6_Files_Mar3/Glycan OMEtif files/`
- **Total files:** 6
- **Samples:** S024, S032, S033, S034, S044, S045

### Metabolite OME-TIFF
- **Location:** `T:/Data/215_Sammy_Data/6_Files_Mar3/Metabolite OMEtif files/`
- **Total files:** 6
- **Samples:** S024, S032, S033, S034, S044, S045

### Peptide OME-TIFF
- **Location:** `T:/Data/215_Sammy_Data/6_Files_Mar3/Peptide OMEtif files/`
- **Total files:** 6
- **Samples:** S024, S032, S033, S034, S044, S045

**Status:** MSI data complete ONLY for 6 ovarian samples. No MSI data for endometrial samples (SO1-SO15, SO19, SO20-SO26, SO30, SO36, SO39, SO40, SO46-SO48, SO50-SO59, SO64, SO67).

---

## STOmics Cellbin Data
- **Total unique chip IDs:** 60 (from 3 locations)
- **Total chips with cellbin.gef:** 35 unique chip IDs

### Distribution by Location
1. **Primary:** `T:/Sammy Data/` (12 chips)
   - A03979E2, A03979G5, C03027C4, C03027F5, C03030F4, C03036D6

2. **Alt Location 1:** `T:/Data/215_Sammy_Data/6_Files_Mar3/` (6 chips)
   - A03979E2 (duplicate), A04100A6, B04101A3, C04139D3, D03451C5, D04165G2

3. **Alt Location 2:** `T:/Data/215_Sammy_Data/Sammy Data/` (42 chips)
   - B04101E6, B04106C3, C03030F4 (duplicate), C03036D6 (duplicate), C03137D4, C03137E3, C03137E6, C03449E3, C03450G3, C04138A5, C04138E4, C04138F3, C04139G2, C04140A3, C04140G4, C04141G2, C04142G4, C04143C6, C04143D3, C04143F3, C04143G2, D03453A6, D03453C2, D04160C4, D04164A1, D04165F6

**Status:** Each cellbin directory contains both `.cellbin.gef` and `.adjusted.cellbin.gef` files in `03.ssDNA_analysis/` subdirectory. Chips duplicated across locations for redundancy.

---

## Data Completeness by Sample Type

### Ovarian Samples (Complete Multimodal)
- **SO24, SO32, SO33, SO34, SO44, SO45** (6 samples)
- **Available:** COMET BS + Visiopharm MLD + GeoJSON + H&E + MSI (all 3 modalities)
- **Status:** Fully complete

### Endometrial Samples - High Coverage
- **SO1, SO4, SO5, SO6, SO8, SO14, SO15, SO22, SO30, SO36, SO39, SO40, SO48, SO50, SO51, SO52, SO53, SO54** (18 samples)
- **Available:** COMET BS + Visiopharm MLD + GeoJSON + H&E
- **Status:** Excellent coverage

### Endometrial Samples - Partial Coverage
- **SO2, SO9, SO10, SO11, SO12, SO26, SO46, SO47** (8 samples)
- **Available:** COMET BS + Visiopharm MLD + GeoJSON (no H&E)
- **Status:** Good COMET/Visiopharm coverage

### Endometrial Samples - MLD/GeoJSON Only
- **SO25** (1 sample)
- **Available:** COMET BS + Visiopharm MLD (Aligned) + GeoJSON (no H&E)
- **Status:** Specialized processing

### Endometrial Samples - Reduced Coverage
- **SO13, SO23, SO55, SO56, SO58, SO59, SO64, SO67** (8 samples)
- **Available:** COMET BS + Visiopharm MLD ± H&E (no GeoJSON)
- **Status:** Missing GeoJSON conversion
- **Note:** SO13 also missing MLD

### Samples with H&E Only
- **SO19, SO57** (2 samples)
- **Status:** Histology only, no COMET/STOmics

---

## Data Gaps and Notes

### Critical Gaps
1. **GeoJSON missing for 9 samples:** SO13, SO23, SO55, SO56, SO58, SO59, SO64, SO67 (22.5% gap)
2. **STOmics cellbin mapping:** No explicit sample-to-chip ID linkage documented. Requires manual matching from project records.
3. **MSI data:** Limited to 6 ovarian samples. No metabolomic/proteomic data for endometrial samples.

### Naming Convention Issues
1. **H&E naming:** Top-level uses `S0XX` format (S024, S032, etc.) while subdirectories use `SO` format (SO1, SO4, etc.)
2. **Duplicate BS images:** Some samples (SO20, SO22, SO25, SO30, SO47, SO59) have multiple versions with suffixes like `(1)`, `(2)`
3. **MLD alignment suffix:** Some MLDs marked `(Aligned)` while BS images use numeric suffixes

### Location Redundancy
- **STOmics chips:** Duplicated across 3 T: drive locations for data redundancy/accessibility
- **Primary location:** `T:/Sammy Data/` has 12 chips
- **Backup locations:** Expanded inventory at `T:/Data/215_Sammy_Data/` (alt 1: 6 chips, alt 2: 42 chips)

---

## Recommendations

1. **GeoJSON Conversion:** Complete GeoJSON exports for remaining 9 samples (SO13, SO23, SO55-SO56, SO58-SO59, SO64, SO67)
2. **Sample-Chip Mapping:** Create formal documentation linking SO sample IDs to STOmics chip IDs
3. **MSI Expansion:** Consider profiling endometrial samples for metabolites/peptides to match ovarian coverage
4. **H&E Standardization:** Consolidate H&E naming convention (S0XX → SO##)
5. **Data Consolidation:** Evaluate consolidating STOmics chips to single primary location to simplify access

---

## File Locations Reference

| Data Type | Primary Location |
|-----------|-----------------|
| COMET BS images | `T:/Sammy Data/BS/` |
| Visiopharm MLD | `T:/Sammy Data/Vis_Analysis_All_Samples/` |
| GeoJSON masks | `T:/Sammy Data/Vis_Analysis_All_Samples/geojson_output/` |
| H&E images | `T:/Sammy Data/HE/` (3 subdirs) |
| MSI Glycans | `T:/Data/215_Sammy_Data/6_Files_Mar3/Glycan OMEtif files/` |
| MSI Metabolites | `T:/Data/215_Sammy_Data/6_Files_Mar3/Metabolite OMEtif files/` |
| MSI Peptides | `T:/Data/215_Sammy_Data/6_Files_Mar3/Peptide OMEtif files/` |
| STOmics Primary | `T:/Sammy Data/{ChipID}/03.ssDNA_analysis/` |
| STOmics Alt 1 | `T:/Data/215_Sammy_Data/6_Files_Mar3/{ChipID}/03.ssDNA_analysis/` |
| STOmics Alt 2 | `T:/Data/215_Sammy_Data/Sammy Data/{ChipID}/03.ssDNA_analysis/` |

---

## Generated Inventory Files

Three master CSV files have been created for reference:
1. **DATA_INVENTORY_MASTER.csv** - Sample-level cross-reference of all data types
2. **STOMICS_CHIPS_INVENTORY.csv** - Complete STOmics chip ID to location mapping
3. **SAMPLE_CHIP_MAPPING.csv** - Detailed sample metadata and availability flags
