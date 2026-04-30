# Data Inventory Index

**Survey Date:** 2026-03-05
**Repository:** MDACC-STOmics-COMET-MSI

This directory now contains comprehensive master lookup files documenting all spatial multi-omics data across the project. Use this index to navigate the inventory files.

## Quick Start

**Just want to find data for a specific sample?**
→ Open: **DATA_INVENTORY_MASTER.csv**

**Want to understand data availability and gaps?**
→ Read: **DATA_COMPLETENESS_REPORT.md**

**Need to locate STOmics chips?**
→ Check: **STOMICS_CHIPS_INVENTORY.csv**

**Want comprehensive file listing?**
→ See: **DETAILED_FILE_LISTING.txt**

---

## Inventory Files (7 total)

### 1. DATA_INVENTORY_MASTER.csv (4.5 KB) ⭐ PRIMARY REFERENCE
**Format:** CSV spreadsheet
**Rows:** 46 samples (SO1-SO67 range)
**Columns:** sample_id | comet_bs_image | visiopharm_mld | geojson_export | he_image | he_location | msi_glycan | msi_metabolite | msi_peptide | stomics_chip_id | stomics_location

**Use this when:** Finding which data types exist for a specific sample

**Key features:**
- One row per sample
- Empty cells = data not available
- Includes both filenames and location hints
- Reference for all sample-level data

**Example usage:**
```
Looking for SO24 data?
Find row SO24, check columns for available files
Result: COMET BS ✓, Visiopharm MLD ✓, GeoJSON ✓, H&E ✓, MSI (all 3) ✓
```

---

### 2. STOMICS_CHIPS_INVENTORY.csv (2.4 KB)
**Format:** CSV spreadsheet
**Rows:** 35 unique STOmics chip IDs
**Columns:** chip_id | primary_location | alt_location_1 | alt_location_2 | notes

**Use this when:** Finding cellbin.gef files for a specific chip

**Key features:**
- Lists all 35 unique chip IDs with cellbin data
- Shows which location(s) have each chip
- Notes duplicates and special cases
- Each chip has 2 files: .cellbin.gef and .adjusted.cellbin.gef

**Example usage:**
```
Looking for A03979E2 cellbin files?
Find chip_id A03979E2, see it's in Primary and Alt Location 1
Navigate to: T:/Sammy Data/A03979E2/03.ssDNA_analysis/
Find: A03979E2.cellbin.gef and A03979E2.adjusted.cellbin.gef
```

---

### 3. SAMPLE_CHIP_MAPPING.csv (3.0 KB)
**Format:** CSV spreadsheet
**Rows:** 46 samples
**Columns:** sample_id | tissue_type | chip_id_primary | chip_id_alt | notes | [6 availability boolean flags]

**Use this when:** Quick overview of sample completeness

**Key features:**
- Boolean flags (TRUE/FALSE) for each data type
- Tissue type (Endometrial vs Ovarian)
- Sample notes (e.g., "Aligned MLD", "Complete multimodal")
- Organized for filtering/sorting

**Example usage:**
```
Find all complete multimodal samples?
Filter where: comet_bs_available=TRUE AND visiopharm_mld_available=TRUE
AND geojson_available=TRUE AND he_available=TRUE AND msi_available=TRUE
Result: SO24, SO32, SO33, SO34, SO44, SO45 (6 samples)
```

---

### 4. DATA_COMPLETENESS_REPORT.md (7.7 KB)
**Format:** Markdown document
**Sections:** 10+ sections including gap analysis and recommendations
**Length:** ~250 lines

**Use this when:** Understanding the big picture, data gaps, and recommendations

**Key sections:**
- Executive Summary
- File Inventory Summary (by type)
- Data Completeness by Sample Type
- Critical Gaps and Analysis
- Recommendations
- File Locations Reference Table

**Example usage:**
```
Want to understand what's missing and why?
Read sections:
- Data Completeness by Sample Type (see coverage breakdown)
- Data Gaps and Notes (understand what's missing)
- Recommendations (actionable next steps)
```

---

### 5. DATA_SURVEY_SUMMARY.txt (8.4 KB)
**Format:** Plain text
**Style:** Structured sections with ASCII formatting
**Content:** Statistics, gaps, and quick reference

**Use this when:** Quick scanning without opening other tools

**Key info:**
- Total file counts per modality
- Detailed breakdown by modality
- Data completeness by sample category
- Critical gaps summary
- File location reference

**Example usage:**
```
Want quick stats without Excel/spreadsheet?
Scan TOTAL INVENTORY section for file counts
Scan CRITICAL DATA GAPS section for blockers
```

---

### 6. DETAILED_FILE_LISTING.txt (11 KB)
**Format:** Plain text with organized sections
**Content:** Complete enumeration of all files
**Layout:** Organized by data type

**Use this when:** Need to see every single filename

**Contents:**
- COMET BS images: 47 files listed
- Visiopharm MLD: 40 files listed
- GeoJSON masks: 31 files listed + list of 9 missing
- H&E images: 32 files (by location)
- MSI data: 18 files (by modality)
- STOmics chips: All 35 chip IDs with locations
- Critical observations section
- File count summary table

**Example usage:**
```
Need to verify a specific file exists?
Find relevant section (e.g., "COMET BS IMAGES")
Scan list for your filename
Or use Ctrl+F to search entire document
```

---

### 7. DATA_INVENTORY_README.md (8.9 KB)
**Format:** Markdown guide
**Purpose:** Navigation and usage instructions
**Content:** How to use all other files

**Includes:**
- File-by-file guide with use cases
- Data location reference with full paths
- Key statistics
- Data gaps summary
- Usage scenarios (5 common questions)
- Important notes on naming conventions

**Example usage:**
```
Confused about file locations?
See "Data Type Locations" section
Find your data type, get full path syntax
```

---

## Statistics at a Glance

| Data Type | Files | Samples | Coverage | Status |
|-----------|-------|---------|----------|--------|
| COMET BS | 47 | 41 | 100% | Complete |
| Visiopharm MLD | 40 | 40 | 97.5% | Complete |
| GeoJSON | 31 | 31 | 77.5% | ⚠ Missing 9 |
| H&E | 32 | 32 | 78% | ⚠ Missing 10 |
| MSI Glycans | 6 | 6 | 14.6% | ⚠ Ovarian only |
| MSI Metabolites | 6 | 6 | 14.6% | ⚠ Ovarian only |
| MSI Peptides | 6 | 6 | 14.6% | ⚠ Ovarian only |
| STOmics Chips | 70* | 35 | ~85% | ⚠ No sample mapping |

*70 = 35 chips × 2 file variants (cellbin + adjusted)

---

## Critical Data Gaps

1. **GeoJSON Conversion (9 missing):** SO13, SO23, SO55, SO56, SO58, SO59, SO64, SO67
2. **Sample-to-Chip Mapping (100% undefined):** No documented SO sample → StOmics chip ID links
3. **MSI Endometrial (35 missing):** Only 6 ovarian samples profiled, no endometrial MSI
4. **H&E Images (10 missing):** SO2, SO9, SO10, SO11, SO12, SO13, SO26, SO46, SO47, SO67

---

## Sample Categories

### Ovarian (Complete Multimodal) - 6 samples
SO24, SO32, SO33, SO34, SO44, SO45
- ✓ All modalities available

### Endometrial (High Coverage) - 18 samples
SO1, SO4, SO5, SO6, SO8, SO14, SO15, SO22, SO30, SO36, SO39, SO40, SO48, SO50, SO51, SO52, SO53, SO54
- ✓ COMET + Visiopharm + GeoJSON + H&E
- ✗ Missing MSI and STOmics mapping

### Endometrial (Medium Coverage) - 8 samples
SO2, SO9, SO10, SO11, SO12, SO26, SO46, SO47
- ✓ COMET + Visiopharm + GeoJSON
- ✗ Missing H&E, MSI, STOmics mapping

### Endometrial (Low Coverage) - 9 samples
SO13, SO23, SO25, SO55, SO56, SO58, SO59, SO64, SO67
- ✓ COMET + Visiopharm
- ✗ Missing GeoJSON, H&E, MSI, STOmics mapping

### Reference Only - 2 samples
SO19, SO57
- ✓ H&E only
- ✗ Missing all other modalities

---

## File Locations (Quick Reference)

```
COMET:              T:/Sammy Data/BS/
Visiopharm MLD:     T:/Sammy Data/Vis_Analysis_All_Samples/
GeoJSON:            T:/Sammy Data/Vis_Analysis_All_Samples/geojson_output/
H&E (3 dirs):       T:/Sammy Data/HE/[, /second set, /third set]/
MSI Glycans:        T:/Data/215_Sammy_Data/6_Files_Mar3/Glycan OMEtif files/
MSI Metabolites:    T:/Data/215_Sammy_Data/6_Files_Mar3/Metabolite OMEtif files/
MSI Peptides:       T:/Data/215_Sammy_Data/6_Files_Mar3/Peptide OMEtif files/
STOmics Primary:    T:/Sammy Data/{ChipID}/03.ssDNA_analysis/
STOmics Alt 1:      T:/Data/215_Sammy_Data/6_Files_Mar3/{ChipID}/03.ssDNA_analysis/
STOmics Alt 2:      T:/Data/215_Sammy_Data/Sammy Data/{ChipID}/03.ssDNA_analysis/
```

---

## Common Use Cases

### "I need all data for sample SO24"
1. Open: **DATA_INVENTORY_MASTER.csv**
2. Find row: SO24
3. Result: All data types available (complete multimodal ovarian sample)

### "Which samples have GeoJSON masks?"
1. Open: **DATA_INVENTORY_MASTER.csv** or **SAMPLE_CHIP_MAPPING.csv**
2. Filter: geojson_available = TRUE
3. Count: 31/41 samples (77.5%)
4. Missing: SO13, SO23, SO55, SO56, SO58, SO59, SO64, SO67

### "Where is the cellbin file for chip C04139D3?"
1. Open: **STOMICS_CHIPS_INVENTORY.csv**
2. Find: C04139D3
3. Check: alt_location_1 = T:/Data/215_Sammy_Data/6_Files_Mar3/
4. Navigate: T:/Data/215_Sammy_Data/6_Files_Mar3/C04139D3/03.ssDNA_analysis/

### "Which samples are complete multimodal?"
1. Open: **SAMPLE_CHIP_MAPPING.csv**
2. Filter: All flags = TRUE
3. Result: SO24, SO32, SO33, SO34, SO44, SO45 (6 ovarian samples)

### "Show me file coverage statistics"
1. Read: **DATA_COMPLETENESS_REPORT.md** section "File Inventory Summary"
2. Or: **DATA_SURVEY_SUMMARY.txt** section "TOTAL INVENTORY"
3. See: Coverage percentages by modality with gaps

---

## Tips for Using These Files

1. **With Excel/Sheets:** CSV files open directly in Excel or Google Sheets for filtering/sorting
2. **With Text Editors:** All TXT/MD files view in any text editor or browser
3. **With Command Line:** Use `grep`, `awk`, or `sort` on CSV files for quick queries
4. **Search Function:** Use Ctrl+F in TXT/MD files or Ctrl+H in spreadsheets for quick lookups
5. **Backup:** These reference files are safe to copy/backup; they reference data, don't contain it

---

## Data Last Updated
**Survey Date:** 2026-03-05

When re-running survey (if new data added):
1. Re-run all find/ls commands from survey protocol
2. Update these reference files
3. Update timestamps
4. Note any changes in critical gaps

---

## Questions or Issues?

- **Data not found where expected?** Check DATA_INVENTORY_README.md for path syntax
- **Want to understand gaps?** Read DATA_COMPLETENESS_REPORT.md recommendations section
- **Need specific file count?** See DETAILED_FILE_LISTING.txt enumeration
- **Confused about sample completeness?** Use SAMPLE_CHIP_MAPPING.csv for quick overview
- **Looking for a specific chip?** Use STOMICS_CHIPS_INVENTORY.csv with chip ID

---

## File Statistics

```
Total inventory files: 7
Total size: ~46 KB
Formats: CSV (3), Markdown (2), Text (2)
Samples referenced: 46 unique (SO1-SO67 range)
Data types documented: 8 major modalities
Files enumerated: 238+ data files
Unique chip IDs: 35
Critical gaps identified: 4 major categories
Recommendations provided: 5 actionable items
```

---

**Generated:** 2026-03-05
**Repository:** T:/0_Organizational/Git/MDACC-STOmics-COMET-MSI/
**Status:** Survey Complete ✓
