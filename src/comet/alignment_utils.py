"""
COMET-STOmics Alignment Pipeline Utilities
===========================================

Pipeline for aligning COMET cell segmentations (from Visiopharm MLD/TIF)
to STOmics spatial transcriptomics data, transferring segmentation masks,
and generating integrated multi-modal datasets.

Workflow:
    1. Load COMET segmentations (GeoJSON from Visiopharm)
    2. Extract COMET DAPI channel for registration
    3. Register COMET DAPI -> STOmics DAPI using VALIS
    4. Warp COMET segmentation polygons to STOmics coordinate space
    5. Load STOmics gene expression (cellbin GEF or h5ad)
    6. Aggregate expression per COMET cell -> AnnData
    7. Map COMET cells <-> STOmics cells
    8. Validate alignment quality
"""

import json
import logging
import shutil
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

logger = logging.getLogger(__name__)


# =============================================================================
# Sample Configuration
# =============================================================================

ENDO_SAMPLES = {
    'SO24': {'chip_id': 'D03451C5', 'aligned': False, 'disease': 'Endo'},
    'SO30': {'chip_id': 'D03453C2', 'aligned': True,  'disease': 'Endo'},
    'SO32': {'chip_id': 'D04160C4', 'aligned': False, 'disease': 'Endo'},
    'SO44': {'chip_id': 'B04101A3', 'aligned': False, 'disease': 'Endo'},
    'SO46': {'chip_id': 'C04139D3', 'aligned': False, 'disease': 'Endo'},
    'SO58': {'chip_id': 'C04139G2', 'aligned': False, 'disease': 'Endo'},
    'SO59': {'chip_id': 'C04143G2', 'aligned': True,  'disease': 'Endo'},
}

MLA_SAMPLES = {
    'SO25': {'chip_id': 'C03450G3', 'aligned': True,  'disease': 'MLA'},
    'SO33': {'chip_id': 'D04165G2', 'aligned': False, 'disease': 'MLA'},
    'SO34': {'chip_id': 'A03979E2', 'aligned': False, 'disease': 'MLA'},
    'SO45': {'chip_id': 'A04100A6', 'aligned': False, 'disease': 'MLA'},
    'SO47': {'chip_id': 'C04138E4', 'aligned': True,  'disease': 'MLA'},
    'SO48': {'chip_id': 'B04106C3', 'aligned': False, 'disease': 'MLA'},
    'SO50': {'chip_id': 'C04138A5', 'aligned': False, 'disease': 'MLA'},
    'SO51': {'chip_id': 'C04140G4', 'aligned': False, 'disease': 'MLA'},
}

ALL_ENDO_MLA = {**ENDO_SAMPLES, **MLA_SAMPLES}

# Ovarian samples with STOmics cellbin data (for pipeline validation)
OVARIAN_WITH_STOMICS = {
    'SO4':  {'chip_id': 'C03027C4', 'aligned': False, 'disease': 'Ovarian'},
    'SO5':  {'chip_id': 'C03027F5', 'aligned': False, 'disease': 'Ovarian'},
    'SO6':  {'chip_id': 'C03036D6', 'aligned': False, 'disease': 'Ovarian'},
    'SO8':  {'chip_id': 'C03030F4', 'aligned': False, 'disease': 'Ovarian'},
    'SO22': {'chip_id': 'A03979G5', 'aligned': True,  'disease': 'Ovarian'},
}

# Full 53-sample mapping (all disease subtypes)
ALL_SAMPLES = {
    # Ovarian
    'SO1':  {'chip_id': 'D03453A6', 'aligned': False, 'disease': 'Ovarian'},
    'SO4':  {'chip_id': 'C03027C4', 'aligned': False, 'disease': 'Ovarian'},
    'SO5':  {'chip_id': 'C03027F5', 'aligned': False, 'disease': 'Ovarian'},
    'SO6':  {'chip_id': 'C03036D6', 'aligned': False, 'disease': 'Ovarian'},
    'SO8':  {'chip_id': 'C03030F4', 'aligned': False, 'disease': 'Ovarian'},
    'SO11': {'chip_id': 'C03137E3', 'aligned': False, 'disease': 'Ovarian'},
    'SO14': {'chip_id': 'C03137D4', 'aligned': False, 'disease': 'Ovarian'},
    'SO15': {'chip_id': 'C03137E6', 'aligned': False, 'disease': 'Ovarian'},
    'SO20': {'chip_id': 'D04164A1', 'aligned': True,  'disease': 'Ovarian'},
    'SO22': {'chip_id': 'A03979G5', 'aligned': True,  'disease': 'Ovarian'},
    'SO26': {'chip_id': 'C03449E3', 'aligned': False, 'disease': 'Ovarian'},
    'SO36': {'chip_id': 'C04139G6', 'aligned': False, 'disease': 'Ovarian'},
    'SO39': {'chip_id': 'C04140E2', 'aligned': False, 'disease': 'Ovarian'},
    'SO40': {'chip_id': 'C04143D3', 'aligned': False, 'disease': 'Ovarian'},
    'SO67': {'chip_id': 'B04106E6', 'aligned': False, 'disease': 'Ovarian'},
    # Endo
    **ENDO_SAMPLES,
    # MLA
    **MLA_SAMPLES,
    # CAH
    'SO52': {'chip_id': 'C04142G4', 'aligned': False, 'disease': 'CAH'},
    'SO53': {'chip_id': 'C04140A3', 'aligned': False, 'disease': 'CAH'},
    'SO54': {'chip_id': 'C04143F3', 'aligned': False, 'disease': 'CAH+Endo'},
    'SO55': {'chip_id': 'C04141G2', 'aligned': False, 'disease': 'CAH'},
    'SO56': {'chip_id': 'C04143C6', 'aligned': False, 'disease': 'CAH'},
    # USC
    'SO64': {'chip_id': 'B04101E6', 'aligned': False, 'disease': 'USC'},
}

# Samples known to have STOmics cellbin GEF on disk (all Endo/MLA + validation Ovarian)
SAMPLES_WITH_CELLBIN = {
    'SO4', 'SO5', 'SO6', 'SO8', 'SO22', 'SO34',  # T:/Sammy Data/{chip}/
    # All 14 remaining Endo/MLA chips at T:/Data/215_Sammy_Data/
    'SO24', 'SO30', 'SO32', 'SO44', 'SO46', 'SO58', 'SO59',  # Endo
    'SO25', 'SO33', 'SO45', 'SO47', 'SO48', 'SO50', 'SO51',  # MLA
}


class DefaultPaths:
    """Default data paths on this machine."""
    BASE = Path('T:/Sammy Data')
    BS_DIR = BASE / 'BS'
    VIS_DIR = BASE / 'Vis_Analysis_All_Samples'
    GEOJSON_DIR = VIS_DIR / 'geojson_output'
    IMAGES_DIR = BASE / 'IMAGES'

    REPO_DIR = Path('T:/0_Organizational/Git/MDACC-STOmics-COMET-MSI')
    OUTPUT_BASE = BASE / 'projects' / 'out'

    # Additional STOmics data locations (chips not under BASE)
    STOMICS_ALT_ROOTS = [
        Path('T:/Data/215_Sammy_Data/6_Files_Mar3'),
        Path('T:/Data/215_Sammy_Data/Sammy Data'),
    ]

    @classmethod
    def comet_bs_path(cls, sample_id):
        """Path to COMET BS OME-TIFF (48GB multi-channel)."""
        return cls.BS_DIR / f'{sample_id}_BS.ome.tiff'

    @classmethod
    def geojson_path(cls, sample_id, aligned=False):
        """Path to pre-exported GeoJSON cell masks."""
        suffix = '(Aligned)' if aligned else ''
        return cls.GEOJSON_DIR / f'{sample_id}{suffix}_BS_masks.geojson'

    @classmethod
    def mld_path(cls, sample_id, aligned=False):
        """Path to Visiopharm MLD file."""
        suffix = '(Aligned)' if aligned else ''
        return cls.VIS_DIR / f'{sample_id}{suffix}_BS.mld'

    @classmethod
    def tif_mask_path(cls, sample_id, aligned=False):
        """Path to Visiopharm TIF segmentation mask."""
        suffix = '(Aligned)' if aligned else ''
        return cls.VIS_DIR / f'{sample_id}{suffix}_BS.tif'

    @classmethod
    def _find_stomics_chip_dir(cls, chip_id):
        """Find the STOmics chip directory across all known locations.

        Searches T:/Sammy Data/{chip_id}/ first, then alternate roots.

        Returns:
            Path to the chip's 03.ssDNA_analysis/ directory, or the
            default (possibly non-existent) path under BASE.
        """
        # Primary location
        primary = cls.BASE / chip_id / '03.ssDNA_analysis'
        if primary.exists():
            return primary

        # Search alternate roots
        for root in cls.STOMICS_ALT_ROOTS:
            alt = root / chip_id / '03.ssDNA_analysis'
            if alt.exists():
                return alt

        # Return primary as default (caller checks .exists())
        return primary

    @classmethod
    def stomics_dapi_path(cls, chip_id):
        """Path to STOmics DAPI image (searches all known locations)."""
        ssdna_dir = cls._find_stomics_chip_dir(chip_id)
        return ssdna_dir / f'ssDNA_{chip_id}_regist.tif'

    @classmethod
    def stomics_cellbin_path(cls, chip_id):
        """Path to STOmics adjusted cellbin GEF (searches all known locations)."""
        ssdna_dir = cls._find_stomics_chip_dir(chip_id)
        return ssdna_dir / f'{chip_id}.adjusted.cellbin.gef'

    @classmethod
    def stomics_raw_cellbin_path(cls, chip_id):
        """Path to STOmics raw cellbin GEF (searches all known locations)."""
        ssdna_dir = cls._find_stomics_chip_dir(chip_id)
        return ssdna_dir / f'{chip_id}.cellbin.gef'

    @classmethod
    def stomics_h5ad_path(cls, chip_id):
        """Path to pre-processed STOmics h5ad (bin200, searches all locations)."""
        ssdna_dir = cls._find_stomics_chip_dir(chip_id)
        return ssdna_dir / f'{chip_id}.bin200_1.0.spatial.cluster.h5ad'

    @classmethod
    def label_geojson_path(cls, sample_id):
        """Path to label GeoJSON (cells with phenotype classifications)."""
        return cls.VIS_DIR / 'geojson' / f'{sample_id}_BS_label.geojson'

    @classmethod
    def roi_geojson_path(cls, sample_id):
        """Path to ROI GeoJSON (tissue subregions)."""
        return cls.VIS_DIR / 'geojson' / f'{sample_id}_BS_roi.geojson'

    @classmethod
    def stomics_microbiome_gef_path(cls, chip_id):
        """Path to microbiome GEF (host + microbial taxa, bin1).

        Searches 02.Microbiome_analysis/ under all known chip directories.
        """
        for root in [cls.BASE] + cls.STOMICS_ALT_ROOTS:
            chip_dir = root / chip_id
            micro_gef = chip_dir / '02.Microbiome_analysis' / f'{chip_id}.host_micro.label.gef'
            if micro_gef.exists():
                return micro_gef
        # Return default path
        return cls.BASE / chip_id / '02.Microbiome_analysis' / f'{chip_id}.host_micro.label.gef'

    @classmethod
    def stomics_tissue_gef_path(cls, chip_id):
        """Path to bin1 tissue GEF (searches StandardWorkflow and ssDNA dirs).

        The tissue.gef contains per-DNB transcript positions and is the input
        for generating cellbin GEFs from custom masks.
        """
        # Search patterns in order of likelihood
        search_dirs = []
        for root in [cls.BASE] + cls.STOMICS_ALT_ROOTS:
            chip_dir = root / chip_id
            search_dirs.extend([
                chip_dir / '01.StandardWorkflow_Result' / 'GeneExpMatrix',
                chip_dir / '03.ssDNA_analysis',
            ])

        for d in search_dirs:
            tissue_gef = d / f'{chip_id}.tissue.gef'
            if tissue_gef.exists():
                return tissue_gef
            raw_gef = d / f'{chip_id}.raw.gef'
            if raw_gef.exists():
                return raw_gef

        # Return default (possibly non-existent) path
        return cls.BASE / chip_id / '01.StandardWorkflow_Result' / 'GeneExpMatrix' / f'{chip_id}.tissue.gef'


def check_sample_data(sample_id, chip_id=None, aligned=None):
    """Check what data files exist for a sample.

    Returns:
        dict: File existence status
    """
    if chip_id is None:
        info = ALL_SAMPLES.get(sample_id, ALL_ENDO_MLA.get(sample_id, OVARIAN_WITH_STOMICS.get(sample_id, {})))
        chip_id = info.get('chip_id')
        if aligned is None:
            aligned = info.get('aligned', False)
    if aligned is None:
        aligned = False

    status = {
        'sample_id': sample_id,
        'chip_id': chip_id,
        'geojson': DefaultPaths.geojson_path(sample_id, aligned).exists(),
        'mld': DefaultPaths.mld_path(sample_id, aligned).exists(),
        'tif_mask': DefaultPaths.tif_mask_path(sample_id, aligned).exists(),
        'comet_bs': DefaultPaths.comet_bs_path(sample_id).exists(),
    }

    if chip_id:
        status['stomics_dapi'] = DefaultPaths.stomics_dapi_path(chip_id).exists()
        status['stomics_cellbin'] = DefaultPaths.stomics_cellbin_path(chip_id).exists()
        status['stomics_h5ad'] = DefaultPaths.stomics_h5ad_path(chip_id).exists()
    else:
        status['stomics_dapi'] = False
        status['stomics_cellbin'] = False
        status['stomics_h5ad'] = False

    return status


# =============================================================================
# GeoJSON Loading & Processing
# =============================================================================

def load_geojson(geojson_path):
    """Load GeoJSON cell segmentation file.

    Args:
        geojson_path: Path to GeoJSON file

    Returns:
        dict: GeoJSON FeatureCollection
    """
    geojson_path = Path(geojson_path)
    logger.info(f"Loading GeoJSON: {geojson_path.name}")

    with open(geojson_path) as f:
        data = json.load(f)

    n_features = len(data.get('features', []))
    logger.info(f"Loaded {n_features} cell features")
    return data


def geojson_centroids(geojson_data):
    """Extract cell centroids from GeoJSON polygon coordinates.

    Computes centroids from actual polygon vertices (not stored properties,
    which may be in a different coordinate system due to x-flipping).

    Returns:
        tuple: (np.ndarray of shape (N, 2), list of cell labels)
    """
    centroids = []
    labels = []

    for feature in geojson_data['features']:
        coords = feature['geometry']['coordinates']
        exterior = np.array(coords[0])
        # Exclude closing point for centroid computation
        cx = np.mean(exterior[:-1, 0])
        cy = np.mean(exterior[:-1, 1])
        centroids.append([cx, cy])
        labels.append(feature['properties'].get('label', None))

    return np.array(centroids), labels


def geojson_to_shapely(geojson_data):
    """Convert GeoJSON features to shapely Polygon objects.

    Returns:
        list of (shapely.Polygon, dict) tuples
    """
    from shapely.geometry import shape

    polygons = []
    for feature in geojson_data['features']:
        try:
            poly = shape(feature['geometry'])
            if not poly.is_valid:
                poly = poly.buffer(0)
            polygons.append((poly, feature['properties']))
        except Exception as e:
            logger.warning(f"Skipping invalid feature: {e}")

    return polygons


# =============================================================================
# Rasterize GeoJSON to Mask
# =============================================================================

def rasterize_geojson_to_mask(geojson_data, mask_shape, output_path=None):
    """Rasterize GeoJSON cell polygons into a labeled uint32 TIFF mask.

    Each polygon is filled with its cell label value. Background is 0.
    Used as input for geftools generateCgef or direct transcript assignment.

    Supports two GeoJSON formats:
        - "masks" format: properties.label (from visiopharm_to_geojson)
        - "label" format: properties.object_index (from mld_to_geojson)
          For label format, uses object_index + 1 as the mask value (0 = background).
          Features with classification.name == "Background" are skipped.

    Also builds and returns a cell_id → phenotype mapping when classification
    data is present in the GeoJSON features.

    Args:
        geojson_data: GeoJSON FeatureCollection dict (warped to STOmics coords)
        mask_shape: (height, width) tuple matching STOmics DAPI dimensions
        output_path: Optional path to save mask as compressed TIFF

    Returns:
        tuple: (mask, phenotype_map)
            mask: np.ndarray uint32 of shape (height, width)
            phenotype_map: dict {cell_id: phenotype_str} or None if no
                           classification data present
    """
    import cv2

    height, width = mask_shape
    mask = np.zeros((height, width), dtype=np.int32)
    features = geojson_data['features']

    logger.info(f"Rasterizing {len(features)} polygons to mask ({height}x{width})")

    n_rasterized = 0
    n_skipped = 0
    n_background = 0
    phenotype_map = {}
    has_classification = False

    for feature in features:
        props = feature['properties']

        # Skip Background features (label format)
        classification = props.get('classification')
        if isinstance(classification, dict):
            class_name = classification.get('name', '')
            if class_name == 'Background':
                n_background += 1
                continue
            has_classification = True
        else:
            class_name = None

        # Determine cell ID: prefer 'label' (masks format), fall back to
        # 'object_index' + 1 (label format, +1 because 0 = background)
        label = props.get('label')
        if label is None:
            obj_idx = props.get('object_index')
            if obj_idx is not None:
                label = obj_idx + 1
            else:
                n_skipped += 1
                continue

        if label == 0:
            n_skipped += 1
            continue

        # Record phenotype mapping
        if class_name is not None:
            phenotype_map[int(label)] = class_name if class_name != 'Cell' else 'Unclassified'

        coords = feature['geometry']['coordinates']
        exterior = np.array(coords[0], dtype=np.float64)

        # Convert to integer pixel coordinates (x, y order for cv2)
        pts = np.round(exterior).astype(np.int32)

        # Clip to mask bounds
        pts[:, 0] = np.clip(pts[:, 0], 0, width - 1)
        pts[:, 1] = np.clip(pts[:, 1], 0, height - 1)

        # cv2.fillPoly works with int32 (not uint32)
        cv2.fillPoly(mask, [pts], int(label))

        # Handle holes (interior rings)
        for hole_coords in coords[1:]:
            hole = np.round(np.array(hole_coords, dtype=np.float64)).astype(np.int32)
            hole[:, 0] = np.clip(hole[:, 0], 0, width - 1)
            hole[:, 1] = np.clip(hole[:, 1], 0, height - 1)
            cv2.fillPoly(mask, [hole], 0)

        n_rasterized += 1

    # Convert to uint32 for output (int32 was needed for cv2.fillPoly)
    mask = mask.astype(np.uint32)

    n_unique = len(np.unique(mask)) - 1  # exclude background
    logger.info(f"Rasterized {n_rasterized} cells, {n_skipped} skipped, "
                f"{n_background} background excluded, "
                f"{n_unique} unique labels in mask")
    if has_classification:
        from collections import Counter
        pheno_counts = Counter(phenotype_map.values())
        logger.info(f"Phenotype distribution ({len(pheno_counts)} classes): "
                    f"{dict(pheno_counts.most_common(5))}...")

    if output_path is not None:
        import tifffile
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            tifffile.imwrite(str(output_path), mask, compression='lzw')
        except KeyError:
            tifffile.imwrite(str(output_path), mask, compression='zlib')
        logger.info(f"Saved mask to {output_path} "
                     f"({output_path.stat().st_size / 1e6:.1f} MB)")

    return mask, (phenotype_map if has_classification else None)


def rasterize_roi_to_mask(roi_geojson_data, mask_shape, output_path=None):
    """Rasterize ROI GeoJSON polygons into a categorical uint8 mask.

    Each pixel is labeled with its ROI class. Overlaps use priority-based
    overwriting: Other (lowest) → TME → Tumor → Vessels (highest).

    Args:
        roi_geojson_data: ROI GeoJSON FeatureCollection dict (warped to
            STOmics coords). Features must have classification.name matching
            one of: "Tumor ROI", "TME ROI", "Vessels ROI", "Other ROI".
        mask_shape: (height, width) tuple matching STOmics DAPI dimensions
        output_path: Optional path to save mask as compressed TIFF

    Returns:
        np.ndarray: uint8 mask of shape (height, width) with values:
            0 = Outside ROI, 1 = Tumor ROI, 2 = TME ROI,
            3 = Vessels ROI, 4 = Other ROI
    """
    import cv2

    height, width = mask_shape
    mask = np.zeros((height, width), dtype=np.uint8)

    roi_name_to_val = {
        'Other ROI': 4,
        'TME ROI': 2,
        'Tumor ROI': 1,
        'Vessels ROI': 3,
    }

    features = roi_geojson_data['features']
    logger.info(f"Rasterizing {len(features)} ROI polygons to mask ({height}x{width})")

    # Group features by class for priority-ordered rasterization
    # Rasterize in order: Other → TME → Tumor → Vessels (last wins)
    raster_order = ['Other ROI', 'TME ROI', 'Tumor ROI', 'Vessels ROI']
    features_by_class = {name: [] for name in raster_order}
    n_skipped = 0

    for feature in features:
        props = feature['properties']
        classification = props.get('classification')
        if isinstance(classification, dict):
            class_name = classification.get('name', '')
        else:
            class_name = props.get('name', '')

        if class_name in features_by_class:
            features_by_class[class_name].append(feature)
        else:
            n_skipped += 1

    for class_name in raster_order:
        class_features = features_by_class[class_name]
        if not class_features:
            continue
        val = roi_name_to_val[class_name]

        for feature in class_features:
            coords = feature['geometry']['coordinates']

            # Handle both Polygon and MultiPolygon
            geom_type = feature['geometry'].get('type', 'Polygon')
            if geom_type == 'MultiPolygon':
                polygon_list = coords
            else:
                polygon_list = [coords]

            for poly_coords in polygon_list:
                exterior = np.array(poly_coords[0], dtype=np.float64)
                pts = np.round(exterior).astype(np.int32)
                pts[:, 0] = np.clip(pts[:, 0], 0, width - 1)
                pts[:, 1] = np.clip(pts[:, 1], 0, height - 1)
                cv2.fillPoly(mask, [pts], int(val))

                # Handle holes
                for hole_coords in poly_coords[1:]:
                    hole = np.round(np.array(hole_coords, dtype=np.float64)).astype(np.int32)
                    hole[:, 0] = np.clip(hole[:, 0], 0, width - 1)
                    hole[:, 1] = np.clip(hole[:, 1], 0, height - 1)
                    cv2.fillPoly(mask, [hole], 0)

        logger.info(f"  {class_name}: {len(class_features)} polygons (val={val})")

    # Summary
    from collections import Counter
    val_counts = Counter(mask.ravel())
    roi_val_to_name = {0: 'Outside', 1: 'Tumor', 2: 'TME', 3: 'Vessels', 4: 'Other'}
    coverage = {roi_val_to_name.get(v, f'val{v}'): c
                for v, c in val_counts.items() if v > 0}
    total_roi_px = sum(c for v, c in val_counts.items() if v > 0)
    logger.info(f"ROI coverage: {total_roi_px / (height * width) * 100:.1f}% of image")
    logger.info(f"ROI pixel distribution: {coverage}")
    if n_skipped:
        logger.warning(f"Skipped {n_skipped} features with unrecognized class")

    if output_path is not None:
        import tifffile
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            tifffile.imwrite(str(output_path), mask, compression='lzw')
        except KeyError:
            tifffile.imwrite(str(output_path), mask, compression='zlib')
        logger.info(f"Saved ROI mask to {output_path} "
                     f"({output_path.stat().st_size / 1e6:.1f} MB)")

    return mask


# =============================================================================
# MLD Loading & Comparison
# =============================================================================

def load_mld_as_geojson(mld_path, image_width=None, image_height=None,
                         pixels_per_mm=4347.8):
    """Load Visiopharm MLD file and convert to GeoJSON via coordinate transform.

    Uses the MLD reader module and coordinate transform module to produce
    GeoJSON in pixel coordinates matching the TIF mask approach.

    Args:
        mld_path: Path to .mld file
        image_width: Image width in pixels (for coordinate transform)
        image_height: Image height in pixels (for coordinate transform)
        pixels_per_mm: MLD mm-to-pixel conversion factor (default: 0.23um/px)

    Returns:
        dict: GeoJSON FeatureCollection
    """
    import importlib.util

    mld_path = Path(mld_path)
    module_dir = Path(__file__).parent

    # Load MLD reader module (has hyphens in name, can't import normally)
    mld_reader_path = module_dir / 'mld-reader-module.py'
    spec = importlib.util.spec_from_file_location('mld_reader', str(mld_reader_path))
    mld_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mld_mod)

    # Load coordinate transform module
    transform_path = module_dir / 'mld_w_coordinate-transform.py'
    spec2 = importlib.util.spec_from_file_location('mld_transform', str(transform_path))
    transform_mod = importlib.util.module_from_spec(spec2)
    # Inject ShapeType into transform module namespace
    transform_mod.ShapeType = mld_mod.ShapeType
    spec2.loader.exec_module(transform_mod)

    # Read MLD file
    reader = mld_mod.MLDReader(str(mld_path))
    data = reader.read()

    logger.info(f"MLD: {data['num_layers']} layers, "
                f"{sum(len(l['objects']) for l in data['layers'])} total objects")

    # Convert to pixel-coordinate GeoJSON
    geojson = transform_mod.export_to_geojson_pix(
        data,
        pixels_per_mm=pixels_per_mm,
        image_width_px=image_width,
        image_height_px=image_height,
        invert_y=True,
        center_origin=True,
    )

    n_features = len(geojson.get('features', []))
    logger.info(f"MLD -> GeoJSON: {n_features} features")
    return geojson


def compare_geojson_sources(geojson_tif, geojson_mld, sample_n=5000):
    """Compare GeoJSON from TIF mask vs MLD extraction.

    Compares polygon count, centroid overlap, and area distributions
    to decide which source is more reliable.

    Args:
        geojson_tif: GeoJSON dict from TIF mask export
        geojson_mld: GeoJSON dict from MLD extraction
        sample_n: Number of cells to use for spatial comparison

    Returns:
        dict: Comparison metrics
    """
    tif_centroids, tif_labels = geojson_centroids(geojson_tif)
    mld_centroids = []
    for f in geojson_mld['features']:
        geom = f['geometry']
        if geom and geom['type'] == 'Polygon':
            coords = np.array(geom['coordinates'][0])
            mld_centroids.append([np.mean(coords[:-1, 0]), np.mean(coords[:-1, 1])])
        elif geom and geom['type'] == 'Point':
            mld_centroids.append(geom['coordinates'])
    mld_centroids = np.array(mld_centroids) if mld_centroids else np.empty((0, 2))

    n_tif = len(geojson_tif['features'])
    n_mld = len(geojson_mld['features'])

    comparison = {
        'n_cells_tif': n_tif,
        'n_cells_mld': n_mld,
        'n_diff': abs(n_tif - n_mld),
        'pct_diff': abs(n_tif - n_mld) / max(n_tif, n_mld, 1) * 100,
    }

    # Coordinate range comparison
    if len(tif_centroids) > 0 and len(mld_centroids) > 0:
        comparison['tif_x_range'] = [float(tif_centroids[:, 0].min()),
                                      float(tif_centroids[:, 0].max())]
        comparison['tif_y_range'] = [float(tif_centroids[:, 1].min()),
                                      float(tif_centroids[:, 1].max())]
        comparison['mld_x_range'] = [float(mld_centroids[:, 0].min()),
                                      float(mld_centroids[:, 0].max())]
        comparison['mld_y_range'] = [float(mld_centroids[:, 1].min()),
                                      float(mld_centroids[:, 1].max())]

        # KD-tree matching: for each TIF centroid, find nearest MLD centroid
        n_query = min(sample_n, len(tif_centroids), len(mld_centroids))
        tif_sample = tif_centroids[:n_query]
        tree = cKDTree(mld_centroids)
        distances, _ = tree.query(tif_sample, k=1)

        comparison['centroid_match_median_dist'] = float(np.median(distances))
        comparison['centroid_match_mean_dist'] = float(np.mean(distances))
        comparison['centroid_match_pct_within_10px'] = float(
            (distances < 10).mean() * 100)
        comparison['centroid_match_pct_within_50px'] = float(
            (distances < 50).mean() * 100)

    # Area comparison (TIF only, MLD may not have area)
    tif_areas = [f['properties'].get('area_px', 0)
                 for f in geojson_tif['features']]
    comparison['tif_area_median'] = float(np.median(tif_areas)) if tif_areas else 0

    # Recommendation
    if comparison['pct_diff'] < 5 and comparison.get('centroid_match_median_dist', 999) < 10:
        comparison['recommendation'] = 'BOTH_EQUIVALENT'
    elif n_tif > n_mld:
        comparison['recommendation'] = 'USE_TIF'
    else:
        comparison['recommendation'] = 'USE_MLD'

    return comparison


# =============================================================================
# COMET DAPI Extraction
# =============================================================================

def extract_comet_dapi(comet_bs_path, output_path, dapi_channel=0):
    """Extract DAPI channel from multi-channel COMET OME-TIFF.

    Reads only the specified channel to avoid loading the full ~48GB image.

    Args:
        comet_bs_path: Path to COMET BS OME-TIFF
        output_path: Path to save single-channel DAPI TIFF
        dapi_channel: Channel index for DAPI (default: 0)

    Returns:
        Path to saved DAPI image
    """
    import tifffile as tiff

    comet_bs_path = Path(comet_bs_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Extracting DAPI (channel {dapi_channel}) from {comet_bs_path.name}")

    with tiff.TiffFile(str(comet_bs_path)) as tf:
        if len(tf.series) > 0:
            dapi = tf.series[0].pages[dapi_channel].asarray()
        else:
            dapi = tf.pages[dapi_channel].asarray()

    logger.info(f"DAPI shape: {dapi.shape}, dtype: {dapi.dtype}")

    tiff.imwrite(str(output_path), dapi)
    logger.info(f"Saved DAPI to {output_path}")

    return output_path


# =============================================================================
# VALIS Registration
# =============================================================================

def setup_registration_inputs(sample_id, chip_id, work_dir,
                               comet_dapi_path=None, stomics_dapi_path=None,
                               comet_bs_path=None, dapi_channel=0):
    """Set up paired images directory for VALIS registration.

    Creates a directory with COMET DAPI and STOmics DAPI images.

    Returns:
        tuple: (src_dir, reference_filename)
    """
    work_dir = Path(work_dir)
    src_dir = work_dir / 'input_images'
    src_dir.mkdir(parents=True, exist_ok=True)

    # Resolve STOmics DAPI
    if stomics_dapi_path is None:
        stomics_dapi_path = DefaultPaths.stomics_dapi_path(chip_id)
    stomics_dapi_path = Path(stomics_dapi_path)

    if not stomics_dapi_path.exists():
        raise FileNotFoundError(f"STOmics DAPI not found: {stomics_dapi_path}")

    # Copy STOmics DAPI to src directory
    stomics_dst = src_dir / stomics_dapi_path.name
    if not stomics_dst.exists():
        shutil.copy2(stomics_dapi_path, stomics_dst)
        logger.info(f"Copied STOmics DAPI to {stomics_dst}")

    # Handle COMET DAPI
    if comet_dapi_path is not None:
        comet_dapi_path = Path(comet_dapi_path)
        comet_dst = src_dir / comet_dapi_path.name
        if not comet_dst.exists():
            shutil.copy2(comet_dapi_path, comet_dst)
    else:
        # Extract DAPI from BS image
        if comet_bs_path is None:
            comet_bs_path = DefaultPaths.comet_bs_path(sample_id)

        comet_dapi_out = src_dir / f'{sample_id}_BS_DAPI.tif'
        if not comet_dapi_out.exists():
            extract_comet_dapi(comet_bs_path, comet_dapi_out, dapi_channel)

    reference_filename = stomics_dapi_path.name
    logger.info(f"Registration inputs ready in {src_dir}")
    logger.info(f"Reference image: {reference_filename}")

    return src_dir, reference_filename


def run_valis_registration(src_dir, dst_dir, reference_img_f):
    """Run VALIS image registration (rigid + non-rigid).

    Registers COMET DAPI to STOmics DAPI, preserving STOmics coordinates.

    Args:
        src_dir: Directory containing paired DAPI images
        dst_dir: Output directory for registration results
        reference_img_f: Filename of the reference image (STOmics DAPI)

    Returns:
        registration.Valis: The registrar object with computed transforms
    """
    from valis import registration

    dst_dir = Path(dst_dir)
    dst_dir.mkdir(parents=True, exist_ok=True)

    logger.info("Starting VALIS registration...")
    logger.info(f"  Source: {src_dir}")
    logger.info(f"  Reference: {reference_img_f}")

    registrar = registration.Valis(
        str(src_dir),
        str(dst_dir),
        align_to_reference=True,
        reference_img_f=reference_img_f
    )

    rigid_registrar, non_rigid_registrar, error_df = registrar.register()

    logger.info("Registration complete")
    if error_df is not None:
        logger.info(f"Registration errors:\n{error_df}")

    return registrar


def get_comet_slide(registrar, slide_name=None):
    """Get the COMET (non-reference) slide from registrar.

    Args:
        registrar: VALIS registrar object
        slide_name: Optional slide name filter

    Returns:
        Slide object for the COMET image
    """
    for name, slide in registrar.slide_dict.items():
        if not slide.is_ref:
            if slide_name is None or slide_name in name:
                return slide

    raise ValueError(
        f"Could not find COMET slide. Available: {list(registrar.slide_dict.keys())}"
    )


def warp_geojson_with_valis(registrar, geojson_path, output_path, slide_name=None):
    """Warp GeoJSON cell segmentations using VALIS transformation.

    Args:
        registrar: VALIS registrar object (after registration)
        geojson_path: Path to input GeoJSON (COMET segmentations)
        output_path: Path to save warped GeoJSON
        slide_name: Name of the COMET slide (auto-detected if None)

    Returns:
        dict: Warped GeoJSON data
    """
    geojson_path = Path(geojson_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    comet_slide = get_comet_slide(registrar, slide_name)
    logger.info(f"Warping GeoJSON using slide: {comet_slide.name}")

    warped_geojson_str = comet_slide.warp_geojson(geojson_f=str(geojson_path))

    if isinstance(warped_geojson_str, str):
        warped_geojson = json.loads(warped_geojson_str)
    else:
        warped_geojson = warped_geojson_str

    with open(output_path, 'w') as f:
        json.dump(warped_geojson, f)

    n_features = len(warped_geojson.get('features', []))
    logger.info(f"Saved warped GeoJSON ({n_features} features) to {output_path}")

    return warped_geojson


def warp_coordinates(registrar, coordinates, slide_name=None):
    """Warp arbitrary (x, y) coordinates using VALIS transformation.

    Args:
        registrar: VALIS registrar object
        coordinates: (N, 2) array of (x, y) in COMET space

    Returns:
        np.ndarray: (N, 2) warped coordinates in STOmics space
    """
    comet_slide = get_comet_slide(registrar, slide_name)
    return comet_slide.warp_xy(coordinates)


def save_registrar(registrar, save_dir):
    """Save VALIS registrar for later reuse.

    Args:
        registrar: VALIS registrar object
        save_dir: Directory to save registration data
    """
    registered_slide_dst = Path(save_dir) / 'registered_slides'
    registered_slide_dst.mkdir(parents=True, exist_ok=True)
    registrar.warp_and_save_slides(str(registered_slide_dst))
    logger.info(f"Saved registered slides to {registered_slide_dst}")


# =============================================================================
# STOmics Data Loading
# =============================================================================

def load_stomics_cellbin_gef(gef_path):
    """Load STOmics cellbin GEF file into AnnData.

    Tries stereopy first, falls back to h5py.

    Args:
        gef_path: Path to .cellbin.gef file

    Returns:
        anndata.AnnData: Cell x gene expression matrix with spatial coordinates
    """
    gef_path = Path(gef_path)
    logger.info(f"Loading STOmics cellbin: {gef_path.name}")

    # Try stereopy
    try:
        import stereo as st
        data = st.io.read_gef(str(gef_path), bin_type='cell_bins')
        adata = data.to_anndata()
        logger.info(f"Loaded via stereopy: {adata.shape}")
        return adata
    except ImportError:
        logger.info("stereopy not available, trying gefpy")
    except Exception as e:
        logger.warning(f"stereopy failed: {e}, trying gefpy")

    # Try gefpy
    try:
        return _load_cellbin_gefpy(gef_path)
    except ImportError:
        logger.info("gefpy not available, trying h5py")
    except Exception as e:
        logger.warning(f"gefpy failed: {e}, trying h5py")

    # Fallback: raw h5py
    return _load_cellbin_h5py(gef_path)


def _load_cellbin_gefpy(gef_path):
    """Load cellbin GEF using gefpy."""
    from gefpy.cgef_reader_cy import CgefR
    import anndata as ad

    cgef = CgefR(str(gef_path))

    # Get cell data
    cell_num = cgef.get_cell_num()
    gene_num = cgef.get_gene_num()

    logger.info(f"gefpy: {cell_num} cells, {gene_num} genes")

    # Read expression into sparse matrix
    from scipy.sparse import csr_matrix
    import h5py

    with h5py.File(str(gef_path), 'r') as f:
        cell_bin = f['cellBin']

        # Gene names
        genes = cell_bin['gene']
        if genes.dtype.names and 'geneName' in genes.dtype.names:
            gene_names = [g['geneName'].decode('utf-8') if isinstance(g['geneName'], bytes)
                          else str(g['geneName']) for g in genes[:]]
        else:
            gene_names = [str(g) for g in genes[:]]

        # Cell info
        cell_data = cell_bin['cell']
        cell_x = cell_data['x'][:]
        cell_y = cell_data['y'][:]

        # Expression data
        exp_data = cell_bin['geneExp']
        exp_cell_ids = exp_data['cellID'][:]
        exp_gene_ids = exp_data['geneID'][:]
        exp_counts = exp_data['count'][:]

    # Build cell ID mapping
    cell_ids = cell_data['cellID'][:]
    cell_id_to_idx = {int(cid): idx for idx, cid in enumerate(cell_ids)}

    row = np.array([cell_id_to_idx.get(int(cid), -1) for cid in exp_cell_ids])
    valid = row >= 0

    X = csr_matrix(
        (exp_counts[valid].astype(np.float32), (row[valid], exp_gene_ids[valid])),
        shape=(cell_num, gene_num)
    )

    obs = pd.DataFrame({
        'cell_id': cell_ids,
        'x': cell_x.astype(np.float64),
        'y': cell_y.astype(np.float64),
    }, index=[str(i) for i in range(cell_num)])

    var = pd.DataFrame(index=gene_names)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm['spatial'] = np.column_stack([cell_x, cell_y]).astype(np.float64)

    logger.info(f"Loaded cellbin via gefpy: {adata.shape}")
    return adata


def _load_cellbin_h5py(gef_path):
    """Load cellbin GEF using raw h5py as fallback."""
    import h5py
    import anndata as ad
    from scipy.sparse import csr_matrix

    with h5py.File(str(gef_path), 'r') as f:
        if 'cellBin' in f:
            cell_bin = f['cellBin']
        elif 'geneExp' in f:
            cell_bin = f['geneExp']
        else:
            raise ValueError(f"Unexpected GEF structure: {list(f.keys())}")

        # Gene names
        genes = cell_bin['gene']
        if genes.dtype.names and 'geneName' in genes.dtype.names:
            gene_names = [g['geneName'].decode('utf-8') if isinstance(g['geneName'], bytes)
                          else str(g['geneName']) for g in genes[:]]
        else:
            raw = genes[:]
            gene_names = [g.decode('utf-8') if isinstance(g, bytes) else str(g) for g in raw]

        n_genes = len(gene_names)

        # Cell data
        cell_data = cell_bin['cell']
        if isinstance(cell_data, h5py.Dataset):
            cell_ids = cell_data['cellID'][:]
            cell_x = cell_data['x'][:].astype(np.float64)
            cell_y = cell_data['y'][:].astype(np.float64)
        else:
            cell_ids = cell_data['cellID'][:]
            cell_x = cell_data['x'][:].astype(np.float64)
            cell_y = cell_data['y'][:].astype(np.float64)

        n_cells = len(cell_ids)

        # Expression data
        exp_data = cell_bin['geneExp']
        if isinstance(exp_data, h5py.Dataset):
            exp_cell_ids = exp_data['cellID'][:]
            exp_gene_ids = exp_data['geneID'][:]
            exp_counts = exp_data['count'][:]
        else:
            exp_cell_ids = exp_data['cellID'][:]
            exp_gene_ids = exp_data['geneID'][:]
            exp_counts = exp_data['count'][:]

    # Build sparse matrix
    cell_id_to_idx = {int(cid): idx for idx, cid in enumerate(cell_ids)}

    row = np.array([cell_id_to_idx.get(int(cid), -1) for cid in exp_cell_ids])
    valid = row >= 0

    X = csr_matrix(
        (exp_counts[valid].astype(np.float32), (row[valid], exp_gene_ids[valid])),
        shape=(n_cells, n_genes)
    )

    obs = pd.DataFrame({
        'cell_id': cell_ids,
        'x': cell_x,
        'y': cell_y,
    }, index=[str(i) for i in range(n_cells)])

    var = pd.DataFrame(index=gene_names)

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obsm['spatial'] = np.column_stack([cell_x, cell_y])

    logger.info(f"Loaded cellbin via h5py: {adata.shape}")
    return adata


def load_stomics_h5ad(h5ad_path):
    """Load pre-processed STOmics h5ad file."""
    import anndata as ad

    h5ad_path = Path(h5ad_path)
    logger.info(f"Loading STOmics h5ad: {h5ad_path.name}")
    adata = ad.read_h5ad(str(h5ad_path))
    logger.info(f"Loaded: {adata.shape}")
    return adata


# =============================================================================
# Cellbin Comparison: Validation, Border Loading, Rasterization
# =============================================================================

def validate_cellbin_comparison_inputs(sample_id, chip_id, output_dir,
                                        warped_geojson_name=None,
                                        cellbin_path=None, dapi_path=None):
    """Preflight validation for COMET vs STOmics cellbin comparison.

    Checks that all required data exists and is usable before running
    the comparison pipeline.

    Args:
        sample_id: Sample ID (e.g., 'SO34')
        chip_id: STOmics chip ID (e.g., 'A03979E2')
        output_dir: Directory where warped GeoJSON was saved
        warped_geojson_name: Override warped GeoJSON filename
            (default: {sample_id}_warped_segmentations.geojson)
        cellbin_path: Explicit path to cellbin GEF (overrides DefaultPaths)
        dapi_path: Explicit path to STOmics DAPI (overrides DefaultPaths)

    Returns:
        dict: Status per check with 'pass', 'message', and overall 'ready'
    """
    import h5py

    output_dir = Path(output_dir)
    status = {
        'sample_id': sample_id,
        'chip_id': chip_id,
        'checks': {},
        'ready': True,
    }

    def _check(name, passed, message):
        status['checks'][name] = {'pass': passed, 'message': message}
        if not passed:
            status['ready'] = False

    # 1. STOmics cellbin GEF exists
    if cellbin_path is None:
        cellbin_path = DefaultPaths.stomics_cellbin_path(chip_id)
    cellbin_path = Path(cellbin_path)
    _check('cellbin_exists',
           cellbin_path.exists(),
           f"Found: {cellbin_path}" if cellbin_path.exists()
           else f"Missing: {cellbin_path}")

    if not cellbin_path.exists():
        # Skip remaining cellbin checks
        _check('cellbin_valid', False, "Skipped (cellbin not found)")
        _check('cellbin_borders', False, "Skipped (cellbin not found)")
    else:
        # 2. Cellbin is valid (>1 cell)
        try:
            with h5py.File(str(cellbin_path), 'r') as f:
                cell_data = f['cellBin']['cell']
                n_cells = len(cell_data)
                valid = n_cells > 1
                _check('cellbin_valid', valid,
                       f"{n_cells:,} cells" if valid
                       else f"Only {n_cells} cell — broken cellbin (super-cell aggregate)")

                # 3. Cellbin has border data
                has_borders = 'cellBorder' in f['cellBin']
                if has_borders:
                    border_shape = f['cellBin']['cellBorder'].shape
                    _check('cellbin_borders', True,
                           f"cellBorder present: shape {border_shape}")
                else:
                    _check('cellbin_borders', False,
                           "cellBorder dataset missing from cellBin group")
        except Exception as e:
            _check('cellbin_valid', False, f"Error reading cellbin: {e}")
            _check('cellbin_borders', False, "Skipped (cellbin read error)")

    # 4. STOmics DAPI exists
    if dapi_path is None:
        dapi_path = DefaultPaths.stomics_dapi_path(chip_id)
    dapi_path = Path(dapi_path)
    _check('stomics_dapi',
           dapi_path.exists(),
           f"Found: {dapi_path}" if dapi_path.exists()
           else f"Missing: {dapi_path}")

    # 5. COMET warped GeoJSON exists
    if warped_geojson_name is None:
        warped_geojson_name = f'{sample_id}_warped_segmentations.geojson'
    warped_path = output_dir / warped_geojson_name
    if warped_path.exists():
        # Quick sanity check: file has features
        try:
            with open(warped_path) as f:
                data = json.load(f)
            n_feat = len(data.get('features', []))
            _check('warped_geojson', n_feat > 0,
                   f"Found: {warped_path.name} ({n_feat:,} features)" if n_feat > 0
                   else f"Found but empty: {warped_path.name}")
        except Exception as e:
            _check('warped_geojson', False, f"Found but unreadable: {e}")
    else:
        # Check if source data exists to guide user
        geojson_exists = DefaultPaths.geojson_path(sample_id).exists()
        hint = (" Source GeoJSON exists — run alignment pipeline first."
                if geojson_exists else " Source GeoJSON also missing.")
        _check('warped_geojson', False,
               f"Missing: {warped_path}{hint}")

    # Log summary
    n_pass = sum(1 for c in status['checks'].values() if c['pass'])
    n_total = len(status['checks'])
    logger.info(f"Preflight validation: {n_pass}/{n_total} checks passed"
                f" — {'READY' if status['ready'] else 'NOT READY'}")
    for name, check in status['checks'].items():
        symbol = 'PASS' if check['pass'] else 'FAIL'
        logger.info(f"  [{symbol}] {name}: {check['message']}")

    return status


def load_stomics_cellbin_borders(gef_path):
    """Load STOmics cellbin cell metadata and border polygons from GEF.

    Reads /cellBin/cell (centroids, areas, IDs) and /cellBin/cellBorder
    (int16 offsets from centroid), reconstructing absolute polygon vertices.

    Args:
        gef_path: Path to .adjusted.cellbin.gef file

    Returns:
        dict with keys:
            'cell_ids': np.ndarray of cell IDs (uint32)
            'centroids': np.ndarray (n_cells, 2) of (x, y) centroids
            'areas': np.ndarray of cell areas
            'gene_counts': np.ndarray of genes per cell
            'dnb_counts': np.ndarray of DNBs per cell
            'polygons': list of np.ndarray (32, 2) absolute (x, y) coordinates
            'n_cells': int
    """
    import h5py

    gef_path = Path(gef_path)
    logger.info(f"Loading cellbin borders: {gef_path.name}")

    with h5py.File(str(gef_path), 'r') as f:
        cell_bin = f['cellBin']
        cell_data = cell_bin['cell'][:]
        border_data = cell_bin['cellBorder'][:]  # (n_cells, 32, 2) int16

    n_cells = len(cell_data)
    logger.info(f"Read {n_cells:,} cells with {border_data.shape[1]}-vertex borders")

    # Extract cell metadata
    cell_ids = cell_data['id'].astype(np.uint32)
    cx = cell_data['x'].astype(np.float64)
    cy = cell_data['y'].astype(np.float64)
    centroids = np.column_stack([cx, cy])
    areas = cell_data['area'].astype(np.int32)
    gene_counts = cell_data['geneCount'].astype(np.int32)
    dnb_counts = cell_data['dnbCount'].astype(np.int32)

    # Reconstruct absolute polygon coordinates from centroid + offset
    # border_data shape: (n_cells, 32, 2) with int16 offsets (dx, dy)
    border_offsets = border_data.astype(np.float64)  # (n_cells, 32, 2)
    centroids_broadcast = centroids[:, np.newaxis, :]  # (n_cells, 1, 2)
    abs_polygons = centroids_broadcast + border_offsets  # (n_cells, 32, 2)

    # Convert to list of per-cell polygon arrays
    polygons = [abs_polygons[i] for i in range(n_cells)]

    # Log coordinate ranges
    all_x = abs_polygons[:, :, 0]
    all_y = abs_polygons[:, :, 1]
    logger.info(f"Polygon coordinate ranges: "
                f"X [{all_x.min():.0f}, {all_x.max():.0f}], "
                f"Y [{all_y.min():.0f}, {all_y.max():.0f}]")
    logger.info(f"Cell area: median={np.median(areas):.0f}, "
                f"mean={np.mean(areas):.0f}, "
                f"gene count: median={np.median(gene_counts):.0f}")

    return {
        'cell_ids': cell_ids,
        'centroids': centroids,
        'areas': areas,
        'gene_counts': gene_counts,
        'dnb_counts': dnb_counts,
        'polygons': polygons,
        'n_cells': n_cells,
    }


def rasterize_cellbin_to_mask(cellbin_data, mask_shape, output_path=None):
    """Rasterize STOmics cellbin polygons into a labeled uint32 mask.

    Args:
        cellbin_data: Dict from load_stomics_cellbin_borders()
        mask_shape: (height, width) matching STOmics DAPI dimensions
        output_path: Optional path to save mask as compressed TIFF

    Returns:
        np.ndarray: uint32 mask of shape (height, width), 0 = background
    """
    import cv2

    height, width = mask_shape
    mask = np.zeros((height, width), dtype=np.int32)

    cell_ids = cellbin_data['cell_ids']
    polygons = cellbin_data['polygons']
    n_cells = cellbin_data['n_cells']

    logger.info(f"Rasterizing {n_cells:,} cellbin polygons to mask "
                f"({height}x{width})")

    n_rasterized = 0
    n_skipped = 0

    for i in range(n_cells):
        poly = polygons[i]  # (32, 2) as float64 (x, y)

        # Skip degenerate polygons (all vertices at centroid)
        if np.allclose(poly, poly[0]):
            n_skipped += 1
            continue

        pts = np.round(poly).astype(np.int32)

        # Clip to mask bounds
        pts[:, 0] = np.clip(pts[:, 0], 0, width - 1)
        pts[:, 1] = np.clip(pts[:, 1], 0, height - 1)

        # Use cell_id as label (cv2.fillPoly needs int32)
        label_val = int(cell_ids[i])
        if label_val == 0:
            label_val = i + 1  # avoid 0 (background)

        cv2.fillPoly(mask, [pts], label_val)
        n_rasterized += 1

    mask = mask.astype(np.uint32)

    n_unique = len(np.unique(mask)) - 1  # exclude background
    fg_frac = (mask > 0).sum() / mask.size * 100
    logger.info(f"Rasterized {n_rasterized:,} cells ({n_skipped} degenerate skipped), "
                f"{n_unique:,} unique labels, {fg_frac:.1f}% foreground")

    if output_path is not None:
        import tifffile
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            tifffile.imwrite(str(output_path), mask, compression='lzw')
        except KeyError:
            tifffile.imwrite(str(output_path), mask, compression='zlib')
        logger.info(f"Saved mask to {output_path} "
                    f"({output_path.stat().st_size / 1e6:.1f} MB)")

    return mask


def export_cellbin_borders_to_geojson(gef_path, output_path, subsample=1):
    """Export STOmics cellbin borders as QuPath-compatible GeoJSON.

    Reconstructs cell polygons from the cellbin GEF and writes them as
    a GeoJSON FeatureCollection that can be loaded in QuPath alongside
    the COMET warped GeoJSON for visual alignment comparison.

    Args:
        gef_path: Path to .adjusted.cellbin.gef file
        output_path: Path to save GeoJSON
        subsample: Export every Nth cell (1 = all, 10 = every 10th).
            Useful for QuPath performance with 500K+ cells.

    Returns:
        dict: GeoJSON FeatureCollection
    """
    gef_path = Path(gef_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    cellbin_data = load_stomics_cellbin_borders(gef_path)

    cell_ids = cellbin_data['cell_ids']
    polygons = cellbin_data['polygons']
    areas = cellbin_data['areas']
    gene_counts = cellbin_data['gene_counts']
    n_cells = cellbin_data['n_cells']

    indices = range(0, n_cells, subsample)

    features = []
    n_exported = 0

    for i in indices:
        poly = polygons[i]

        # Skip degenerate polygons
        if np.allclose(poly, poly[0]):
            continue

        # Close the polygon ring
        coords = poly.tolist()
        if coords[0] != coords[-1]:
            coords.append(coords[0])

        feature = {
            'type': 'Feature',
            'geometry': {
                'type': 'Polygon',
                'coordinates': [coords],
            },
            'properties': {
                'label': int(cell_ids[i]),
                'area_px': int(areas[i]),
                'geneCount': int(gene_counts[i]),
                'source': 'stomics_cellbin',
            },
        }
        features.append(feature)
        n_exported += 1

    geojson = {
        'type': 'FeatureCollection',
        'features': features,
    }

    with open(output_path, 'w') as f:
        json.dump(geojson, f)

    file_size = output_path.stat().st_size / 1e6
    logger.info(f"Exported {n_exported:,} cellbin polygons to {output_path.name} "
                f"({file_size:.1f} MB, subsample={subsample})")

    return geojson


# =============================================================================
# Gene Expression Aggregation
# =============================================================================

def aggregate_expression_per_comet_cell(warped_geojson, stomics_adata,
                                         method='nearest', max_distance=50):
    """Aggregate STOmics gene expression within each warped COMET cell.

    For each COMET cell (warped to STOmics space), finds the STOmics data
    points within or near the cell and aggregates their gene expression.

    Args:
        warped_geojson: Warped GeoJSON dict (COMET cells in STOmics coords)
        stomics_adata: AnnData with STOmics expression and spatial coordinates
        method: 'nearest' (KD-tree, fast) or 'polygon' (point-in-polygon, precise)
        max_distance: Max distance for nearest neighbor assignment (pixels)

    Returns:
        anndata.AnnData: COMET cells x genes with aggregated expression
    """
    import anndata as ad

    # Get STOmics coordinates
    if 'spatial' in stomics_adata.obsm:
        stomics_coords = stomics_adata.obsm['spatial']
    else:
        stomics_coords = np.column_stack([
            stomics_adata.obs['x'].values,
            stomics_adata.obs['y'].values
        ])

    # Get expression matrix as dense
    X = stomics_adata.X
    if hasattr(X, 'toarray'):
        X = X.toarray()
    X = np.array(X, dtype=np.float32)

    features = warped_geojson['features']
    n_cells = len(features)
    n_genes = X.shape[1]

    logger.info(f"Aggregating expression for {n_cells} COMET cells "
                f"from {len(stomics_coords)} STOmics points, {n_genes} genes")

    if method == 'nearest':
        return _aggregate_by_nearest(features, stomics_coords, X,
                                     stomics_adata, n_cells, n_genes, max_distance)
    elif method == 'polygon':
        return _aggregate_by_polygon(features, stomics_coords, X,
                                     stomics_adata, n_cells, n_genes)
    else:
        raise ValueError(f"Unknown method: {method}. Use 'nearest' or 'polygon'.")


def _aggregate_by_nearest(features, stomics_coords, X, stomics_adata,
                           n_cells, n_genes, max_distance):
    """Aggregate by nearest STOmics point to COMET cell centroid (fast)."""
    import anndata as ad

    # Compute COMET centroids from warped polygons
    centroids = []
    cell_labels = []
    cell_areas = []

    for feature in features:
        coords = feature['geometry']['coordinates']
        exterior = np.array(coords[0])
        cx = np.mean(exterior[:-1, 0])
        cy = np.mean(exterior[:-1, 1])
        centroids.append([cx, cy])
        cell_labels.append(feature['properties'].get('label', 0))
        cell_areas.append(feature['properties'].get('area_px', 0))

    centroids = np.array(centroids)

    # KD-tree nearest neighbor
    tree = cKDTree(stomics_coords)
    distances, indices = tree.query(centroids, k=1)

    # Assign expression values
    agg_X = X[indices].copy()

    # Zero out cells beyond max_distance
    beyond = distances > max_distance
    agg_X[beyond] = 0

    obs = pd.DataFrame({
        'cell_label': cell_labels,
        'centroid_x': centroids[:, 0],
        'centroid_y': centroids[:, 1],
        'area_px': cell_areas,
        'nn_distance': distances,
        'within_threshold': ~beyond,
    }, index=[str(i) for i in range(n_cells)])

    var = stomics_adata.var.copy()

    result = ad.AnnData(X=agg_X, obs=obs, var=var)
    result.obsm['spatial'] = centroids

    n_within = (~beyond).sum()
    logger.info(f"Nearest-neighbor assignment: {n_within}/{n_cells} cells "
                f"within {max_distance}px")
    logger.info(f"Distance stats: min={distances.min():.1f}, "
                f"median={np.median(distances):.1f}, max={distances.max():.1f}")

    return result


def _aggregate_by_polygon(features, stomics_coords, X, stomics_adata,
                           n_cells, n_genes):
    """Aggregate by point-in-polygon containment (precise, slower).

    Uses shapely STRtree for efficient spatial queries.
    For large datasets (>50K cells), consider using 'nearest' instead.
    """
    import anndata as ad
    from shapely.geometry import shape, Point
    from shapely import STRtree

    # Build spatial index for STOmics points
    stomics_points = [Point(x, y) for x, y in stomics_coords]
    tree = STRtree(stomics_points)

    # Pre-build index mapping
    point_indices = np.arange(len(stomics_points))

    agg_X = np.zeros((n_cells, n_genes), dtype=np.float32)
    cell_labels = []
    cell_centroids = []
    cell_areas = []
    n_stomics_per_cell = []

    for i, feature in enumerate(features):
        if i % 10000 == 0 and i > 0:
            logger.info(f"  Processing cell {i}/{n_cells}...")

        try:
            poly = shape(feature['geometry'])
            if not poly.is_valid:
                poly = poly.buffer(0)
        except Exception:
            cell_labels.append(feature['properties'].get('label', i))
            cell_centroids.append([0, 0])
            cell_areas.append(0)
            n_stomics_per_cell.append(0)
            continue

        # Query spatial index: find points inside polygon
        inside_idx = tree.query(poly, predicate='contains')

        if len(inside_idx) > 0:
            agg_X[i] = np.mean(X[inside_idx], axis=0)

        cell_labels.append(feature['properties'].get('label', i))
        centroid = poly.centroid
        cell_centroids.append([centroid.x, centroid.y])
        cell_areas.append(poly.area)
        n_stomics_per_cell.append(len(inside_idx))

    obs = pd.DataFrame({
        'cell_label': cell_labels,
        'centroid_x': [c[0] for c in cell_centroids],
        'centroid_y': [c[1] for c in cell_centroids],
        'area_px': cell_areas,
        'n_stomics_points': n_stomics_per_cell,
    }, index=[str(i) for i in range(n_cells)])

    var = stomics_adata.var.copy()

    result = ad.AnnData(X=agg_X, obs=obs, var=var)
    result.obsm['spatial'] = np.array(cell_centroids)

    n_with_data = sum(1 for n in n_stomics_per_cell if n > 0)
    logger.info(f"Polygon aggregation: {n_with_data}/{n_cells} cells have expression data")

    return result


# =============================================================================
# COMET <-> STOmics Cell Mapping
# =============================================================================

def map_comet_to_stomics_cells(comet_centroids, stomics_centroids,
                                max_distance=50):
    """Map COMET cells to nearest STOmics cells using KD-tree.

    Args:
        comet_centroids: (N, 2) COMET cell centroids (in STOmics space)
        stomics_centroids: (M, 2) STOmics cell centroids
        max_distance: Maximum distance for valid mapping (pixels)

    Returns:
        pd.DataFrame: Mapping table
    """
    tree = cKDTree(stomics_centroids)
    distances, indices = tree.query(comet_centroids, k=1)

    mapping = pd.DataFrame({
        'comet_idx': range(len(comet_centroids)),
        'stomics_idx': indices,
        'distance': distances,
        'mapped': distances <= max_distance,
    })

    n_mapped = mapping['mapped'].sum()
    logger.info(f"Cell mapping: {n_mapped}/{len(comet_centroids)} COMET cells "
                f"mapped within {max_distance}px")

    return mapping


def compare_segmentation_cells(comet_centroids, stomics_centroids,
                                max_distance=30):
    """Compare two cell segmentations via bidirectional centroid matching.

    Classifies every cell from both segmentations into one of:
    - 1:1 match: mutual nearest neighbors within threshold
    - fragmented: one COMET cell is nearest for multiple STOmics cells
    - merged: multiple COMET cells share the same nearest STOmics cell
    - comet_only: COMET cell with no STOmics neighbor within threshold
    - stomics_only: STOmics cell with no COMET neighbor within threshold

    Args:
        comet_centroids: (N, 2) array of COMET cell centroids
        stomics_centroids: (M, 2) array of STOmics cell centroids
        max_distance: Maximum distance (px) for a valid match

    Returns:
        dict with:
            'summary': counts per category
            'matches_1to1': DataFrame of mutual nearest-neighbor pairs
            'comet_classification': array of category per COMET cell
            'stomics_classification': array of category per STOmics cell
            'comet_to_stomics': DataFrame (comet_idx, stomics_idx, distance)
            'stomics_to_comet': DataFrame (stomics_idx, comet_idx, distance)
    """
    comet_centroids = np.asarray(comet_centroids)
    stomics_centroids = np.asarray(stomics_centroids)
    n_comet = len(comet_centroids)
    n_stomics = len(stomics_centroids)

    logger.info(f"Comparing segmentations: {n_comet:,} COMET vs "
                f"{n_stomics:,} STOmics cells (max_distance={max_distance}px)")

    # Build KD-trees
    comet_tree = cKDTree(comet_centroids)
    stomics_tree = cKDTree(stomics_centroids)

    # Forward: each COMET cell's nearest STOmics cell
    c2s_dist, c2s_idx = stomics_tree.query(comet_centroids, k=1)
    # Reverse: each STOmics cell's nearest COMET cell
    s2c_dist, s2c_idx = comet_tree.query(stomics_centroids, k=1)

    # Classify COMET cells
    comet_class = np.full(n_comet, 'comet_only', dtype='U20')
    stomics_class = np.full(n_stomics, 'stomics_only', dtype='U20')

    # Step 1: Find mutual nearest neighbors (1:1 matches)
    for ci in range(n_comet):
        if c2s_dist[ci] > max_distance:
            continue
        si = c2s_idx[ci]
        # Check if this STOmics cell's nearest COMET cell is ci (mutual)
        if s2c_idx[si] == ci and s2c_dist[si] <= max_distance:
            comet_class[ci] = '1to1'
            stomics_class[si] = '1to1'

    # Step 2: Classify remaining within-threshold cells
    # COMET cells within threshold but not 1:1
    for ci in range(n_comet):
        if comet_class[ci] != 'comet_only':
            continue
        if c2s_dist[ci] <= max_distance:
            # This COMET cell has a nearby STOmics cell but isn't its mutual NN
            # → multiple COMET cells share one STOmics cell = merged
            comet_class[ci] = 'merged'

    # STOmics cells within threshold but not 1:1
    for si in range(n_stomics):
        if stomics_class[si] != 'stomics_only':
            continue
        if s2c_dist[si] <= max_distance:
            # This STOmics cell has a nearby COMET cell but isn't its mutual NN
            # → one COMET cell maps to multiple STOmics cells = fragmented
            stomics_class[si] = 'fragmented'

    # Build summary
    from collections import Counter
    comet_counts = Counter(comet_class)
    stomics_counts = Counter(stomics_class)

    n_1to1 = comet_counts.get('1to1', 0)
    summary = {
        'n_comet': n_comet,
        'n_stomics': n_stomics,
        'n_1to1': n_1to1,
        'n_comet_merged': comet_counts.get('merged', 0),
        'n_comet_only': comet_counts.get('comet_only', 0),
        'n_stomics_fragmented': stomics_counts.get('fragmented', 0),
        'n_stomics_only': stomics_counts.get('stomics_only', 0),
        'pct_comet_matched': (n_comet - comet_counts.get('comet_only', 0)) / n_comet * 100,
        'pct_stomics_matched': (n_stomics - stomics_counts.get('stomics_only', 0)) / n_stomics * 100,
    }

    # Build 1:1 match DataFrame
    match_mask = comet_class == '1to1'
    match_comet_idx = np.where(match_mask)[0]
    match_stomics_idx = c2s_idx[match_mask]
    match_distances = c2s_dist[match_mask]

    matches_1to1 = pd.DataFrame({
        'comet_idx': match_comet_idx,
        'stomics_idx': match_stomics_idx,
        'distance': match_distances,
    })

    # Forward/reverse mapping DataFrames
    c2s_df = pd.DataFrame({
        'comet_idx': range(n_comet),
        'stomics_idx': c2s_idx,
        'distance': c2s_dist,
        'classification': comet_class,
    })
    s2c_df = pd.DataFrame({
        'stomics_idx': range(n_stomics),
        'comet_idx': s2c_idx,
        'distance': s2c_dist,
        'classification': stomics_class,
    })

    logger.info(f"Match results: {n_1to1:,} 1:1 matches, "
                f"{summary['n_comet_merged']:,} COMET merged, "
                f"{summary['n_comet_only']:,} COMET-only, "
                f"{summary['n_stomics_fragmented']:,} STOmics fragmented, "
                f"{summary['n_stomics_only']:,} STOmics-only")
    logger.info(f"Match rate: {summary['pct_comet_matched']:.1f}% COMET, "
                f"{summary['pct_stomics_matched']:.1f}% STOmics")
    if len(matches_1to1) > 0:
        logger.info(f"1:1 match distances: median={matches_1to1['distance'].median():.1f}px, "
                    f"mean={matches_1to1['distance'].mean():.1f}px")

    return {
        'summary': summary,
        'matches_1to1': matches_1to1,
        'comet_classification': comet_class,
        'stomics_classification': stomics_class,
        'comet_to_stomics': c2s_df,
        'stomics_to_comet': s2c_df,
    }


# =============================================================================
# Cellbin Comparison: Pixel Overlap & Expression Concordance
# =============================================================================

def compute_mask_overlap_metrics(comet_mask, stomics_mask, matches_1to1=None):
    """Compute pixel-level overlap metrics between two cell segmentation masks.

    Args:
        comet_mask: uint32 mask (H, W) — COMET cell labels, 0=background
        stomics_mask: uint32 mask (H, W) — STOmics cell labels, 0=background
        matches_1to1: Optional DataFrame with 'comet_idx' and 'stomics_idx'
            columns from compare_segmentation_cells() for per-cell IoU

    Returns:
        dict with overlap metrics
    """
    assert comet_mask.shape == stomics_mask.shape, \
        f"Mask shapes differ: {comet_mask.shape} vs {stomics_mask.shape}"

    comet_fg = comet_mask > 0
    stomics_fg = stomics_mask > 0

    total_px = comet_mask.size
    both_bg = (~comet_fg) & (~stomics_fg)
    both_fg = comet_fg & stomics_fg
    comet_only = comet_fg & (~stomics_fg)
    stomics_only = (~comet_fg) & stomics_fg

    # Pixel agreement: fraction where both agree (both fg or both bg)
    pixel_agreement = (both_bg.sum() + both_fg.sum()) / total_px

    # Foreground Jaccard: intersection / union of foreground regions
    fg_union = comet_fg | stomics_fg
    fg_jaccard = both_fg.sum() / fg_union.sum() if fg_union.any() else 0.0

    metrics = {
        'total_pixels': int(total_px),
        'comet_fg_pixels': int(comet_fg.sum()),
        'stomics_fg_pixels': int(stomics_fg.sum()),
        'both_fg_pixels': int(both_fg.sum()),
        'comet_only_pixels': int(comet_only.sum()),
        'stomics_only_pixels': int(stomics_only.sum()),
        'comet_coverage_pct': float(comet_fg.sum() / total_px * 100),
        'stomics_coverage_pct': float(stomics_fg.sum() / total_px * 100),
        'pixel_agreement': float(pixel_agreement),
        'foreground_jaccard': float(fg_jaccard),
    }

    logger.info(f"Mask overlap: COMET {metrics['comet_coverage_pct']:.1f}% fg, "
                f"STOmics {metrics['stomics_coverage_pct']:.1f}% fg, "
                f"Jaccard {fg_jaccard:.3f}, agreement {pixel_agreement:.3f}")

    # Per-cell IoU for 1:1 matched pairs (if provided)
    if matches_1to1 is not None and len(matches_1to1) > 0:
        comet_labels_in_mask = np.unique(comet_mask)
        stomics_labels_in_mask = np.unique(stomics_mask)

        # Sample up to 5000 pairs for performance
        n_sample = min(len(matches_1to1), 5000)
        sampled = matches_1to1.sample(n=n_sample, random_state=42) \
            if len(matches_1to1) > n_sample else matches_1to1

        ious = []
        for _, row in sampled.iterrows():
            # We need actual cell label values, not indices
            # The match DataFrame has indices into the centroid arrays
            # For pixel IoU we need to find which label values correspond
            # This is expensive for 500K cells, so we compute a simpler
            # metric: for each matched pair's centroid region, check overlap
            ci, si = int(row['comet_idx']), int(row['stomics_idx'])
            # Skip — we'll compute a centroid-neighborhood overlap instead
            pass

        # Simpler approach: compute overlap fraction in matched-pair regions
        # For each 1:1 pair, check if the pixel at the STOmics centroid
        # is also foreground in the COMET mask (and vice versa)
        # This is a fast proxy for spatial agreement
        logger.info(f"Per-cell IoU deferred — using foreground Jaccard as proxy")
        metrics['n_matched_pairs'] = len(matches_1to1)

    return metrics


def compare_matched_cell_expression(comet_adata, stomics_adata, matches_1to1,
                                     min_cells=10):
    """Compare gene expression between matched COMET and STOmics cells.

    For each 1:1 matched cell pair, compares expression vectors to assess
    concordance between the two segmentation approaches.

    Args:
        comet_adata: AnnData from aggregate_transcripts_by_mask (COMET mask)
        stomics_adata: AnnData from load_stomics_cellbin_gef (STOmics mask)
        matches_1to1: DataFrame with 'comet_idx' and 'stomics_idx' columns
        min_cells: Minimum cells expressing a gene for per-gene correlation

    Returns:
        dict with:
            'per_cell': DataFrame with per-cell concordance metrics
            'per_gene': DataFrame with per-gene correlation across matched cells
            'summary': dict with aggregate metrics
    """
    from scipy.stats import pearsonr, spearmanr
    from scipy.sparse import issparse

    # Find common genes
    comet_genes = set(comet_adata.var_names)
    stomics_genes = set(stomics_adata.var_names)
    common_genes = sorted(comet_genes & stomics_genes)
    n_common = len(common_genes)

    logger.info(f"Expression comparison: {len(comet_genes):,} COMET genes, "
                f"{len(stomics_genes):,} STOmics genes, "
                f"{n_common:,} common")

    if n_common == 0:
        logger.warning("No common genes between COMET and STOmics AnnData")
        return {'per_cell': pd.DataFrame(), 'per_gene': pd.DataFrame(),
                'summary': {'n_common_genes': 0}}

    # Subset to common genes
    comet_sub = comet_adata[:, common_genes]
    stomics_sub = stomics_adata[:, common_genes]

    # Get expression matrices
    comet_X = comet_sub.X.toarray() if issparse(comet_sub.X) else np.asarray(comet_sub.X)
    stomics_X = stomics_sub.X.toarray() if issparse(stomics_sub.X) else np.asarray(stomics_sub.X)

    # Per-cell metrics for matched pairs
    n_pairs = len(matches_1to1)
    per_cell_records = []

    # Sample for performance if too many pairs
    n_sample = min(n_pairs, 10000)
    sampled = matches_1to1.sample(n=n_sample, random_state=42) \
        if n_pairs > n_sample else matches_1to1

    for _, row in sampled.iterrows():
        ci = int(row['comet_idx'])
        si = int(row['stomics_idx'])

        if ci >= comet_X.shape[0] or si >= stomics_X.shape[0]:
            continue

        c_expr = comet_X[ci]
        s_expr = stomics_X[si]

        c_total = c_expr.sum()
        s_total = s_expr.sum()
        c_ngenes = (c_expr > 0).sum()
        s_ngenes = (s_expr > 0).sum()
        shared = ((c_expr > 0) & (s_expr > 0)).sum()

        # Pearson correlation (only if both have variation)
        if c_expr.std() > 0 and s_expr.std() > 0:
            r, p = pearsonr(c_expr, s_expr)
        else:
            r, p = 0.0, 1.0

        per_cell_records.append({
            'comet_idx': ci,
            'stomics_idx': si,
            'distance': row['distance'],
            'comet_total_counts': float(c_total),
            'stomics_total_counts': float(s_total),
            'comet_n_genes': int(c_ngenes),
            'stomics_n_genes': int(s_ngenes),
            'shared_genes': int(shared),
            'pearson_r': float(r),
            'pearson_p': float(p),
        })

    per_cell_df = pd.DataFrame(per_cell_records)

    # Per-gene correlation across all matched cells
    per_gene_records = []
    comet_matched = comet_X[sampled['comet_idx'].values]
    stomics_matched = stomics_X[sampled['stomics_idx'].values]

    for gi, gene in enumerate(common_genes):
        c_vals = comet_matched[:, gi]
        s_vals = stomics_matched[:, gi]

        # Only compute for genes expressed in enough cells
        n_expr = ((c_vals > 0) | (s_vals > 0)).sum()
        if n_expr < min_cells:
            continue

        if c_vals.std() > 0 and s_vals.std() > 0:
            r, p = pearsonr(c_vals, s_vals)
        else:
            r, p = 0.0, 1.0

        per_gene_records.append({
            'gene': gene,
            'n_expressing': int(n_expr),
            'comet_mean': float(c_vals.mean()),
            'stomics_mean': float(s_vals.mean()),
            'pearson_r': float(r),
            'pearson_p': float(p),
        })

    per_gene_df = pd.DataFrame(per_gene_records)

    # Aggregate summary
    summary = {
        'n_common_genes': n_common,
        'n_pairs_compared': len(per_cell_df),
        'n_genes_correlated': len(per_gene_df),
    }
    if len(per_cell_df) > 0:
        summary['median_pearson_r'] = float(per_cell_df['pearson_r'].median())
        summary['mean_pearson_r'] = float(per_cell_df['pearson_r'].mean())
        summary['pct_r_above_0.5'] = float((per_cell_df['pearson_r'] > 0.5).mean() * 100)
        summary['median_shared_genes'] = float(per_cell_df['shared_genes'].median())
        summary['median_comet_counts'] = float(per_cell_df['comet_total_counts'].median())
        summary['median_stomics_counts'] = float(per_cell_df['stomics_total_counts'].median())

    if len(per_gene_df) > 0:
        summary['median_gene_pearson_r'] = float(per_gene_df['pearson_r'].median())
        summary['n_genes_r_above_0.5'] = int((per_gene_df['pearson_r'] > 0.5).sum())

    logger.info(f"Expression concordance: {summary.get('median_pearson_r', 'N/A')} "
                f"median per-cell r, {summary.get('n_genes_r_above_0.5', 0)} genes "
                f"with r>0.5")

    return {
        'per_cell': per_cell_df,
        'per_gene': per_gene_df,
        'summary': summary,
    }


# =============================================================================
# Validation
# =============================================================================

def compute_alignment_metrics(warped_comet_centroids, stomics_centroids):
    """Compute alignment quality metrics.

    Returns:
        dict with distance statistics and quality assessment
    """
    tree = cKDTree(stomics_centroids)
    distances, _ = tree.query(warped_comet_centroids, k=1)

    metrics = {
        'n_comet_cells': len(warped_comet_centroids),
        'n_stomics_cells': len(stomics_centroids),
        'median_distance': float(np.median(distances)),
        'mean_distance': float(np.mean(distances)),
        'std_distance': float(np.std(distances)),
        'min_distance': float(np.min(distances)),
        'max_distance': float(np.max(distances)),
        'pct_within_30px': float((distances < 30).mean() * 100),
        'pct_within_50px': float((distances < 50).mean() * 100),
        'pct_within_100px': float((distances < 100).mean() * 100),
    }

    if metrics['median_distance'] < 30:
        metrics['quality'] = 'EXCELLENT'
    elif metrics['median_distance'] < 50:
        metrics['quality'] = 'GOOD'
    elif metrics['median_distance'] < 100:
        metrics['quality'] = 'ACCEPTABLE'
    else:
        metrics['quality'] = 'POOR'

    return metrics


def plot_alignment_validation(warped_geojson, stomics_adata=None,
                               stomics_dapi_path=None, output_path=None,
                               sample_n=5000, title=None):
    """Generate alignment validation plots.

    Creates up to 3 subplots:
    1. Warped COMET cell polygons
    2. COMET vs STOmics cell overlay (if stomics_adata provided)
    3. Cells on STOmics DAPI (if stomics_dapi_path provided)
    """
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon
    from matplotlib.collections import PatchCollection

    features = warped_geojson['features']

    # Subsample for performance
    if len(features) > sample_n:
        indices = np.random.choice(len(features), sample_n, replace=False)
        features = [features[i] for i in sorted(indices)]

    n_plots = 1
    if stomics_adata is not None:
        n_plots += 1
    if stomics_dapi_path is not None:
        n_plots += 1

    fig, axes = plt.subplots(1, n_plots, figsize=(7 * n_plots, 7))
    if n_plots == 1:
        axes = [axes]

    ax_idx = 0

    # Plot 1: Warped COMET polygons
    ax = axes[ax_idx]
    patches = []
    for feature in features:
        coords = feature['geometry']['coordinates']
        exterior = np.array(coords[0])
        patches.append(MplPolygon(exterior, closed=True))

    pc = PatchCollection(patches, alpha=0.3, edgecolor='red',
                          facecolor='lightcoral', linewidth=0.5)
    ax.add_collection(pc)

    all_coords = np.vstack([np.array(f['geometry']['coordinates'][0])
                             for f in features])
    ax.set_xlim(all_coords[:, 0].min() - 100, all_coords[:, 0].max() + 100)
    ax.set_ylim(all_coords[:, 1].min() - 100, all_coords[:, 1].max() + 100)
    ax.set_aspect('equal')
    ax.set_title(f'Warped COMET Cells (n={len(features)})')
    ax.set_xlabel('X (STOmics pixels)')
    ax.set_ylabel('Y (STOmics pixels)')
    ax_idx += 1

    # Plot 2: Overlay with STOmics cells
    if stomics_adata is not None:
        ax = axes[ax_idx]

        if 'spatial' in stomics_adata.obsm:
            st_coords = stomics_adata.obsm['spatial']
        else:
            st_coords = np.column_stack([
                stomics_adata.obs['x'], stomics_adata.obs['y']
            ])

        ax.scatter(st_coords[:, 0], st_coords[:, 1],
                   s=0.5, alpha=0.3, c='blue', label='STOmics cells')

        # COMET centroids
        comet_centroids = []
        for feature in features:
            coords = np.array(feature['geometry']['coordinates'][0])
            cx, cy = np.mean(coords[:-1, 0]), np.mean(coords[:-1, 1])
            comet_centroids.append([cx, cy])
        comet_centroids = np.array(comet_centroids)

        ax.scatter(comet_centroids[:, 0], comet_centroids[:, 1],
                   s=1, alpha=0.5, c='red', label='COMET cells')
        ax.set_aspect('equal')
        ax.legend()
        ax.set_title('COMET vs STOmics Cell Overlay')
        ax.set_xlabel('X (STOmics pixels)')
        ax_idx += 1

    # Plot 3: DAPI overlay
    if stomics_dapi_path is not None:
        import tifffile as tiff
        ax = axes[ax_idx]
        dapi = tiff.imread(str(stomics_dapi_path))
        ax.imshow(dapi, cmap='gray', alpha=0.7)

        patches2 = []
        for feature in features[:2000]:
            coords = feature['geometry']['coordinates']
            exterior = np.array(coords[0])
            patches2.append(MplPolygon(exterior, closed=True))

        pc2 = PatchCollection(patches2, alpha=0.2, edgecolor='red',
                               facecolor='none', linewidth=0.3)
        ax.add_collection(pc2)
        ax.set_title('COMET Cells on STOmics DAPI')
        ax_idx += 1

    if title:
        fig.suptitle(title, fontsize=14, fontweight='bold')

    plt.tight_layout()

    if output_path:
        plt.savefig(str(output_path), dpi=150, bbox_inches='tight')
        logger.info(f"Saved validation plot to {output_path}")

    return fig


def plot_distance_distribution(distances, output_path=None,
                                title='Alignment Distance Distribution'):
    """Plot distance distribution histogram with quality thresholds."""
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 5))

    ax.hist(distances, bins=100, alpha=0.7, color='steelblue', edgecolor='white')
    ax.axvline(np.median(distances), color='red', linestyle='--',
               label=f'Median: {np.median(distances):.1f}px')
    ax.axvline(30, color='green', linestyle=':', alpha=0.7, label='30px (excellent)')
    ax.axvline(50, color='orange', linestyle=':', alpha=0.7, label='50px (good)')

    ax.set_xlabel('Distance to nearest point (pixels)')
    ax.set_ylabel('Count')
    ax.set_title(title)
    ax.legend()

    if output_path:
        plt.savefig(str(output_path), dpi=150, bbox_inches='tight')

    return fig


# =============================================================================
# Protein-Gene Correlation Validation
# =============================================================================

# Known COMET protein channels that map to gene names
PROTEIN_GENE_MAP = {
    'CD4': 'CD4',
    'CD8': 'CD8A',
    'CD20': 'MS4A1',
    'CD68': 'CD68',
    'CD163': 'CD163',
    'CD45': 'PTPRC',
    'CD45RO': 'PTPRC',
    'CD31': 'PECAM1',
    'CD10': 'MME',
    'CD133': 'PROM1',
    'CD11C': 'ITGAX',
    'CD56': 'NCAM1',
    'CD66B': 'CEACAM8',
    'CD86': 'CD86',
    'FOXp3': 'FOXP3',
    'TOX1': 'TOX',
    'GZMB': 'GZMB',
    'Periostin': 'POSTN',
    'COL1A1': 'COL1A1',
    'KERATIN8': 'KRT8',
    'aSMA': 'ACTA2',
}


def validate_protein_gene_correlation(integrated_adata, comet_protein_df,
                                       protein_gene_map=None):
    """Validate alignment by correlating COMET protein intensity with STOmics gene expression.

    High correlations (e.g., CD4 protein vs CD4 gene) indicate good alignment.

    Args:
        integrated_adata: AnnData from aggregate_expression_per_comet_cell
        comet_protein_df: DataFrame with COMET protein intensities
            (columns = protein names, index = cell labels)
        protein_gene_map: Dict mapping protein names to gene names
            (defaults to PROTEIN_GENE_MAP)

    Returns:
        pd.DataFrame: Correlation results per protein-gene pair
    """
    from scipy.stats import spearmanr

    if protein_gene_map is None:
        protein_gene_map = PROTEIN_GENE_MAP

    results = []
    gene_names = set(integrated_adata.var_names)

    for protein, gene in protein_gene_map.items():
        if protein not in comet_protein_df.columns:
            continue
        if gene not in gene_names:
            continue

        # Match cells by label
        protein_vals = comet_protein_df[protein].values
        gene_idx = list(integrated_adata.var_names).index(gene)
        gene_vals = integrated_adata.X[:, gene_idx]
        if hasattr(gene_vals, 'toarray'):
            gene_vals = gene_vals.toarray().ravel()
        gene_vals = np.asarray(gene_vals, dtype=float).ravel()

        # Only use cells within threshold
        if 'within_threshold' in integrated_adata.obs.columns:
            mask = integrated_adata.obs['within_threshold'].values
            protein_vals = protein_vals[mask]
            gene_vals = gene_vals[mask]

        # Remove zeros/nans
        valid = (~np.isnan(protein_vals)) & (~np.isnan(gene_vals))
        valid &= (protein_vals > 0) | (gene_vals > 0)

        if valid.sum() < 10:
            continue

        rho, pval = spearmanr(protein_vals[valid], gene_vals[valid])
        results.append({
            'protein': protein,
            'gene': gene,
            'spearman_rho': rho,
            'p_value': pval,
            'n_cells': int(valid.sum()),
        })

    results_df = pd.DataFrame(results)
    if len(results_df) > 0:
        results_df = results_df.sort_values('spearman_rho', ascending=False)
        logger.info(f"Protein-gene correlations ({len(results_df)} pairs):")
        for _, row in results_df.iterrows():
            sig = '***' if row['p_value'] < 0.001 else ('**' if row['p_value'] < 0.01 else '*' if row['p_value'] < 0.05 else '')
            logger.info(f"  {row['protein']:>12s} vs {row['gene']:<10s}: "
                        f"rho={row['spearman_rho']:.3f}{sig} (n={row['n_cells']})")

    return results_df


# =============================================================================
# GeoJSON Export from TIF Masks
# =============================================================================

def export_geojson_from_tif(tif_path, output_path, image_width=None,
                             simplify_tol=1.5, flip_x=True):
    """Export cell segmentation GeoJSON from Visiopharm TIF mask.

    Replicates the logic from visiopharm_to_geojson.ipynb.

    Args:
        tif_path: Path to Visiopharm segmentation TIF
        output_path: Path to save GeoJSON
        image_width: Image width for x-flipping (auto-detected if None)
        simplify_tol: Polygon simplification tolerance
        flip_x: Whether to flip x-coordinates

    Returns:
        dict: GeoJSON FeatureCollection
    """
    import tifffile as tiff
    from skimage.measure import regionprops, find_contours, approximate_polygon, label
    from tqdm import tqdm

    tif_path = Path(tif_path)
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logger.info(f"Loading segmentation mask: {tif_path.name}")
    v_img = tiff.imread(str(tif_path))

    if image_width is None:
        image_width = v_img.shape[2] if len(v_img.shape) > 2 else v_img.shape[1]

    single_cell = label(v_img[0, ...] > 0)
    props = regionprops(single_cell)

    logger.info(f"Found {len(props)} cells, image_width={image_width}")

    def _signed_area(ring_xy):
        x, y = ring_xy[:, 0], ring_xy[:, 1]
        return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)

    def _ensure_closed(r):
        if not np.allclose(r[0], r[-1]):
            r = np.vstack([r, r[0]])
        return r

    features = []
    skipped = 0
    for prop in tqdm(props, desc="Converting cells"):
        # Skip cells too small for contour extraction
        if prop.image.shape[0] < 2 or prop.image.shape[1] < 2:
            skipped += 1
            continue
        cnts = find_contours(prop.image.astype(float), 0.5)
        if not cnts:
            continue

        y0, x0, y1, x1 = prop.bbox
        rings = []

        for cnt in cnts:
            cnt[:, 0] += y0
            cnt[:, 1] += x0
            ring_xy = np.c_[cnt[:, 1], cnt[:, 0]]
            ring_xy = approximate_polygon(ring_xy, tolerance=simplify_tol)
            ring_xy = _ensure_closed(ring_xy)

            if flip_x:
                ring_xy[:, 0] = (image_width - 1) - ring_xy[:, 0]

            rings.append(ring_xy)

        areas = [abs(_signed_area(r)) for r in rings]
        ext_idx = int(np.argmax(areas))
        exterior = rings[ext_idx]
        holes = [rings[i] for i in range(len(rings)) if i != ext_idx]

        if _signed_area(exterior) < 0:
            exterior = exterior[::-1]
        for j, h in enumerate(holes):
            if _signed_area(h) > 0:
                holes[j] = h[::-1]

        rings_xy = [exterior.tolist()] + [h.tolist() for h in holes]

        features.append({
            "type": "Feature",
            "geometry": {"type": "Polygon", "coordinates": rings_xy},
            "properties": {
                "label": int(prop.label),
                "area_px": int(prop.area),
                "centroid_x": float(prop.centroid[1]),
                "centroid_y": float(prop.centroid[0]),
                "bbox": [int(v) for v in prop.bbox],
            }
        })

    if skipped:
        logger.info(f"Skipped {skipped} cells smaller than 2x2 pixels")

    fc = {"type": "FeatureCollection", "features": features}

    with open(output_path, 'w') as f:
        json.dump(fc, f)

    logger.info(f"Saved GeoJSON: {output_path} ({len(features)} cells)")
    return fc


# =============================================================================
# Cellbin GEF Generation from COMET Mask (Path A: SAW-native)
# =============================================================================

def generate_cellbin_with_comet_mask(tissue_gef_path, mask_path, output_dir,
                                      method='FAST', expand_distance=0):
    """Generate cellbin GEF from bin1 tissue GEF + COMET rasterized mask.

    Uses stereopy cell_correct (preferred) or geftools CLI as fallback.
    The output is a SAW-compatible cellbin GEF.

    Note: geftools re-numbers cell IDs via connectedComponentsWithStats.
    Original COMET labels are NOT preserved. Use centroid-based spatial
    join to re-map back to COMET cells.

    Args:
        tissue_gef_path: Path to bin1 tissue.gef or raw.gef
        mask_path: Path to rasterized uint32 TIFF cell mask
        output_dir: Directory for output cellbin GEF files
        method: cell_correct method ('FAST' for full-cell masks, 'EDM' for nuclear)
        expand_distance: Border expansion in pixels (0 for FAST, 10+ for EDM)

    Returns:
        dict with output paths and cell count
    """
    tissue_gef_path = Path(tissue_gef_path)
    mask_path = Path(mask_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if not tissue_gef_path.exists():
        raise FileNotFoundError(f"Tissue GEF not found: {tissue_gef_path}")
    if not mask_path.exists():
        raise FileNotFoundError(f"Mask not found: {mask_path}")

    result = {
        'tissue_gef': str(tissue_gef_path),
        'mask': str(mask_path),
        'method': method,
    }

    # Try stereopy
    try:
        import stereo as st
        logger.info(f"Generating cellbin GEF via stereopy cell_correct "
                     f"(method={method})...")

        st.tools.cell_correct(
            bgef_path=str(tissue_gef_path),
            mask_path=str(mask_path),
            out_dir=str(output_dir),
            method=method,
            distance=expand_distance,
        )

        # Find output files
        adjusted = list(output_dir.glob('*.adjusted.cellbin.gef'))
        raw_cb = list(output_dir.glob('*.raw.cellbin.gef'))

        if adjusted:
            result['adjusted_cellbin_gef'] = str(adjusted[0])
            logger.info(f"Adjusted cellbin GEF: {adjusted[0]}")
        if raw_cb:
            result['raw_cellbin_gef'] = str(raw_cb[0])
            logger.info(f"Raw cellbin GEF: {raw_cb[0]}")

        result['status'] = 'success'
        result['backend'] = 'stereopy'
        return result

    except ImportError:
        logger.info("stereopy not available, trying geftools CLI")
    except Exception as e:
        logger.warning(f"stereopy cell_correct failed: {e}, trying geftools CLI")

    # Fallback: geftools CLI
    import subprocess

    chip_id = tissue_gef_path.stem.split('.')[0]
    output_gef = output_dir / f'{chip_id}.comet.cellbin.gef'

    # Search for geftools binary
    geftools_paths = [
        Path(__file__).parent.parent / 'stomics' / 'saw-8.2.2' / 'lib' / 'geftools' / 'bin' / 'geftools',
        Path(shutil.which('geftools') or ''),
    ]

    geftools_bin = None
    for p in geftools_paths:
        if p.exists():
            geftools_bin = str(p)
            break

    if geftools_bin is None:
        logger.warning("geftools binary not found. Path A (SAW-native) unavailable. "
                        "Use Path B (Python direct) instead.")
        result['status'] = 'failed_no_geftools'
        return result

    try:
        cmd = [
            geftools_bin, 'cellCut', 'cgef',
            '-i', str(tissue_gef_path),
            '-m', str(mask_path),
            '-o', str(output_gef),
        ]
        logger.info(f"Running geftools: {' '.join(cmd)}")
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=600)

        if proc.returncode == 0:
            result['cellbin_gef'] = str(output_gef)
            result['status'] = 'success'
            result['backend'] = 'geftools_cli'
            logger.info(f"Generated cellbin GEF: {output_gef}")
        else:
            logger.error(f"geftools failed: {proc.stderr}")
            result['status'] = 'failed_geftools'
            result['stderr'] = proc.stderr

    except subprocess.TimeoutExpired:
        logger.error("geftools timed out after 600s")
        result['status'] = 'failed_timeout'

    return result


# =============================================================================
# Direct Transcript-to-Cell Assignment (Path B: Python)
# =============================================================================

def aggregate_transcripts_by_mask(tissue_gef_path, mask, geojson_data=None,
                                   phenotype_map=None, roi_mask=None):
    """Assign raw bin1 transcripts to COMET cells via mask lookup.

    For each DNB position in the tissue GEF, looks up the cell assignment
    from the rasterized mask. Aggregates gene expression per cell into a
    sparse matrix, preserving original COMET cell labels.

    Args:
        tissue_gef_path: Path to bin1 tissue.gef
        mask: uint32 mask array (height x width) or path to TIFF
        geojson_data: Optional warped GeoJSON for centroid/area metadata
        phenotype_map: Optional dict {cell_id: phenotype_str} from
            rasterize_geojson_to_mask. Populates obs['comet_phenotype'].
        roi_mask: Optional uint8 ROI mask (same shape as cell mask).
            Populates obs['roi_class'] and captures extracellular transcripts.

    Returns:
        tuple: (adata, extracellular_df)
            adata: anndata.AnnData with COMET cells x genes
            extracellular_df: pd.DataFrame with per-ROI extracellular gene
                counts, or None if roi_mask not provided
    """
    import h5py
    import anndata as ad
    from scipy.sparse import coo_matrix, csr_matrix

    tissue_gef_path = Path(tissue_gef_path)

    # Load mask if path provided
    if isinstance(mask, (str, Path)):
        import tifffile
        mask = tifffile.imread(str(mask))
    mask = np.asarray(mask, dtype=np.uint32)
    mask_h, mask_w = mask.shape

    logger.info(f"Loading bin1 GEF: {tissue_gef_path.name}")
    logger.info(f"Mask shape: {mask_h}x{mask_w}")

    with h5py.File(str(tissue_gef_path), 'r') as f:
        bin1 = f['geneExp']['bin1']
        gene_data = bin1['gene'][:]
        expression = bin1['expression']

        # Extract gene names and offset/count info
        gene_names = [g['gene'].decode('utf-8') if isinstance(g['gene'], bytes)
                      else str(g['gene']) for g in gene_data]
        offsets = gene_data['offset']
        counts = gene_data['count']
        n_genes = len(gene_names)

        logger.info(f"Found {n_genes} genes, {expression.shape[0]} total DNB records")

        # ROI class names for extracellular capture
        roi_class_names = {1: 'Tumor ROI', 2: 'TME ROI', 3: 'Vessels ROI',
                           4: 'Other ROI'}
        collect_extracellular = roi_mask is not None

        # Process in gene-by-gene chunks for memory efficiency
        # Collect COO triplets: (cell_id, gene_idx, count)
        row_ids = []  # cell IDs
        col_ids = []  # gene indices
        vals = []     # counts

        # Extracellular accumulators: {roi_class_int: {gene_idx: total_count}}
        if collect_extracellular:
            extra_counts = {k: np.zeros(n_genes, dtype=np.int64)
                           for k in roi_class_names}

        total_assigned = 0
        total_outside = 0
        total_extracellular = 0

        for gene_idx in range(n_genes):
            if gene_idx % 5000 == 0 and gene_idx > 0:
                logger.info(f"  Processing gene {gene_idx}/{n_genes}...")

            offset = int(offsets[gene_idx])
            count = int(counts[gene_idx])
            if count == 0:
                continue

            # Read expression records for this gene
            records = expression[offset:offset + count]
            x = records['x']
            y = records['y']
            cnts = records['count']

            # Clip to mask bounds
            valid = (x >= 0) & (x < mask_w) & (y >= 0) & (y < mask_h)
            x_valid = x[valid]
            y_valid = y[valid]
            cnts_valid = cnts[valid]

            # Vectorized mask lookup
            cell_ids = mask[y_valid, x_valid]

            # Keep only assigned (non-background) DNBs
            assigned = cell_ids > 0
            if assigned.any():
                row_ids.append(cell_ids[assigned])
                col_ids.append(np.full(assigned.sum(), gene_idx, dtype=np.int32))
                vals.append(cnts_valid[assigned].astype(np.float32))
                total_assigned += assigned.sum()

            # Extracellular: not in any cell but within an ROI
            if collect_extracellular:
                not_assigned = ~assigned
                if not_assigned.any():
                    roi_ids = roi_mask[y_valid[not_assigned],
                                       x_valid[not_assigned]]
                    for roi_val in roi_class_names:
                        in_roi = roi_ids == roi_val
                        if in_roi.any():
                            extra_counts[roi_val][gene_idx] += int(
                                cnts_valid[not_assigned][in_roi].sum())
                            total_extracellular += int(in_roi.sum())

            total_outside += (~assigned).sum() + (~valid).sum()

    logger.info(f"Assigned {total_assigned} DNBs to cells, "
                f"{total_outside} outside cells/bounds")
    if collect_extracellular:
        logger.info(f"Captured {total_extracellular} extracellular DNBs "
                    f"within ROI regions")

    # Concatenate all triplets
    if not row_ids:
        logger.warning("No DNBs were assigned to any cell!")
        empty = ad.AnnData()
        return (empty, None) if collect_extracellular else (empty, None)

    all_cell_ids = np.concatenate(row_ids)
    all_gene_ids = np.concatenate(col_ids)
    all_counts = np.concatenate(vals)

    # Build cell_id -> row_index mapping (preserve original COMET labels)
    unique_cells = np.unique(all_cell_ids)
    cell_id_to_idx = {int(cid): idx for idx, cid in enumerate(unique_cells)}
    n_cells = len(unique_cells)

    logger.info(f"Building sparse matrix: {n_cells} cells x {n_genes} genes")

    row_indices = np.array([cell_id_to_idx[int(cid)] for cid in all_cell_ids])

    # Sum duplicate entries (same cell, same gene) via COO -> CSR
    X = coo_matrix(
        (all_counts, (row_indices, all_gene_ids)),
        shape=(n_cells, n_genes)
    ).tocsr()

    # Build obs metadata
    obs_data = {
        'cell_label': unique_cells.astype(int),
        'n_transcripts': np.array(X.sum(axis=1)).ravel().astype(int),
        'n_genes': np.array((X > 0).sum(axis=1)).ravel().astype(int),
    }

    # Add centroid and area from GeoJSON if available
    # Supports both "masks" format (properties.label) and "label" format
    # (properties.object_index)
    if geojson_data is not None:
        label_to_props = {}
        for feature in geojson_data['features']:
            props = feature['properties']
            label = props.get('label')
            if label is None:
                obj_idx = props.get('object_index')
                if obj_idx is not None:
                    label = obj_idx + 1
            if label is not None:
                coords = np.array(feature['geometry']['coordinates'][0])
                cx = np.mean(coords[:-1, 0])
                cy = np.mean(coords[:-1, 1])
                area = props.get('area_px', 0)
                label_to_props[int(label)] = (cx, cy, area)

        centroids_x = []
        centroids_y = []
        areas = []
        for cid in unique_cells:
            lprops = label_to_props.get(int(cid), (np.nan, np.nan, 0))
            centroids_x.append(lprops[0])
            centroids_y.append(lprops[1])
            areas.append(lprops[2])

        obs_data['centroid_x'] = centroids_x
        obs_data['centroid_y'] = centroids_y
        obs_data['area_px'] = areas

    # Add COMET phenotype from phenotype_map
    if phenotype_map:
        obs_data['comet_phenotype'] = [
            phenotype_map.get(int(cid), 'Unknown')
            for cid in unique_cells
        ]
        from collections import Counter
        pheno_dist = Counter(obs_data['comet_phenotype'])
        logger.info(f"Phenotype assignment: {len(pheno_dist)} classes, "
                    f"top 5: {dict(pheno_dist.most_common(5))}")

    # Assign cells to ROI regions based on centroid position
    if roi_mask is not None and 'centroid_x' in obs_data:
        roi_class_names_full = {0: 'Outside ROI', 1: 'Tumor ROI',
                                2: 'TME ROI', 3: 'Vessels ROI',
                                4: 'Other ROI'}
        roi_labels = []
        roi_h, roi_w = roi_mask.shape
        for cx, cy in zip(obs_data['centroid_x'], obs_data['centroid_y']):
            if np.isnan(cx) or np.isnan(cy):
                roi_labels.append('Unknown')
            else:
                rx, ry = int(round(cx)), int(round(cy))
                if 0 <= rx < roi_w and 0 <= ry < roi_h:
                    roi_val = int(roi_mask[ry, rx])
                    roi_labels.append(roi_class_names_full.get(roi_val,
                                                               'Outside ROI'))
                else:
                    roi_labels.append('Outside ROI')
        obs_data['roi_class'] = roi_labels
        roi_dist = Counter(roi_labels)
        logger.info(f"ROI assignment: {dict(roi_dist)}")

    obs = pd.DataFrame(obs_data, index=[str(i) for i in range(n_cells)])
    var = pd.DataFrame(index=gene_names)

    adata = ad.AnnData(X=X, obs=obs, var=var)

    # Add spatial coordinates
    if 'centroid_x' in obs.columns:
        adata.obsm['spatial'] = np.column_stack([
            obs['centroid_x'].values.astype(np.float64),
            obs['centroid_y'].values.astype(np.float64)
        ])

    adata.uns['method'] = 'direct_mask_lookup'
    adata.uns['total_dnbs_assigned'] = int(total_assigned)
    adata.uns['total_dnbs_outside'] = int(total_outside)

    logger.info(f"Created AnnData: {adata.shape}")
    logger.info(f"  Median transcripts/cell: {np.median(obs_data['n_transcripts']):.0f}")
    logger.info(f"  Median genes/cell: {np.median(obs_data['n_genes']):.0f}")

    # Build extracellular DataFrame
    extracellular_df = None
    if collect_extracellular:
        extra_data = {}
        for roi_val, roi_name in roi_class_names.items():
            total = int(extra_counts[roi_val].sum())
            if total > 0:
                extra_data[roi_name] = extra_counts[roi_val]
        if extra_data:
            extracellular_df = pd.DataFrame(extra_data, index=gene_names)
            logger.info(f"Extracellular transcripts by ROI:")
            for col in extracellular_df.columns:
                total = int(extracellular_df[col].sum())
                n_genes_detected = int((extracellular_df[col] > 0).sum())
                logger.info(f"  {col}: {total:,} transcripts, "
                           f"{n_genes_detected} genes")
        adata.uns['total_extracellular'] = int(total_extracellular)

    return adata, extracellular_df


# Microbial taxa prefixes used in STOmics microbiome GEF
MICROBIAL_PREFIXES = ('g__', 's__', 'f__', 'o__', 'c__', 'p__')


def integrate_microbiome(adata, microbiome_gef_path, cell_mask, roi_mask=None):
    """Assign microbial taxa from microbiome GEF to cells in an existing AnnData.

    Runs aggregate_transcripts_by_mask on the microbiome GEF using the same
    cell mask, then separates microbial taxa from human genes and stores the
    microbial per-cell matrix in adata.obsm['microbiome'].

    Args:
        adata: AnnData from aggregate_transcripts_by_mask (human transcripts)
        microbiome_gef_path: Path to {chip}.host_micro.label.gef
        cell_mask: uint32 cell mask (same as used for human transcripts)
        roi_mask: Optional uint8 ROI mask for extracellular microbiome capture

    Returns:
        tuple: (adata, micro_extracellular_df)
            adata: Input AnnData with .obsm['microbiome'] added and
                obs['n_microbial_reads'] populated
            micro_extracellular_df: Extracellular microbiome DataFrame or None
    """
    microbiome_gef_path = Path(microbiome_gef_path)
    if not microbiome_gef_path.exists():
        logger.warning(f"Microbiome GEF not found: {microbiome_gef_path}")
        return adata, None

    logger.info(f"Processing microbiome GEF: {microbiome_gef_path.name}")

    micro_adata, micro_extra = aggregate_transcripts_by_mask(
        microbiome_gef_path, cell_mask, roi_mask=roi_mask
    )

    if micro_adata.n_obs == 0:
        logger.warning("No cells received microbiome reads")
        return adata, micro_extra

    # Identify microbial taxa by prefix
    gene_names = list(micro_adata.var_names)
    is_microbial = np.array([any(g.startswith(p) for p in MICROBIAL_PREFIXES)
                             for g in gene_names])
    micro_idx = np.where(is_microbial)[0]
    n_micro = len(micro_idx)
    n_human = len(gene_names) - n_micro
    logger.info(f"Microbiome GEF: {n_human} human genes, {n_micro} microbial taxa")

    if n_micro == 0:
        logger.warning("No microbial taxa found in microbiome GEF")
        return adata, micro_extra

    # Extract microbial matrix for cells that match the main adata
    micro_taxa_names = [gene_names[i] for i in micro_idx]

    # Build cell_label → row_index mapping for micro_adata
    micro_labels = micro_adata.obs['cell_label'].values
    micro_label_to_idx = {int(lbl): idx for idx, lbl in enumerate(micro_labels)}

    # For each cell in adata, look up its microbial counts
    main_labels = adata.obs['cell_label'].values
    n_main_cells = len(main_labels)
    micro_matrix = np.zeros((n_main_cells, n_micro), dtype=np.float32)
    n_matched = 0

    micro_X = micro_adata.X[:, micro_idx]
    for i, lbl in enumerate(main_labels):
        micro_row = micro_label_to_idx.get(int(lbl))
        if micro_row is not None:
            row_data = micro_X[micro_row, :].toarray().ravel() \
                if hasattr(micro_X, 'toarray') else micro_X[micro_row, :]
            micro_matrix[i, :] = row_data
            n_matched += 1

    adata.obsm['microbiome'] = micro_matrix
    adata.uns['microbiome_taxa'] = micro_taxa_names

    # Summary stats
    total_per_cell = micro_matrix.sum(axis=1)
    adata.obs['n_microbial_reads'] = total_per_cell.astype(int)
    cells_with_microbes = (total_per_cell > 0).sum()

    logger.info(f"Microbiome integration: {n_matched}/{n_main_cells} cells matched")
    logger.info(f"  Cells with microbial reads: {cells_with_microbes}")
    logger.info(f"  Total microbial reads: {int(total_per_cell.sum())}")
    if cells_with_microbes > 0:
        positive = total_per_cell[total_per_cell > 0]
        logger.info(f"  Median reads (positive cells): {np.median(positive):.0f}")

    # Filter extracellular to microbial-only
    micro_extracellular_df = None
    if micro_extra is not None:
        micro_extra_taxa = micro_extra.loc[micro_extra.index.isin(micro_taxa_names)]
        if micro_extra_taxa.sum().sum() > 0:
            micro_extracellular_df = micro_extra_taxa
            logger.info(f"Extracellular microbiome by ROI:")
            for col in micro_extracellular_df.columns:
                total = int(micro_extracellular_df[col].sum())
                logger.info(f"  {col}: {total:,} microbial reads")

    return adata, micro_extracellular_df


# =============================================================================
# COMET Phenotype Annotation
# =============================================================================

# Phenotype classification tree (CD45+ hierarchy from COMET READ_ME)
PHENOTYPE_MARKERS = {
    'Treg': {'positive': ['CD45', 'CD4', 'FOXp3'], 'negative': []},
    'CD4+ T cell': {'positive': ['CD45', 'CD4'], 'negative': ['FOXp3']},
    'CD8+ T cell': {'positive': ['CD45', 'CD8'], 'negative': []},
    'Macrophage': {'positive': ['CD45', 'CD68'], 'negative': []},
    'M2 Macrophage': {'positive': ['CD45', 'CD163'], 'negative': []},
    'B cell': {'positive': ['CD45', 'CD20'], 'negative': []},
    'NK cell': {'positive': ['CD45', 'CD56'], 'negative': []},
    'DC': {'positive': ['CD45', 'CD11C'], 'negative': []},
    'Neutrophil': {'positive': ['CD45', 'CD66B'], 'negative': []},
    'Immune (other)': {'positive': ['CD45'], 'negative': []},
    'Epithelial': {'positive': ['KERATIN8'], 'negative': ['CD45']},
    'Fibroblast': {'positive': ['COL1A1'], 'negative': ['CD45']},
    'Endothelial': {'positive': ['CD31'], 'negative': ['CD45']},
}


def annotate_with_comet_phenotypes(adata, comet_protein_path, warped_geojson=None,
                                    match_method='label', match_threshold=30,
                                    intensity_threshold_pctile=75):
    """Merge COMET protein intensities and phenotype classifications into AnnData.

    For Path B (direct mask lookup): joins on cell_label (exact match).
    For Path A (SAW-native): joins via centroid spatial KD-tree.

    Args:
        adata: AnnData from aggregate_transcripts_by_mask or loaded cellbin GEF
        comet_protein_path: Path to {SO}_BS_protein_combined.parquet
        warped_geojson: Warped GeoJSON (needed for Path A centroid mapping)
        match_method: 'label' (Path B, exact join) or 'spatial' (Path A, KD-tree)
        match_threshold: Max distance in pixels for spatial matching (Path A)
        intensity_threshold_pctile: Percentile threshold for positive marker call

    Returns:
        anndata.AnnData with phenotype annotations in .obs
    """
    comet_protein_path = Path(comet_protein_path)
    protein_df = pd.read_parquet(str(comet_protein_path))

    logger.info(f"Loaded COMET protein data: {protein_df.shape}")
    logger.info(f"Protein columns: {list(protein_df.columns)}")

    # Identify protein channels (exclude coordinate columns)
    coord_cols = {'centroid_x', 'centroid_y', 'x', 'y', 'cell_id', 'label'}
    protein_cols = [c for c in protein_df.columns if c not in coord_cols]

    if match_method == 'label':
        # Path B: direct label match
        if 'cell_label' not in adata.obs.columns:
            raise ValueError("AnnData must have 'cell_label' in .obs for label matching")

        # Build mapping: protein_df index = cell label (1-based from regionprops)
        protein_df = protein_df.copy()
        protein_df['_label'] = protein_df.index + 1  # regionprops labels are 1-based

        adata_labels = adata.obs['cell_label'].astype(int).values
        label_to_protein_idx = {int(lab): idx for idx, lab in enumerate(protein_df['_label'])}

        matched_indices = []
        matched_adata_indices = []
        for i, label in enumerate(adata_labels):
            pidx = label_to_protein_idx.get(label)
            if pidx is not None:
                matched_indices.append(pidx)
                matched_adata_indices.append(i)

        n_matched = len(matched_indices)
        logger.info(f"Label matching: {n_matched}/{len(adata_labels)} cells matched")

    elif match_method == 'spatial':
        # Path A: KD-tree centroid matching
        if warped_geojson is None:
            raise ValueError("warped_geojson required for spatial matching")

        # Get COMET centroids from warped GeoJSON
        comet_centroids, comet_labels = geojson_centroids(warped_geojson)

        # Get cellbin centroids from AnnData
        if 'spatial' in adata.obsm:
            cellbin_centroids = adata.obsm['spatial']
        elif 'x' in adata.obs.columns and 'y' in adata.obs.columns:
            cellbin_centroids = np.column_stack([
                adata.obs['x'].values, adata.obs['y'].values
            ])
        else:
            raise ValueError("AnnData must have spatial coordinates")

        # Match cellbin cells to COMET cells
        tree = cKDTree(comet_centroids)
        distances, indices = tree.query(cellbin_centroids, k=1)

        within = distances <= match_threshold
        matched_adata_indices = np.where(within)[0].tolist()
        matched_comet_indices = indices[within].tolist()

        # Map from COMET GeoJSON index to protein_df index (same ordering)
        matched_indices = matched_comet_indices

        n_matched = len(matched_indices)
        logger.info(f"Spatial matching: {n_matched}/{len(cellbin_centroids)} cells "
                     f"matched within {match_threshold}px")
        logger.info(f"Match distance: median={np.median(distances[within]):.1f}px")

    else:
        raise ValueError(f"Unknown match_method: {match_method}")

    # Add protein intensities to .obs
    for col in protein_cols:
        adata.obs[f'protein_{col}'] = np.nan
        if n_matched > 0:
            values = protein_df[col].values[matched_indices]
            adata.obs.loc[adata.obs.index[matched_adata_indices], f'protein_{col}'] = values

    # Store protein matrix in .obsm
    protein_matrix = np.full((len(adata), len(protein_cols)), np.nan)
    if n_matched > 0:
        protein_matrix[matched_adata_indices] = protein_df[protein_cols].values[matched_indices]
    adata.obsm['comet_protein'] = protein_matrix

    # Classify phenotypes
    _classify_phenotypes(adata, protein_cols, intensity_threshold_pctile)

    return adata


def _classify_phenotypes(adata, protein_cols, threshold_pctile=75):
    """Assign phenotype labels based on COMET protein marker intensities.

    Uses a hierarchical decision tree: most specific match wins.
    """
    # Compute per-marker thresholds from non-NaN values
    thresholds = {}
    for col in protein_cols:
        vals = adata.obs[f'protein_{col}'].dropna()
        if len(vals) > 0:
            thresholds[col] = np.percentile(vals, threshold_pctile)

    phenotypes = ['Unclassified'] * len(adata)

    # Apply classification tree (order matters: most specific first)
    classification_order = [
        'Treg', 'CD4+ T cell', 'CD8+ T cell', 'M2 Macrophage',
        'Macrophage', 'B cell', 'NK cell', 'DC', 'Neutrophil',
        'Immune (other)', 'Epithelial', 'Fibroblast', 'Endothelial',
    ]

    for pheno in classification_order:
        if pheno not in PHENOTYPE_MARKERS:
            continue
        rules = PHENOTYPE_MARKERS[pheno]

        # Check if all required markers have thresholds
        pos_markers = rules['positive']
        neg_markers = rules['negative']

        pos_available = all(m in thresholds for m in pos_markers)
        neg_available = all(m in thresholds for m in neg_markers)

        if not pos_available:
            continue

        for i in range(len(adata)):
            if phenotypes[i] != 'Unclassified':
                continue

            # Check positive markers
            all_pos = True
            for marker in pos_markers:
                val = adata.obs.iloc[i].get(f'protein_{marker}', np.nan)
                if pd.isna(val) or val < thresholds[marker]:
                    all_pos = False
                    break

            if not all_pos:
                continue

            # Check negative markers
            all_neg_clear = True
            if neg_available:
                for marker in neg_markers:
                    val = adata.obs.iloc[i].get(f'protein_{marker}', np.nan)
                    if not pd.isna(val) and val >= thresholds[marker]:
                        all_neg_clear = False
                        break

            if all_neg_clear:
                phenotypes[i] = pheno

    adata.obs['phenotype'] = phenotypes

    # Log distribution
    counts = pd.Series(phenotypes).value_counts()
    logger.info(f"Phenotype classification:")
    for pheno, count in counts.items():
        logger.info(f"  {pheno}: {count} ({count/len(adata)*100:.1f}%)")


# =============================================================================
# Full Pipeline Runner
# =============================================================================

def run_alignment_pipeline(sample_id, chip_id, work_dir,
                            geojson_path=None, comet_bs_path=None,
                            stomics_dapi_path=None, stomics_data_path=None,
                            skip_registration=False, registrar=None,
                            aggregation_method='nearest', max_distance=50,
                            save_registered_slides=False,
                            method='valis',
                            qupath_transforms_path=None,
                            ecc_warp_mode='euclidean',
                            ecc_downsample=4):
    """Run the complete COMET -> STOmics alignment pipeline for one sample.

    Steps:
        1. Check data availability
        2. Load COMET GeoJSON segmentations
        3. Run registration (COMET DAPI -> STOmics DAPI)
        4. Warp segmentations to STOmics space
        5. Load STOmics gene expression
        6. Aggregate expression per COMET cell -> AnnData
        7. Map COMET <-> STOmics cells
        8. Validate alignment

    Args:
        sample_id: e.g., 'SO34'
        chip_id: e.g., 'A03979E2'
        work_dir: Working directory for outputs
        geojson_path: Path to COMET GeoJSON (auto-resolved if None)
        comet_bs_path: Path to COMET BS OME-TIFF
        stomics_dapi_path: Path to STOmics DAPI
        stomics_data_path: Path to STOmics cellbin GEF or h5ad
        skip_registration: Skip registration (use pre-computed registrar)
        registrar: Pre-computed VALIS registrar
        aggregation_method: 'nearest' (fast) or 'polygon' (precise)
        max_distance: Max distance for nearest-neighbor methods (pixels)
        save_registered_slides: Save warped slide images
        method: Registration method:
            'valis' — VALIS rigid + non-rigid (default)
            'qupath_affine' — QuPath manual affine only
            'qupath_ecc' — QuPath affine + OpenCV ECC refinement
            'qupath_valis' — QuPath affine pre-warp + VALIS refinement
        qupath_transforms_path: Path to config/qupath_transforms.json
            (required for qupath_affine and qupath_ecc methods)
        ecc_warp_mode: ECC warp mode: 'euclidean' (default) or 'affine'
        ecc_downsample: Downscale factor for ECC on large images (default 4)

    Returns:
        dict: Pipeline results (AnnData, mapping, metrics, paths)
    """
    work_dir = Path(work_dir)
    sample_dir = work_dir / f'{sample_id}_{chip_id}'
    sample_dir.mkdir(parents=True, exist_ok=True)

    results = {'sample_id': sample_id, 'chip_id': chip_id, 'status': 'started'}

    # --- Resolve paths ---
    info = ALL_SAMPLES.get(sample_id, ALL_ENDO_MLA.get(sample_id, OVARIAN_WITH_STOMICS.get(sample_id, {})))
    aligned = info.get('aligned', False)

    if geojson_path is None:
        geojson_path = DefaultPaths.geojson_path(sample_id, aligned)
    geojson_path = Path(geojson_path)

    if comet_bs_path is None:
        comet_bs_path = DefaultPaths.comet_bs_path(sample_id)
    comet_bs_path = Path(comet_bs_path)

    if stomics_dapi_path is None:
        stomics_dapi_path = DefaultPaths.stomics_dapi_path(chip_id)
    stomics_dapi_path = Path(stomics_dapi_path)

    if stomics_data_path is None:
        stomics_data_path = DefaultPaths.stomics_cellbin_path(chip_id)
    stomics_data_path = Path(stomics_data_path)

    # Check files
    logger.info(f"\n{'='*60}")
    logger.info(f"Pipeline: {sample_id} ({chip_id})")
    logger.info(f"{'='*60}")

    file_status = {
        'geojson': geojson_path.exists(),
        'comet_bs': comet_bs_path.exists(),
        'stomics_dapi': stomics_dapi_path.exists(),
        'stomics_data': stomics_data_path.exists(),
    }

    for name, exists in file_status.items():
        logger.info(f"  {name}: {'OK' if exists else 'MISSING'}")

    results['file_status'] = file_status

    # --- Step 1: Load GeoJSON ---
    if not file_status['geojson']:
        logger.error(f"GeoJSON not found: {geojson_path}")
        results['status'] = 'failed_no_geojson'
        return results

    logger.info("\nStep 1: Loading COMET segmentations...")
    geojson_data = load_geojson(geojson_path)
    results['n_comet_cells'] = len(geojson_data['features'])

    # --- Step 2 & 3: Registration + Warp ---
    valid_methods = ('valis', 'qupath_affine', 'qupath_ecc', 'qupath_valis')
    if method not in valid_methods:
        raise ValueError(f"Unknown method '{method}'. Use one of: {valid_methods}")

    warped_geojson = None
    warped_path = sample_dir / f'{sample_id}_warped_segmentations.geojson'

    if method == 'valis':
        # --- VALIS registration (existing behavior) ---
        if not skip_registration and registrar is None:
            if not file_status['stomics_dapi']:
                logger.error("STOmics DAPI not found - cannot register")
                results['status'] = 'failed_no_stomics_dapi'
                return results

            logger.info("\nStep 2: VALIS Registration...")
            src_dir, ref_img = setup_registration_inputs(
                sample_id, chip_id, sample_dir,
                comet_bs_path=comet_bs_path,
                stomics_dapi_path=stomics_dapi_path
            )

            reg_dir = sample_dir / 'registration_output'
            registrar = run_valis_registration(src_dir, reg_dir, ref_img)
            results['registrar'] = registrar

            if save_registered_slides:
                save_registrar(registrar, sample_dir)
        elif registrar is not None:
            results['registrar'] = registrar

        if registrar is not None:
            logger.info("\nStep 3: Warping segmentations (VALIS)...")
            warped_geojson = warp_geojson_with_valis(registrar, geojson_path, warped_path)
            results['warped_geojson'] = warped_geojson
            results['warped_geojson_path'] = warped_path
        else:
            logger.warning("No registrar - skipping warp")
            results['status'] = 'failed_no_registrar'
            return results

    elif method in ('qupath_affine', 'qupath_ecc'):
        # --- QuPath-based registration ---
        from comet.registration import (
            register_qupath_affine, register_qupath_ecc,
            apply_affine_to_geojson, plot_registration_overlay
        )

        if qupath_transforms_path is None:
            raise ValueError(
                f"qupath_transforms_path is required for method='{method}'"
            )

        if method == 'qupath_affine':
            logger.info("\nStep 2: QuPath Affine Registration...")
            forward_matrix = register_qupath_affine(
                sample_id, qupath_transforms_path, modality='comet'
            )
            results['registration'] = {
                'method': 'qupath_affine',
                'affine_matrix': forward_matrix.tolist(),
            }

        else:  # qupath_ecc
            logger.info("\nStep 2: QuPath + ECC Registration...")
            if not file_status['stomics_dapi']:
                logger.error("STOmics DAPI not found - cannot run ECC")
                results['status'] = 'failed_no_stomics_dapi'
                return results

            # Extract COMET DAPI for ECC comparison
            comet_dapi_path = sample_dir / f'{sample_id}_dapi.tif'
            if not comet_dapi_path.exists() and file_status['comet_bs']:
                extract_comet_dapi(comet_bs_path, comet_dapi_path)

            forward_matrix, ecc_score, ecc_converged = register_qupath_ecc(
                sample_id, qupath_transforms_path,
                source_img_path=str(comet_dapi_path),
                target_img_path=str(stomics_dapi_path),
                modality='comet',
                warp_mode=ecc_warp_mode,
                downsample=ecc_downsample
            )
            results['registration'] = {
                'method': 'qupath_ecc',
                'affine_matrix': forward_matrix.tolist(),
                'ecc_score': ecc_score,
                'ecc_converged': ecc_converged,
                'ecc_warp_mode': ecc_warp_mode,
            }

            # Save verification overlay
            overlay_path = sample_dir / f'{sample_id}_registration_overlay.png'
            try:
                import cv2
                comet_dapi = cv2.imread(str(comet_dapi_path), cv2.IMREAD_GRAYSCALE)
                stomics_dapi = cv2.imread(str(stomics_dapi_path), cv2.IMREAD_GRAYSCALE)
                if comet_dapi is not None and stomics_dapi is not None:
                    plot_registration_overlay(
                        comet_dapi, stomics_dapi, forward_matrix,
                        overlay_path,
                        title=f'{sample_id} QuPath+ECC (score={ecc_score:.4f})'
                    )
            except Exception as e:
                logger.warning(f"Could not generate overlay: {e}")

        # Step 3: Warp GeoJSON with affine
        logger.info("\nStep 3: Warping segmentations (affine)...")
        warped_geojson = apply_affine_to_geojson(
            geojson_data, forward_matrix, output_path=warped_path
        )
        results['warped_geojson'] = warped_geojson
        results['warped_geojson_path'] = warped_path

    elif method == 'qupath_valis':
        # --- QuPath pre-warp + VALIS refinement ---
        from comet.registration import (
            register_qupath_valis, apply_affine_to_geojson,
            register_qupath_affine
        )

        if qupath_transforms_path is None:
            raise ValueError(
                "qupath_transforms_path is required for method='qupath_valis'"
            )
        if not file_status['stomics_dapi']:
            logger.error("STOmics DAPI not found - cannot run VALIS")
            results['status'] = 'failed_no_stomics_dapi'
            return results

        logger.info("\nStep 2: QuPath + VALIS Registration...")

        # Extract COMET DAPI if needed
        comet_dapi_path = sample_dir / f'{sample_id}_dapi.tif'
        if not comet_dapi_path.exists() and file_status['comet_bs']:
            extract_comet_dapi(comet_bs_path, comet_dapi_path)

        valis_work_dir = sample_dir / 'qupath_valis'
        registrar, qupath_matrix = register_qupath_valis(
            sample_id, qupath_transforms_path,
            source_img_path=str(comet_dapi_path),
            target_img_path=str(stomics_dapi_path),
            work_dir=str(valis_work_dir),
            modality='comet'
        )
        results['registrar'] = registrar
        results['registration'] = {
            'method': 'qupath_valis',
            'qupath_affine_matrix': qupath_matrix.tolist(),
        }

        # Step 3: Warp GeoJSON — first apply QuPath affine, then VALIS non-rigid
        logger.info("\nStep 3: Warping segmentations (QuPath affine + VALIS)...")

        # Apply QuPath affine to GeoJSON first
        affine_warped_path = sample_dir / f'{sample_id}_affine_warped.geojson'
        affine_warped_geojson = apply_affine_to_geojson(
            geojson_data, qupath_matrix, output_path=affine_warped_path
        )

        # Then apply VALIS non-rigid warp
        warped_geojson = warp_geojson_with_valis(
            registrar, affine_warped_path, warped_path
        )
        results['warped_geojson'] = warped_geojson
        results['warped_geojson_path'] = warped_path

    # --- Step 4: Load STOmics data ---
    stomics_adata = None
    if file_status['stomics_data']:
        logger.info("\nStep 4: Loading STOmics data...")
        suffix = stomics_data_path.suffix
        if suffix == '.gef':
            stomics_adata = load_stomics_cellbin_gef(stomics_data_path)
        elif suffix == '.h5ad':
            stomics_adata = load_stomics_h5ad(stomics_data_path)
        results['stomics_adata'] = stomics_adata
    else:
        logger.warning("STOmics data not found - skipping expression aggregation")

    # --- Step 5: Aggregate expression ---
    if stomics_adata is not None:
        logger.info("\nStep 5: Aggregating expression per COMET cell...")
        integrated_adata = aggregate_expression_per_comet_cell(
            warped_geojson, stomics_adata,
            method=aggregation_method, max_distance=max_distance
        )

        h5ad_path = sample_dir / f'{sample_id}_{chip_id}_comet_stomics_integrated.h5ad'
        integrated_adata.uns['sample_id'] = sample_id
        integrated_adata.uns['chip_id'] = chip_id
        integrated_adata.uns['disease'] = info.get('disease', 'Unknown')
        integrated_adata.uns['aggregation_method'] = aggregation_method
        integrated_adata.write_h5ad(str(h5ad_path))
        results['integrated_adata'] = integrated_adata
        results['integrated_h5ad_path'] = h5ad_path
        logger.info(f"Saved: {h5ad_path}")

    # --- Step 6: Cell mapping ---
    if stomics_adata is not None:
        logger.info("\nStep 6: Mapping COMET <-> STOmics cells...")
        comet_centroids, comet_labels = geojson_centroids(warped_geojson)

        if 'spatial' in stomics_adata.obsm:
            st_centroids = stomics_adata.obsm['spatial']
        else:
            st_centroids = np.column_stack([
                stomics_adata.obs['x'], stomics_adata.obs['y']
            ])

        cell_mapping = map_comet_to_stomics_cells(
            comet_centroids, st_centroids, max_distance
        )
        cell_mapping['comet_label'] = comet_labels

        mapping_path = sample_dir / f'{sample_id}_{chip_id}_cell_mapping.csv'
        cell_mapping.to_csv(mapping_path, index=False)
        results['cell_mapping'] = cell_mapping
        results['cell_mapping_path'] = mapping_path

    # --- Step 7: Validation ---
    logger.info("\nStep 7: Validation...")
    if stomics_adata is not None:
        comet_centroids, _ = geojson_centroids(warped_geojson)

        if 'spatial' in stomics_adata.obsm:
            st_centroids = stomics_adata.obsm['spatial']
        else:
            st_centroids = np.column_stack([
                stomics_adata.obs['x'], stomics_adata.obs['y']
            ])

        metrics = compute_alignment_metrics(comet_centroids, st_centroids)
        metrics['registration_method'] = method
        if 'registration' in results:
            metrics['registration_details'] = results['registration']
        results['alignment_metrics'] = metrics

        # Save metrics
        metrics_path = sample_dir / f'{sample_id}_{chip_id}_alignment_metrics.json'
        with open(metrics_path, 'w') as f:
            json.dump(metrics, f, indent=2)

        logger.info(f"\nAlignment Quality: {metrics['quality']}")
        logger.info(f"  Median distance: {metrics['median_distance']:.1f}px")
        logger.info(f"  Within 30px: {metrics['pct_within_30px']:.1f}%")
        logger.info(f"  Within 50px: {metrics['pct_within_50px']:.1f}%")

    # --- Create Validation Manifest ---
    logger.info("\nCreating validation manifest for signoff...")
    try:
        from comet.alignment_validation import AlignmentValidator

        files_used = {
            'geojson_source': str(geojson_path),
            'comet_bs': str(comet_bs_path),
            'stomics_dapi': str(stomics_dapi_path),
            'stomics_data': str(stomics_data_path),
        }
        if 'warped_geojson_path' in results:
            files_used['warped_geojson'] = str(results['warped_geojson_path'])
        if 'integrated_h5ad_path' in results:
            files_used['integrated_h5ad'] = str(results['integrated_h5ad_path'])
        if 'cell_mapping_path' in results:
            files_used['cell_mapping'] = str(results['cell_mapping_path'])
        reg_dir = sample_dir / 'registration_output'
        if reg_dir.exists():
            files_used['registration_dir'] = str(reg_dir)

        validator = AlignmentValidator(sample_dir)
        validator.create_manifest(
            files_used=files_used,
            registration_method=method if method else 'valis_qupath_affine',
            extra_metadata={
                'aggregation_method': aggregation_method,
                'max_distance': max_distance,
                'alignment_metrics': results.get('alignment_metrics', {}),
            },
        )
        results['validator'] = validator
        results['manifest_path'] = sample_dir / 'validation_manifest.json'
        logger.info(f"  Manifest created: {results['manifest_path']}")
        logger.info(f"  Status: PENDING — run QC review and sign_off() before downstream analysis")
    except ImportError:
        logger.warning("alignment_validation module not available — skipping manifest creation")
    except Exception as e:
        logger.warning(f"Could not create validation manifest: {e}")

    results['status'] = 'complete'
    logger.info(f"\n{'='*60}")
    logger.info(f"Pipeline complete for {sample_id}")
    logger.info(f"{'='*60}\n")

    return results
