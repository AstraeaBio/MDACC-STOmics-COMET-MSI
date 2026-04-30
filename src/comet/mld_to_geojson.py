"""
Convert Visiopharm MLD files to GeoJSON format for COMET-STOmics-MSI project.

Exports both tissue-level (ROI) and cell-level (Label) layers from MLD files,
transforming from Visiopharm mm coordinates to COMET top-left pixel coordinates.
Output GeoJSONs are QuPath-compatible for review and alignment work.

Adapted from the Maitra project solution:
    T:/0_Organizational/Git/maitra_sow124/co_registration/mld_to_geojson.py

Usage:
    # Single sample
    python comet/mld_to_geojson.py --sample SO1_BS

    # All samples
    python comet/mld_to_geojson.py --all

    # With verification overlays
    python comet/mld_to_geojson.py --sample SO1_BS --verify

    # Export rasterized masks
    python comet/mld_to_geojson.py --sample SO1_BS --masks

    # Both ROI + Label layers (default is ROI only; Label is large)
    python comet/mld_to_geojson.py --sample SO1_BS --layers ROI Label
"""

import os
import sys
import io
import json
import uuid
import argparse
from datetime import datetime
from collections import Counter
import numpy as np
import cv2
import tifffile
import xml.etree.ElementTree as ET
from shapely.geometry import shape, mapping

# Import MLDReader from the local module
_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
from importlib.machinery import SourceFileLoader as _SFL
_mld_module = _SFL("mld_reader", os.path.join(_SCRIPT_DIR, "mld-reader-module.py")).load_module()
MLDReader = _mld_module.MLDReader
ShapeType = _mld_module.ShapeType


# ──────────────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────────────

MLD_DIR = r"T:/Sammy Data/Vis_Analysis_All_Samples"
OUT_DIR = os.path.expanduser("~/ext_hd_sammy/projects/out/out_comet_mld")

# Label maps from Visiopharm TSV export:
#   T:/Sammy Data/Vis_Analysis/Vis_Analysis_All_Samples/40SamplesProject_Analysis_objects.tsv
#   Columns: "Object Info - LabelType/LabelName" and "Object Info - ROIType/ROIName"

# ROI layer: tissue-level regions (object_type -> (name, [R,G,B]))
LABELS_ROI = {
    0:  ("Background",      [180, 180, 180]),
    3:  ("Tumor ROI",       [255, 0, 0]),
    4:  ("TME ROI",         [0, 200, 0]),
    5:  ("Vessels ROI",     [0, 100, 255]),
    6:  ("Other ROI",       [255, 165, 0]),
}

# Label layer: cell-level phenotypes (object_type -> (name, [R,G,B]))
LABELS_CELL = {
    0:  ("Background",                  [180, 180, 180]),
    2:  ("Cell",                        [200, 200, 200]),
    5:  ("IMMUNE CELLS",                [255, 255, 0]),
    6:  ("CD4 T CELLS",                 [0, 200, 0]),
    7:  ("T REG CELLS",                 [0, 128, 0]),
    8:  ("MEMORY HELPER T CELLS",       [100, 200, 100]),
    9:  ("CD8 T CELLS",                 [0, 0, 255]),
    10: ("ACTIVATED CD8 T CELLS",       [50, 50, 255]),
    11: ("EXHAUSTED CD8 T CELLS",       [100, 100, 255]),
    12: ("MEMORY CD8 T CELLS",          [150, 150, 255]),
    13: ("B CELLS",                     [255, 128, 0]),
    14: ("DENDRITIC CELLS",             [200, 0, 255]),
    15: ("DENDRITIC CELLS CD45RO+",     [170, 0, 220]),
    16: ("MDSC",                        [128, 0, 200]),
    17: ("MACROPHAGES",                 [255, 0, 255]),
    18: ("M1 MACROPHAGES",              [255, 50, 200]),
    19: ("M2 MACROPHAGES",              [200, 50, 255]),
    20: ("TUMOR CELLS CD56+",           [255, 0, 0]),
    22: ("NEUTROPHILS",                 [255, 200, 0]),
    23: ("IMMUNOSUPPRESSIVE NEUTROPHILS", [200, 160, 0]),
    24: ("NK CELLS",                    [0, 255, 128]),
    25: ("ACTIVATED NK CELLS",          [0, 200, 100]),
    26: ("NKT CELLS",                   [0, 255, 200]),
    28: ("BLOOD VESSELS",               [128, 0, 0]),
    29: ("BLOOD VESSELS aSMA+",         [180, 0, 0]),
    30: ("TUMOR CELLS",                 [220, 50, 50]),
    31: ("CAFs",                        [160, 82, 45]),
    32: ("FIBROBLASTS",                 [210, 180, 140]),
    33: ("CAFs COL1A1+",                [139, 69, 19]),
    35: ("FIBROBLASTS COL1A1+ CD10+",   [188, 143, 143]),
    36: ("CAFs COL1A1+ CD10+",          [160, 120, 90]),
    38: ("FIBROBLASTS PERIOSTIN+",      [222, 184, 135]),
    39: ("CAFs PERIOSTIN+",             [205, 133, 63]),
}

# Per-layer label maps
LABELS_SAMMY = {
    "ROI": LABELS_ROI,
    "Label": LABELS_CELL,
}


# ──────────────────────────────────────────────────────────────────────────────
# OME metadata extraction
# ──────────────────────────────────────────────────────────────────────────────

def get_comet_image_params(tif_path):
    """
    Parse OME metadata from a COMET TIF to get image dimensions and pixel size.

    Returns dict with:
        size_x, size_y: image dimensions in pixels
        physical_size_x_um: pixel size in micrometers
        pixels_per_mm: conversion factor from mm to pixels
        position_x_um, position_y_um: stage position in micrometers (if available)
    """
    with tifffile.TiffFile(tif_path) as tf:
        if not tf.ome_metadata:
            raise ValueError(f"No OME metadata in {tif_path}")

        root = ET.fromstring(tf.ome_metadata)
        ns_uri = root.tag.split("}")[0].strip("{") if "}" in root.tag else ""
        ns = {"ome": ns_uri} if ns_uri else {}

        pixels = root.find(".//ome:Pixels", ns) if ns else root.find(".//Pixels")
        if pixels is None:
            raise ValueError(f"No Pixels element in OME metadata of {tif_path}")

        size_x = int(pixels.get("SizeX"))
        size_y = int(pixels.get("SizeY"))
        phys_x = float(pixels.get("PhysicalSizeX", "0.23"))

        params = {
            "size_x": size_x,
            "size_y": size_y,
            "physical_size_x_um": phys_x,
            "pixels_per_mm": 1000.0 / phys_x,
        }

        plane = root.find(".//ome:Plane", ns) if ns else root.find(".//Plane")
        if plane is not None:
            pos_x = plane.get("PositionX")
            pos_y = plane.get("PositionY")
            if pos_x is not None:
                params["position_x_um"] = float(pos_x)
            if pos_y is not None:
                params["position_y_um"] = float(pos_y)

        return params


# ──────────────────────────────────────────────────────────────────────────────
# Coordinate transform
# ──────────────────────────────────────────────────────────────────────────────

def make_mm_to_px_converter(pixels_per_mm, image_width_px, image_height_px,
                            stage_x_mm=None, stage_y_mm=None):
    """
    Create coordinate converter: MLD mm (center-origin, Y-up) -> top-left px.

    With stage position:
        x_px = (x_mm - stage_x_mm) * ppm
        y_px = (stage_y_mm - y_mm) * ppm

    Without stage position (center-origin fallback):
        x_px = x_mm * ppm + width/2
        y_px = height/2 - y_mm * ppm
    """
    if stage_x_mm is not None and stage_y_mm is not None:
        def convert(point):
            x_mm, y_mm = point
            x_px = (x_mm - stage_x_mm) * pixels_per_mm
            y_px = (stage_y_mm - y_mm) * pixels_per_mm
            return [x_px, y_px]
    else:
        def convert(point):
            x_mm, y_mm = point
            x_px = x_mm * pixels_per_mm + (image_width_px / 2.0)
            y_px = (image_height_px / 2.0) - y_mm * pixels_per_mm
            return [x_px, y_px]
    return convert


# ──────────────────────────────────────────────────────────────────────────────
# Shape conversion
# ──────────────────────────────────────────────────────────────────────────────

def _convert_shape_to_geometry(obj, converter):
    """Convert an MLD shape object to a GeoJSON geometry."""
    shape_type = obj["shape"]

    if shape_type in [ShapeType.POLYGON, ShapeType.POLYGON2]:
        points = obj["shape_data"]["points"]
        if points and points[0] != points[-1]:
            points = points + [points[0]]
        coordinates = [[converter(p) for p in points]]
        return {"type": "Polygon", "coordinates": coordinates}

    elif shape_type == ShapeType.RECTANGLE:
        sd = obj["shape_data"]
        if sd.get("width", 0) > 1e10 or sd.get("height", 0) > 1e10:
            return None  # bounding box sentinel
        ox, oy = sd["origin"]
        hw, hh = sd["width"] / 2, sd["height"] / 2
        corners = [
            converter([ox - hw, oy - hh]), converter([ox + hw, oy - hh]),
            converter([ox + hw, oy + hh]), converter([ox - hw, oy + hh]),
            converter([ox - hw, oy - hh]),
        ]
        return {"type": "Polygon", "coordinates": [corners]}

    elif shape_type == ShapeType.POLYLINE:
        points = obj["shape_data"]["points"]
        return {"type": "LineString", "coordinates": [converter(p) for p in points]}

    elif shape_type == ShapeType.CIRCLE:
        return {"type": "Point", "coordinates": converter(obj["shape_data"]["origin"])}

    elif shape_type == ShapeType.TEXT_POINT:
        return {"type": "Point", "coordinates": converter(obj["shape_data"]["origin"])}

    return None


def _simplify_ring(ring, tolerance):
    """Simplify a polygon ring using Douglas-Peucker algorithm."""
    if len(ring) < 5 or tolerance <= 0:
        return ring
    pts = np.array(ring, dtype=np.float32)
    pts_cv = pts.reshape(-1, 1, 2)
    simplified = cv2.approxPolyDP(pts_cv, tolerance, closed=True)
    coords = simplified.reshape(-1, 2)
    result = [[float(c[0]), float(c[1])] for c in coords]
    # Ensure ring is closed
    if result[0] != result[-1]:
        result.append(result[0])
    # Need at least 4 points for a valid polygon ring
    if len(result) < 4:
        return ring
    return result


# ──────────────────────────────────────────────────────────────────────────────
# MLD reading
# ──────────────────────────────────────────────────────────────────────────────

def read_mld_quiet(mld_path):
    """Read an MLD file, suppressing MLDReader's print output."""
    old_stdout = sys.stdout
    old_argv = sys.argv
    sys.stdout = io.StringIO()
    sys.argv = ["mld_reader"]
    try:
        reader = MLDReader(mld_path)
        data = reader.read()
    finally:
        sys.stdout = old_stdout
        sys.argv = old_argv
    return data


# ──────────────────────────────────────────────────────────────────────────────
# GeoJSON export
# ──────────────────────────────────────────────────────────────────────────────

def mld_to_geojson(mld_path, tif_path, output_path, layers=None,
                   label_map=None, max_objects=None, simplify_tolerance=0,
                   detection_layers=None, exclude_types=None):
    """
    Convert MLD layers to GeoJSON in COMET pixel coordinates.

    Args:
        mld_path: Path to the .mld file
        tif_path: Path to the paired .tif (for OME metadata)
        output_path: Path for the output GeoJSON (pixel coords)
        layers: List of layer names to export (default: ["ROI"])
        label_map: Dict mapping object_type -> (name, [R,G,B])
        max_objects: Optional limit on objects per layer (for testing)
        simplify_tolerance: Douglas-Peucker tolerance in pixels (0 = no simplification).
                           1.5 is a good default for reducing file size while preserving shape.
        detection_layers: Set of layer names to export as QuPath "detection" objects
                         (default: None -> all layers exported as "annotation")
        exclude_types: Dict mapping layer name -> set of object_type ints to skip
                      (e.g., {"ROI": {0}} to skip Background)

    Returns:
        The GeoJSON dict
    """
    if layers is None:
        layers = ["ROI"]
    if detection_layers is None:
        detection_layers = set()
    if exclude_types is None:
        exclude_types = {}

    mld_data = read_mld_quiet(mld_path)

    params = get_comet_image_params(tif_path)
    ppm = params["pixels_per_mm"]
    w, h = params["size_x"], params["size_y"]

    has_stage = "position_x_um" in params and "position_y_um" in params
    if has_stage:
        stage_x_mm = params["position_x_um"] / 1000.0
        stage_y_mm = params["position_y_um"] / 1000.0
        print(f"  Image: {w}x{h} px, {params['physical_size_x_um']:.4f} um/px, ppm={ppm:.2f}")
        print(f"  Stage position: ({stage_x_mm:.4f}, {stage_y_mm:.4f}) mm")
        converter = make_mm_to_px_converter(ppm, w, h, stage_x_mm, stage_y_mm)
    else:
        print(f"  Image: {w}x{h} px, {params['physical_size_x_um']:.4f} um/px, ppm={ppm:.2f}")
        print(f"  No stage position -> center-origin transform")
        converter = make_mm_to_px_converter(ppm, w, h)

    if simplify_tolerance > 0:
        print(f"  Simplification: tolerance={simplify_tolerance} px")

    features = []
    skipped_bbox = 0
    skipped_excluded = 0
    skipped_unfixable = 0
    total_fixed = 0
    fixed_indices = []       # (layer_name, obj_idx, object_type)
    unfixable_indices = []   # (layer_name, obj_idx, object_type)
    unmapped_types = set()
    total_pts_before = 0
    total_pts_after = 0
    type_counts = Counter()  # (layer_name, object_type) -> count
    layer_stats = {}         # layer_name -> {total, exported, excluded, ...}

    for layer in mld_data["layers"]:
        lname = layer["name"]
        if lname not in layers:
            continue

        n_objects = len(layer["objects"])
        limit = min(n_objects, max_objects) if max_objects else n_objects
        print(f"  Layer '{lname}': {n_objects} objects" +
              (f" (exporting first {limit})" if limit < n_objects else ""))

        layer_excludes = exclude_types.get(lname, set())
        layer_exported = 0
        layer_excluded = 0

        for obj_idx, obj in enumerate(layer["objects"][:limit]):
            otype = int(obj["type"])
            if otype in layer_excludes:
                skipped_excluded += 1
                layer_excluded += 1
                continue
            geom = _convert_shape_to_geometry(obj, converter)
            if geom is None:
                skipped_bbox += 1
                continue

            # Simplify polygon rings if requested
            if geom["type"] == "Polygon":
                for ring in geom["coordinates"]:
                    total_pts_before += len(ring)
                if simplify_tolerance > 0:
                    geom["coordinates"] = [
                        _simplify_ring(ring, simplify_tolerance)
                        for ring in geom["coordinates"]
                    ]
                for ring in geom["coordinates"]:
                    total_pts_after += len(ring)

            # Fix invalid geometries (self-intersections) via Shapely
            if geom["type"] == "Polygon":
                shp = shape(geom)
                if not shp.is_valid:
                    shp = shp.buffer(0)
                    if shp.is_valid:
                        geom = json.loads(json.dumps(mapping(shp)))
                        total_fixed += 1
                        fixed_indices.append((lname, obj_idx, otype))
                    else:
                        skipped_unfixable += 1
                        unfixable_indices.append((lname, obj_idx, otype))
                        continue

            is_detection = lname in detection_layers

            props = {
                "objectType": "detection" if is_detection else "annotation",
                "object_type": otype,
                "object_index": obj_idx,
            }

            # Resolve label map: support per-layer dict or flat dict
            lm = label_map
            if isinstance(label_map, dict) and lname in label_map and isinstance(label_map[lname], dict):
                lm = label_map[lname]

            if lm and otype in lm:
                name, color = lm[otype]
                props["classification"] = {
                    "name": name,
                    "color": color,
                }
            else:
                unmapped_types.add((lname, otype))

            type_counts[(lname, otype)] += 1
            layer_exported += 1

            features.append({
                "type": "Feature",
                "id": f"{lname}_{obj_idx}",
                "geometry": geom,
                "properties": props,
            })

        layer_stats[lname] = {
            "total_objects": n_objects,
            "limit_applied": limit if limit < n_objects else None,
            "exported": layer_exported,
            "excluded_by_type": layer_excluded,
        }

    geojson = {
        "type": "FeatureCollection",
        "features": features,
        "metadata": {
            "source_mld": os.path.basename(mld_path),
            "source_tif": os.path.basename(tif_path),
            "coordinate_system": "COMET_pixel",
            "image_width_px": w,
            "image_height_px": h,
            "pixels_per_mm": ppm,
            "physical_size_x_um": params["physical_size_x_um"],
            "has_stage_position": has_stage,
            "layers_exported": layers,
        },
    }

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(geojson, f, indent=2)

    total_skipped = skipped_bbox + skipped_excluded + skipped_unfixable
    print(f"  -> {len(features)} features saved to {output_path}")
    if total_skipped:
        parts = []
        if skipped_bbox:
            parts.append(f"{skipped_bbox} bbox/unsupported")
        if skipped_excluded:
            parts.append(f"{skipped_excluded} excluded by type")
        if skipped_unfixable:
            parts.append(f"{skipped_unfixable} unfixable geometry")
        print(f"  Skipped {total_skipped} objects ({', '.join(parts)})")
    if total_fixed:
        print(f"  Fixed {total_fixed} invalid geometries (self-intersections)")
    if total_pts_before > 0 and simplify_tolerance > 0:
        pct = (1 - total_pts_after / total_pts_before) * 100
        print(f"  Simplified: {total_pts_before:,} -> {total_pts_after:,} points ({pct:.1f}% reduction)")

    # Build export log
    output_size = os.path.getsize(output_path)
    export_log = {
        "timestamp": datetime.now().isoformat(),
        "input": {
            "mld_path": os.path.abspath(mld_path),
            "tif_path": os.path.abspath(tif_path),
        },
        "output": {
            "geojson_path": os.path.abspath(output_path),
            "file_size_mb": round(output_size / 1e6, 2),
        },
        "image_params": {
            "width_px": w,
            "height_px": h,
            "physical_size_um": params["physical_size_x_um"],
            "pixels_per_mm": round(ppm, 2),
            "stage_position_mm": [round(stage_x_mm, 4), round(stage_y_mm, 4)] if has_stage else None,
            "coordinate_transform": "stage_position" if has_stage else "center_origin",
        },
        "export_settings": {
            "layers": layers,
            "simplify_tolerance_px": simplify_tolerance,
            "detection_layers": list(detection_layers),
            "exclude_types": {k: list(v) for k, v in exclude_types.items()} if exclude_types else {},
            "max_objects": max_objects,
        },
        "results": {
            "total_features_exported": len(features),
            "skipped_bbox_unsupported": skipped_bbox,
            "skipped_excluded_types": skipped_excluded,
            "skipped_unfixable_geometry": skipped_unfixable,
            "fixed_self_intersections": total_fixed,
        },
        "layer_details": {},
        "type_distribution": {},
    }

    # Per-layer stats
    for lname, stats in layer_stats.items():
        export_log["layer_details"][lname] = stats

    # Type distribution with names
    for (lname, otype), cnt in sorted(type_counts.items()):
        lm = label_map
        if isinstance(label_map, dict) and lname in label_map and isinstance(label_map[lname], dict):
            lm = label_map[lname]
        type_name = lm[otype][0] if (lm and otype in lm) else f"type_{otype}"
        export_log["type_distribution"][f"{lname}:{otype}:{type_name}"] = cnt

    # Fixed geometry details
    if fixed_indices:
        export_log["fixed_geometries"] = [
            {"layer": l, "object_index": i, "object_type": t}
            for l, i, t in fixed_indices
        ]

    # Unfixable geometry details
    if unfixable_indices:
        export_log["unfixable_geometries"] = [
            {"layer": l, "object_index": i, "object_type": t}
            for l, i, t in unfixable_indices
        ]

    # Unmapped types
    if unmapped_types:
        export_log["unmapped_types"] = [
            {"layer": l, "object_type": t} for l, t in sorted(unmapped_types)
        ]

    # Simplification stats
    if total_pts_before > 0 and simplify_tolerance > 0:
        export_log["simplification"] = {
            "points_before": total_pts_before,
            "points_after": total_pts_after,
            "reduction_pct": round((1 - total_pts_after / total_pts_before) * 100, 1),
        }

    # Write log file alongside GeoJSON
    log_path = output_path.replace(".geojson", "_export_log.json")
    with open(log_path, "w") as f:
        json.dump(export_log, f, indent=2)
    print(f"  Log -> {log_path}")

    return geojson


# ──────────────────────────────────────────────────────────────────────────────
# Verification overlay
# ──────────────────────────────────────────────────────────────────────────────

def verify_overlay(tif_path, geojson_data, output_png, layer_filter=None):
    """
    Plot GeoJSON polygons overlaid on the TIF label mask for visual QC.

    Args:
        tif_path: Path to the COMET .tif (uses layer 1 = ROI mask as background)
        geojson_data: GeoJSON dict from mld_to_geojson()
        output_png: Output path for the verification PNG
        layer_filter: Optional layer name to filter features (e.g., "ROI")
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon as MplPolygon, Patch

    img = tifffile.imread(tif_path)
    if img.ndim == 3 and img.shape[0] >= 2:
        tissue = img[1]  # ROI mask layer
    else:
        tissue = img if img.ndim == 2 else img[0]

    # Downsample for display
    target_w = 2000
    factor = max(1, tissue.shape[1] // target_w)
    tissue_ds = tissue[::factor, ::factor]

    fig, ax = plt.subplots(1, 1, figsize=(14, 14))
    ax.imshow(tissue_ds, cmap="gray", aspect="equal",
              extent=[0, tissue.shape[1], tissue.shape[0], 0])

    features = geojson_data["features"]
    if layer_filter:
        features = [f for f in features if f["properties"].get("layer") == layer_filter]

    cmap = plt.cm.tab20
    all_types = sorted(set(f["properties"]["object_type"] for f in features))
    type_to_color = {t: cmap(i / max(len(all_types), 1)) for i, t in enumerate(all_types)}

    for feature in features:
        geom = feature["geometry"]
        otype = feature["properties"]["object_type"]
        color = type_to_color[otype]

        if geom["type"] == "Polygon":
            for ring in geom["coordinates"]:
                coords = np.array(ring)
                poly = MplPolygon(coords, closed=True, facecolor=color,
                                  edgecolor="white", linewidth=0.3, alpha=0.5)
                ax.add_patch(poly)

    sample = os.path.basename(tif_path).replace(".tif", "")
    layer_str = f" [{layer_filter}]" if layer_filter else ""
    ax.set_title(f"{sample}{layer_str}: {len(features)} features", fontsize=13)
    ax.set_xlabel("X (pixels)")
    ax.set_ylabel("Y (pixels)")

    legend_elements = [Patch(facecolor=type_to_color[t], edgecolor="gray",
                             label=f"Type {t}") for t in all_types]
    if len(legend_elements) <= 20:
        ax.legend(handles=legend_elements, loc="upper right", fontsize=8)

    fig.tight_layout()
    os.makedirs(os.path.dirname(output_png), exist_ok=True)
    fig.savefig(output_png, dpi=150, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    print(f"  Verification plot -> {output_png}")


# ──────────────────────────────────────────────────────────────────────────────
# Rasterization: MLD polygons -> label masks
# ──────────────────────────────────────────────────────────────────────────────

def rasterize_mld_to_mask(mld_path, tif_path, output_path, layers=None,
                          label_map=None):
    """
    Rasterize MLD polygons into a uint8 label mask at image resolution.

    Pixel values correspond to object_type values.
    Background types rendered first, specific types on top.

    Args:
        mld_path: Path to the .mld file
        tif_path: Path to the paired .tif (for OME metadata)
        output_path: Path for the output TIFF mask
        layers: Layer names to rasterize (default: ["ROI"])
        label_map: Optional label map (determines render order)

    Returns:
        dict with mask_shape, unique_values, n_polygons, output_path
    """
    if layers is None:
        layers = ["ROI"]

    mld_data = read_mld_quiet(mld_path)
    params = get_comet_image_params(tif_path)
    ppm = params["pixels_per_mm"]
    w, h = params["size_x"], params["size_y"]

    has_stage = "position_x_um" in params and "position_y_um" in params
    if has_stage:
        stage_x_mm = params["position_x_um"] / 1000.0
        stage_y_mm = params["position_y_um"] / 1000.0
        converter = make_mm_to_px_converter(ppm, w, h, stage_x_mm, stage_y_mm)
    else:
        converter = make_mm_to_px_converter(ppm, w, h)

    # Collect polygons grouped by object_type
    polys_by_type = {}
    n_polys = 0
    for layer in mld_data["layers"]:
        if layer["name"] not in layers:
            continue
        for obj in layer["objects"]:
            shape_type = obj["shape"]
            if shape_type not in [ShapeType.POLYGON, ShapeType.POLYGON2]:
                # Skip non-polygon shapes for rasterization
                if shape_type == ShapeType.RECTANGLE:
                    sd = obj["shape_data"]
                    if sd.get("width", 0) > 1e10 or sd.get("height", 0) > 1e10:
                        continue
                    ox, oy = sd["origin"]
                    hw, hh = sd["width"] / 2, sd["height"] / 2
                    corners = [
                        [ox - hw, oy - hh], [ox + hw, oy - hh],
                        [ox + hw, oy + hh], [ox - hw, oy + hh],
                    ]
                    pts = np.array([converter(c) for c in corners], dtype=np.int32)
                else:
                    continue
            else:
                points = obj["shape_data"]["points"]
                if not points:
                    continue
                if points[0] != points[-1]:
                    points = points + [points[0]]
                pts = np.array([converter(p) for p in points], dtype=np.float64)
                pts = pts[:, :2].astype(np.int32)

            otype = int(obj["type"])
            polys_by_type.setdefault(otype, []).append(pts)
            n_polys += 1

    # Render: lower type numbers first, higher on top
    mask = np.zeros((h, w), dtype=np.uint8)
    for otype in sorted(polys_by_type.keys()):
        for pts in polys_by_type[otype]:
            cv2.fillPoly(mask, [pts], color=int(otype))

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    tifffile.imwrite(output_path, mask, compression="lzw")

    uvals = np.unique(mask)
    print(f"  Mask: {h}x{w}, {n_polys} polygons, unique={uvals.tolist()} -> {output_path}")

    return {
        "mask_shape": (h, w),
        "unique_values": uvals.tolist(),
        "n_polygons": n_polys,
        "output_path": output_path,
    }


# ──────────────────────────────────────────────────────────────────────────────
# Sample discovery
# ──────────────────────────────────────────────────────────────────────────────

def discover_samples(mld_dir=MLD_DIR):
    """
    Discover MLD + TIF pairs in the data directory.

    Returns list of dicts with keys: sample_id, mld_path, tif_path
    """
    import glob
    samples = []
    for mld_path in sorted(glob.glob(os.path.join(mld_dir, "*.mld"))):
        stem = os.path.splitext(os.path.basename(mld_path))[0]
        tif_path = os.path.join(mld_dir, f"{stem}.tif")
        if os.path.exists(tif_path):
            samples.append({
                "sample_id": stem,
                "mld_path": mld_path,
                "tif_path": tif_path,
            })
        else:
            print(f"  WARNING: no TIF for {stem}")
    return samples


def process_sample(sample, out_dir, layers, label_map, verify=False,
                   masks=False, max_objects=None, simplify_tolerance=0,
                   detection_layers=None, exclude_types=None):
    """Process a single sample: export GeoJSON and optionally masks/overlays."""
    sid = sample["sample_id"]
    safe_name = sid.replace(" ", "_")
    print(f"\n{'='*60}")
    print(f"Sample: {sid}")
    print(f"{'='*60}")

    # Build layer-specific output names
    layer_tag = "_".join(layers).lower()

    # GeoJSON export
    geojson_path = os.path.join(out_dir, "geojson", f"{safe_name}_{layer_tag}.geojson")
    geojson_data = mld_to_geojson(
        sample["mld_path"], sample["tif_path"], geojson_path,
        layers=layers, label_map=label_map, max_objects=max_objects,
        simplify_tolerance=simplify_tolerance,
        detection_layers=detection_layers,
        exclude_types=exclude_types,
    )

    # Verification overlay
    if verify:
        for lyr in layers:
            png_path = os.path.join(out_dir, "verify", f"{safe_name}_{lyr.lower()}_verify.png")
            verify_overlay(sample["tif_path"], geojson_data, png_path, layer_filter=lyr)

    # Rasterized masks
    if masks:
        for lyr in layers:
            mask_path = os.path.join(out_dir, "masks", f"{safe_name}_{lyr.lower()}_mask.tif")
            rasterize_mld_to_mask(
                sample["mld_path"], sample["tif_path"], mask_path,
                layers=[lyr], label_map=label_map,
            )

    return geojson_data


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Export Visiopharm MLD files to GeoJSON (COMET pixel coordinates)"
    )
    parser.add_argument("--sample", type=str, default=None,
                        help="Sample ID to process (e.g., SO1_BS). Default: none.")
    parser.add_argument("--all", action="store_true",
                        help="Process all samples in MLD_DIR")
    parser.add_argument("--layers", nargs="+", default=["ROI"],
                        help="Layer names to export (default: ROI). Use 'ROI Label' for both.")
    parser.add_argument("--verify", action="store_true",
                        help="Generate verification overlay PNGs")
    parser.add_argument("--masks", action="store_true",
                        help="Generate rasterized label mask TIFFs")
    parser.add_argument("--simplify", type=float, default=1.5,
                        help="Douglas-Peucker simplification tolerance in pixels (default: 1.5, 0=off)")
    parser.add_argument("--detection-layers", nargs="+", default=[],
                        help="Layers to export as QuPath 'detection' objects (default: none)")
    parser.add_argument("--exclude-types", type=str, default=None,
                        help="Exclude object types per layer, e.g. 'ROI:0' or 'ROI:0,Label:2'")
    parser.add_argument("--max-objects", type=int, default=None,
                        help="Limit objects per layer (for testing)")
    parser.add_argument("--mld-dir", type=str, default=MLD_DIR,
                        help=f"MLD directory (default: {MLD_DIR})")
    parser.add_argument("--out-dir", type=str, default=OUT_DIR,
                        help=f"Output directory (default: {OUT_DIR})")
    args = parser.parse_args()

    samples = discover_samples(args.mld_dir)
    print(f"Found {len(samples)} MLD+TIF pairs in {args.mld_dir}")

    if args.sample:
        matches = [s for s in samples if s["sample_id"] == args.sample]
        if not matches:
            # Try partial match
            matches = [s for s in samples if args.sample in s["sample_id"]]
        if not matches:
            print(f"ERROR: sample '{args.sample}' not found. Available:")
            for s in samples:
                print(f"  {s['sample_id']}")
            sys.exit(1)
        targets = matches
    elif args.all:
        targets = samples
    else:
        parser.print_help()
        print("\nSpecify --sample <ID> or --all")
        sys.exit(1)

    # Parse --exclude-types (e.g. "ROI:0" or "ROI:0,Label:2")
    exclude_types = {}
    if args.exclude_types:
        for part in args.exclude_types.split(","):
            layer_name, type_id = part.strip().split(":")
            exclude_types.setdefault(layer_name.strip(), set()).add(int(type_id.strip()))
        print(f"Excluding types: {exclude_types}")

    total_ok, total_fail = 0, 0
    for sample in targets:
        try:
            process_sample(
                sample, args.out_dir, args.layers, LABELS_SAMMY,
                verify=args.verify, masks=args.masks,
                max_objects=args.max_objects,
                simplify_tolerance=args.simplify,
                detection_layers=set(args.detection_layers),
                exclude_types=exclude_types,
            )
            total_ok += 1
        except Exception as e:
            print(f"  ERROR: {e}")
            import traceback
            traceback.print_exc()
            total_fail += 1

    print(f"\n{'='*60}")
    print(f"Done: {total_ok} succeeded, {total_fail} failed")
    print(f"Output: {args.out_dir}")


if __name__ == "__main__":
    main()
