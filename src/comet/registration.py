"""
Multi-method image registration for spatial multi-omics alignment.

Methods:
    - QuPath affine: Manual 2x3 affine from QuPath interactive alignment
    - QuPath + ECC: QuPath affine followed by OpenCV ECC refinement
    - VALIS: Handled separately in alignment_utils.py (unchanged)

Reference coordinate space: STOmics DAPI (where gene expression coordinates live).

Adapted from Chiba project ECC registration:
    T:/Analysis/116_MDACC_Chiba/chiba/out_feb2025_imc_registration/
    script04_ecc_refined/register_imc_gbm_ecc_refined.py
"""

import json
import logging
from copy import deepcopy
from pathlib import Path

import cv2
import numpy as np

logger = logging.getLogger(__name__)


# ── Config I/O ──────────────────────────────────────────────────────────────


def load_qupath_transforms(json_path, sample_id=None):
    """Load QuPath affine transforms from JSON config file.

    Automatically inverts matrices where invert=true (default when omitted).

    Args:
        json_path: Path to qupath_transforms.json
        sample_id: If provided, return only this sample's transforms.
                   Otherwise return all samples.

    Returns:
        dict: Sample transforms with matrices as np.ndarray (2x3).
              If sample_id given, returns that sample's dict directly.
    """
    with open(json_path, 'r') as f:
        config = json.load(f)

    # Remove comment keys
    config = {k: v for k, v in config.items() if not k.startswith('_')}

    # Process each sample's transforms
    for sid, sample_cfg in config.items():
        for modality, tform in sample_cfg.get('transforms', {}).items():
            matrix = np.array(tform['matrix'], dtype=np.float64)
            if matrix.shape != (2, 3):
                raise ValueError(
                    f"Matrix for {sid}/{modality} has shape {matrix.shape}, "
                    f"expected (2, 3)"
                )

            # Invert by default (QuPath gives reference→source direction)
            should_invert = tform.get('invert', True)
            if should_invert:
                tform['matrix_raw'] = matrix.copy()
                matrix = invert_affine(matrix)

            tform['matrix'] = matrix

    if sample_id is not None:
        if sample_id not in config:
            raise KeyError(
                f"Sample '{sample_id}' not found in {json_path}. "
                f"Available: {list(config.keys())}"
            )
        return config[sample_id]

    return config


def invert_affine(matrix_2x3):
    """Invert a 2x3 affine matrix to get the forward transform.

    For M = [[a, b, tx], [c, d, ty]]:
        det = a*d - b*c
        M_inv = [[ d/det, -b/det, (b*ty - d*tx)/det],
                 [-c/det,  a/det, (c*tx - a*ty)/det]]

    Args:
        matrix_2x3: np.ndarray of shape (2, 3)

    Returns:
        np.ndarray: Inverted 2x3 affine matrix
    """
    M = np.array(matrix_2x3, dtype=np.float64)
    a, b, tx = M[0]
    c, d, ty = M[1]

    det = a * d - b * c
    if abs(det) < 1e-12:
        raise ValueError(f"Affine matrix is singular (det={det})")

    M_inv = np.array([
        [ d / det, -b / det, (b * ty - d * tx) / det],
        [-c / det,  a / det, (c * tx - a * ty) / det]
    ], dtype=np.float64)

    return M_inv


# ── Affine Operations ───────────────────────────────────────────────────────


def apply_affine_to_coordinates(coords_xy, matrix_2x3):
    """Warp (N, 2) coordinate array through a 2x3 affine transform.

    Args:
        coords_xy: np.ndarray of shape (N, 2) with [x, y] columns
        matrix_2x3: np.ndarray of shape (2, 3)

    Returns:
        np.ndarray: Warped coordinates of shape (N, 2)
    """
    M = np.array(matrix_2x3, dtype=np.float64)
    coords = np.array(coords_xy, dtype=np.float64)

    # coords_warped = coords @ M[:2,:2].T + M[:, 2]
    linear = M[:, :2]  # (2, 2)
    translation = M[:, 2]  # (2,)
    warped = coords @ linear.T + translation

    return warped


def apply_affine_to_geojson(geojson, matrix_2x3, output_path=None):
    """Transform all polygon vertices in a GeoJSON FeatureCollection.

    Handles Polygon and MultiPolygon geometries. Preserves all properties.

    Args:
        geojson: dict (GeoJSON FeatureCollection) or path to file
        matrix_2x3: np.ndarray of shape (2, 3)
        output_path: If provided, write warped GeoJSON to this path.

    Returns:
        dict: Warped GeoJSON FeatureCollection
    """
    if isinstance(geojson, (str, Path)):
        with open(geojson, 'r') as f:
            geojson = json.load(f)

    warped = deepcopy(geojson)
    M = np.array(matrix_2x3, dtype=np.float64)
    n_features = len(warped.get('features', []))

    for feature in warped.get('features', []):
        geom = feature.get('geometry', {})
        geom_type = geom.get('type', '')

        if geom_type == 'Polygon':
            geom['coordinates'] = _warp_polygon_coords(
                geom['coordinates'], M
            )
        elif geom_type == 'MultiPolygon':
            geom['coordinates'] = [
                _warp_polygon_coords(poly, M)
                for poly in geom['coordinates']
            ]

    logger.info(f"Warped {n_features} features with affine transform")

    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(warped, f, indent=2)
        logger.info(f"Saved warped GeoJSON to {output_path}")

    return warped


def _warp_polygon_coords(rings, matrix_2x3):
    """Warp polygon coordinate rings through affine transform.

    Args:
        rings: list of rings, each ring is list of [x, y] pairs
        matrix_2x3: np.ndarray (2, 3)

    Returns:
        list: Warped rings
    """
    warped_rings = []
    for ring in rings:
        coords = np.array(ring, dtype=np.float64)
        warped = apply_affine_to_coordinates(coords, matrix_2x3)
        warped_rings.append(warped.tolist())
    return warped_rings


def warp_image_affine(image, matrix_2x3, target_shape):
    """Warp a raster image using cv2.warpAffine.

    Args:
        image: np.ndarray (grayscale or color)
        matrix_2x3: np.ndarray of shape (2, 3)
        target_shape: (height, width) tuple for output size

    Returns:
        np.ndarray: Warped image
    """
    M = np.array(matrix_2x3, dtype=np.float64)
    h, w = target_shape[:2]
    warped = cv2.warpAffine(
        image, M, (w, h),
        flags=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_CONSTANT,
        borderValue=0
    )
    return warped


# ── ECC Refinement ──────────────────────────────────────────────────────────


def refine_with_ecc(source_img, target_img, init_matrix=None,
                    warp_mode='euclidean', max_iterations=200,
                    epsilon=1e-6, gaussian_ksize=5, downsample=None):
    """OpenCV Enhanced Correlation Coefficient refinement.

    Adapted from Chiba project (register_imc_gbm_ecc_refined.py).

    OpenCV convention: templateImage=target_img, inputImage=source_img.
    The returned matrix warps source pixels to align with target.

    Both images must be grayscale uint8.

    Args:
        source_img: Grayscale uint8 source (e.g., COMET DAPI after affine
                    pre-warp into STOmics space)
        target_img: Grayscale uint8 target (e.g., STOmics DAPI)
        init_matrix: Initial 2x3 warp (default: identity). Typically identity
                     because QuPath affine is already baked into the
                     pre-warped source image.
        warp_mode: 'euclidean' (rotation+translation, 3 params, default)
                   or 'affine' (full 6 params)
        max_iterations: ECC convergence max iterations
        epsilon: ECC convergence tolerance
        gaussian_ksize: Gaussian blur kernel size (applied to both images)
        downsample: Optional int downscale factor (e.g., 4). Both images
                    are downsampled before ECC and translation components
                    are scaled back up in the result.

    Returns:
        tuple: (refined_matrix_2x3, ecc_score, converged)
            If ECC fails to converge, returns (init_matrix, 0.0, False).
    """
    # Select warp mode
    if warp_mode == 'euclidean':
        motion_type = cv2.MOTION_EUCLIDEAN
    elif warp_mode == 'affine':
        motion_type = cv2.MOTION_AFFINE
    else:
        raise ValueError(f"Unknown warp_mode: {warp_mode}. "
                         f"Use 'euclidean' or 'affine'.")

    # Default init matrix
    if init_matrix is None:
        if motion_type == cv2.MOTION_EUCLIDEAN:
            init_matrix = np.eye(2, 3, dtype=np.float32)
        else:
            init_matrix = np.eye(2, 3, dtype=np.float32)

    warp_matrix = np.array(init_matrix, dtype=np.float32)

    # Ensure grayscale uint8
    src = _ensure_gray_uint8(source_img)
    tgt = _ensure_gray_uint8(target_img)

    # Optional downsampling for large images
    scale_factor = 1
    if downsample is not None and downsample > 1:
        scale_factor = downsample
        src = cv2.resize(src, None, fx=1.0/scale_factor, fy=1.0/scale_factor,
                         interpolation=cv2.INTER_AREA)
        tgt = cv2.resize(tgt, None, fx=1.0/scale_factor, fy=1.0/scale_factor,
                         interpolation=cv2.INTER_AREA)
        # Scale translation components of init matrix
        warp_matrix[0, 2] /= scale_factor
        warp_matrix[1, 2] /= scale_factor
        logger.info(f"Downsampled images by {scale_factor}x for ECC: "
                    f"{src.shape[1]}x{src.shape[0]}")

    # Apply Gaussian blur
    ksize = (gaussian_ksize, gaussian_ksize)
    src_blur = cv2.GaussianBlur(src, ksize, 0)
    tgt_blur = cv2.GaussianBlur(tgt, ksize, 0)

    # ECC convergence criteria
    criteria = (
        cv2.TERM_CRITERIA_EPS | cv2.TERM_CRITERIA_COUNT,
        max_iterations,
        epsilon
    )

    # Run ECC
    try:
        ecc_score, warp_matrix = cv2.findTransformECC(
            tgt_blur,      # templateImage (target)
            src_blur,       # inputImage (source to warp)
            warp_matrix,
            motion_type,
            criteria
        )
        converged = True

        # Extract angle for logging
        if motion_type == cv2.MOTION_EUCLIDEAN:
            angle_deg = np.degrees(np.arctan2(warp_matrix[1, 0],
                                              warp_matrix[0, 0]))
            logger.info(f"ECC converged: score={ecc_score:.4f}, "
                        f"angle={angle_deg:.2f}deg, "
                        f"tx={warp_matrix[0,2]:.1f}, "
                        f"ty={warp_matrix[1,2]:.1f}")
        else:
            logger.info(f"ECC converged: score={ecc_score:.4f}")

    except cv2.error as e:
        logger.warning(f"ECC did not converge: {e}")
        ecc_score = 0.0
        converged = False
        # warp_matrix remains as init_matrix

    # Scale translation back up if downsampled
    result_matrix = warp_matrix.astype(np.float64)
    if scale_factor > 1:
        result_matrix[0, 2] *= scale_factor
        result_matrix[1, 2] *= scale_factor

    return result_matrix, float(ecc_score), converged


def compose_affine_transforms(first, second):
    """Compose two 2x3 affine transforms: result = second(first(x)).

    Converts to 3x3 homogeneous, multiplies, returns 2x3.

    Example: For QuPath+ECC workflow:
        compose_affine_transforms(qupath_forward, ecc_correction)
    yields the single matrix that replaces both steps.

    Args:
        first: np.ndarray (2, 3) — applied first
        second: np.ndarray (2, 3) — applied second

    Returns:
        np.ndarray: Composed 2x3 affine matrix
    """
    # Convert to 3x3 homogeneous
    M1 = np.vstack([np.array(first, dtype=np.float64), [0, 0, 1]])
    M2 = np.vstack([np.array(second, dtype=np.float64), [0, 0, 1]])

    composed = M2 @ M1
    return composed[:2, :]


def _ensure_gray_uint8(img):
    """Convert image to grayscale uint8 if needed."""
    if img.ndim == 3:
        img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    if img.dtype != np.uint8:
        # Normalize to 0-255
        img_min, img_max = img.min(), img.max()
        if img_max > img_min:
            img = ((img - img_min) / (img_max - img_min) * 255).astype(np.uint8)
        else:
            img = np.zeros_like(img, dtype=np.uint8)
    return img


# ── Pipeline Helpers ────────────────────────────────────────────────────────


def register_qupath_affine(sample_id, json_path, modality='comet'):
    """Load QuPath affine transform, return ready-to-use forward matrix.

    Args:
        sample_id: e.g., 'SO34'
        json_path: Path to qupath_transforms.json
        modality: Transform key, e.g., 'comet', 'msi_glycan'

    Returns:
        np.ndarray: Forward 2x3 affine matrix (already inverted if needed)
    """
    sample_cfg = load_qupath_transforms(json_path, sample_id)
    transforms = sample_cfg.get('transforms', {})

    if modality not in transforms:
        available = list(transforms.keys())
        raise KeyError(
            f"Modality '{modality}' not found for {sample_id}. "
            f"Available: {available}"
        )

    matrix = transforms[modality]['matrix']
    logger.info(f"Loaded QuPath affine for {sample_id}/{modality}")
    return matrix


def register_qupath_ecc(sample_id, json_path, source_img_path,
                        target_img_path, modality='comet',
                        warp_mode='euclidean', downsample=4,
                        **ecc_kwargs):
    """Full QuPath + ECC registration pipeline.

    Steps:
        1. Load QuPath affine from JSON (auto-inverts)
        2. Load source and target images as grayscale uint8
        3. Warp source image with QuPath affine
        4. Run ECC refinement (warped source vs target, starting from identity)
        5. Compose QuPath affine + ECC correction
        6. Return final matrix + metrics

    If ECC fails to converge, falls back to QuPath affine only (logs warning).

    Args:
        sample_id: e.g., 'SO34'
        json_path: Path to qupath_transforms.json
        source_img_path: Path to source image (e.g., COMET DAPI)
        target_img_path: Path to target image (e.g., STOmics DAPI)
        modality: Transform key, e.g., 'comet'
        warp_mode: 'euclidean' (default) or 'affine'
        downsample: Downscale factor for ECC (default 4, recommended for
                    large images >10K pixels)
        **ecc_kwargs: Additional kwargs passed to refine_with_ecc
                      (max_iterations, epsilon, gaussian_ksize)

    Returns:
        tuple: (composed_matrix_2x3, ecc_score, converged)
    """
    # Step 1: Load QuPath affine
    qupath_matrix = register_qupath_affine(sample_id, json_path, modality)

    # Step 2: Load images — use memory-mapped reads for large TIFs that
    # fail with tifffile.imread on Windows (OSError: Invalid argument)
    import tifffile

    def _load_tif_safe(path):
        """Load a TIFF safely, handling large files that fail with imread."""
        try:
            img = tifffile.imread(str(path))
        except OSError:
            # Fall back to page-level reading with memory mapping
            logger.info(f"imread failed for {Path(path).name}, using mmap fallback")
            with tifffile.TiffFile(str(path)) as tif:
                page = tif.pages[0]
                img = page.asarray(out='memmap')
                img = np.array(img)  # copy from memmap to regular array
        if img is None:
            raise FileNotFoundError(f"Could not load image: {path}")
        # Handle multi-channel
        if img.ndim == 3:
            if img.shape[0] <= 4:  # (C, H, W) — take first channel
                img = img[0]
            else:  # (H, W, C)
                img = cv2.cvtColor(img, cv2.COLOR_RGB2GRAY)
        # Convert to uint8
        if img.dtype != np.uint8:
            p1 = np.percentile(img[img > 0], 1) if img.max() > 0 else 0
            p99 = np.percentile(img[img > 0], 99) if img.max() > 0 else 1
            img = np.clip((img.astype(np.float32) - p1) / max(p99 - p1, 1) * 255,
                          0, 255).astype(np.uint8)
        return img

    source_img = _load_tif_safe(source_img_path)
    target_img = _load_tif_safe(target_img_path)

    logger.info(f"Source image: {source_img.shape}, Target image: {target_img.shape}")

    # Step 3: Warp source with QuPath affine
    target_shape = target_img.shape[:2]
    # Use strip-based warping for images exceeding OpenCV's 32K pixel limit
    if max(source_img.shape[:2] + target_shape) > 32000:
        logger.info("Using strip-based warp for large image "
                     f"({source_img.shape} -> {target_shape})")
        source_warped = warp_image_large(source_img, qupath_matrix, target_shape)
    else:
        source_warped = warp_image_affine(source_img, qupath_matrix, target_shape)

    # Step 4: ECC refinement (starting from identity — affine already baked in)
    # Get per-transform warp_mode override if present
    sample_cfg = load_qupath_transforms(json_path, sample_id)
    tform = sample_cfg['transforms'].get(modality, {})
    warp_mode_override = tform.get('ecc_warp_mode', warp_mode)

    ecc_matrix, ecc_score, converged = refine_with_ecc(
        source_warped, target_img,
        warp_mode=warp_mode_override,
        downsample=downsample,
        **ecc_kwargs
    )

    # Step 5: Compose transforms
    if converged:
        composed = compose_affine_transforms(qupath_matrix, ecc_matrix)
        logger.info(f"QuPath+ECC registration complete for "
                    f"{sample_id}/{modality}: ECC score={ecc_score:.4f}")
    else:
        composed = qupath_matrix
        logger.warning(f"ECC did not converge for {sample_id}/{modality}. "
                       f"Falling back to QuPath affine only.")

    return composed, ecc_score, converged


def warp_image_large(image, matrix_2x3, target_shape, strip_height=2048):
    """Warp a raster image using manual inverse mapping.

    Handles images larger than OpenCV's 32K pixel limit (SHRT_MAX)
    by processing in row strips with numpy coordinate mapping.

    Args:
        image: np.ndarray (grayscale, any dtype)
        matrix_2x3: np.ndarray (2, 3) forward transform
        target_shape: (height, width) tuple
        strip_height: Rows to process at once (controls memory usage)

    Returns:
        np.ndarray: Warped image of target_shape, same dtype as input
    """
    M_fwd = np.array(matrix_2x3, dtype=np.float64)
    M_inv = invert_affine(M_fwd)

    out_h, out_w = target_shape[:2]
    src_h, src_w = image.shape[:2]
    output = np.zeros((out_h, out_w), dtype=image.dtype)

    for y0 in range(0, out_h, strip_height):
        y1 = min(y0 + strip_height, out_h)
        ys = np.arange(y0, y1)
        xs = np.arange(out_w)
        xg, yg = np.meshgrid(xs, ys)

        src_x = M_inv[0, 0] * xg + M_inv[0, 1] * yg + M_inv[0, 2]
        src_y = M_inv[1, 0] * xg + M_inv[1, 1] * yg + M_inv[1, 2]

        src_xi = np.round(src_x).astype(np.int32)
        src_yi = np.round(src_y).astype(np.int32)

        valid = (src_xi >= 0) & (src_xi < src_w) & \
                (src_yi >= 0) & (src_yi < src_h)

        strip = np.zeros((y1 - y0, out_w), dtype=image.dtype)
        strip[valid] = image[src_yi[valid], src_xi[valid]]
        output[y0:y1, :] = strip

    return output


def register_qupath_valis(sample_id, json_path, source_img_path,
                          target_img_path, work_dir, modality='comet'):
    """Full QuPath + VALIS registration pipeline.

    Pre-warps the source image with QuPath affine into target space,
    then runs VALIS rigid + non-rigid registration for fine alignment.

    Requires VALIS (valis-wsi) to be installed. Use conda env
    'maitra_spatial' on A16 laptop.

    Steps:
        1. Load QuPath affine from JSON (auto-inverts)
        2. Pre-warp source image into target space using QuPath affine
        3. Run VALIS on pre-warped source vs target
        4. Return VALIS registrar for GeoJSON warping

    Args:
        sample_id: e.g., 'SO1'
        json_path: Path to qupath_transforms.json
        source_img_path: Path to source image (e.g., COMET DAPI)
        target_img_path: Path to target image (e.g., STOmics DAPI)
        work_dir: Working directory for VALIS outputs
        modality: Transform key, e.g., 'comet'

    Returns:
        tuple: (registrar, qupath_matrix)
            registrar: VALIS Valis object with computed transforms
            qupath_matrix: The QuPath forward affine (for composing later)
    """
    try:
        import tifffile
    except ImportError:
        raise ImportError("tifffile required for register_qupath_valis")

    try:
        from valis import registration
    except ImportError:
        raise ImportError(
            "VALIS not installed. Use conda env 'maitra_spatial' or "
            "install with: pip install valis-wsi"
        )

    work_dir = Path(work_dir)

    # Step 1: Load QuPath affine
    qupath_matrix = register_qupath_affine(sample_id, json_path, modality)
    logger.info(f"Loaded QuPath affine for {sample_id}/{modality}")

    # Step 2: Pre-warp source into target space
    input_dir = work_dir / 'input_images'
    input_dir.mkdir(parents=True, exist_ok=True)

    prewarp_path = input_dir / f'{sample_id}_prewarped.tif'
    target_filename = Path(target_img_path).name
    target_dst = input_dir / target_filename

    # Copy target if not already there
    if not target_dst.exists():
        import shutil
        shutil.copy2(str(target_img_path), str(target_dst))

    # Pre-warp source if not already done
    if not prewarp_path.exists():
        logger.info("Pre-warping source image with QuPath affine...")
        source_img = tifffile.imread(str(source_img_path))
        target_img = tifffile.imread(str(target_img_path))
        target_shape = target_img.shape[:2]

        prewarped = warp_image_large(source_img, qupath_matrix, target_shape)
        tifffile.imwrite(str(prewarp_path), prewarped)
        logger.info(f"Saved pre-warped image: {prewarp_path}")
        del source_img, target_img, prewarped
    else:
        logger.info(f"Using cached pre-warped image: {prewarp_path}")

    # Step 3: Run VALIS
    logger.info("Running VALIS registration on pre-warped images...")
    reg_dir = str(work_dir / 'valis_output')
    registrar = registration.Valis(
        str(input_dir),
        reg_dir,
        align_to_reference=True,
        reference_img_f=target_filename,
    )

    rigid_registrar, non_rigid_registrar, error_df = registrar.register()
    logger.info(f"VALIS registration complete")
    logger.info(f"Error summary:\n{error_df}")

    return registrar, qupath_matrix


# ── Verification ────────────────────────────────────────────────────────────


def plot_registration_overlay(source_img, target_img, matrix_2x3,
                              output_path, title=None):
    """Save verification PNG with checkerboard and alpha blend overlays.

    Args:
        source_img: Grayscale or color source image
        target_img: Grayscale or color target image
        matrix_2x3: Forward affine transform (source → target space)
        output_path: Path to save PNG
        title: Optional figure title
    """
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available, skipping overlay plot")
        return

    M = np.array(matrix_2x3, dtype=np.float64)
    target_shape = target_img.shape[:2]

    # Warp source to target space
    src_gray = _ensure_gray_uint8(source_img)
    tgt_gray = _ensure_gray_uint8(target_img)
    warped = warp_image_affine(src_gray, M, target_shape)

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # Panel 1: Target
    axes[0].imshow(tgt_gray, cmap='gray')
    axes[0].set_title('Target (STOmics DAPI)')
    axes[0].axis('off')

    # Panel 2: Checkerboard
    block = 64
    h, w = target_shape
    checker = np.zeros((h, w), dtype=np.uint8)
    for i in range(0, h, block):
        for j in range(0, w, block):
            bi, bj = i // block, j // block
            i2 = min(i + block, h)
            j2 = min(j + block, w)
            if (bi + bj) % 2 == 0:
                checker[i:i2, j:j2] = tgt_gray[i:i2, j:j2]
            else:
                checker[i:i2, j:j2] = warped[i:i2, j:j2]
    axes[1].imshow(checker, cmap='gray')
    axes[1].set_title('Checkerboard')
    axes[1].axis('off')

    # Panel 3: Alpha blend (magenta/green)
    blend = np.zeros((h, w, 3), dtype=np.uint8)
    blend[..., 0] = warped      # Red channel = source
    blend[..., 1] = tgt_gray    # Green channel = target
    blend[..., 2] = warped      # Blue channel = source (magenta)
    axes[2].imshow(blend)
    axes[2].set_title('Blend (magenta=source, green=target)')
    axes[2].axis('off')

    if title:
        fig.suptitle(title, fontsize=14)

    plt.tight_layout()
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close(fig)
    logger.info(f"Saved registration overlay to {output_path}")


def compute_registration_metrics(source_img, target_img, matrix_2x3):
    """Compute registration quality metrics after warping.

    Args:
        source_img: Grayscale source image
        target_img: Grayscale target image
        matrix_2x3: Forward affine transform

    Returns:
        dict: {ncc_score, mutual_information, overlap_pct}
    """
    target_shape = target_img.shape[:2]
    src_gray = _ensure_gray_uint8(source_img)
    tgt_gray = _ensure_gray_uint8(target_img)
    warped = warp_image_affine(src_gray, matrix_2x3, target_shape)

    # Tissue masks (threshold > 0)
    src_mask = warped > 10
    tgt_mask = tgt_gray > 10

    # Overlap percentage
    intersection = np.logical_and(src_mask, tgt_mask).sum()
    union = np.logical_or(src_mask, tgt_mask).sum()
    overlap_pct = (intersection / union * 100) if union > 0 else 0.0

    # NCC on overlapping region
    overlap_mask = np.logical_and(src_mask, tgt_mask)
    if overlap_mask.sum() > 100:
        s = warped[overlap_mask].astype(np.float64)
        t = tgt_gray[overlap_mask].astype(np.float64)
        s_norm = s - s.mean()
        t_norm = t - t.mean()
        denom = np.sqrt((s_norm ** 2).sum() * (t_norm ** 2).sum())
        ncc = float((s_norm * t_norm).sum() / denom) if denom > 0 else 0.0
    else:
        ncc = 0.0

    return {
        'ncc_score': ncc,
        'overlap_pct': overlap_pct,
        'overlap_pixels': int(intersection),
    }
