#!/usr/bin/env python3
"""MSI Intensity Extraction CLI Tool.

Extracts pixel intensities from multi-channel OME-TIFF files commonly used
in Mass Spectrometry Imaging (MSI). Supports batch processing with parallel
execution for both images and channels.

Usage:
    python script01_intensity_extraction.py \
        --input-dir /path/to/tiffs \
        --output-dir /path/to/output \
        --modality glycans \
        --pattern "*.ome.tif" \
        --min-intensity 1.0 \
        --workers-images 4 \
        --workers-channels 8

Output formats:
    --output-format single-file (default, recommended):
        Creates one parquet file per image with all channels as columns.
        Output: output_dir/sample_id.parquet
        Columns: x, y, <channel_name_1>, <channel_name_2>, ...
        Preserves original channel names exactly (no filesystem restrictions).

    --output-format per-channel (legacy):
        Creates per-channel parquet files with columns: x, y, intensity
        Output: output_dir/sample_id/channel_name.parquet
        Channel names are sanitized for filesystem compatibility.
"""

import argparse
import logging
import sys
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Optional, Tuple
import time

import numpy as np
import pandas as pd
import tifffile

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S'
)
logger = logging.getLogger(__name__)


def read_ome_tiff_metadata(filepath: Path) -> dict:
    """Read metadata from an OME-TIFF file without loading full data."""
    with tifffile.TiffFile(filepath) as tif:
        series = tif.series[0]
        shape = series.shape
        dtype = series.dtype

        # Try to extract channel names from OME metadata
        channel_names = []
        if tif.ome_metadata:
            try:
                import xml.etree.ElementTree as ET
                root = ET.fromstring(tif.ome_metadata)
                # Try common OME namespace variants
                for ns_uri in ['http://www.openmicroscopy.org/Schemas/OME/2016-06',
                              'http://www.openmicroscopy.org/Schemas/OME/2015-01']:
                    ns = {'ome': ns_uri}
                    channels = root.findall('.//ome:Channel', ns)
                    if channels:
                        channel_names = [ch.get('Name', f'channel_{i}')
                                        for i, ch in enumerate(channels)]
                        break
            except Exception as e:
                logger.debug(f"Could not parse OME metadata: {e}")

        return {
            'shape': shape,
            'dtype': str(dtype),
            'n_channels': shape[0] if len(shape) == 3 else 1,
            'height': shape[1] if len(shape) == 3 else shape[0],
            'width': shape[2] if len(shape) == 3 else shape[1],
            'channel_names': channel_names if channel_names else None,
        }


def extract_single_channel(
    filepath: Path,
    channel_idx: int,
    min_intensity: float,
    channel_name: str,
) -> Tuple[str, pd.DataFrame, int]:
    """Extract intensities from a single channel.

    Returns:
        Tuple of (channel_name, DataFrame, n_pixels)
    """
    with tifffile.TiffFile(filepath) as tif:
        # Read just this channel
        if len(tif.series[0].shape) == 3:
            channel_data = tif.asarray(key=channel_idx)
        else:
            channel_data = tif.asarray()

    # Find pixels above threshold
    mask = channel_data > min_intensity
    y_coords, x_coords = np.where(mask)
    intensities = channel_data[mask]

    df = pd.DataFrame({
        'x': x_coords.astype(np.int32),
        'y': y_coords.astype(np.int32),
        'intensity': intensities.astype(np.float32),
    })

    return channel_name, df, len(df)


def extract_all_channels_combined(
    filepath: Path,
    min_intensity: float,
    channel_names: List[str],
) -> Tuple[pd.DataFrame, int]:
    """Extract all channels and combine into a single wide DataFrame.

    Each unique (x, y) coordinate gets one row, with each channel as a column.
    Only coordinates where at least one channel exceeds min_intensity are kept.

    Returns:
        Tuple of (combined DataFrame, total non-null values)
    """
    with tifffile.TiffFile(filepath) as tif:
        # Load full image data
        data = tif.asarray()

    if len(data.shape) == 2:
        # Single channel image
        data = data[np.newaxis, ...]

    n_channels, height, width = data.shape

    # Create mask of any pixel above threshold in any channel
    combined_mask = np.any(data > min_intensity, axis=0)
    y_coords, x_coords = np.where(combined_mask)

    if len(x_coords) == 0:
        # No pixels above threshold
        cols = {'x': np.array([], dtype=np.int32), 'y': np.array([], dtype=np.int32)}
        for name in channel_names:
            cols[name] = np.array([], dtype=np.float32)
        return pd.DataFrame(cols), 0

    # Build DataFrame with x, y and all channel values
    result = {
        'x': x_coords.astype(np.int32),
        'y': y_coords.astype(np.int32),
    }

    total_values = 0
    for i, name in enumerate(channel_names):
        channel_values = data[i, y_coords, x_coords].astype(np.float32)
        # Set values below threshold to NaN (sparse representation)
        channel_values[channel_values <= min_intensity] = np.nan
        result[name] = channel_values
        total_values += np.sum(~np.isnan(channel_values))

    df = pd.DataFrame(result)
    return df, int(total_values)


def process_single_image(
    filepath: Path,
    output_dir: Path,
    min_intensity: float,
    workers_channels: int,
    output_format: str = 'single-file',
) -> Tuple[str, int, int, float]:
    """Process all channels from a single image file.

    Args:
        filepath: Path to OME-TIFF file
        output_dir: Output directory
        min_intensity: Minimum intensity threshold
        workers_channels: Number of parallel workers for channels (per-channel mode only)
        output_format: 'single-file' (one parquet per image) or 'per-channel' (legacy)

    Returns:
        Tuple of (sample_id, n_channels, total_pixels, elapsed_time)
    """
    start_time = time.time()

    # Get sample ID from filename
    sample_id = filepath.stem.replace('.ome', '')

    # Read metadata
    metadata = read_ome_tiff_metadata(filepath)
    n_channels = metadata['n_channels']

    # Determine channel names
    if metadata['channel_names']:
        channel_names = metadata['channel_names']
    else:
        channel_names = [f'channel_{i}' for i in range(n_channels)]

    logger.info(f"Processing {filepath.name}: {n_channels} channels")

    if output_format == 'single-file':
        # New format: single parquet with all channels as columns
        output_dir.mkdir(parents=True, exist_ok=True)

        df, total_values = extract_all_channels_combined(
            filepath, min_intensity, channel_names
        )

        output_path = output_dir / f"{sample_id}.parquet"
        df.to_parquet(output_path, index=False, compression='snappy')

        elapsed = time.time() - start_time
        logger.info(f"  Completed {sample_id}: {len(df):,} pixels, {total_values:,} values in {elapsed:.1f}s")

        return sample_id, n_channels, total_values, elapsed

    else:
        # Legacy format: per-channel parquet files
        sample_dir = output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        total_pixels = 0
        results = []

        if workers_channels <= 1:
            # Sequential channel processing
            for i, name in enumerate(channel_names):
                ch_name, df, n_px = extract_single_channel(
                    filepath, i, min_intensity, name
                )
                results.append((ch_name, df, n_px))
        else:
            # Parallel channel processing
            with ProcessPoolExecutor(max_workers=workers_channels) as executor:
                futures = {
                    executor.submit(
                        extract_single_channel,
                        filepath, i, min_intensity, name
                    ): i
                    for i, name in enumerate(channel_names)
                }

                for future in as_completed(futures):
                    try:
                        result = future.result()
                        results.append(result)
                    except Exception as e:
                        idx = futures[future]
                        logger.error(f"Error processing channel {idx}: {e}")

        # Save results
        for ch_name, df, n_px in results:
            # Clean channel name for filename
            safe_name = "".join(c if c.isalnum() or c in '-_' else '_' for c in str(ch_name))
            output_path = sample_dir / f"{safe_name}.parquet"
            df.to_parquet(output_path, index=False, compression='snappy')
            total_pixels += n_px

        elapsed = time.time() - start_time
        logger.info(f"  Completed {sample_id}: {total_pixels:,} pixels in {elapsed:.1f}s")

        return sample_id, n_channels, total_pixels, elapsed


def batch_extract(
    input_dir: Path,
    output_dir: Path,
    pattern: str,
    min_intensity: float,
    workers_images: int,
    workers_channels: int,
    output_format: str = 'single-file',
) -> List[dict]:
    """Extract intensities from all matching files in a directory."""
    # Find input files
    input_files = list(input_dir.glob(pattern))

    if not input_files:
        logger.error(f"No files found matching pattern '{pattern}' in {input_dir}")
        return []

    logger.info(f"Found {len(input_files)} files to process")
    logger.info(f"Output directory: {output_dir}")
    logger.info(f"Output format: {output_format}")
    logger.info(f"Min intensity threshold: {min_intensity}")
    if output_format == 'per-channel':
        logger.info(f"Workers: {workers_images} images, {workers_channels} channels")
    else:
        logger.info(f"Workers: {workers_images} images")

    output_dir.mkdir(parents=True, exist_ok=True)
    results = []
    total_start = time.time()

    if workers_images <= 1:
        # Sequential image processing
        for filepath in input_files:
            try:
                sample_id, n_ch, n_px, elapsed = process_single_image(
                    filepath, output_dir, min_intensity, workers_channels, output_format
                )
                results.append({
                    'sample_id': sample_id,
                    'n_channels': n_ch,
                    'n_pixels': n_px,
                    'elapsed_time': elapsed,
                    'status': 'success',
                })
            except Exception as e:
                logger.error(f"Error processing {filepath}: {e}")
                results.append({
                    'sample_id': filepath.stem,
                    'status': 'error',
                    'error': str(e),
                })
    else:
        # Parallel image processing
        with ProcessPoolExecutor(max_workers=workers_images) as executor:
            futures = {
                executor.submit(
                    process_single_image,
                    filepath, output_dir, min_intensity, workers_channels, output_format
                ): filepath
                for filepath in input_files
            }

            for future in as_completed(futures):
                filepath = futures[future]
                try:
                    sample_id, n_ch, n_px, elapsed = future.result()
                    results.append({
                        'sample_id': sample_id,
                        'n_channels': n_ch,
                        'n_pixels': n_px,
                        'elapsed_time': elapsed,
                        'status': 'success',
                    })
                except Exception as e:
                    logger.error(f"Error processing {filepath}: {e}")
                    results.append({
                        'sample_id': filepath.stem,
                        'status': 'error',
                        'error': str(e),
                    })

    total_elapsed = time.time() - total_start

    # Summary
    successful = [r for r in results if r['status'] == 'success']
    total_pixels = sum(r['n_pixels'] for r in successful)
    total_channels = sum(r['n_channels'] for r in successful)

    logger.info("=" * 60)
    logger.info("EXTRACTION COMPLETE")
    logger.info(f"  Files processed: {len(successful)}/{len(input_files)}")
    logger.info(f"  Total channels: {total_channels}")
    logger.info(f"  Total pixels: {total_pixels:,}")
    logger.info(f"  Total time: {total_elapsed:.1f}s")
    logger.info("=" * 60)

    # Save summary
    summary_df = pd.DataFrame(results)
    summary_path = output_dir / "extraction_summary.csv"
    summary_df.to_csv(summary_path, index=False)
    logger.info(f"Summary saved to {summary_path}")

    return results


def main():
    parser = argparse.ArgumentParser(
        description="Extract pixel intensities from MSI OME-TIFF files.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Extract glycan intensities (default: single-file format)
  python script01_intensity_extraction.py \\
      --input-dir /data/msi/glycans \\
      --output-dir /out/glycans_parquet \\
      --modality glycans

  # Extract with parallel processing
  python script01_intensity_extraction.py \\
      --input-dir /data/msi/metabolites \\
      --output-dir /out/metab_parquet \\
      --modality metabolites \\
      --min-intensity 0.5 \\
      --workers-images 4

  # Legacy per-channel format (for backward compatibility)
  python script01_intensity_extraction.py \\
      --input-dir /data/msi/peptides \\
      --output-dir /out/peptides_parquet \\
      --modality peptides \\
      --output-format per-channel \\
      --workers-images 4 \\
      --workers-channels 8
        """
    )

    parser.add_argument(
        '--input-dir', '-i',
        type=Path,
        required=True,
        help='Directory containing input OME-TIFF files'
    )
    parser.add_argument(
        '--output-dir', '-o',
        type=Path,
        required=True,
        help='Directory for output parquet files'
    )
    parser.add_argument(
        '--modality', '-m',
        type=str,
        choices=['glycans', 'metabolites', 'peptides'],
        help='Data modality (for logging/organization)'
    )
    parser.add_argument(
        '--pattern', '-p',
        type=str,
        default='*.ome.tif',
        help='Glob pattern for input files (default: *.ome.tif)'
    )
    parser.add_argument(
        '--min-intensity',
        type=float,
        default=1.0,
        help='Minimum intensity threshold for pixel inclusion (default: 1.0)'
    )
    parser.add_argument(
        '--workers-images',
        type=int,
        default=1,
        help='Number of parallel workers for processing images (default: 1)'
    )
    parser.add_argument(
        '--workers-channels',
        type=int,
        default=1,
        help='Number of parallel workers for processing channels (per-channel mode only, default: 1)'
    )
    parser.add_argument(
        '--output-format', '-f',
        type=str,
        choices=['single-file', 'per-channel'],
        default='single-file',
        help='Output format: single-file (one parquet per image, recommended) or per-channel (legacy). Default: single-file'
    )
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose logging'
    )

    args = parser.parse_args()

    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    # Validate input directory
    if not args.input_dir.exists():
        logger.error(f"Input directory does not exist: {args.input_dir}")
        sys.exit(1)

    if args.modality:
        logger.info(f"Modality: {args.modality}")

    # Run extraction
    results = batch_extract(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        pattern=args.pattern,
        min_intensity=args.min_intensity,
        workers_images=args.workers_images,
        workers_channels=args.workers_channels,
        output_format=args.output_format,
    )

    # Exit with error code if any failed
    if any(r['status'] == 'error' for r in results):
        sys.exit(1)


if __name__ == '__main__':
    main()
