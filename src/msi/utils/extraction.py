"""TIFF reading and intensity extraction utilities for MSI data.

This module provides functions for extracting pixel intensities from
multi-channel OME-TIFF files commonly used in Mass Spectrometry Imaging.
"""

import logging
from pathlib import Path
from typing import Optional, Tuple, Union
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
import pandas as pd
import tifffile

logger = logging.getLogger(__name__)


def read_ome_tiff(
    filepath: Union[str, Path],
    channel: Optional[int] = None,
) -> Tuple[np.ndarray, dict]:
    """Read an OME-TIFF file and return image array and metadata.

    Parameters
    ----------
    filepath : str or Path
        Path to the OME-TIFF file.
    channel : int, optional
        Specific channel to read. If None, reads all channels.

    Returns
    -------
    image : np.ndarray
        Image array. Shape is (C, Y, X) for multi-channel or (Y, X) for single channel.
    metadata : dict
        Dictionary containing image metadata (shape, dtype, channel names if available).
    """
    filepath = Path(filepath)

    with tifffile.TiffFile(filepath) as tif:
        # Get image shape info
        series = tif.series[0]
        shape = series.shape
        dtype = series.dtype

        # Try to extract channel names from OME metadata
        channel_names = []
        if tif.ome_metadata:
            try:
                import xml.etree.ElementTree as ET
                root = ET.fromstring(tif.ome_metadata)
                ns = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}
                channels = root.findall('.//ome:Channel', ns)
                channel_names = [ch.get('Name', f'channel_{i}') for i, ch in enumerate(channels)]
            except Exception as e:
                logger.debug(f"Could not parse OME metadata: {e}")

        # Read image data
        if channel is not None:
            # Read specific channel
            image = tif.asarray(key=channel)
        else:
            image = tif.asarray()

    metadata = {
        'shape': shape,
        'dtype': str(dtype),
        'channel_names': channel_names if channel_names else None,
        'n_channels': shape[0] if len(shape) == 3 else 1,
    }

    return image, metadata


def extract_channel_intensities(
    image: np.ndarray,
    channel_idx: int,
    min_intensity: float = 0.0,
    channel_name: Optional[str] = None,
) -> pd.DataFrame:
    """Extract non-zero pixel intensities from a single channel.

    Parameters
    ----------
    image : np.ndarray
        Full image array with shape (C, Y, X) or single channel (Y, X).
    channel_idx : int
        Channel index to extract (ignored if image is 2D).
    min_intensity : float, default=0.0
        Minimum intensity threshold. Pixels below this are excluded.
    channel_name : str, optional
        Name for the intensity column. Defaults to 'intensity'.

    Returns
    -------
    df : pd.DataFrame
        DataFrame with columns ['x', 'y', 'intensity'] or ['x', 'y', channel_name].
        Only pixels above min_intensity are included.
    """
    # Get the channel data
    if image.ndim == 3:
        channel_data = image[channel_idx, :, :]
    else:
        channel_data = image

    # Find pixels above threshold
    mask = channel_data > min_intensity
    y_coords, x_coords = np.where(mask)
    intensities = channel_data[mask]

    # Create dataframe
    col_name = channel_name if channel_name else 'intensity'
    df = pd.DataFrame({
        'x': x_coords,
        'y': y_coords,
        col_name: intensities,
    })

    return df


def _process_single_channel(args: tuple) -> Tuple[int, pd.DataFrame]:
    """Worker function for parallel channel processing.

    Parameters
    ----------
    args : tuple
        (filepath, channel_idx, min_intensity, channel_name)

    Returns
    -------
    channel_idx : int
        The channel index that was processed.
    df : pd.DataFrame
        Extracted intensities for this channel.
    """
    filepath, channel_idx, min_intensity, channel_name = args

    # Read just this channel
    image, _ = read_ome_tiff(filepath, channel=channel_idx)

    # Extract intensities
    df = extract_channel_intensities(
        image,
        channel_idx=0,  # Already single channel
        min_intensity=min_intensity,
        channel_name=channel_name,
    )

    return channel_idx, df


def extract_all_channels(
    filepath: Union[str, Path],
    min_intensity: float = 0.0,
    workers: int = 1,
    channel_names: Optional[list] = None,
) -> dict:
    """Extract intensities from all channels in an OME-TIFF file.

    Parameters
    ----------
    filepath : str or Path
        Path to the OME-TIFF file.
    min_intensity : float, default=0.0
        Minimum intensity threshold for pixel inclusion.
    workers : int, default=1
        Number of parallel workers for channel processing.
    channel_names : list, optional
        List of channel names. If None, uses names from metadata or indices.

    Returns
    -------
    results : dict
        Dictionary mapping channel names/indices to DataFrames.
    """
    filepath = Path(filepath)

    # Get metadata to determine channels
    _, metadata = read_ome_tiff(filepath)
    n_channels = metadata['n_channels']

    # Determine channel names
    if channel_names is None:
        channel_names = metadata.get('channel_names')
    if channel_names is None:
        channel_names = [f'channel_{i}' for i in range(n_channels)]

    # Prepare arguments for each channel
    args_list = [
        (filepath, i, min_intensity, channel_names[i])
        for i in range(n_channels)
    ]

    results = {}

    if workers == 1:
        # Sequential processing
        for args in args_list:
            idx, df = _process_single_channel(args)
            results[channel_names[idx]] = df
            logger.info(f"Extracted channel {idx}: {len(df)} pixels")
    else:
        # Parallel processing
        with ProcessPoolExecutor(max_workers=workers) as executor:
            futures = {executor.submit(_process_single_channel, args): args[1]
                      for args in args_list}

            for future in as_completed(futures):
                idx, df = future.result()
                results[channel_names[idx]] = df
                logger.info(f"Extracted channel {idx}: {len(df)} pixels")

    return results


def batch_extract_images(
    input_dir: Union[str, Path],
    output_dir: Union[str, Path],
    pattern: str = "*.ome.tif",
    min_intensity: float = 0.0,
    workers_images: int = 1,
    workers_channels: int = 1,
) -> list:
    """Extract intensities from multiple TIFF files in a directory.

    Parameters
    ----------
    input_dir : str or Path
        Directory containing input TIFF files.
    output_dir : str or Path
        Directory to save output parquet files.
    pattern : str, default="*.ome.tif"
        Glob pattern to match input files.
    min_intensity : float, default=0.0
        Minimum intensity threshold.
    workers_images : int, default=1
        Number of parallel workers for processing images.
    workers_channels : int, default=1
        Number of parallel workers for processing channels within each image.

    Returns
    -------
    output_files : list
        List of paths to generated parquet files.
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find all matching files
    input_files = list(input_dir.glob(pattern))
    logger.info(f"Found {len(input_files)} files matching pattern '{pattern}'")

    output_files = []

    def process_single_image(filepath: Path) -> list:
        """Process a single image file."""
        sample_id = filepath.stem.replace('.ome', '')
        sample_dir = output_dir / sample_id
        sample_dir.mkdir(parents=True, exist_ok=True)

        logger.info(f"Processing {filepath.name}...")

        # Extract all channels
        results = extract_all_channels(
            filepath,
            min_intensity=min_intensity,
            workers=workers_channels,
        )

        # Save each channel as parquet
        saved_files = []
        for channel_name, df in results.items():
            # Clean channel name for filename
            safe_name = "".join(c if c.isalnum() or c in '-_' else '_' for c in str(channel_name))
            output_path = sample_dir / f"{safe_name}.parquet"
            df.to_parquet(output_path, index=False)
            saved_files.append(output_path)
            logger.info(f"  Saved {output_path.name}: {len(df)} pixels")

        return saved_files

    if workers_images == 1:
        # Sequential image processing
        for filepath in input_files:
            files = process_single_image(filepath)
            output_files.extend(files)
    else:
        # Parallel image processing
        with ProcessPoolExecutor(max_workers=workers_images) as executor:
            futures = {executor.submit(process_single_image, fp): fp for fp in input_files}

            for future in as_completed(futures):
                filepath = futures[future]
                try:
                    files = future.result()
                    output_files.extend(files)
                except Exception as e:
                    logger.error(f"Error processing {filepath}: {e}")

    logger.info(f"Extraction complete. Generated {len(output_files)} parquet files.")
    return output_files


def get_image_info(filepath: Union[str, Path]) -> dict:
    """Get information about an OME-TIFF file without loading full data.

    Parameters
    ----------
    filepath : str or Path
        Path to the OME-TIFF file.

    Returns
    -------
    info : dict
        Dictionary with image information including shape, dtype, channels.
    """
    filepath = Path(filepath)

    with tifffile.TiffFile(filepath) as tif:
        series = tif.series[0]

        info = {
            'filepath': str(filepath),
            'filename': filepath.name,
            'shape': series.shape,
            'dtype': str(series.dtype),
            'n_dimensions': len(series.shape),
        }

        if len(series.shape) == 3:
            info['n_channels'] = series.shape[0]
            info['height'] = series.shape[1]
            info['width'] = series.shape[2]
        elif len(series.shape) == 2:
            info['n_channels'] = 1
            info['height'] = series.shape[0]
            info['width'] = series.shape[1]

        # Check for OME metadata
        info['has_ome_metadata'] = tif.ome_metadata is not None

    return info
