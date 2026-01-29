import os
import pandas as pd
import tifffile as tiff
from skimage.measure import label, regionprops
from multiprocessing import Pool
from tqdm import tqdm
from pathlib import Path


home_path = os.path.expanduser("~")
data_path = Path(home_path) / 'ext_hd_sammy' / 'data' / 'COMET'

image_path = data_path / 'original'
visio_out = data_path / 'visio_output' / '6_Samples_Trial' / '6_Samples_mld_xml'

save_dir = Path(home_path) / 'ext_hd_sammy' / 'projects' / 'out' /'out_comet'
os.makedirs(Path(save_dir), exist_ok=True)

print(image_path)
print(visio_out)


### markers

markers = [
    'DAPI',
    'CY5',
    'TRITC',
    'CD4',
    'CD31',
    'CD8',
    'CD45RO',
    'CD163',
    'FOXp3',
    'TOX1',
    'CD20',
    'CD86',
    'CD68',
    'CD10',
    'CD133',
    'CD11C',
    'GZMB',
    'Periostin',
    'CD66B',
    'COL1A1',
    'CD56',
    'KERATIN8',
    'CD45',
    'aSMA']

def create_batches(cell_data, batch_size):
    """Yield batches of regionprops objects."""
    for i in range(0, len(cell_data), batch_size):
        yield cell_data[i:i + batch_size]


def centroid_protein_intensity_extraction(args):
    """Extract protein intensity for a batch of cells."""
    batch, orig_img, single_cell, markers, save_dir, visio_object = args
    protein_df = pd.DataFrame(columns=markers + ['centroid_x', 'centroid_y'])

    for prop in batch:
        centroid_y, centroid_x = prop.centroid
        lbl = prop.label
        mask = single_cell == lbl

        protein_expression = []
        for channel in range(orig_img.shape[0]):
            channel_img = orig_img[channel, ...]
            mask_mtx = channel_img[mask]
            expression = mask_mtx.mean() if mask_mtx.size > 0 else 0
            protein_expression.append(expression)

        protein_expression.extend([centroid_x, centroid_y])
        protein_df.loc[lbl] = protein_expression

    # Write per-batch output
    batch_lbls = [prop.label for prop in batch]
    batch_id = f"{min(batch_lbls)}-{max(batch_lbls)}"
    base_name = visio_object.split('.tif')[0]
    batch_parquet_path = os.path.join(save_dir, f"{base_name}_batch{batch_id}_protein.parquet")
    protein_df.to_parquet(batch_parquet_path, index=False)

    return batch_parquet_path


def process_visio_object(visio_object, visio_out, orig_images, image_path, markers, save_dir, batch_size):
    """Process one Visio object image and extract protein intensity in batches."""
    base_name = visio_object.split('.tif')[0]
    print(base_name)

    # Read segmentation mask image
    v_img = tiff.imread(os.path.join(visio_out, visio_object))
    single_cell = label(v_img[0, ...] > 0)
    props = regionprops(single_cell)

    # Match original image
    try:
        orig_img = next(
            tiff.imread(os.path.join(image_path, img))
            for img in orig_images if base_name in img
        )
    except StopIteration:
        print(f"[Warning] Original image not found for: {visio_object}")
        return

    # Create batches and prepare arguments
    batches = list(create_batches(props, batch_size=batch_size))
    args_list = [(batch, orig_img, single_cell, markers, save_dir, visio_object)
                 for batch in batches]

    # Run multiprocessing
    with Pool(processes=4, maxtasksperchild=1) as pool:
        list(tqdm(pool.map(centroid_protein_intensity_extraction, args_list),
                  total=len(args_list), desc=f"Processing batches for {base_name}"))

    # Combine all batch files into one final output
    all_batch_files = [os.path.join(save_dir, f) for f in os.listdir(save_dir)
                       if f.startswith(base_name) and f.endswith('_protein.parquet')]

    if all_batch_files:
        combined_df = pd.concat([pd.read_parquet(f) for f in all_batch_files], axis=0)
        combined_path = os.path.join(save_dir, f"{base_name}_protein_combined.parquet")
        combined_df.to_parquet(combined_path, index=False)

        # Delete temp batch parquet files
        for f in all_batch_files:
            os.remove(f)

if __name__ == "__main__":
    # === USER CONFIGURATION ===
    visio_out = visio_out
    image_path = image_path
    save_dir = save_dir
    markers = markers  # Modify this list to your marker names
    batch_size = 150
    num_processes = 12  # or use os.cpu_count() // 2

    # ==========================
    os.makedirs(save_dir, exist_ok=True)

    visio_objects = sorted([f for f in os.listdir(visio_out) if f.endswith('.tif')])
    orig_images = sorted([f for f in os.listdir(image_path) if f.endswith('.tiff')])

    for visio_object in tqdm(visio_objects, desc="Processing visio objects"):
        process_visio_object(
            visio_object=visio_object,
            visio_out=visio_out,
            orig_images=orig_images,
            image_path=image_path,
            markers=markers,
            save_dir=save_dir,
            batch_size=batch_size
        )