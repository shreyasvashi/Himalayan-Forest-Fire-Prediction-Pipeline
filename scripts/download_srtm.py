import boto3
from botocore import UNSIGNED
from botocore.config import Config
import rasterio
from rasterio.merge import merge
from rasterio.warp import calculate_default_transform, reproject, Resampling
import os
import glob

os.makedirs("data/srtm_tiles", exist_ok=True)

s3 = boto3.client(
    "s3",
    region_name="eu-central-1",
    config=Config(signature_version=UNSIGNED),
)

BUCKET = "copernicus-dem-30m"

lat_range = range(27, 31)
lon_range = range(79, 88)

for lat in lat_range:
    for lon in lon_range:
        folder = f"Copernicus_DSM_COG_10_N{lat:02d}_00_E{lon:03d}_00_DEM"
        s3_key = f"{folder}/{folder}.tif"
        local_path = f"data/srtm_tiles/N{lat:02d}E{lon:03d}.tif"

        if os.path.exists(local_path):
            print(f"Already have N{lat:02d}E{lon:03d}, skipping")
            continue

        try:
            print(f"Downloading N{lat:02d}E{lon:03d} ...")
            s3.download_file(BUCKET, s3_key, local_path)
            print(f"Saved N{lat:02d}E{lon:03d}")
        except Exception as e:
            print(f"Skipping N{lat:02d}E{lon:03d}: {e}")

tile_files = sorted(glob.glob("data/srtm_tiles/*.tif"))
if not tile_files:
    raise RuntimeError("No tiles downloaded. Check internet connection.")

print(f"Merging {len(tile_files)} tiles ...")
datasets = [rasterio.open(f) for f in tile_files]
merged_array, merged_transform = merge(datasets)
merged_profile = datasets[0].profile.copy()
merged_profile.update({
    "height": merged_array.shape[1],
    "width": merged_array.shape[2],
    "transform": merged_transform,
})

with rasterio.open("data/srtm_merged.tif", "w", **merged_profile) as f:
    f.write(merged_array)

for ds in datasets:
    ds.close()

print("Reprojecting to UTM Zone 44N (EPSG:32644) ...")
with rasterio.open("data/srtm_merged.tif") as src:
    dst_crs = "EPSG:32644"
    dst_transform, dst_width, dst_height = calculate_default_transform(
        src.crs, dst_crs, src.width, src.height, *src.bounds
    )
    dst_profile = src.meta.copy()
    dst_profile.update({
        "crs": dst_crs,
        "transform": dst_transform,
        "width": dst_width,
        "height": dst_height,
        "dtype": "float32",
    })
    with rasterio.open("data/srtm_himalaya_30m.tif", "w", **dst_profile) as dst:
        reproject(
            source=rasterio.band(src, 1),
            destination=rasterio.band(dst, 1),
            src_transform=src.transform,
            src_crs=src.crs,
            dst_transform=dst_transform,
            dst_crs=dst_crs,
            resampling=Resampling.bilinear,
        )

print("Done: data/srtm_himalaya_30m.tif")