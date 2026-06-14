import rasterio
from rasterio.enums import Resampling
input_path = "data/srtm_himalaya_30m.tif"
output_path = "data/srtm_himalaya_250m.tif"

downsample_factor = 8  # 30m x 8 = 240m, close enough to 250m

with rasterio.open(input_path) as src:
    new_width = src.width // downsample_factor
    new_height = src.height // downsample_factor

    data = src.read(
        out_shape=(1, new_height, new_width),
        resampling=Resampling.bilinear,
    )

    new_transform = src.transform * src.transform.scale(
        src.width / new_width,
        src.height / new_height,
    )

    profile = src.profile.copy()
    profile.update({
        "width": new_width,
        "height": new_height,
        "transform": new_transform,
        "dtype": "float32",
    })

    with rasterio.open(output_path, "w", **profile) as dst:
        dst.write(data)

print(f"Done: {output_path}")
print(f"New size: {new_width} x {new_height} pixels")