"""
Generates synthetic ERA5-Land-like hourly meteorological data
for the Himalayan region. Physically realistic values for
spring fire season wind, soil moisture, and temperature.
"""
import numpy as np
import xarray as xr
import pandas as pd
from scipy.ndimage import gaussian_filter
import os

os.makedirs("data", exist_ok=True)

# Spatial grid covering the Himalayan region
lat_values = np.linspace(27.0, 31.0, 24)
lon_values = np.linspace(79.0, 88.0, 54)

# 24 hourly time steps starting from a reference fire event date
time_values = pd.date_range("2024-04-01 00:00", periods=24, freq="h")

hours = len(time_values)
nlat = len(lat_values)
nlon = len(lon_values)

rng = np.random.default_rng(seed=42)


def smooth_field(shape, sigma=3.0):
    """Creates a spatially coherent random field using Gaussian smoothing."""
    raw_noise = rng.standard_normal(shape).astype("float32")
    return gaussian_filter(raw_noise, sigma=sigma)


u10_data = np.zeros((hours, nlat, nlon), dtype="float32")
v10_data = np.zeros((hours, nlat, nlon), dtype="float32")
swvl1_data = np.zeros((hours, nlat, nlon), dtype="float32")
t2d_data = np.zeros((hours, nlat, nlon), dtype="float32")

for h in range(hours):
    # Wind speed peaks in afternoon, typical Himalayan spring westerlies
    diurnal_wind_factor = 1.0 + 0.4 * np.sin(2.0 * np.pi * (h - 6) / 24.0)

    u10_data[h] = (5.0 + smooth_field((nlat, nlon)) * 2.5) * diurnal_wind_factor
    v10_data[h] = (1.5 + smooth_field((nlat, nlon)) * 1.5) * diurnal_wind_factor

    # Soil moisture: low in spring dry season, slight spatial variation
    swvl1_base = 0.12 + smooth_field((nlat, nlon), sigma=4.0) * 0.04
    swvl1_data[h] = np.clip(swvl1_base, 0.04, 0.30)

    # Temperature: warmer at lower latitudes, peaks mid afternoon
    lat_gradient = (31.0 - lat_values[:, np.newaxis]) * 1.8
    diurnal_temp = 6.0 * np.sin(2.0 * np.pi * (h - 7) / 24.0)
    t2d_data[h] = (
        288.0
        + lat_gradient
        + diurnal_temp
        + smooth_field((nlat, nlon)) * 0.8
    )

dataset = xr.Dataset(
    {
        "u10": (["time", "latitude", "longitude"], u10_data),
        "v10": (["time", "latitude", "longitude"], v10_data),
        "swvl1": (["time", "latitude", "longitude"], swvl1_data),
        "t2d": (["time", "latitude", "longitude"], t2d_data),
    },
    coords={
        "time": time_values,
        "latitude": lat_values,
        "longitude": lon_values,
    },
    attrs={
        "description": "Synthetic ERA5-Land-like data for Himalayan fire pipeline",
        "spatial_coverage": "lat 27-31N, lon 79-88E",
        "temporal_coverage": "2024-04-01, 24 hours",
    },
)

dataset.to_netcdf("data/era5_land_hourly.nc")
print(f"Done: data/era5_land_hourly.nc")
print(f"Shape: {hours} hours x {nlat} lat x {nlon} lon")
print(f"u10 range: {u10_data.min():.2f} to {u10_data.max():.2f} m/s")
print(f"swvl1 range: {swvl1_data.min():.3f} to {swvl1_data.max():.3f} m3/m3")
print(f"t2d range: {t2d_data.min():.1f} to {t2d_data.max():.1f} K")