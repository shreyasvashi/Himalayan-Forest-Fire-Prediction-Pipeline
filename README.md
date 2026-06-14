# Himalayan Forest Fire Prediction Pipeline

## Overview

This project implements an hourly forest fire spread prediction system for Himalayan terrain. It combines a Cellular Automata fire spread model with Markov Chain Monte Carlo calibration, and produces probabilistic threat level maps that can be used for emergency evacuation planning.

The pipeline integrates three real world data sources:

1. NASA SRTM 30 meter Digital Elevation Model (GeoTIFF) for slope and aspect calculation.
2. ECMWF ERA5-Land hourly NetCDF data for wind vectors and fuel moisture estimation.
3. NASA FIRMS VIIRS 375 meter active fire detections (CSV) for initializing and validating the burning state of the grid.

The output of the pipeline is a set of GeoTIFF rasters showing burn probability and a five level threat classification, plus a JSON summary describing the calibrated model parameters and the estimated area at each threat level.

## How the Pipeline Works

### 1. Topography Processing

The Digital Elevation Model is read in fixed size blocks (default 512 by 512 pixels) with a one pixel halo on each side. For each block, slope in degrees and aspect in degrees (measured clockwise from north) are computed using a Sobel gradient filter. This block based approach keeps memory usage low even for very large mountain regions.

### 2. Meteorology Processing

ERA5-Land data is opened with xarray using time based chunking. For each simulation hour, the pipeline extracts the 10 meter wind components (u10 and v10), converts them into wind speed and wind direction, and estimates a fuel moisture coefficient from soil moisture (swvl1) and 2 meter air temperature (t2d). A higher fuel moisture coefficient represents wetter fuel and reduces the chance of fire spreading.

### 3. Active Fire Initialization

FIRMS VIIRS active fire points are loaded as a GeoDataFrame, filtered by detection confidence, reprojected to match the DEM coordinate reference system, and rasterized onto the model grid. Cells that contain a confident active fire detection are marked as burning at the start of the simulation.

### 4. Cellular Automata Fire Spread Model

Each grid cell can be in one of three states: unburned, burning, or burned out. At every hourly time step, each burning cell can ignite its eight neighboring cells. The probability of ignition for each neighbor depends on three factors:

- Slope alignment: fire spreads faster uphill, based on the aspect of the cell and the steepness of the slope.
- Wind alignment and speed: fire spreads faster when the wind blows toward the neighboring cell, scaled by wind speed.
- Fuel moisture: wetter fuel reduces the probability of ignition.

After a burning cell has had a chance to ignite its neighbors for one hour, it transitions to the burned out state.

### 5. MCMC Calibration

The Cellular Automata model has four free parameters:

- base_spread_probability: the baseline chance that a burning cell ignites a neighbor under neutral conditions.
- slope_weight: how strongly slope alignment influences spread probability.
- wind_weight: how strongly wind alignment and speed influence spread probability.
- moisture_weight: how strongly fuel moisture reduces spread probability.

The pipeline uses the emcee ensemble sampler to explore the space of these four parameters. For each candidate parameter set, the model is run forward over a short calibration window using real ERA5-Land meteorology, and the resulting simulated burn extent is compared to the observed FIRMS burned area using an intersection over union score. The log likelihood is built from this score, and the posterior mean of the sampled chain (after discarding a burn in period) becomes the calibrated parameter set used for the full simulation.

### 6. Forward Simulation

Using the calibrated parameters, the model runs forward hour by hour for the requested simulation horizon (24 hours by default). At each hour, the model reads the corresponding ERA5-Land wind and moisture fields, advances the Cellular Automata by one step, and records the per cell ignition probability. The pipeline keeps a running maximum of these hourly probabilities, producing a single cumulative probability grid that represents the highest risk observed for each cell across the entire simulation period.

### 7. Threat Map and Evacuation Summary

The cumulative probability grid is classified into five threat levels:

- Safe: probability below 0.05
- Low: probability between 0.05 and 0.2
- Moderate: probability between 0.2 and 0.5
- High: probability between 0.5 and 0.8
- Critical: probability of 0.8 or higher, or any cell that is already burning or burned out

The pipeline writes two GeoTIFF rasters (the threat classification and the raw probability grid) and a JSON summary listing the cell count and area in square kilometers for each threat level, along with the total number of cells recommended for evacuation (high and critical combined).

## Project Structure

```
project_root/
  data/
    srtm_himalaya_30m.tif        SRTM 30 meter DEM, GeoTIFF
    era5_land_hourly.nc          ERA5-Land hourly NetCDF
    firms_viirs_active_fire.csv  FIRMS VIIRS active fire detections
  output/
    slope.tif                    Generated slope raster
    aspect.tif                   Generated aspect raster
    threat_map.tif               Five level threat classification raster
    probability_map.tif          Cumulative ignition probability raster
    run_summary.json             Calibrated parameters and evacuation summary
  fire_prediction_pipeline.py    Main pipeline script
  README.md                      This file
```

## Requirements

Python 3.10 or later is recommended.

Install the required packages with pip:

```
pip install numpy rasterio geopandas xarray netCDF4 scipy emcee
```

Notes on dependencies:

- rasterio requires the GDAL system library. On Ubuntu or Debian, install it first with `sudo apt-get install gdal-bin libgdal-dev` before installing rasterio with pip.
- geopandas depends on shapely, fiona, and pyproj. If pip installation fails, consider using conda or mamba with the conda-forge channel, which provides prebuilt binaries for all of these packages.
- xarray needs netCDF4 (or h5netcdf) installed to open ERA5-Land NetCDF files.
- emcee is a pure Python package and installs without additional system dependencies.

## Data Preparation

### SRTM Digital Elevation Model

1. Download a 30 meter SRTM tile covering your area of interest, for example from the USGS Earth Explorer service or from OpenTopography.
2. If your area of interest spans multiple tiles, merge them into a single GeoTIFF using a tool such as `gdal_merge.py` or `gdalwarp`.
3. Reproject the merged DEM to a projected coordinate reference system with meter based units (for example a local UTM zone), since the slope and aspect calculations in this pipeline assume the cell size is in meters.
4. Save the final file as `data/srtm_himalaya_30m.tif`.

### ERA5-Land Meteorology

1. Use the Climate Data Store (CDS) API from the Copernicus Climate Change Service to request ERA5-Land hourly data for your region and time period of interest.
2. Request at least the following variables: 10 meter u wind component (u10), 10 meter v wind component (v10), volumetric soil water layer 1 (swvl1), and 2 meter temperature (t2d).
3. Download the result as a NetCDF file and save it as `data/era5_land_hourly.nc`.
4. Confirm that the time dimension of this file covers at least the calibration window plus the full simulation horizon (default 24 hours), and that the spatial grid overlaps the DEM area.

### FIRMS VIIRS Active Fire Data

1. Visit the NASA FIRMS active fire data portal and request VIIRS 375 meter active fire detections for your region and date range, in CSV format.
2. Confirm the CSV contains at minimum the columns `latitude`, `longitude`, and `confidence`.
3. Save the file as `data/firms_viirs_active_fire.csv`.

## Configuration

All configuration values are defined in the `PipelineConfig` dataclass at the top of `fire_prediction_pipeline.py`. The most important fields to review before running the pipeline are:

- `dem_path`, `era5_path`, `firms_csv_path`: file paths to your input data.
- `output_dir`: directory where output rasters and the summary file will be written. This directory must exist before running the pipeline.
- `block_size`: size in pixels of the blocks used for slope and aspect processing. Reduce this value if you encounter memory issues on very large DEM files.
- `cell_size_m`: the size of each DEM pixel in meters. This must match the resolution of your reprojected DEM.
- `hours_to_simulate`: the number of hourly steps to simulate forward from the initial fire state.
- `mcmc_walkers`, `mcmc_steps`, `mcmc_burn_in`: control the MCMC calibration. More walkers and steps give a more thorough exploration of the parameter space but increase runtime.

## Running the Pipeline

1. Create the output directory if it does not already exist:

```
mkdir -p output
```

2. Edit the file paths and configuration values inside the `main` function of `fire_prediction_pipeline.py` to match your data, or replace the values when constructing `PipelineConfig`.

3. Run the script:

```
python fire_prediction_pipeline.py
```

4. The pipeline will perform the following steps in order, logging progress to the console:

   - Compute slope and aspect from the DEM and write them to `output/slope.tif` and `output/aspect.tif`.
   - Load ERA5-Land meteorology and FIRMS active fire points.
   - Rasterize the initial burning state from FIRMS detections.
   - Run MCMC calibration of the Cellular Automata parameters using a short calibration window.
   - Run the full hourly forward simulation using the calibrated parameters.
   - Classify the cumulative probability grid into threat levels.
   - Write `output/threat_map.tif`, `output/probability_map.tif`, and `output/run_summary.json`.

## Interpreting the Output

- `output/probability_map.tif` contains floating point values between 0 and 1 for each cell, representing the maximum hourly ignition probability observed during the simulation. This can be used to produce continuous risk surfaces or heat maps.
- `output/threat_map.tif` contains integer values from 0 to 4, corresponding to the safe, low, moderate, high, and critical threat classes described earlier. This raster is suitable for direct use in emergency planning systems, for example by highlighting high and critical zones for evacuation routing.
- `output/run_summary.json` contains the calibrated Cellular Automata parameters and a breakdown of the area (in square kilometers) and cell count for each threat level, along with the total number of cells recommended for evacuation.

## Performance and Memory Considerations

- Slope and aspect computation processes the DEM in blocks, so memory usage does not scale with the full size of the DEM. Reduce `block_size` if memory is constrained.
- The Cellular Automata step operates on full grid arrays using vectorized numpy operations across eight neighbor directions. For very large grids, consider tiling the simulation domain and processing each tile independently, then stitching the results at tile boundaries.
- MCMC calibration runs the Cellular Automata model once per walker per step, so the calibration window length and grid size directly affect runtime. Using a smaller calibration window (a few hours) and a coarser resolution copy of the grid for calibration can significantly speed up this stage while still producing a usable parameter set for the full resolution forward simulation.
- ERA5-Land data is opened with time based chunking through xarray, so only the meteorology for the current hour is loaded into memory at any given time.

## Extending the Pipeline

- Additional fuel type or vegetation density layers can be incorporated by adding a new grid to the Cellular Automata step function and including an additional weight parameter in the MCMC calibration.
- The threat level thresholds defined in `THREAT_THRESHOLDS` can be adjusted to match local emergency management guidelines.
- The neighbor offset list used in the Cellular Automata step can be extended to a larger neighborhood (for example a two cell radius) to model longer range ember transport, at the cost of additional computation per step.

## Limitations

- The Cellular Automata model assumes a uniform fuel type across the landscape. In reality, fuel type and load vary significantly and would improve prediction accuracy if included as an additional input layer.
- The fuel moisture estimate is a simplified proxy derived from soil moisture and air temperature. A dedicated fuel moisture model or field measurements would provide more accurate results.
- The MCMC calibration in this pipeline uses a single observed burned area snapshot for validation. Calibrating against multiple historical fire events would produce more robust parameter estimates.
