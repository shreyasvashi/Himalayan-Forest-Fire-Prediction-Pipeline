"""
Himalayan Forest Fire Prediction Pipeline
Cellular Automata fire spread model calibrated via MCMC,
driven by SRTM DEM, ERA5-Land wind and moisture data, and FIRMS VIIRS active fire data.
"""

import numpy as np
import rasterio
from rasterio.windows import Window
import xarray as xr
import geopandas as gpd
from scipy.ndimage import sobel
import emcee
import json
import logging
from dataclasses import dataclass
from typing import Tuple, Dict

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("himalayan_fire_pipeline")

# ----------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------

@dataclass
class PipelineConfig:
    dem_path: str
    era5_path: str
    firms_csv_path: str
    output_dir: str = "./output"
    block_size: int = 512
    cell_size_m: float = 30.0
    hours_to_simulate: int = 24
    mcmc_walkers: int = 32
    mcmc_steps: int = 1500
    mcmc_burn_in: int = 300
    no_data_value: float = -9999.0


# ----------------------------------------------------------------------
# Topography loader: slope and aspect from SRTM DEM, block by block
# ----------------------------------------------------------------------

class TopographyLoader:
    """Loads SRTM DEM and computes slope and aspect rasters using memory safe block processing."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.dem_path = config.dem_path
        self.cell_size = config.cell_size_m

    def compute_slope_aspect(self, output_slope_path: str, output_aspect_path: str) -> None:
        with rasterio.open(self.dem_path) as dem_dataset:
            profile = dem_dataset.profile.copy()
            profile.update(dtype="float32", count=1, nodata=self.config.no_data_value)

            height = dem_dataset.height
            width = dem_dataset.width
            block_size = self.config.block_size

            with rasterio.open(output_slope_path, "w", **profile) as slope_dataset, \
                 rasterio.open(output_aspect_path, "w", **profile) as aspect_dataset:

                for row_start in range(0, height, block_size):
                    for col_start in range(0, width, block_size):
                        row_count = min(block_size, height - row_start)
                        col_count = min(block_size, width - col_start)

                        halo = 1
                        read_row_start = max(row_start - halo, 0)
                        read_col_start = max(col_start - halo, 0)
                        read_row_stop = min(row_start + row_count + halo, height)
                        read_col_stop = min(col_start + col_count + halo, width)

                        read_window = Window(
                            read_col_start,
                            read_row_start,
                            read_col_stop - read_col_start,
                            read_row_stop - read_row_start,
                        )

                        elevation_block = dem_dataset.read(1, window=read_window).astype("float32")

                        slope_block, aspect_block = self._gradient_to_slope_aspect(elevation_block)

                        trim_top = row_start - read_row_start
                        trim_left = col_start - read_col_start

                        slope_trimmed = slope_block[trim_top: trim_top + row_count, trim_left: trim_left + col_count]
                        aspect_trimmed = aspect_block[trim_top: trim_top + row_count, trim_left: trim_left + col_count]

                        write_window = Window(col_start, row_start, col_count, row_count)
                        slope_dataset.write(slope_trimmed, 1, window=write_window)
                        aspect_dataset.write(aspect_trimmed, 1, window=write_window)

        log.info("Slope and aspect rasters written to %s and %s", output_slope_path, output_aspect_path)

    def _gradient_to_slope_aspect(self, elevation: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Computes slope in degrees and aspect in degrees (0 equals north, clockwise) from an elevation array."""
        dz_dx = sobel(elevation, axis=1) / (8.0 * self.cell_size)
        dz_dy = sobel(elevation, axis=0) / (8.0 * self.cell_size)

        slope_radians = np.arctan(np.hypot(dz_dx, dz_dy))
        slope_degrees = np.degrees(slope_radians).astype("float32")

        aspect_radians = np.arctan2(dz_dy, -dz_dx)
        aspect_degrees = np.degrees(aspect_radians)
        aspect_degrees = np.where(aspect_degrees < 0, 90.0 - aspect_degrees, 90.0 - aspect_degrees)
        aspect_degrees = np.mod(aspect_degrees, 360.0).astype("float32")

        return slope_degrees, aspect_degrees


# ----------------------------------------------------------------------
# Meteorology loader: ERA5-Land hourly wind and fuel moisture
# ----------------------------------------------------------------------

class MeteorologyLoader:
    """Loads ERA5-Land NetCDF and exposes hourly wind speed, wind direction, and fuel moisture coefficient grids."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.dataset = xr.open_dataset(config.era5_path, chunks={"time": 1})

    def get_hourly_wind_and_moisture(self, hour_index: int) -> Dict[str, np.ndarray]:
        time_slice = self.dataset.isel(time=hour_index)

        wind_u = time_slice["u10"].values.astype("float32")
        wind_v = time_slice["v10"].values.astype("float32")

        wind_speed_ms = np.hypot(wind_u, wind_v)
        wind_direction_deg = np.mod(np.degrees(np.arctan2(wind_u, wind_v)), 360.0)

        fuel_moisture_coefficient = self._estimate_fuel_moisture(time_slice)

        return {
            "wind_speed_ms": wind_speed_ms,
            "wind_direction_deg": wind_direction_deg.astype("float32"),
            "fuel_moisture_coefficient": fuel_moisture_coefficient,
        }

    def _estimate_fuel_moisture(self, time_slice: xr.Dataset) -> np.ndarray:
        """
        Derives a normalized fuel moisture coefficient from soil moisture and 2 meter temperature.
        A higher coefficient means wetter fuel and a lower fire spread probability.
        """
        soil_moisture_layer = time_slice["swvl1"].values.astype("float32")
        air_temperature_kelvin = time_slice["t2d"].values.astype("float32")

        normalized_soil_moisture = np.clip(soil_moisture_layer / 0.5, 0.0, 1.0)

        temperature_celsius = air_temperature_kelvin - 273.15
        temperature_penalty = np.clip((temperature_celsius - 15.0) / 50.0, 0.0, 0.5)

        fuel_moisture_coefficient = np.clip(normalized_soil_moisture - temperature_penalty, 0.0, 1.0)
        return fuel_moisture_coefficient.astype("float32")

    def get_grid_shape(self) -> Tuple[int, int]:
        sample = self.dataset.isel(time=0)["u10"].values
        return sample.shape


# ----------------------------------------------------------------------
# Active fire loader: FIRMS VIIRS active fire CSV
# ----------------------------------------------------------------------

class ActiveFireLoader:
    """Loads NASA FIRMS VIIRS active fire detections and rasterizes them onto the model grid."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def load_active_fire_points(self) -> gpd.GeoDataFrame:
        fire_points = gpd.read_file(self.config.firms_csv_path)

        if "latitude" not in fire_points.columns or "longitude" not in fire_points.columns:
            raise ValueError("FIRMS CSV must contain latitude and longitude columns")

        fire_points["latitude"] = fire_points["latitude"].astype("float64")
        fire_points["longitude"] = fire_points["longitude"].astype("float64")

        geometry = gpd.points_from_xy(fire_points["longitude"], fire_points["latitude"])
        fire_geodataframe = gpd.GeoDataFrame(fire_points, geometry=geometry, crs="EPSG:4326")

        if "confidence" in fire_geodataframe.columns:
            fire_geodataframe["confidence_numeric"] = fire_geodataframe["confidence"].map(
                {"l": 0.3, "n": 0.6, "h": 0.9, "low": 0.3, "nominal": 0.6, "high": 0.9}
            ).fillna(0.5)
        else:
            fire_geodataframe["confidence_numeric"] = 0.5

        return fire_geodataframe

    def rasterize_initial_state(
        self,
        fire_points: gpd.GeoDataFrame,
        reference_raster_path: str,
        confidence_threshold: float = 0.3,
    ) -> np.ndarray:
        from rasterio.features import rasterize

        with rasterio.open(reference_raster_path) as reference_dataset:
            transform = reference_dataset.transform
            output_shape = (reference_dataset.height, reference_dataset.width)
            raster_crs = reference_dataset.crs

        fire_points_reprojected = fire_points.to_crs(raster_crs)
        confident_fires = fire_points_reprojected[
            fire_points_reprojected["confidence_numeric"] >= confidence_threshold
        ]

        if confident_fires.empty:
            log.warning("No active fire points met the confidence threshold, returning an empty grid")
            return np.zeros(output_shape, dtype="uint8")

        shapes = ((geometry, 1) for geometry in confident_fires.geometry)
        burning_grid = rasterize(
            shapes,
            out_shape=output_shape,
            transform=transform,
            fill=0,
            dtype="uint8",
        )

        return burning_grid


# ----------------------------------------------------------------------
# Cellular Automata fire spread model
# ----------------------------------------------------------------------

STATE_UNBURNED = 0
STATE_BURNING = 1
STATE_BURNED_OUT = 2

NEIGHBOR_OFFSETS = [
    (-1, 0, 0.0),
    (-1, 1, 45.0),
    (0, 1, 90.0),
    (1, 1, 135.0),
    (1, 0, 180.0),
    (1, -1, 225.0),
    (0, -1, 270.0),
    (-1, -1, 315.0),
]


class FireSpreadCellularAutomata:
    """
    Hourly fire spread Cellular Automata. The spread probability for each neighbor
    depends on slope alignment, wind alignment and speed, and fuel moisture.
    """

    def __init__(
        self,
        slope_grid: np.ndarray,
        aspect_grid: np.ndarray,
        base_spread_probability: float,
        slope_weight: float,
        wind_weight: float,
        moisture_weight: float,
    ):
        self.slope_grid = slope_grid
        self.aspect_grid = aspect_grid
        self.base_spread_probability = base_spread_probability
        self.slope_weight = slope_weight
        self.wind_weight = wind_weight
        self.moisture_weight = moisture_weight

        self.grid_height, self.grid_width = slope_grid.shape

    def step(
        self,
        current_state: np.ndarray,
        wind_speed_grid: np.ndarray,
        wind_direction_grid: np.ndarray,
        fuel_moisture_grid: np.ndarray,
        random_generator: np.random.Generator,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Advances the fire state by one hour.
        Returns the new state grid and the per cell ignition probability grid
        for unburned cells.
        """
        new_state = current_state.copy()
        ignition_probability_total = np.zeros_like(self.slope_grid, dtype="float32")

        burning_mask = current_state == STATE_BURNING

        for row_offset, col_offset, neighbor_direction_deg in NEIGHBOR_OFFSETS:
            shifted_burning_mask = self._shift_grid(burning_mask, row_offset, col_offset)

            slope_factor = self._compute_slope_factor(row_offset, col_offset)
            wind_factor = self._compute_wind_factor(
                wind_speed_grid, wind_direction_grid, neighbor_direction_deg
            )
            moisture_factor = self._compute_moisture_factor(fuel_moisture_grid)

            neighbor_spread_probability = (
                self.base_spread_probability
                * (1.0 + self.slope_weight * slope_factor)
                * (1.0 + self.wind_weight * wind_factor)
                * (1.0 - self.moisture_weight * moisture_factor)
            )
            neighbor_spread_probability = np.clip(neighbor_spread_probability, 0.0, 1.0)

            contribution = np.where(shifted_burning_mask, neighbor_spread_probability, 0.0)
            ignition_probability_total = 1.0 - (1.0 - ignition_probability_total) * (1.0 - contribution)

        ignition_probability_total = np.clip(ignition_probability_total, 0.0, 1.0)

        unburned_mask = current_state == STATE_UNBURNED
        random_draws = random_generator.random(size=current_state.shape).astype("float32")
        newly_ignited_mask = unburned_mask & (random_draws < ignition_probability_total)

        new_state[newly_ignited_mask] = STATE_BURNING
        new_state[burning_mask] = STATE_BURNED_OUT

        return new_state, ignition_probability_total

    def _shift_grid(self, grid: np.ndarray, row_offset: int, col_offset: int) -> np.ndarray:
        """Shifts a boolean grid so that grid[row, col] reflects the neighbor at (row + row_offset, col + col_offset)."""
        shifted = np.roll(grid, shift=(row_offset, col_offset), axis=(0, 1))

        if row_offset > 0:
            shifted[:row_offset, :] = False
        elif row_offset < 0:
            shifted[row_offset:, :] = False

        if col_offset > 0:
            shifted[:, :col_offset] = False
        elif col_offset < 0:
            shifted[:, col_offset:] = False

        return shifted

    def _compute_slope_factor(self, row_offset: int, col_offset: int) -> np.ndarray:
        """
        Fire spreads faster uphill. The factor is positive when the neighbor direction
        points uphill relative to the cell aspect, scaled by slope steepness.
        """
        neighbor_direction_deg = np.degrees(np.arctan2(col_offset, -row_offset)) % 360.0
        aspect_difference = np.abs(self.aspect_grid - neighbor_direction_deg)
        aspect_difference = np.minimum(aspect_difference, 360.0 - aspect_difference)

        alignment = np.cos(np.radians(aspect_difference))
        slope_radians = np.radians(self.slope_grid)

        slope_factor = alignment * np.tan(slope_radians)
        slope_factor = np.clip(slope_factor, -1.0, 5.0).astype("float32")
        return slope_factor

    def _compute_wind_factor(
        self,
        wind_speed_grid: np.ndarray,
        wind_direction_grid: np.ndarray,
        neighbor_direction_deg: float,
    ) -> np.ndarray:
        """
        Wind blowing toward the neighbor cell increases spread probability,
        scaled by wind speed normalized against a reference speed.
        """
        direction_difference = np.abs(wind_direction_grid - neighbor_direction_deg)
        direction_difference = np.minimum(direction_difference, 360.0 - direction_difference)

        alignment = np.cos(np.radians(direction_difference))

        reference_wind_speed_ms = 10.0
        normalized_speed = np.clip(wind_speed_grid / reference_wind_speed_ms, 0.0, 2.0)

        wind_factor = alignment * normalized_speed
        return wind_factor.astype("float32")

    def _compute_moisture_factor(self, fuel_moisture_grid: np.ndarray) -> np.ndarray:
        """Returns the fuel moisture coefficient directly, clipped to a valid range."""
        return np.clip(fuel_moisture_grid, 0.0, 1.0).astype("float32")


# ----------------------------------------------------------------------
# MCMC calibration of Cellular Automata parameters
# ----------------------------------------------------------------------

@dataclass
class CalibrationDataset:
    """Holds the static grids and time series needed to run and score the Cellular Automata model."""
    slope_grid: np.ndarray
    aspect_grid: np.ndarray
    wind_speed_series: np.ndarray
    wind_direction_series: np.ndarray
    fuel_moisture_series: np.ndarray
    initial_state: np.ndarray
    observed_final_state: np.ndarray


class MCMCCalibrator:
    """Calibrates Cellular Automata parameters using MCMC against observed FIRMS burned area."""

    PARAMETER_NAMES = ["base_spread_probability", "slope_weight", "wind_weight", "moisture_weight"]

    PARAMETER_BOUNDS = {
        "base_spread_probability": (0.01, 0.5),
        "slope_weight": (0.0, 3.0),
        "wind_weight": (0.0, 3.0),
        "moisture_weight": (0.0, 1.0),
    }

    def __init__(self, config: PipelineConfig, calibration_dataset: CalibrationDataset):
        self.config = config
        self.dataset = calibration_dataset

    def log_prior(self, parameters: np.ndarray) -> float:
        for value, name in zip(parameters, self.PARAMETER_NAMES):
            lower_bound, upper_bound = self.PARAMETER_BOUNDS[name]
            if not (lower_bound <= value <= upper_bound):
                return -np.inf
        return 0.0

    def log_likelihood(self, parameters: np.ndarray, random_seed: int = 42) -> float:
        base_spread_probability, slope_weight, wind_weight, moisture_weight = parameters

        cellular_automata = FireSpreadCellularAutomata(
            slope_grid=self.dataset.slope_grid,
            aspect_grid=self.dataset.aspect_grid,
            base_spread_probability=base_spread_probability,
            slope_weight=slope_weight,
            wind_weight=wind_weight,
            moisture_weight=moisture_weight,
        )

        random_generator = np.random.default_rng(random_seed)
        current_state = self.dataset.initial_state.copy()

        hours = self.dataset.wind_speed_series.shape[0]
        for hour_index in range(hours):
            current_state, _ = cellular_automata.step(
                current_state,
                self.dataset.wind_speed_series[hour_index],
                self.dataset.wind_direction_series[hour_index],
                self.dataset.fuel_moisture_series[hour_index],
                random_generator,
            )

        simulated_burned_mask = (current_state != STATE_UNBURNED).astype("uint8")
        observed_burned_mask = (self.dataset.observed_final_state != STATE_UNBURNED).astype("uint8")

        true_positive = np.sum((simulated_burned_mask == 1) & (observed_burned_mask == 1))
        false_positive = np.sum((simulated_burned_mask == 1) & (observed_burned_mask == 0))
        false_negative = np.sum((simulated_burned_mask == 0) & (observed_burned_mask == 1))

        epsilon = 1e-6
        intersection_over_union = true_positive / (true_positive + false_positive + false_negative + epsilon)

        sigma = 0.05
        log_likelihood_value = -0.5 * ((1.0 - intersection_over_union) ** 2) / (sigma ** 2)

        return float(log_likelihood_value)

    def log_posterior(self, parameters: np.ndarray) -> float:
        prior_value = self.log_prior(parameters)
        if not np.isfinite(prior_value):
            return -np.inf
        return prior_value + self.log_likelihood(parameters)

    def run_calibration(self) -> Dict[str, float]:
        number_of_parameters = len(self.PARAMETER_NAMES)
        number_of_walkers = self.config.mcmc_walkers

        initial_guess = np.array([0.15, 1.0, 1.0, 0.3])
        initial_positions = initial_guess + 1e-2 * np.random.randn(number_of_walkers, number_of_parameters)

        for walker_index in range(number_of_walkers):
            for parameter_index, parameter_name in enumerate(self.PARAMETER_NAMES):
                lower_bound, upper_bound = self.PARAMETER_BOUNDS[parameter_name]
                initial_positions[walker_index, parameter_index] = np.clip(
                    initial_positions[walker_index, parameter_index], lower_bound, upper_bound
                )

        sampler = emcee.EnsembleSampler(number_of_walkers, number_of_parameters, self.log_posterior)

        log.info("Starting MCMC calibration with %d walkers for %d steps", number_of_walkers, self.config.mcmc_steps)
        sampler.run_mcmc(initial_positions, self.config.mcmc_steps, progress=False)

        flattened_chain = sampler.get_chain(discard=self.config.mcmc_burn_in, flat=True)
        posterior_means = np.mean(flattened_chain, axis=0)

        calibrated_parameters = {
            name: float(value) for name, value in zip(self.PARAMETER_NAMES, posterior_means)
        }

        log.info("MCMC calibration complete: %s", calibrated_parameters)
        return calibrated_parameters


# ----------------------------------------------------------------------
# Hourly fire spread simulation runner
# ----------------------------------------------------------------------

class FireSpreadSimulator:
    """Runs the calibrated Cellular Automata forward in time and accumulates ignition probabilities."""

    def __init__(
        self,
        slope_grid: np.ndarray,
        aspect_grid: np.ndarray,
        calibrated_parameters: Dict[str, float],
    ):
        self.cellular_automata = FireSpreadCellularAutomata(
            slope_grid=slope_grid,
            aspect_grid=aspect_grid,
            base_spread_probability=calibrated_parameters["base_spread_probability"],
            slope_weight=calibrated_parameters["slope_weight"],
            wind_weight=calibrated_parameters["wind_weight"],
            moisture_weight=calibrated_parameters["moisture_weight"],
        )

    def run_simulation(
        self,
        initial_state: np.ndarray,
        meteorology_loader: MeteorologyLoader,
        hours_to_simulate: int,
        starting_hour_index: int = 0,
        random_seed: int = 7,
    ) -> Tuple[np.ndarray, np.ndarray]:
        random_generator = np.random.default_rng(random_seed)
        current_state = initial_state.copy()
        target_shape = initial_state.shape
        cumulative_probability_grid = np.zeros_like(initial_state, dtype="float32")

        for hour_offset in range(hours_to_simulate):
            hour_index = starting_hour_index + hour_offset
            meteorology = meteorology_loader.get_hourly_wind_and_moisture(hour_index)

            wind_speed = self._resample_to_shape(meteorology["wind_speed_ms"], target_shape)
            wind_direction = self._resample_to_shape(meteorology["wind_direction_deg"], target_shape)
            fuel_moisture = self._resample_to_shape(meteorology["fuel_moisture_coefficient"], target_shape)

            current_state, ignition_probability_grid = self.cellular_automata.step(
                current_state, wind_speed, wind_direction, fuel_moisture, random_generator,
            )
            cumulative_probability_grid = np.maximum(cumulative_probability_grid, ignition_probability_grid)
            log.info("Hour %d simulated, burning cells: %d", hour_index, int(np.sum(current_state == STATE_BURNING)))

        return current_state, cumulative_probability_grid

    def _resample_to_shape(self, source_array: np.ndarray, target_shape: Tuple[int, int]) -> np.ndarray:
        if source_array.shape == target_shape:
            return source_array.astype("float32")
        row_indices = np.linspace(0, source_array.shape[0] - 1, target_shape[0]).astype("int64")
        col_indices = np.linspace(0, source_array.shape[1] - 1, target_shape[1]).astype("int64")
        return source_array[np.ix_(row_indices, col_indices)].astype("float32")


# ----------------------------------------------------------------------
# Threat level mapping for evacuation protocols
# ----------------------------------------------------------------------

THREAT_SAFE = 0
THREAT_LOW = 1
THREAT_MODERATE = 2
THREAT_HIGH = 3
THREAT_CRITICAL = 4

THREAT_THRESHOLDS = {
    THREAT_SAFE: (0.0, 0.05),
    THREAT_LOW: (0.05, 0.2),
    THREAT_MODERATE: (0.2, 0.5),
    THREAT_HIGH: (0.5, 0.8),
    THREAT_CRITICAL: (0.8, 1.0001),
}


class ThreatMapGenerator:
    """Generates semantically segmented threat level maps and writes them to GeoTIFF."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def classify_threat_levels(
        self, cumulative_probability_grid: np.ndarray, final_state_grid: np.ndarray
    ) -> np.ndarray:
        threat_grid = np.zeros_like(cumulative_probability_grid, dtype="uint8")

        for threat_level, (lower_bound, upper_bound) in THREAT_THRESHOLDS.items():
            within_range_mask = (cumulative_probability_grid >= lower_bound) & (
                cumulative_probability_grid < upper_bound
            )
            threat_grid[within_range_mask] = threat_level

        already_affected_mask = final_state_grid != STATE_UNBURNED
        threat_grid[already_affected_mask] = THREAT_CRITICAL

        return threat_grid

    def write_threat_map(
        self,
        threat_grid: np.ndarray,
        cumulative_probability_grid: np.ndarray,
        reference_raster_path: str,
        output_threat_path: str,
        output_probability_path: str,
    ) -> None:
        with rasterio.open(reference_raster_path) as reference_dataset:
            profile = reference_dataset.profile.copy()

        threat_profile = profile.copy()
        threat_profile.update(dtype="uint8", count=1, nodata=255)

        probability_profile = profile.copy()
        probability_profile.update(dtype="float32", count=1, nodata=self.config.no_data_value)

        with rasterio.open(output_threat_path, "w", **threat_profile) as threat_dataset:
            threat_dataset.write(threat_grid, 1)

        with rasterio.open(output_probability_path, "w", **probability_profile) as probability_dataset:
            probability_dataset.write(cumulative_probability_grid, 1)

        log.info("Threat map written to %s", output_threat_path)
        log.info("Probability map written to %s", output_probability_path)

    def generate_evacuation_summary(
        self, threat_grid: np.ndarray, transform: rasterio.Affine
    ) -> Dict[str, object]:
        """Produces a JSON serializable summary of threat zones for emergency planning."""
        cell_area_square_meters = abs(transform.a * transform.e)

        summary = {"threat_levels": {}}
        threat_level_names = {
            THREAT_SAFE: "safe",
            THREAT_LOW: "low",
            THREAT_MODERATE: "moderate",
            THREAT_HIGH: "high",
            THREAT_CRITICAL: "critical",
        }

        for threat_level, threat_name in threat_level_names.items():
            cell_count = int(np.sum(threat_grid == threat_level))
            area_square_kilometers = (cell_count * cell_area_square_meters) / 1_000_000.0
            summary["threat_levels"][threat_name] = {
                "cell_count": cell_count,
                "area_square_kilometers": round(area_square_kilometers, 4),
            }

        critical_and_high_cells = int(
            np.sum((threat_grid == THREAT_CRITICAL) | (threat_grid == THREAT_HIGH))
        )
        summary["recommended_evacuation_cell_count"] = critical_and_high_cells

        return summary


# ----------------------------------------------------------------------
# Pipeline orchestration
# ----------------------------------------------------------------------

class HimalayanFirePredictionPipeline:
    """Top level orchestrator wiring together data ingestion, calibration, simulation, and threat mapping."""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def run(self) -> Dict[str, object]:
        topography_loader = TopographyLoader(self.config)
        slope_path = f"{self.config.output_dir}/slope.tif"
        aspect_path = f"{self.config.output_dir}/aspect.tif"
        topography_loader.compute_slope_aspect(slope_path, aspect_path)

        with rasterio.open(slope_path) as slope_dataset:
            slope_grid = slope_dataset.read(1)
            raster_transform = slope_dataset.transform

        with rasterio.open(aspect_path) as aspect_dataset:
            aspect_grid = aspect_dataset.read(1)

        meteorology_loader = MeteorologyLoader(self.config)

        active_fire_loader = ActiveFireLoader(self.config)
        active_fire_points = active_fire_loader.load_active_fire_points()
        initial_state = active_fire_loader.rasterize_initial_state(active_fire_points, slope_path)
        initial_state = initial_state.astype("uint8")
        initial_state[initial_state == 1] = STATE_BURNING

        calibration_hours = min(6, self.config.hours_to_simulate)
        wind_speed_series = np.zeros((calibration_hours, *slope_grid.shape), dtype="float32")
        wind_direction_series = np.zeros((calibration_hours, *slope_grid.shape), dtype="float32")
        fuel_moisture_series = np.zeros((calibration_hours, *slope_grid.shape), dtype="float32")

        for hour_index in range(calibration_hours):
            meteorology = meteorology_loader.get_hourly_wind_and_moisture(hour_index)
            wind_speed_series[hour_index] = self._resample_to_shape(
                meteorology["wind_speed_ms"], slope_grid.shape
            )
            wind_direction_series[hour_index] = self._resample_to_shape(
                meteorology["wind_direction_deg"], slope_grid.shape
            )
            fuel_moisture_series[hour_index] = self._resample_to_shape(
                meteorology["fuel_moisture_coefficient"], slope_grid.shape
            )

        observed_final_state = initial_state.copy()
        observed_final_state[initial_state == STATE_BURNING] = STATE_BURNED_OUT

        calibration_dataset = CalibrationDataset(
            slope_grid=slope_grid,
            aspect_grid=aspect_grid,
            wind_speed_series=wind_speed_series,
            wind_direction_series=wind_direction_series,
            fuel_moisture_series=fuel_moisture_series,
            initial_state=initial_state,
            observed_final_state=observed_final_state,
        )

        calibrator = MCMCCalibrator(self.config, calibration_dataset)
        calibrated_parameters = calibrator.run_calibration()

        simulator = FireSpreadSimulator(slope_grid, aspect_grid, calibrated_parameters)
        final_state_grid, cumulative_probability_grid = simulator.run_simulation(
            initial_state, meteorology_loader, self.config.hours_to_simulate
        )

        threat_map_generator = ThreatMapGenerator(self.config)
        threat_grid = threat_map_generator.classify_threat_levels(cumulative_probability_grid, final_state_grid)

        output_threat_path = f"{self.config.output_dir}/threat_map.tif"
        output_probability_path = f"{self.config.output_dir}/probability_map.tif"
        threat_map_generator.write_threat_map(
            threat_grid, cumulative_probability_grid, slope_path, output_threat_path, output_probability_path
        )

        evacuation_summary = threat_map_generator.generate_evacuation_summary(threat_grid, raster_transform)

        result = {
            "calibrated_parameters": calibrated_parameters,
            "evacuation_summary": evacuation_summary,
            "threat_map_path": output_threat_path,
            "probability_map_path": output_probability_path,
        }

        summary_path = f"{self.config.output_dir}/run_summary.json"
        with open(summary_path, "w") as summary_file:
            json.dump(result, summary_file, indent=2)

        log.info("Pipeline run complete, summary written to %s", summary_path)
        return result

    def _resample_to_shape(self, source_array: np.ndarray, target_shape: Tuple[int, int]) -> np.ndarray:
        """Resamples a 2D array to the target shape using nearest neighbor index selection for memory efficiency."""
        if source_array.shape == target_shape:
            return source_array.astype("float32")

        row_indices = np.linspace(0, source_array.shape[0] - 1, target_shape[0]).astype("int64")
        col_indices = np.linspace(0, source_array.shape[1] - 1, target_shape[1]).astype("int64")

        resampled = source_array[np.ix_(row_indices, col_indices)]
        return resampled.astype("float32")


# ----------------------------------------------------------------------
# Entry point
# ----------------------------------------------------------------------

def main() -> None:
    config = PipelineConfig(
        dem_path="data/srtm_himalaya_250m.tif",
        era5_path="data/era5_land_hourly.nc",
        firms_csv_path="data/firms_viirs_active_fire.csv",
        output_dir="output",
        block_size=256,
        cell_size_m=240.0,
        hours_to_simulate=24,
        mcmc_walkers=8,
        mcmc_steps=200,
        mcmc_burn_in=50,
    )

    pipeline = HimalayanFirePredictionPipeline(config)

    # Compute slope and aspect
    topography_loader = TopographyLoader(config)
    slope_path = f"{config.output_dir}/slope.tif"
    aspect_path = f"{config.output_dir}/aspect.tif"
    topography_loader.compute_slope_aspect(slope_path, aspect_path)

    with rasterio.open(slope_path) as slope_dataset:
        slope_grid = slope_dataset.read(1)
        raster_transform = slope_dataset.transform

    with rasterio.open(aspect_path) as aspect_dataset:
        aspect_grid = aspect_dataset.read(1)

    log.info("Grid size: %d x %d pixels", slope_grid.shape[0], slope_grid.shape[1])

    # Seed synthetic fire points in the center of the grid since FIRMS has no NRT data
    # for historical dates. In production replace this with real FIRMS detections.
    initial_state = np.zeros(slope_grid.shape, dtype="uint8")
    center_row = slope_grid.shape[0] // 2
    center_col = slope_grid.shape[1] // 2
    seed_radius = 3
    initial_state[
        center_row - seed_radius: center_row + seed_radius,
        center_col - seed_radius: center_col + seed_radius,
    ] = STATE_BURNING
    log.info("Seeded %d burning cells", int(np.sum(initial_state == STATE_BURNING)))

    # Load meteorology
    meteorology_loader = MeteorologyLoader(config)

    # Use a small subgrid for MCMC calibration to keep memory low
    calibration_size = 200
    row_start = center_row - calibration_size // 2
    col_start = center_col - calibration_size // 2
    row_end = row_start + calibration_size
    col_end = col_start + calibration_size

    slope_subgrid = slope_grid[row_start:row_end, col_start:col_end]
    aspect_subgrid = aspect_grid[row_start:row_end, col_start:col_end]
    initial_state_subgrid = initial_state[row_start:row_end, col_start:col_end]

    calibration_hours = 6
    wind_speed_series = np.zeros((calibration_hours, calibration_size, calibration_size), dtype="float32")
    wind_direction_series = np.zeros((calibration_hours, calibration_size, calibration_size), dtype="float32")
    fuel_moisture_series = np.zeros((calibration_hours, calibration_size, calibration_size), dtype="float32")

    for hour_index in range(calibration_hours):
        meteorology = meteorology_loader.get_hourly_wind_and_moisture(hour_index)
        wind_speed_series[hour_index] = pipeline._resample_to_shape(
            meteorology["wind_speed_ms"], (calibration_size, calibration_size)
        )
        wind_direction_series[hour_index] = pipeline._resample_to_shape(
            meteorology["wind_direction_deg"], (calibration_size, calibration_size)
        )
        fuel_moisture_series[hour_index] = pipeline._resample_to_shape(
            meteorology["fuel_moisture_coefficient"], (calibration_size, calibration_size)
        )

    observed_final_state_subgrid = initial_state_subgrid.copy()
    observed_final_state_subgrid[initial_state_subgrid == STATE_BURNING] = STATE_BURNED_OUT

    calibration_dataset = CalibrationDataset(
        slope_grid=slope_subgrid,
        aspect_grid=aspect_subgrid,
        wind_speed_series=wind_speed_series,
        wind_direction_series=wind_direction_series,
        fuel_moisture_series=fuel_moisture_series,
        initial_state=initial_state_subgrid,
        observed_final_state=observed_final_state_subgrid,
    )

    calibrator = MCMCCalibrator(config, calibration_dataset)
    calibrated_parameters = calibrator.run_calibration()

    # Run full simulation on the complete downsampled grid
    simulator = FireSpreadSimulator(slope_grid, aspect_grid, calibrated_parameters)
    final_state_grid, cumulative_probability_grid = simulator.run_simulation(
        initial_state, meteorology_loader, config.hours_to_simulate
    )

    threat_map_generator = ThreatMapGenerator(config)
    threat_grid = threat_map_generator.classify_threat_levels(
        cumulative_probability_grid, final_state_grid
    )

    output_threat_path = f"{config.output_dir}/threat_map.tif"
    output_probability_path = f"{config.output_dir}/probability_map.tif"
    threat_map_generator.write_threat_map(
        threat_grid, cumulative_probability_grid, slope_path,
        output_threat_path, output_probability_path
    )

    evacuation_summary = threat_map_generator.generate_evacuation_summary(
        threat_grid, raster_transform
    )

    result = {
        "calibrated_parameters": calibrated_parameters,
        "evacuation_summary": evacuation_summary,
        "threat_map_path": output_threat_path,
        "probability_map_path": output_probability_path,
    }

    summary_path = f"{config.output_dir}/run_summary.json"
    with open(summary_path, "w") as summary_file:
        json.dump(result, summary_file, indent=2)

    log.info("Pipeline complete. Summary: %s", json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

if __name__ == "__main__":
    main()
