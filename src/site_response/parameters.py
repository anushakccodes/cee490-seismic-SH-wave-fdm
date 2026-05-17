"""Parameters and Stage 1 checks for the 1D SH-wave site-response model."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class SimulationConfig:
    """Numerical and physical settings for the finite-difference model."""
    model_depth_m: float = 400.0
    grid_spacing_m: float = 2.0
    time_step_s: float = 0.00125
    total_time_s: float = 1.8
    dominant_frequency_hz: float = 4.0
    maximum_interpreted_frequency_hz: float = 12.0
    soil_layer_thickness_m: float = 50.0
    soil_density_kg_m3: float = 1800.0
    soil_shear_velocity_m_s: float = 250.0
    rock_density_kg_m3: float = 2200.0
    rock_shear_velocity_m_s: float = 800.0
    absorbing_layer_thickness_m: float = 80.0
    absorbing_strength: float = 3.5
    project_root: Path = Path(__file__).resolve().parents[2]

    @property
    def number_of_depth_nodes(self) -> int:
        return int(round(self.model_depth_m / self.grid_spacing_m)) + 1

    @property
    def number_of_time_steps(self) -> int:
        return int(round(self.total_time_s / self.time_step_s)) + 1

    @property
    def source_peak_time_s(self) -> float:
        return 1.5 / self.dominant_frequency_hz

    @property
    def maximum_shear_velocity_m_s(self) -> float:
        return max(self.soil_shear_velocity_m_s, self.rock_shear_velocity_m_s)

    @property
    def minimum_shear_velocity_m_s(self) -> float:
        return min(self.soil_shear_velocity_m_s, self.rock_shear_velocity_m_s)

    @property
    def figures_dir(self) -> Path:
        return self.project_root / "figures"

    @property
    def results_dir(self) -> Path:
        return self.project_root / "results"


def default_config() -> SimulationConfig:
    """Return the baseline configuration for Stages 1--3."""
    return SimulationConfig()


def make_depth_axis(config: SimulationConfig) -> np.ndarray:
    """Create the depth axis with z = 0 at the free surface."""
    return np.linspace(0.0, config.model_depth_m, config.number_of_depth_nodes)


def make_time_axis(config: SimulationConfig) -> np.ndarray:
    """Create the simulation time axis."""
    return np.linspace(0.0, config.total_time_s, config.number_of_time_steps)


def compute_parameter_checks(config: SimulationConfig) -> dict[str, float | bool]:
    """Compute CFL, grid-spacing, frequency-resolution, and Nyquist checks."""
    courant_number = config.maximum_shear_velocity_m_s * config.time_step_s / config.grid_spacing_m
    shortest_wavelength_m = config.minimum_shear_velocity_m_s / config.maximum_interpreted_frequency_hz
    points_per_shortest_wavelength = shortest_wavelength_m / config.grid_spacing_m
    return {
        "courant_number": courant_number,
        "cfl_pass": courant_number <= 1.0,
        "cfl_preferred_pass": courant_number <= 0.8,
        "shortest_wavelength_m": shortest_wavelength_m,
        "points_per_shortest_wavelength": points_per_shortest_wavelength,
        "grid_resolution_pass": points_per_shortest_wavelength >= 10.0,
        "frequency_resolution_hz": 1.0 / config.total_time_s,
        "nyquist_frequency_hz": 1.0 / (2.0 * config.time_step_s),
    }


def format_parameter_summary(config: SimulationConfig) -> str:
    """Return a readable Stage 1 parameter summary."""
    checks = compute_parameter_checks(config)
    soil_mu = config.soil_density_kg_m3 * config.soil_shear_velocity_m_s**2
    rock_mu = config.rock_density_kg_m3 * config.rock_shear_velocity_m_s**2
    return (
        "CEE 490 finite-difference site-response model\n"
        + "=" * 58 + "\n"
        + f"Depth domain:                  {config.model_depth_m:.1f} m\n"
        + f"Grid spacing dz:               {config.grid_spacing_m:.3f} m\n"
        + f"Time step dt:                  {config.time_step_s:.6f} s\n"
        + f"Total duration T:              {config.total_time_s:.3f} s\n"
        + f"Ricker dominant frequency f0:  {config.dominant_frequency_hz:.3f} Hz\n"
        + f"Soil Vs, rho, mu:              {config.soil_shear_velocity_m_s:.1f} m/s, {config.soil_density_kg_m3:.1f} kg/m^3, {soil_mu:.3e} Pa\n"
        + f"Rock Vs, rho, mu:              {config.rock_shear_velocity_m_s:.1f} m/s, {config.rock_density_kg_m3:.1f} kg/m^3, {rock_mu:.3e} Pa\n"
        + f"CFL number C:                  {checks['courant_number']:.3f} ({'PASS' if checks['cfl_pass'] else 'FAIL'})\n"
        + f"Points / shortest wavelength:  {checks['points_per_shortest_wavelength']:.2f} ({'PASS' if checks['grid_resolution_pass'] else 'CHECK'})\n"
        + f"Frequency resolution df:       {checks['frequency_resolution_hz']:.3f} Hz\n"
        + f"Nyquist frequency:             {checks['nyquist_frequency_hz']:.1f} Hz\n"
    )


def ensure_output_directories(config: SimulationConfig) -> None:
    """Create output folders used by scripts and notebooks."""
    config.figures_dir.mkdir(parents=True, exist_ok=True)
    config.results_dir.mkdir(parents=True, exist_ok=True)
