"""2D material models for SH-wave finite-difference simulations."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class Simulation2DConfig:
    """Numerical and physical settings for the 2D SH-wave FD model."""
    model_width_m: float = 800.0
    model_depth_m: float = 400.0
    dx_m: float = 4.0
    dz_m: float = 4.0
    time_step_s: float = 0.0015
    total_time_s: float = 2.0
    dominant_frequency_hz: float = 4.0
    soil_layer_thickness_m: float = 50.0
    soil_density_kg_m3: float = 1800.0
    soil_shear_velocity_m_s: float = 250.0
    rock_density_kg_m3: float = 2200.0
    rock_shear_velocity_m_s: float = 800.0
    absorbing_layer_thickness_m: float = 80.0
    absorbing_strength: float = 3.5
    project_root: Path = Path(__file__).resolve().parents[2]

    @property
    def nx(self) -> int:
        return int(round(self.model_width_m / self.dx_m)) + 1

    @property
    def nz(self) -> int:
        return int(round(self.model_depth_m / self.dz_m)) + 1

    @property
    def number_of_time_steps(self) -> int:
        return int(round(self.total_time_s / self.time_step_s)) + 1

    @property
    def source_peak_time_s(self) -> float:
        return 1.5 / self.dominant_frequency_hz

    @property
    def figures_2d_dir(self) -> Path:
        return self.project_root / "figures_2d"

    @property
    def results_2d_dir(self) -> Path:
        return self.project_root / "results_2d"


@dataclass(frozen=True)
class Material2D:
    """2D arrays of density, shear velocity, and shear modulus."""
    x_m: np.ndarray
    z_m: np.ndarray
    density_kg_m3: np.ndarray
    shear_velocity_m_s: np.ndarray
    shear_modulus_pa: np.ndarray
    label: str
    interface_depth_m: float | None = None

    @property
    def impedance_kg_m2_s(self) -> np.ndarray:
        return self.density_kg_m3 * self.shear_velocity_m_s


def make_2d_axes(config: Simulation2DConfig) -> tuple[np.ndarray, np.ndarray]:
    """Return (x_m, z_m) 1D coordinate arrays."""
    x_m = np.linspace(0.0, config.model_width_m, config.nx)
    z_m = np.linspace(0.0, config.model_depth_m, config.nz)
    return x_m, z_m


def _shear_modulus(density: np.ndarray, vs: np.ndarray) -> np.ndarray:
    return density * vs**2


def homogeneous_2d_model(config: Simulation2DConfig) -> Material2D:
    """Return a uniform rock model for solver verification."""
    x_m, z_m = make_2d_axes(config)
    density = np.full((config.nz, config.nx), config.rock_density_kg_m3, dtype=float)
    vs = np.full((config.nz, config.nx), config.rock_shear_velocity_m_s, dtype=float)
    return Material2D(
        x_m=x_m, z_m=z_m,
        density_kg_m3=density,
        shear_velocity_m_s=vs,
        shear_modulus_pa=_shear_modulus(density, vs),
        label="Homogeneous rock",
    )


def layered_2d_model(config: Simulation2DConfig) -> Material2D:
    """Return a horizontally layered soil-over-rock model."""
    x_m, z_m = make_2d_axes(config)
    # Build full (nz, nx) boolean mask; z increases downward
    in_soil = (z_m[:, np.newaxis] <= config.soil_layer_thickness_m) * np.ones((1, config.nx), dtype=bool)
    density = np.where(in_soil, config.soil_density_kg_m3, config.rock_density_kg_m3).astype(float)
    vs = np.where(in_soil, config.soil_shear_velocity_m_s, config.rock_shear_velocity_m_s).astype(float)
    return Material2D(
        x_m=x_m, z_m=z_m,
        density_kg_m3=density,
        shear_velocity_m_s=vs,
        shear_modulus_pa=_shear_modulus(density, vs),
        label="Layered soil over rock",
        interface_depth_m=config.soil_layer_thickness_m,
    )


def basin_2d_model(
    config: Simulation2DConfig,
    base_depth_m: float = 20.0,
    extra_depth_m: float = 60.0,
    basin_half_width_m: float = 200.0,
) -> Material2D:
    """Return a Gaussian-shaped soft basin embedded in rock.

    Basin depth varies as: base_depth + extra_depth * exp(-((x - x_center)/half_width)^2)
    """
    x_m, z_m = make_2d_axes(config)
    x_center = config.model_width_m / 2.0
    basin_depth = base_depth_m + extra_depth_m * np.exp(-((x_m - x_center) / basin_half_width_m) ** 2)

    # Build full (nz, nx) boolean mask: True where z_m[iz] <= basin_depth[ix]
    in_soil = z_m[:, np.newaxis] <= basin_depth[np.newaxis, :]

    density = np.where(in_soil, config.soil_density_kg_m3, config.rock_density_kg_m3).astype(float)
    vs = np.where(in_soil, config.soil_shear_velocity_m_s, config.rock_shear_velocity_m_s).astype(float)
    return Material2D(
        x_m=x_m, z_m=z_m,
        density_kg_m3=density,
        shear_velocity_m_s=vs,
        shear_modulus_pa=_shear_modulus(density, vs),
        label="Gaussian basin",
    )


def stress_grid_mu_x(material: Material2D) -> np.ndarray:
    """Harmonic-average shear modulus at x-staggered stress locations.

    Returns array of shape (nz, nx-1).
    """
    mu = material.shear_modulus_pa
    mu_left = mu[:, :-1]
    mu_right = mu[:, 1:]
    return 2.0 * mu_left * mu_right / (mu_left + mu_right)


def stress_grid_mu_z(material: Material2D) -> np.ndarray:
    """Harmonic-average shear modulus at z-staggered stress locations.

    Returns array of shape (nz-1, nx).
    """
    mu = material.shear_modulus_pa
    mu_top = mu[:-1, :]
    mu_bot = mu[1:, :]
    return 2.0 * mu_top * mu_bot / (mu_top + mu_bot)


def cfl_check_2d(config: Simulation2DConfig) -> dict[str, float | bool]:
    """Return 2D CFL diagnostics."""
    vs_max = config.rock_shear_velocity_m_s
    dt = config.time_step_s
    dx = config.dx_m
    dz = config.dz_m
    courant = vs_max * dt * np.sqrt(1.0 / dx**2 + 1.0 / dz**2)
    return {
        "vs_max_m_s": vs_max,
        "courant_2d": courant,
        "cfl_pass": bool(courant <= 1.0),
        "cfl_preferred_pass": bool(courant <= 0.6),
        "dt_max_s": dx / (vs_max * np.sqrt(2.0)),
    }


def format_2d_parameter_summary(config: Simulation2DConfig) -> str:
    """Return a readable parameter summary for the 2D model."""
    cfl = cfl_check_2d(config)
    soil_mu = config.soil_density_kg_m3 * config.soil_shear_velocity_m_s**2
    rock_mu = config.rock_density_kg_m3 * config.rock_shear_velocity_m_s**2
    return (
        "CEE 490  2D SH-wave finite-difference model\n"
        + "=" * 58 + "\n"
        + f"Domain:            {config.model_width_m:.0f} m wide x {config.model_depth_m:.0f} m deep\n"
        + f"Grid spacing dx/dz: {config.dx_m:.2f} m / {config.dz_m:.2f} m\n"
        + f"Grid nodes (nx,nz): {config.nx} x {config.nz}\n"
        + f"Time step dt:       {config.time_step_s:.6f} s\n"
        + f"Total duration T:   {config.total_time_s:.3f} s\n"
        + f"Number of steps:    {config.number_of_time_steps}\n"
        + f"Ricker f0:          {config.dominant_frequency_hz:.2f} Hz\n"
        + f"Soil Vs, rho, mu:   {config.soil_shear_velocity_m_s:.0f} m/s, {config.soil_density_kg_m3:.0f} kg/m^3, {soil_mu:.3e} Pa\n"
        + f"Rock Vs, rho, mu:   {config.rock_shear_velocity_m_s:.0f} m/s, {config.rock_density_kg_m3:.0f} kg/m^3, {rock_mu:.3e} Pa\n"
        + f"2D CFL number:      {cfl['courant_2d']:.4f} ({'PASS' if cfl['cfl_pass'] else 'FAIL'})"
        + (" [preferred]" if cfl['cfl_preferred_pass'] else " [marginal — consider reducing dt]") + "\n"
        + f"Max stable dt:      {cfl['dt_max_s']:.6f} s\n"
    )


def ensure_2d_output_directories(config: Simulation2DConfig) -> None:
    """Create figures_2d/ and results_2d/ directories."""
    config.figures_2d_dir.mkdir(parents=True, exist_ok=True)
    config.results_2d_dir.mkdir(parents=True, exist_ok=True)
