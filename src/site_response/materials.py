"""Material models for homogeneous and layered 1D SH-wave simulations."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .parameters import SimulationConfig, make_depth_axis


@dataclass(frozen=True)
class MaterialProfile:
    """Depth-dependent density, shear velocity, and shear modulus."""
    depth_m: np.ndarray
    density_kg_m3: np.ndarray
    shear_velocity_m_s: np.ndarray
    shear_modulus_pa: np.ndarray
    label: str
    interface_depth_m: float | None = None

    @property
    def impedance_kg_m2_s(self) -> np.ndarray:
        """Return shear impedance Z = rho * Vs."""
        return self.density_kg_m3 * self.shear_velocity_m_s


def compute_shear_modulus(density_kg_m3: np.ndarray, shear_velocity_m_s: np.ndarray) -> np.ndarray:
    """Compute shear modulus using mu = rho * Vs**2."""
    return density_kg_m3 * shear_velocity_m_s**2


def homogeneous_rock_model(config: SimulationConfig) -> MaterialProfile:
    """Return the homogeneous rock model used for solver verification."""
    depth_m = make_depth_axis(config)
    density = np.full_like(depth_m, config.rock_density_kg_m3, dtype=float)
    shear_velocity = np.full_like(depth_m, config.rock_shear_velocity_m_s, dtype=float)
    return MaterialProfile(depth_m, density, shear_velocity, compute_shear_modulus(density, shear_velocity), "Homogeneous rock reference")


def layered_soil_over_rock_model(config: SimulationConfig) -> MaterialProfile:
    """Return a soft soil layer over a stiffer rock halfspace."""
    depth_m = make_depth_axis(config)
    in_soil = depth_m <= config.soil_layer_thickness_m
    density = np.where(in_soil, config.soil_density_kg_m3, config.rock_density_kg_m3).astype(float)
    shear_velocity = np.where(in_soil, config.soil_shear_velocity_m_s, config.rock_shear_velocity_m_s).astype(float)
    return MaterialProfile(
        depth_m,
        density,
        shear_velocity,
        compute_shear_modulus(density, shear_velocity),
        "Soft soil layer over rock",
        config.soil_layer_thickness_m,
    )


def stress_grid_shear_modulus(material: MaterialProfile) -> np.ndarray:
    """Return harmonic-average shear modulus at stress-grid locations."""
    mu_left = material.shear_modulus_pa[:-1]
    mu_right = material.shear_modulus_pa[1:]
    return 2.0 * mu_left * mu_right / (mu_left + mu_right)


def material_summary(material: MaterialProfile) -> str:
    """Return a readable summary of the material profile."""
    return (
        f"{material.label}\n"
        f"Depth range:     {material.depth_m.min():.1f} to {material.depth_m.max():.1f} m\n"
        f"Vs range:        {material.shear_velocity_m_s.min():.1f} to {material.shear_velocity_m_s.max():.1f} m/s\n"
        f"Density range:   {material.density_kg_m3.min():.1f} to {material.density_kg_m3.max():.1f} kg/m^3\n"
        f"Impedance range: {material.impedance_kg_m2_s.min():.3e} to {material.impedance_kg_m2_s.max():.3e} kg/(m^2 s)\n"
    )
