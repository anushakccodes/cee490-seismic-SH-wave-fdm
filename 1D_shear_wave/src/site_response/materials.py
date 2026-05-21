"""Material models for homogeneous and layered 1D SH-wave simulations."""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np

from .parameters import SimulationConfig, make_depth_axis


@dataclass(frozen=True)
class MaterialProfile:
    depth_m: np.ndarray
    density_kg_m3: np.ndarray
    shear_velocity_m_s: np.ndarray
    shear_modulus_pa: np.ndarray
    label: str
    interface_depth_m: float | None = None

    @property
    def impedance_kg_m2_s(self):
        """Shear impedance Z = rho * Vs."""
        return self.density_kg_m3 * self.shear_velocity_m_s


def compute_shear_modulus(density_kg_m3, shear_velocity_m_s):
    """Shear modulus mu = rho * Vs^2."""
    return density_kg_m3 * shear_velocity_m_s**2


def homogeneous_rock_model(config):
    depth_m = make_depth_axis(config)
    density = np.full_like(depth_m, config.rock_density_kg_m3, dtype=float)
    shear_velocity = np.full_like(depth_m, config.rock_shear_velocity_m_s, dtype=float)
    return MaterialProfile(depth_m, density, shear_velocity, compute_shear_modulus(density, shear_velocity), "Homogeneous rock reference")


def layered_soil_over_rock_model(config):
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


def stress_grid_shear_modulus(material):
    """Harmonic average of mu across each stress-grid interface."""
    mu_left = material.shear_modulus_pa[:-1]
    mu_right = material.shear_modulus_pa[1:]
    return 2.0 * mu_left * mu_right / (mu_left + mu_right)


def material_summary(material):
    lines = [
        material.label,
        f"Total depth range:  {material.depth_m.min():.1f} to {material.depth_m.max():.1f} m",
    ]
    if material.interface_depth_m is not None:
        h = material.interface_depth_m
        rock_bottom = material.depth_m.max()
        lines.append(f"  Soil layer:       0.0 to {h:.1f} m  (thickness = {h:.1f} m)")
        lines.append(f"  Rock halfspace:   {h:.1f} to {rock_bottom:.1f} m  (thickness = {rock_bottom - h:.1f} m)")
    lines += [
        f"Vs range:           {material.shear_velocity_m_s.min():.1f} to {material.shear_velocity_m_s.max():.1f} m/s",
        f"Density range:      {material.density_kg_m3.min():.1f} to {material.density_kg_m3.max():.1f} kg/m^3",
        f"Impedance range:    {material.impedance_kg_m2_s.min():.3e} to {material.impedance_kg_m2_s.max():.3e} kg/(m^2 s)",
    ]
    return "\n".join(lines) + "\n"