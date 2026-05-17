"""Explicit 1D velocity-stress finite-difference solver for SH waves."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np

from .materials import MaterialProfile, stress_grid_shear_modulus
from .parameters import SimulationConfig, make_time_axis
from .source import build_source_time_function


@dataclass(frozen=True)
class SimulationResult:
    """Numerical result from a 1D velocity-stress run."""
    time_s: np.ndarray
    depth_m: np.ndarray
    stored_time_s: np.ndarray
    velocity_wavefield_m_s: np.ndarray
    surface_velocity_m_s: np.ndarray
    receiver_velocity_m_s: np.ndarray
    source_time_function: np.ndarray
    source_depth_m: float
    receiver_depth_m: float
    source_index: int
    receiver_index: int
    material: MaterialProfile
    config: SimulationConfig


def nearest_grid_index(depth_m: np.ndarray, requested_depth_m: float) -> int:
    """Return the grid index closest to a requested depth."""
    return int(np.argmin(np.abs(depth_m - requested_depth_m)))


def build_bottom_damping_mask(depth_m: np.ndarray, damping_thickness_m: float, damping_strength: float) -> np.ndarray:
    """Create a smooth damping mask near the finite bottom boundary."""
    mask = np.ones_like(depth_m, dtype=float)
    damping_start_m = depth_m.max() - damping_thickness_m
    damping_nodes = depth_m > damping_start_m
    if np.any(damping_nodes):
        eta = (depth_m[damping_nodes] - damping_start_m) / damping_thickness_m
        mask[damping_nodes] = np.exp(-damping_strength * eta**2)
    return mask


def run_velocity_stress_solver(
    config: SimulationConfig,
    material: MaterialProfile,
    source_depth_m: float = 40.0,
    receiver_depth_m: float = 260.0,
    store_every: int = 2,
    source_scale: float = 2.0e-4,
) -> SimulationResult:
    """Run the staggered 1D SH-wave solver.

    The solver advances rho * dv/dt = d tau/dz and d tau/dt = mu * dv/dz.
    Velocity is stored at depth nodes; shear stress is stored between nodes.
    """
    if store_every < 1:
        raise ValueError("store_every must be at least 1.")

    time_s, source_values = build_source_time_function(config)
    depth_m = material.depth_m
    density = material.density_kg_m3
    mu_stress = stress_grid_shear_modulus(material)

    velocity = np.zeros_like(depth_m, dtype=float)
    shear_stress = np.zeros(len(depth_m) - 1, dtype=float)

    source_index = nearest_grid_index(depth_m, source_depth_m)
    receiver_index = nearest_grid_index(depth_m, receiver_depth_m)
    actual_source_depth_m = float(depth_m[source_index])
    actual_receiver_depth_m = float(depth_m[receiver_index])

    velocity_damping = build_bottom_damping_mask(depth_m, config.absorbing_layer_thickness_m, config.absorbing_strength)
    stress_depth_m = 0.5 * (depth_m[:-1] + depth_m[1:])
    stress_damping = build_bottom_damping_mask(stress_depth_m, config.absorbing_layer_thickness_m, config.absorbing_strength)

    stored_velocity_profiles: list[np.ndarray] = []
    stored_time_values: list[float] = []
    surface_velocity = np.zeros_like(time_s)
    receiver_velocity = np.zeros_like(time_s)

    dz = config.grid_spacing_m
    dt = config.time_step_s

    for step_index, current_time_s in enumerate(time_s):
        # Velocity update with stress-free top boundary tau(z=0,t)=0.
        velocity[0] += (dt / density[0]) * (shear_stress[0] - 0.0) / dz
        velocity[1:-1] += (dt / density[1:-1]) * (shear_stress[1:] - shear_stress[:-1]) / dz
        velocity[-1] += (dt / density[-1]) * (0.0 - shear_stress[-1]) / dz

        # Compact horizontal body-force source at one velocity node.
        velocity[source_index] += source_scale * source_values[step_index]

        # Stress update between velocity nodes.
        shear_stress += mu_stress * dt * np.diff(velocity) / dz

        # Approximate absorbing treatment near the finite bottom boundary.
        velocity *= velocity_damping
        shear_stress *= stress_damping

        surface_velocity[step_index] = velocity[0]
        receiver_velocity[step_index] = velocity[receiver_index]

        if step_index % store_every == 0:
            stored_velocity_profiles.append(velocity.copy())
            stored_time_values.append(float(current_time_s))

    return SimulationResult(
        time_s=time_s,
        depth_m=depth_m,
        stored_time_s=np.asarray(stored_time_values),
        velocity_wavefield_m_s=np.asarray(stored_velocity_profiles),
        surface_velocity_m_s=surface_velocity,
        receiver_velocity_m_s=receiver_velocity,
        source_time_function=source_values,
        source_depth_m=actual_source_depth_m,
        receiver_depth_m=actual_receiver_depth_m,
        source_index=source_index,
        receiver_index=receiver_index,
        material=material,
        config=config,
    )
