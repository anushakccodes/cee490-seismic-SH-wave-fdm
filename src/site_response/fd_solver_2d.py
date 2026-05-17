"""Explicit 2D velocity-stress finite-difference solver for SH waves."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .materials_2d import (
    Material2D,
    Simulation2DConfig,
    stress_grid_mu_x,
    stress_grid_mu_z,
)
from .source import ricker_wavelet


@dataclass(frozen=True)
class Simulation2DResult:
    """Numerical result from a 2D SH velocity-stress run."""
    time_s: np.ndarray          # (nt,)
    x_m: np.ndarray             # (nx,)
    z_m: np.ndarray             # (nz,)
    stored_time_s: np.ndarray   # (n_snap,)
    velocity_snapshots_m_s: np.ndarray  # (n_snap, nz, nx)
    receiver_records_m_s: dict[str, np.ndarray]   # name -> (nt,)
    receiver_positions: dict[str, tuple[float, float]]  # name -> (x_m, z_m)
    source_time_function: np.ndarray  # (nt,)
    source_x_m: float
    source_z_m: float
    material: Material2D
    config: Simulation2DConfig


def _build_damping_masks(config: Simulation2DConfig) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Build 2D Gaussian damping masks for velocity and both stress grids.

    Absorbs on left, right, and bottom edges. Top is free surface — no damping.
    Returns: vel_mask (nz, nx), stress_xy_mask (nz, nx-1), stress_zy_mask (nz-1, nx)
    """
    nz, nx = config.nz, config.nx
    thickness = config.absorbing_layer_thickness_m
    strength = config.absorbing_strength
    dx, dz = config.dx_m, config.dz_m
    Lx = config.model_width_m
    Lz = config.model_depth_m

    x_m = np.linspace(0.0, Lx, nx)
    z_m = np.linspace(0.0, Lz, nz)

    # Velocity grid mask (nz, nx)
    mx = np.ones(nx, dtype=float)
    mz = np.ones(nz, dtype=float)

    # left absorbing
    left_end = thickness
    in_left = x_m < left_end
    eta = (left_end - x_m[in_left]) / thickness
    mx[in_left] = np.exp(-strength * eta**2)

    # right absorbing
    right_start = Lx - thickness
    in_right = x_m > right_start
    eta = (x_m[in_right] - right_start) / thickness
    mx[in_right] = np.exp(-strength * eta**2)

    # bottom absorbing (no damping at top = free surface)
    bot_start = Lz - thickness
    in_bot = z_m > bot_start
    eta = (z_m[in_bot] - bot_start) / thickness
    mz[in_bot] = np.exp(-strength * eta**2)

    vel_mask = mz[:, np.newaxis] * mx[np.newaxis, :]  # (nz, nx)

    # stress_xy at (nz, nx-1): average x-mask between adjacent x-nodes
    x_sx = 0.5 * (x_m[:-1] + x_m[1:])
    mxs = np.ones(nx - 1, dtype=float)
    in_left_s = x_sx < left_end
    eta = (left_end - x_sx[in_left_s]) / thickness
    mxs[in_left_s] = np.exp(-strength * eta**2)
    in_right_s = x_sx > right_start
    eta = (x_sx[in_right_s] - right_start) / thickness
    mxs[in_right_s] = np.exp(-strength * eta**2)
    stress_xy_mask = mz[:, np.newaxis] * mxs[np.newaxis, :]  # (nz, nx-1)

    # stress_zy at (nz-1, nx): average z-mask between adjacent z-nodes
    z_sz = 0.5 * (z_m[:-1] + z_m[1:])
    mzs = np.ones(nz - 1, dtype=float)
    in_bot_s = z_sz > bot_start
    eta = (z_sz[in_bot_s] - bot_start) / thickness
    mzs[in_bot_s] = np.exp(-strength * eta**2)
    stress_zy_mask = mzs[:, np.newaxis] * mx[np.newaxis, :]  # (nz-1, nx)

    return vel_mask, stress_xy_mask, stress_zy_mask


def _nearest_index(arr: np.ndarray, value: float) -> int:
    return int(np.argmin(np.abs(arr - value)))


def run_2d_solver(
    config: Simulation2DConfig,
    material: Material2D,
    source_x_m: float | None = None,
    source_z_m: float = 200.0,
    receivers: dict[str, tuple[float, float]] | None = None,
    store_every: int = 10,
    source_scale: float = 1.0e-4,
    plane_wave_source: bool = False,
) -> Simulation2DResult:
    """Run the 2D staggered-grid SH velocity-stress solver.

    Fields:
      velocity_y  (nz, nx)    — out-of-plane particle velocity
      stress_xy   (nz, nx-1)  — shear stress, x-staggered
      stress_zy   (nz-1, nx)  — shear stress, z-staggered

    Free surface: tau_zy at z=0 is excluded (z-staggered grid starts between
    row 0 and row 1), so the top boundary is naturally stress-free.

    Parameters
    ----------
    plane_wave_source : bool
        If True, apply source over all x-nodes at source_z_m (1D-like input).
        If False (default), apply source at single (source_x_m, source_z_m).
    """
    if source_x_m is None:
        source_x_m = config.model_width_m / 2.0

    nz, nx = config.nz, config.nx
    dx, dz, dt = config.dx_m, config.dz_m, config.time_step_s
    x_m = material.x_m
    z_m = material.z_m
    density = material.density_kg_m3   # (nz, nx)
    mu_x = stress_grid_mu_x(material)  # (nz, nx-1)
    mu_z = stress_grid_mu_z(material)  # (nz-1, nx)

    # Time axis and source
    nt = config.number_of_time_steps
    time_s = np.linspace(0.0, config.total_time_s, nt)
    source_values = ricker_wavelet(time_s, config.dominant_frequency_hz, config.source_peak_time_s)

    # Source/receiver indices
    src_zi = _nearest_index(z_m, source_z_m)
    src_xi = _nearest_index(x_m, source_x_m)

    # Build receiver index map
    if receivers is None:
        # Default: surface array + one internal verification receiver
        surface_xs = np.arange(100.0, config.model_width_m, 100.0)
        receivers = {f"surf_{int(xr):d}": (float(xr), 0.0) for xr in surface_xs}
        receivers["internal"] = (source_x_m + 100.0, source_z_m)

    rec_indices: dict[str, tuple[int, int]] = {}
    for name, (rx, rz) in receivers.items():
        rec_indices[name] = (_nearest_index(z_m, rz), _nearest_index(x_m, rx))

    # Damping masks
    vel_mask, sxy_mask, szy_mask = _build_damping_masks(config)

    # Field arrays
    vy = np.zeros((nz, nx), dtype=float)
    sxy = np.zeros((nz, nx - 1), dtype=float)
    szy = np.zeros((nz - 1, nx), dtype=float)

    # Pre-compute dt/rho for velocity update
    dt_over_rho = dt / density  # (nz, nx)

    # Storage
    n_snap = nt // store_every + 1
    snapshots: list[np.ndarray] = []
    snap_times: list[float] = []
    rec_records: dict[str, list[float]] = {name: [] for name in receivers}

    for step in range(nt):
        # ------------------------------------------------------------------ #
        # 1. Velocity update
        #    rho dv/dt = d tau_xy/dx + d tau_zy/dz
        #
        # Interior x (columns 1..nx-2): both x-stress gradient terms available
        # Left col (x=0) and right col (x=nx-1): handled by boundary condition
        # Top row (z=0): free surface — tau_zy(z=0) = 0
        # ------------------------------------------------------------------ #

        # x-stress gradient: (sxy[:, 1:] - sxy[:, :-1]) / dx  -> (nz, nx-2)
        dsx_dx = (sxy[:, 1:] - sxy[:, :-1]) / dx  # (nz, nx-2) for cols 1..nx-2

        # z-stress gradient: (szy[1:,:] - szy[:-1,:]) / dz  -> (nz-2, nx) for rows 1..nz-2
        # For the free surface row (iz=0), tau_zy just above is zero:
        #   d tau_zy/dz at row 0 = (szy[0,:] - 0) / dz  (one-sided)
        # For the bottom row (iz=nz-1), tau_zy below is absorbed:
        #   d tau_zy/dz at row nz-1 = (0 - szy[-1,:]) / dz

        dsz_dz_interior = (szy[1:, :] - szy[:-1, :]) / dz  # (nz-2, nx) rows 1..nz-2
        dsz_dz_top = szy[0, :] / dz                          # (nx,)  one-sided free surface
        dsz_dz_bot = -szy[-1, :] / dz                        # (nx,)

        # Update interior rows and columns
        vy[0, 1:-1] += dt_over_rho[0, 1:-1] * (dsx_dx[0, :] + dsz_dz_top[1:-1])
        vy[1:-1, 1:-1] += dt_over_rho[1:-1, 1:-1] * (dsx_dx[1:-1, :] + dsz_dz_interior[:, 1:-1])
        vy[-1, 1:-1] += dt_over_rho[-1, 1:-1] * (dsx_dx[-1, :] + dsz_dz_bot[1:-1])

        # Left boundary (ix=0): no x-stress gradient term to the left
        vy[0, 0] += dt_over_rho[0, 0] * (sxy[0, 0] / dx + dsz_dz_top[0])
        vy[1:-1, 0] += dt_over_rho[1:-1, 0] * (sxy[1:-1, 0] / dx + dsz_dz_interior[:, 0])
        vy[-1, 0] += dt_over_rho[-1, 0] * (sxy[-1, 0] / dx + dsz_dz_bot[0])

        # Right boundary (ix=nx-1): no x-stress gradient term to the right
        vy[0, -1] += dt_over_rho[0, -1] * (-sxy[0, -1] / dx + dsz_dz_top[-1])
        vy[1:-1, -1] += dt_over_rho[1:-1, -1] * (-sxy[1:-1, -1] / dx + dsz_dz_interior[:, -1])
        vy[-1, -1] += dt_over_rho[-1, -1] * (-sxy[-1, -1] / dx + dsz_dz_bot[-1])

        # ------------------------------------------------------------------ #
        # 2. Source injection
        # ------------------------------------------------------------------ #
        if plane_wave_source:
            vy[src_zi, :] += source_scale * source_values[step]
        else:
            vy[src_zi, src_xi] += source_scale * source_values[step]

        # ------------------------------------------------------------------ #
        # 3. Stress updates
        #    d tau_xy/dt = mu * dv/dx   -> (nz, nx-1)
        #    d tau_zy/dt = mu * dv/dz   -> (nz-1, nx)
        # ------------------------------------------------------------------ #
        sxy += mu_x * dt * np.diff(vy, axis=1) / dx
        szy += mu_z * dt * np.diff(vy, axis=0) / dz

        # ------------------------------------------------------------------ #
        # 4. Absorbing boundary damping
        # ------------------------------------------------------------------ #
        vy *= vel_mask
        sxy *= sxy_mask
        szy *= szy_mask

        # ------------------------------------------------------------------ #
        # 5. Record receivers and snapshots
        # ------------------------------------------------------------------ #
        for name, (iz, ix) in rec_indices.items():
            rec_records[name].append(float(vy[iz, ix]))

        if step % store_every == 0:
            snapshots.append(vy.copy())
            snap_times.append(float(time_s[step]))

    receiver_arrays = {name: np.asarray(vals) for name, vals in rec_records.items()}

    return Simulation2DResult(
        time_s=time_s,
        x_m=x_m,
        z_m=z_m,
        stored_time_s=np.asarray(snap_times),
        velocity_snapshots_m_s=np.asarray(snapshots),
        receiver_records_m_s=receiver_arrays,
        receiver_positions=receivers,
        source_time_function=source_values,
        source_x_m=float(x_m[src_xi]),
        source_z_m=float(z_m[src_zi]),
        material=material,
        config=config,
    )
