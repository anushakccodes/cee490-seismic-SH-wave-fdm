"""Explicit 2D velocity-stress finite-difference solver for SH waves."""
from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
import numpy as np

from .materials_2d import (Material2D, Simulation2DConfig, stress_grid_mu_x, stress_grid_mu_z)
from .source import ricker_wavelet

@dataclass(frozen=True)
class Simulation2DResult:
    """Numerical result from a 2D SH velocity-stress run."""
    time_s: np.ndarray
    x_m: np.ndarray
    z_m: np.ndarray
    stored_time_s: np.ndarray
    velocity_snapshots_m_s: np.ndarray
    receiver_records_m_s: dict
    receiver_positions: dict
    source_time_function: np.ndarray
    source_x_m: float
    source_z_m: float
    material: Material2D
    config: Simulation2DConfig


def _build_damping_masks(config):
    """Gaussian absorbing masks for left, right, and bottom boundaries; top is free surface."""
    nz, nx = config.nz, config.nx
    thickness = config.absorbing_layer_thickness_m
    strength = config.absorbing_strength
    dx, dz = config.dx_m, config.dz_m
    Lx = config.model_width_m
    Lz = config.model_depth_m

    x_m = np.linspace(0.0, Lx, nx)
    z_m = np.linspace(0.0, Lz, nz)

    mx = np.ones(nx, dtype=float)
    mz = np.ones(nz, dtype=float)

    left_end = thickness
    in_left = x_m < left_end
    eta = (left_end - x_m[in_left]) / thickness
    mx[in_left] = np.exp(-strength * eta**2)

    right_start = Lx - thickness
    in_right = x_m > right_start
    eta = (x_m[in_right] - right_start) / thickness
    mx[in_right] = np.exp(-strength * eta**2)

    bot_start = Lz - thickness
    in_bot = z_m > bot_start
    eta = (z_m[in_bot] - bot_start) / thickness
    mz[in_bot] = np.exp(-strength * eta**2)

    vel_mask = mz[:, np.newaxis] * mx[np.newaxis, :]

    x_sx = 0.5 * (x_m[:-1] + x_m[1:])
    mxs = np.ones(nx - 1, dtype=float)
    in_left_s = x_sx < left_end
    eta = (left_end - x_sx[in_left_s]) / thickness
    mxs[in_left_s] = np.exp(-strength * eta**2)
    in_right_s = x_sx > right_start
    eta = (x_sx[in_right_s] - right_start) / thickness
    mxs[in_right_s] = np.exp(-strength * eta**2)
    stress_xy_mask = mz[:, np.newaxis] * mxs[np.newaxis, :]

    z_sz = 0.5 * (z_m[:-1] + z_m[1:])
    mzs = np.ones(nz - 1, dtype=float)
    in_bot_s = z_sz > bot_start
    eta = (z_sz[in_bot_s] - bot_start) / thickness
    mzs[in_bot_s] = np.exp(-strength * eta**2)
    stress_zy_mask = mzs[:, np.newaxis] * mx[np.newaxis, :]

    return vel_mask, stress_xy_mask, stress_zy_mask


def _nearest_index(arr, value):
    return int(np.argmin(np.abs(arr - value)))


def run_2d_solver(
    config,
    material,
    source_x_m=None,
    source_z_m=200.0,
    receivers=None,
    store_every=10,
    source_amplitude=1.0e-4,
    plane_wave_source=False,
    progress_every=None,
):
    """Staggered-grid leapfrog solver for the 2D SH velocity-stress system.

    Fields: vy (nz, nx), tau_xy (nz, nx-1), tau_zy (nz-1, nx).
    The z-staggered grid has no row above z=0, enforcing the stress-free surface.
    """

    if source_x_m is None:
        source_x_m = config.model_width_m / 2.0

    nz, nx = config.nz, config.nx
    dx, dz, dt = config.dx_m, config.dz_m, config.time_step_s
    x_m = material.x_m
    z_m = material.z_m
    density = material.density_kg_m3
    mu_x = stress_grid_mu_x(material)
    mu_z = stress_grid_mu_z(material)

    nt = config.number_of_time_steps
    time_s = np.linspace(0.0, config.total_time_s, nt)
    source_values = ricker_wavelet(time_s, config.dominant_frequency_hz, config.source_peak_time_s)

    src_zi = _nearest_index(z_m, source_z_m)
    src_xi = _nearest_index(x_m, source_x_m)

    if receivers is None:
        surface_xs = np.arange(100.0, config.model_width_m, 100.0)
        receivers = {f"surf_{int(xr):d}": (float(xr), 0.0) for xr in surface_xs}
        receivers["internal"] = (source_x_m + 100.0, source_z_m)

    rec_indices = {}
    for name, (rx, rz) in receivers.items():
        rec_indices[name] = (_nearest_index(z_m, rz), _nearest_index(x_m, rx))

    vel_mask, sxy_mask, szy_mask = _build_damping_masks(config)

    vy = np.zeros((nz, nx), dtype=float)
    sxy = np.zeros((nz, nx - 1), dtype=float)
    szy = np.zeros((nz - 1, nx), dtype=float)

    dt_over_rho = dt / density

    if progress_every is None:
        progress_every = max(1, nt // 10)

    snapshots = []
    snap_times = []
    rec_records = {name: [] for name in receivers}

    for step in range(nt):
        if step % progress_every == 0:
            pct = 100.0 * step / nt
            t_now = step * dt
            peak_vy = float(np.abs(vy).max())
            print(f'  step {step:5d}/{nt}  ({pct:5.1f}%)  t = {t_now:.3f} s  |vy|_max = {peak_vy:.3e} m/s')
        dsx_dx = (sxy[:, 1:] - sxy[:, :-1]) / dx

        # one-sided at free surface (tau_zy=0 above z=0) and at bottom
        dsz_dz_interior = (szy[1:, :] - szy[:-1, :]) / dz
        dsz_dz_top = szy[0, :] / dz
        dsz_dz_bot = -szy[-1, :] / dz

        vy[0, 1:-1]    += dt_over_rho[0, 1:-1]    * (dsx_dx[0, :]    + dsz_dz_top[1:-1])
        vy[1:-1, 1:-1] += dt_over_rho[1:-1, 1:-1] * (dsx_dx[1:-1, :] + dsz_dz_interior[:, 1:-1])
        vy[-1, 1:-1]   += dt_over_rho[-1, 1:-1]   * (dsx_dx[-1, :]   + dsz_dz_bot[1:-1])

        vy[0, 0]    += dt_over_rho[0, 0]    * (sxy[0, 0]    / dx + dsz_dz_top[0])
        vy[1:-1, 0] += dt_over_rho[1:-1, 0] * (sxy[1:-1, 0] / dx + dsz_dz_interior[:, 0])
        vy[-1, 0]   += dt_over_rho[-1, 0]   * (sxy[-1, 0]   / dx + dsz_dz_bot[0])

        vy[0, -1]    += dt_over_rho[0, -1]    * (-sxy[0, -1]    / dx + dsz_dz_top[-1])
        vy[1:-1, -1] += dt_over_rho[1:-1, -1] * (-sxy[1:-1, -1] / dx + dsz_dz_interior[:, -1])
        vy[-1, -1]   += dt_over_rho[-1, -1]   * (-sxy[-1, -1]   / dx + dsz_dz_bot[-1])

        if plane_wave_source:
            vy[src_zi, :] += source_amplitude * source_values[step]
        else:
            vy[src_zi, src_xi] += source_amplitude * source_values[step]

        sxy += mu_x * dt * np.diff(vy, axis=1) / dx
        szy += mu_z * dt * np.diff(vy, axis=0) / dz

        vy  *= vel_mask
        sxy *= sxy_mask
        szy *= szy_mask

        for name, (iz, ix) in rec_indices.items():
            rec_records[name].append(float(vy[iz, ix]))

        if step % store_every == 0:
            snapshots.append(vy.copy())
            snap_times.append(float(time_s[step]))

    print(f'  step {nt}/{nt}  (100.0%)  t = {time_s[-1]:.3f} s  |vy|_max = {float(np.abs(vy).max()):.3e} m/s  — done.')

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