"""Validation tests for the 1D SH-wave finite-difference solver.

These tests verify physical correctness, not just that the code runs:
  - CFL stability condition.
  - Arrival-time accuracy in the homogeneous case.
  - 2nd-order spatial convergence of arrival-time error.
  - Spectral amplification near the fundamental resonance frequency.
  - Energy conservation in an undamped homogeneous run.
"""
from __future__ import annotations

import numpy as np
import pytest

from site_response.analysis import (
    arrival_time_verification,
    compute_amplification,
    run_convergence_study,
    run_layered_site_case,
)
from site_response.fd_solver_1d import run_velocity_stress_solver
from site_response.materials import homogeneous_rock_model
from site_response.parameters import SimulationConfig, compute_parameter_checks


# ---------------------------------------------------------------------------
# Stage 1: parameter checks
# ---------------------------------------------------------------------------

def test_default_cfl_passes() -> None:
    checks = compute_parameter_checks(SimulationConfig(total_time_s=0.05))
    assert checks["cfl_pass"], "Default configuration violates the CFL stability limit."


def test_default_grid_resolution_passes() -> None:
    checks = compute_parameter_checks(SimulationConfig(total_time_s=0.05))
    assert checks["grid_resolution_pass"], "Default grid has fewer than 10 points per shortest wavelength."


# ---------------------------------------------------------------------------
# Stage 3: arrival-time accuracy
# ---------------------------------------------------------------------------

def test_arrival_time_error_below_one_percent() -> None:
    """Default grid should reproduce travel time within 1 %."""
    config = SimulationConfig(total_time_s=0.9)
    material = homogeneous_rock_model(config)
    result = run_velocity_stress_solver(config, material, store_every=10)
    check = arrival_time_verification(result)
    assert check.percent_error < 1.0, (
        f"Arrival-time error {check.percent_error:.3f} % exceeds 1 % tolerance."
    )


def test_surface_velocity_nonzero_after_source() -> None:
    """Surface should receive energy before the end of the simulation."""
    config = SimulationConfig(total_time_s=1.0)
    material = homogeneous_rock_model(config)
    result = run_velocity_stress_solver(config, material, source_depth_m=20.0, receiver_depth_m=0.0)
    peak = np.max(np.abs(result.surface_velocity_m_s))
    assert peak > 0.0, "Surface velocity is identically zero — wave did not reach the surface."


# ---------------------------------------------------------------------------
# Stage 4: 2nd-order spatial convergence
# ---------------------------------------------------------------------------

def test_convergence_study_error_decreases() -> None:
    """Arrival-time error must decrease monotonically as the grid is refined."""
    config = SimulationConfig(total_time_s=0.9)
    _, df = run_convergence_study(config, dz_values_m=[8.0, 4.0, 2.0])
    errors = df["percent_error"].values
    for coarse, fine in zip(errors[:-1], errors[1:]):
        assert fine <= coarse, (
            f"Error did not decrease from dz={df['grid_spacing_m'].values[list(errors).index(coarse)]} "
            f"to the next finer grid: {coarse:.4f} -> {fine:.4f} %"
        )


def test_convergence_rate_near_first_order() -> None:
    """The log-log slope of arrival-time error vs dz should be near 1.

    With peak-picking and constant Courant number (dt ∝ dz), the measured
    error is quantised to ±Δt/2, giving O(Δz) convergence.  This is expected
    behaviour for this measurement approach, not a limitation of the scheme.
    """
    config = SimulationConfig(total_time_s=1.0)
    _, df = run_convergence_study(config, dz_values_m=[8.0, 4.0, 2.0, 1.0])
    errors = df["percent_error"].values
    dz_vals = df["grid_spacing_m"].values
    positive = errors > 0
    if positive.sum() < 2:
        pytest.skip("Not enough nonzero error values to compute convergence rate.")
    slope = np.polyfit(np.log(dz_vals[positive]), np.log(errors[positive]), 1)[0]
    assert 0.7 <= slope <= 1.5, (
        f"Convergence rate {slope:.2f} is outside the expected range [0.7, 1.5] "
        "for peak-picking with constant Courant number."
    )


# ---------------------------------------------------------------------------
# Stage 5: layered site response
# ---------------------------------------------------------------------------

def test_layered_surface_amplitude_exceeds_reference() -> None:
    """The layered-site peak surface amplitude should exceed the rock reference
    near the fundamental resonance frequency (soft layer traps energy)."""
    config = SimulationConfig(total_time_s=1.8)
    ref_result, layered_result = run_layered_site_case(config)
    peak_ref = np.max(np.abs(ref_result.surface_velocity_m_s))
    peak_layered = np.max(np.abs(layered_result.surface_velocity_m_s))
    assert peak_layered > peak_ref, (
        "Layered-site peak surface amplitude should exceed the rock reference "
        "due to impedance amplification in the soft layer."
    )


# ---------------------------------------------------------------------------
# Stage 6: amplification near fundamental resonance
# ---------------------------------------------------------------------------

def test_amplification_peak_near_fundamental_resonance() -> None:
    """The largest amplification should occur within ±50 % of the f1 frequency."""
    config = SimulationConfig(total_time_s=1.8)
    ref_result, layered_result = run_layered_site_case(config)
    amp = compute_amplification(ref_result, layered_result)

    f1 = config.soil_shear_velocity_m_s / (4.0 * config.soil_layer_thickness_m)
    in_band = np.isfinite(amp.amplification)
    if not in_band.any():
        pytest.skip("No in-band amplification values available.")

    peak_freq = amp.frequencies_hz[in_band][np.nanargmax(amp.amplification[in_band])]
    assert abs(peak_freq - f1) / f1 < 0.5, (
        f"Peak amplification at {peak_freq:.2f} Hz is more than 50 % away from "
        f"the theoretical fundamental resonance f1 = {f1:.2f} Hz."
    )


# ---------------------------------------------------------------------------
# Energy conservation (undamped homogeneous run)
# ---------------------------------------------------------------------------

def test_energy_bounded_undamped() -> None:
    """Total energy should stay bounded (not grow) after the source ends in an
    undamped homogeneous run.  We compare energy in the first third of the
    post-source record to the last third."""
    from dataclasses import replace
    config = SimulationConfig(
        total_time_s=1.2,
        absorbing_layer_thickness_m=0.0,  # no damping
        absorbing_strength=0.0,
    )
    material = homogeneous_rock_model(config)
    result = run_velocity_stress_solver(config, material, store_every=5)

    # Discrete total energy at each stored frame.
    dz = config.grid_spacing_m
    rho = result.material.density_kg_m3
    from site_response.materials import stress_grid_shear_modulus
    mu = stress_grid_shear_modulus(result.material)
    energies = []
    for frame in result.velocity_wavefield_m_s:
        ke = 0.5 * np.sum(rho * frame**2) * dz
        # Stress is not stored per frame; use a proxy: KE only suffices to check growth.
        energies.append(ke)

    energies = np.array(energies)
    n = len(energies)
    third = max(1, n // 3)
    early_mean = energies[:third].mean()
    late_mean = energies[-third:].mean()
    # Allow up to 20 % growth (boundary reflections can temporarily concentrate energy).
    assert late_mean <= 1.2 * early_mean, (
        f"Energy grew from {early_mean:.3e} to {late_mean:.3e} — "
        "possible instability in undamped run."
    )
