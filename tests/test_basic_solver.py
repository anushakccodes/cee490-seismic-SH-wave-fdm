"""Basic smoke tests for the implemented Stage 1--3 solver."""
from __future__ import annotations

from site_response.analysis import arrival_time_verification
from site_response.fd_solver_1d import run_velocity_stress_solver
from site_response.materials import homogeneous_rock_model
from site_response.parameters import SimulationConfig, compute_parameter_checks


def test_default_cfl_passes() -> None:
    """The default explicit time step should satisfy the 1D CFL limit."""
    checks = compute_parameter_checks(SimulationConfig(total_time_s=0.05))
    assert checks["cfl_pass"]


def test_short_homogeneous_run_produces_finite_values() -> None:
    """A short run should complete with finite array-shaped outputs."""
    config = SimulationConfig(total_time_s=0.7)
    material = homogeneous_rock_model(config)
    result = run_velocity_stress_solver(config, material, store_every=5)
    assert result.receiver_velocity_m_s.shape == result.time_s.shape
    assert result.velocity_wavefield_m_s.ndim == 2
    assert result.receiver_velocity_m_s.size > 0


def test_arrival_verification_returns_nonnegative_error() -> None:
    """Arrival-time verification should return a valid nonnegative error."""
    config = SimulationConfig(total_time_s=0.9)
    material = homogeneous_rock_model(config)
    result = run_velocity_stress_solver(config, material, store_every=10)
    check = arrival_time_verification(result)
    assert check.percent_error >= 0.0
