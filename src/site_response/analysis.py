"""Verification and numerical-analysis utilities."""
from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .fd_solver_1d import SimulationResult
from .materials import homogeneous_rock_model
from .parameters import SimulationConfig


@dataclass(frozen=True)
class ArrivalTimeCheck:
    """Arrival-time comparison for a homogeneous verification run."""
    source_depth_m: float
    receiver_depth_m: float
    travel_distance_m: float
    shear_velocity_m_s: float
    source_peak_time_s: float
    theoretical_peak_arrival_s: float
    numerical_peak_arrival_s: float
    absolute_error_s: float
    percent_error: float

    def as_dataframe(self) -> pd.DataFrame:
        """Return the check as a one-row pandas table."""
        return pd.DataFrame([self.__dict__])


def theoretical_travel_time(distance_m: float, shear_velocity_m_s: float) -> float:
    """Compute t = distance / Vs for a homogeneous medium."""
    return distance_m / shear_velocity_m_s


def pick_peak_arrival_time(
    time_s: np.ndarray,
    signal: np.ndarray,
    expected_peak_time_s: float,
    search_half_width_s: float = 0.08,
) -> float:
    """Pick the largest absolute signal peak near an expected arrival time."""
    window = (time_s >= expected_peak_time_s - search_half_width_s) & (time_s <= expected_peak_time_s + search_half_width_s)
    if not np.any(window):
        raise ValueError("Arrival search window does not overlap the time axis.")
    local_time = time_s[window]
    local_signal = np.abs(signal[window])
    return float(local_time[int(np.argmax(local_signal))])


def arrival_time_verification(result: SimulationResult) -> ArrivalTimeCheck:
    """Compare numerical and theoretical arrival time in a homogeneous model."""
    unique_vs = np.unique(np.round(result.material.shear_velocity_m_s, decimals=10))
    if len(unique_vs) != 1:
        raise ValueError("Arrival-time verification requires a homogeneous material profile.")
    shear_velocity_m_s = float(unique_vs[0])
    travel_distance_m = abs(result.receiver_depth_m - result.source_depth_m)
    theoretical_peak_arrival_s = result.config.source_peak_time_s + theoretical_travel_time(travel_distance_m, shear_velocity_m_s)
    numerical_peak_arrival_s = pick_peak_arrival_time(result.time_s, result.receiver_velocity_m_s, theoretical_peak_arrival_s)
    absolute_error_s = abs(numerical_peak_arrival_s - theoretical_peak_arrival_s)
    percent_error = absolute_error_s / theoretical_peak_arrival_s * 100.0
    return ArrivalTimeCheck(
        source_depth_m=result.source_depth_m,
        receiver_depth_m=result.receiver_depth_m,
        travel_distance_m=travel_distance_m,
        shear_velocity_m_s=shear_velocity_m_s,
        source_peak_time_s=result.config.source_peak_time_s,
        theoretical_peak_arrival_s=theoretical_peak_arrival_s,
        numerical_peak_arrival_s=numerical_peak_arrival_s,
        absolute_error_s=absolute_error_s,
        percent_error=percent_error,
    )


def run_homogeneous_verification_case(config: SimulationConfig) -> tuple[SimulationResult, ArrivalTimeCheck]:
    """Run the default homogeneous model and return result plus verification metrics."""
    from .fd_solver_1d import run_velocity_stress_solver

    material = homogeneous_rock_model(config)
    result = run_velocity_stress_solver(config=config, material=material)
    return result, arrival_time_verification(result)
