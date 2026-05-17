"""Verification and numerical-analysis utilities."""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
import pandas as pd

from .fd_solver_1d import SimulationResult
from .materials import homogeneous_rock_model, layered_soil_over_rock_model
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
        return pd.DataFrame([self.__dict__])


@dataclass(frozen=True)
class ConvergenceRow:
    grid_spacing_m: float
    courant_number: float
    theoretical_peak_arrival_s: float
    numerical_peak_arrival_s: float
    percent_error: float


@dataclass(frozen=True)
class AmplificationResult:
    """Spectral amplification result for a layered site."""
    frequencies_hz: np.ndarray
    ref_spectrum: np.ndarray
    layered_spectrum: np.ndarray
    amplification: np.ndarray
    resonance_frequencies_hz: np.ndarray

    def resonance_table(self) -> pd.DataFrame:
        rows = [{"mode_n": n + 1, "resonance_frequency_hz": float(f)} for n, f in enumerate(self.resonance_frequencies_hz)]
        return pd.DataFrame(rows)


def theoretical_travel_time(distance_m: float, shear_velocity_m_s: float) -> float:
    return distance_m / shear_velocity_m_s


def pick_peak_arrival_time(
    time_s: np.ndarray,
    signal: np.ndarray,
    expected_peak_time_s: float,
    search_half_width_s: float = 0.08,
) -> float:
    window = (time_s >= expected_peak_time_s - search_half_width_s) & (time_s <= expected_peak_time_s + search_half_width_s)
    if not np.any(window):
        raise ValueError("Arrival search window does not overlap the time axis.")
    local_time = time_s[window]
    local_signal = np.abs(signal[window])
    return float(local_time[int(np.argmax(local_signal))])


def arrival_time_verification(result: SimulationResult) -> ArrivalTimeCheck:
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
    from .fd_solver_1d import run_velocity_stress_solver
    material = homogeneous_rock_model(config)
    result = run_velocity_stress_solver(config=config, material=material)
    return result, arrival_time_verification(result)


# ---------------------------------------------------------------------------
# Stage 4: Grid-refinement convergence study
# ---------------------------------------------------------------------------

def run_convergence_study(
    config: SimulationConfig,
    dz_values_m: list[float] | None = None,
) -> tuple[list[ConvergenceRow], pd.DataFrame]:
    """Run the homogeneous model at several grid spacings and collect arrival-time errors.

    For each dz, dt is scaled to maintain a constant Courant number of 0.5.
    """
    from .fd_solver_1d import run_velocity_stress_solver

    if dz_values_m is None:
        dz_values_m = [8.0, 4.0, 2.0, 1.0]

    base_courant = 0.5
    rows: list[ConvergenceRow] = []

    for dz in dz_values_m:
        dt = base_courant * dz / config.maximum_shear_velocity_m_s
        refined_config = replace(config, grid_spacing_m=dz, time_step_s=dt)
        material = homogeneous_rock_model(refined_config)
        result = run_velocity_stress_solver(config=refined_config, material=material)
        check = arrival_time_verification(result)
        courant = config.maximum_shear_velocity_m_s * dt / dz
        rows.append(ConvergenceRow(
            grid_spacing_m=dz,
            courant_number=courant,
            theoretical_peak_arrival_s=check.theoretical_peak_arrival_s,
            numerical_peak_arrival_s=check.numerical_peak_arrival_s,
            percent_error=check.percent_error,
        ))

    df = pd.DataFrame([r.__dict__ for r in rows])
    return rows, df


# ---------------------------------------------------------------------------
# Stage 5: Layered site-response model and amplification
# ---------------------------------------------------------------------------

def surface_fourier_spectrum(
    time_s: np.ndarray,
    surface_velocity: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return one-sided Fourier amplitude spectrum of the surface velocity."""
    dt = float(time_s[1] - time_s[0])
    frequencies_hz = np.fft.rfftfreq(len(surface_velocity), d=dt)
    amplitude = np.abs(np.fft.rfft(surface_velocity))
    return frequencies_hz, amplitude


def compute_amplification(
    ref_result: SimulationResult,
    layered_result: SimulationResult,
    epsilon: float = 1e-20,
) -> AmplificationResult:
    """Compute the spectral amplification A(f) = |V_layered| / (|V_ref| + eps)."""
    freq_ref, spec_ref = surface_fourier_spectrum(ref_result.time_s, ref_result.surface_velocity_m_s)
    freq_lay, spec_lay = surface_fourier_spectrum(layered_result.time_s, layered_result.surface_velocity_m_s)

    # Interpolate onto the reference frequency axis if lengths differ.
    if len(freq_ref) != len(freq_lay):
        spec_lay = np.interp(freq_ref, freq_lay, spec_lay)
        freq_lay = freq_ref

    amplification = spec_lay / (spec_ref + epsilon)
    resonance_freqs = _quarter_wavelength_frequencies(layered_result.config, n_modes=5)

    return AmplificationResult(
        frequencies_hz=freq_ref,
        ref_spectrum=spec_ref,
        layered_spectrum=spec_lay,
        amplification=amplification,
        resonance_frequencies_hz=resonance_freqs,
    )


def _quarter_wavelength_frequencies(config: SimulationConfig, n_modes: int = 5) -> np.ndarray:
    """Compute quarter-wavelength resonance frequencies for the soft layer."""
    H = config.soil_layer_thickness_m
    Vs_soil = config.soil_shear_velocity_m_s
    return np.array([(2 * n - 1) * Vs_soil / (4.0 * H) for n in range(1, n_modes + 1)])


def run_layered_site_case(config: SimulationConfig) -> tuple[SimulationResult, SimulationResult]:
    """Run both the rock reference and layered-site simulations.

    The source is placed below the soft layer (in the rock halfspace) to drive
    the column from below, which is the standard engineering arrangement.
    """
    from .fd_solver_1d import run_velocity_stress_solver

    source_depth_m = config.soil_layer_thickness_m + 20.0
    receiver_depth_m = 0.0

    ref_material = homogeneous_rock_model(config)
    ref_result = run_velocity_stress_solver(
        config=config,
        material=ref_material,
        source_depth_m=source_depth_m,
        receiver_depth_m=receiver_depth_m,
    )

    layered_material = layered_soil_over_rock_model(config)
    layered_result = run_velocity_stress_solver(
        config=config,
        material=layered_material,
        source_depth_m=source_depth_m,
        receiver_depth_m=receiver_depth_m,
    )

    return ref_result, layered_result
