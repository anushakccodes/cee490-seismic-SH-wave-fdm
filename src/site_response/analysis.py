"""Verification and numerical-analysis utilities."""
from __future__ import annotations

from dataclasses import dataclass, replace
import numpy as np
import pandas as pd

from .fd_solver_1d import SimulationResult
from .materials import homogeneous_rock_model, layered_soil_over_rock_model
from .parameters import SimulationConfig


# ---------------------------------------------------------------------------
# Stage 3: Homogeneous arrival-time verification
# ---------------------------------------------------------------------------

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


def theoretical_travel_time(distance_m: float, shear_velocity_m_s: float) -> float:
    return distance_m / shear_velocity_m_s


def pick_peak_arrival_time(
    time_s: np.ndarray,
    signal: np.ndarray,
    expected_peak_time_s: float,
    search_half_width_s: float = 0.15,
) -> float:
    """Pick the largest absolute signal peak near an expected arrival time.

    search_half_width_s is intentionally wider than one grid cell so that
    the picker works correctly even for coarse grids in the convergence study.
    """
    window = (time_s >= expected_peak_time_s - search_half_width_s) & (time_s <= expected_peak_time_s + search_half_width_s)
    if not np.any(window):
        # Fall back to the full record if the window misses.
        window = np.ones(len(time_s), dtype=bool)
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

@dataclass(frozen=True)
class ConvergenceRow:
    grid_spacing_m: float
    courant_number: float
    theoretical_peak_arrival_s: float
    numerical_peak_arrival_s: float
    percent_error: float


def run_convergence_study(
    config: SimulationConfig,
    dz_values_m: list[float] | None = None,
) -> tuple[list[ConvergenceRow], pd.DataFrame]:
    """Run the homogeneous model at several grid spacings and measure arrival-time error.

    Setup:
        - Source at 200 m, receiver at the free surface (z = 0).  The direct
          upgoing wave arrives cleanly with no free-surface-reflection
          contamination within the measurement window.
        - The Courant number is held constant at 0.5 (dt ∝ dz), so both
          spatial and temporal errors decrease together with grid refinement.
        - Source (200 m) is a multiple of all tested spacings.

    Convergence note:
        The arrival-time error is measured by peak-picking, so it is quantised
        to ±Δt/2.  With dt ∝ dz the measured error decreases as O(Δz) — each
        grid halving halves the error — even though the underlying scheme is
        2nd-order in space.  The O(Δz) slope is expected and is not a sign of
        a 1st-order scheme; it is an artefact of the peak-picking measurement.
        The log-log plot should show a clean straight line with slope ≈ 1.
    """
    from .fd_solver_1d import run_velocity_stress_solver

    if dz_values_m is None:
        dz_values_m = [8.0, 4.0, 2.0, 1.0]

    source_depth_m = 200.0   # multiple of 8, 4, 2, 1
    receiver_depth_m = 0.0   # free surface
    vs = config.rock_shear_velocity_m_s

    rows: list[ConvergenceRow] = []
    for dz in dz_values_m:
        courant = 0.5  # constant by construction
        dt = courant * dz / vs
        refined_config = replace(
            config,
            grid_spacing_m=dz,
            time_step_s=dt,
        )
        material = homogeneous_rock_model(refined_config)
        result = run_velocity_stress_solver(
            config=refined_config,
            material=material,
            source_depth_m=source_depth_m,
            receiver_depth_m=receiver_depth_m,
        )
        check = arrival_time_verification(result)
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
# Stage 5: Layered site-response model
# ---------------------------------------------------------------------------

def run_layered_site_case(config: SimulationConfig) -> tuple[SimulationResult, SimulationResult]:
    """Run rock reference and layered-site simulations with matched source geometry.

    Source placement:
    - In both runs the source node is placed in the rock halfspace, one grid
      cell below the soil-rock interface (z = H + dz).  This ensures the
      identical source excites both models from the same depth.
    - The receiver is at the free surface (z = 0) for both runs.

    Transfer-function definition:
    - The rock-reference run uses a purely homogeneous rock column (no soil
      layer).  Dividing the layered surface spectrum by the reference surface
      spectrum gives a numerical estimate of the site amplification function —
      the ratio of surface motion in the layered site to surface motion in the
      equivalent rock site for the same input.
    - The reference run thus plays the role of the "input motion at the base
      of the soil column" expressed back at the surface through the rock alone.
    """
    from .fd_solver_1d import run_velocity_stress_solver

    source_depth_m = config.soil_layer_thickness_m + config.grid_spacing_m
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


# ---------------------------------------------------------------------------
# Stage 6: Site amplification analysis
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class AmplificationResult:
    """Spectral amplification result for a layered site."""
    frequencies_hz: np.ndarray
    ref_spectrum: np.ndarray
    layered_spectrum: np.ndarray
    amplification: np.ndarray
    resonance_frequencies_hz: np.ndarray
    # frequency band over which the amplification is considered reliable
    f_min_hz: float
    f_max_hz: float

    def resonance_table(self) -> pd.DataFrame:
        rows = [{"mode_n": n + 1, "resonance_frequency_hz": float(f)} for n, f in enumerate(self.resonance_frequencies_hz)]
        return pd.DataFrame(rows)


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
    f_min_hz: float | None = None,
    f_max_hz: float | None = None,
    smoothing_octaves: float = 0.1,
) -> AmplificationResult:
    """Compute spectral amplification A(f) = |V_layered(f)| / |V_ref(f)|.

    Stability measures applied here:
    1. Band-limiting: the ratio is only evaluated where the source has
       meaningful energy (between f_min and f_max derived from the source
       bandwidth).  Outside that band the amplification is set to NaN so it
       is not plotted or interpreted.
    2. Reference-spectrum floor: the denominator is floored at 1 % of its
       in-band peak before dividing, preventing division by values that are
       near zero only because of numerical noise.
    3. Optional octave-band smoothing on the raw spectra before dividing —
       this reduces spectral variance and makes resonance peaks easier to
       identify.
    """
    config = ref_result.config
    freq, spec_ref = surface_fourier_spectrum(ref_result.time_s, ref_result.surface_velocity_m_s)
    _, spec_lay = surface_fourier_spectrum(layered_result.time_s, layered_result.surface_velocity_m_s)

    # Default band: 0.5 f0 to 2.5 * f_max_interpreted
    if f_min_hz is None:
        f_min_hz = max(0.1, 0.5 * config.dominant_frequency_hz)
    if f_max_hz is None:
        f_max_hz = 2.0 * config.maximum_interpreted_frequency_hz

    # Smooth spectra with a running geometric-mean window (octave smoothing).
    if smoothing_octaves > 0.0:
        spec_ref = _octave_smooth(freq, spec_ref, smoothing_octaves)
        spec_lay = _octave_smooth(freq, spec_lay, smoothing_octaves)

    # Floor the reference spectrum at 1 % of its in-band peak.
    in_band = (freq >= f_min_hz) & (freq <= f_max_hz)
    ref_peak_in_band = spec_ref[in_band].max() if in_band.any() else spec_ref.max()
    floor = 0.01 * max(ref_peak_in_band, 1e-30)
    spec_ref_floored = np.maximum(spec_ref, floor)

    amplification = spec_lay / spec_ref_floored

    # Mask outside the usable band with NaN so downstream plots skip it.
    amplification = np.where(in_band, amplification, np.nan)

    resonance_freqs = _quarter_wavelength_frequencies(config, n_modes=5)

    return AmplificationResult(
        frequencies_hz=freq,
        ref_spectrum=spec_ref,
        layered_spectrum=spec_lay,
        amplification=amplification,
        resonance_frequencies_hz=resonance_freqs,
        f_min_hz=float(f_min_hz),
        f_max_hz=float(f_max_hz),
    )


def _octave_smooth(freq: np.ndarray, spectrum: np.ndarray, octaves: float) -> np.ndarray:
    """Apply running geometric-mean smoothing over ±octaves/2 around each frequency."""
    smoothed = spectrum.copy()
    for i, f0 in enumerate(freq):
        if f0 <= 0.0:
            continue
        factor = 2.0 ** (octaves / 2.0)
        mask = (freq >= f0 / factor) & (freq <= f0 * factor)
        if mask.any():
            smoothed[i] = np.exp(np.mean(np.log(spectrum[mask] + 1e-40)))
    return smoothed


def _quarter_wavelength_frequencies(config: SimulationConfig, n_modes: int = 5) -> np.ndarray:
    """Compute quarter-wavelength resonance frequencies for the soft layer."""
    H = config.soil_layer_thickness_m
    Vs_soil = config.soil_shear_velocity_m_s
    return np.array([(2 * n - 1) * Vs_soil / (4.0 * H) for n in range(1, n_modes + 1)])
