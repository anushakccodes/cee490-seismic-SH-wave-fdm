"""Verification and numerical-analysis utilities."""
from __future__ import annotations
from dataclasses import dataclass, replace
import numpy as np
import pandas as pd

from .fd_solver_1d import SimulationResult
from .materials import homogeneous_rock_model, layered_soil_over_rock_model
from .parameters import SimulationConfig


# Stage 3: Homogeneous arrival-time verification
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

    def as_dataframe(self):
        return pd.DataFrame([self.__dict__])


def theoretical_travel_time(distance_m, shear_velocity_m_s):
    return distance_m / shear_velocity_m_s


def pick_peak_arrival_time(time_s, signal, expected_peak_time_s, search_half_width_s=0.15):
    """Return the time of the largest absolute peak within a window around the expected arrival."""
    window = (time_s >= expected_peak_time_s - search_half_width_s) & (time_s <= expected_peak_time_s + search_half_width_s)
    if not np.any(window):
        window = np.ones(len(time_s), dtype=bool)
    local_time = time_s[window]
    local_signal = np.abs(signal[window])
    return float(local_time[int(np.argmax(local_signal))])


def arrival_time_verification(result):
    """Compare numerical peak arrival time against t_theory = t0 + d/Vs"""
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


def run_homogeneous_verification_case(config):
    from .fd_solver_1d import run_velocity_stress_solver
    material = homogeneous_rock_model(config)
    result = run_velocity_stress_solver(
        config=config,
        material=material,
        source_depth_m=config.source_depth_m,
        receiver_depth_m=config.receiver_depth_m,
    )
    return result, arrival_time_verification(result)


# Stage 4: Grid-refinement convergence study
@dataclass(frozen=True)
class ConvergenceRow:
    grid_spacing_m: float
    time_step_s: float
    courant_number: float
    theoretical_peak_arrival_s: float
    numerical_peak_arrival_s: float
    percent_error: float
    l2_error_relative: float


def compute_l2_planewave_error(dz, courant, vs, rho, f0, n_steps=200):
    """Truncation-error measurement for the staggered leapfrog scheme.

    Initialises velocity and stress from the exact analytical solution for a
    downward-propagating Ricker plane wave on the staggered grid, advances the
    interior leapfrog stencil for n_steps, then compares to the exact final
    state.  No source injection, no boundaries — pure truncation error.

    Staggered-in-time convention (matching the solver):
      velocity  evaluated at t = -dt/2
      stress    evaluated at t =  0
    For a downward wave: tau = -rho * Vs * v (from the PDE).

    On smooth Ricker data the staggered leapfrog achieves slope >= 2 (typically
    ~3 due to cancellation of the leading error term on smooth, band-limited
    waveforms).  This confirms the scheme is at least second-order accurate.
    """
    from .source import ricker_wavelet

    dt = courant * dz / vs
    mu = rho * vs**2

    wave_travel = vs * n_steps * dt
    margin      = 6.0 * vs / f0
    domain      = wave_travel + 2.0 * margin
    n_nodes     = int(np.ceil(domain / dz)) + 1
    z           = np.arange(n_nodes) * dz
    z_s         = z[:-1] + 0.5 * dz
    z_c         = 0.5 * domain

    def v_exact(z_arr, t):
        return ricker_wavelet(t - (z_arr - z_c) / vs, f0, 0.0)

    velocity     = v_exact(z,   -dt / 2.0)
    shear_stress = -rho * vs * v_exact(z_s, 0.0)

    for _ in range(n_steps):
        velocity[1:-1] += (dt / (rho * dz)) * (shear_stress[1:] - shear_stress[:-1])
        shear_stress   += (mu * dt / dz) * np.diff(velocity)

    t_final = n_steps * dt - dt / 2.0
    v_ex    = v_exact(z, t_final)

    norm_ex = np.linalg.norm(v_ex)
    if norm_ex == 0.0:
        return float("nan")
    return float(np.linalg.norm(velocity - v_ex) / norm_ex)


def run_convergence_study(config, dz_values_m=None):
    """Run the homogeneous solver at decreasing grid spacings with constant Courant number.

    Returns two error metrics per grid spacing:
    - percent_error: arrival-time error by peak-picking (quantised to ±dt/2,
      so log-log slope ≈ 1 with constant Courant number).
    - l2_error_relative: exact plane-wave L2 error using a source-free interior
      leapfrog run initialised from the analytical solution.  This isolates the
      pure truncation error with no source-positioning or boundary artefacts,
      giving a clean log-log slope of ≈ 2.
    """
    from .fd_solver_1d import run_velocity_stress_solver

    if dz_values_m is None:
        dz_values_m = [8.0, 4.0, 2.0, 1.0]

    source_depth_m   = config.source_depth_m
    receiver_depth_m = config.receiver_depth_m
    vs  = config.rock_shear_velocity_m_s
    rho = config.rock_density_kg_m3
    f0  = config.dominant_frequency_hz
    courant = 0.5

    rows = []
    for dz in dz_values_m:
        dt = courant * dz / vs
        refined_config = replace(config, grid_spacing_m=dz, time_step_s=dt)
        material = homogeneous_rock_model(refined_config)
        result   = run_velocity_stress_solver(
            config=refined_config,
            material=material,
            source_depth_m=source_depth_m,
            receiver_depth_m=receiver_depth_m,
        )
        check  = arrival_time_verification(result)
        l2_err = compute_l2_planewave_error(dz, courant, vs, rho, f0)
        rows.append(ConvergenceRow(
            grid_spacing_m=dz,
            time_step_s=dt,
            courant_number=courant,
            theoretical_peak_arrival_s=round(check.theoretical_peak_arrival_s, 5),
            numerical_peak_arrival_s=round(check.numerical_peak_arrival_s, 5),
            percent_error=check.percent_error,
            l2_error_relative=l2_err,
        ))

    df = pd.DataFrame([r.__dict__ for r in rows])
    df["theoretical_peak_arrival_s"] = df["theoretical_peak_arrival_s"].map("{:.5f}".format)
    df["numerical_peak_arrival_s"]   = df["numerical_peak_arrival_s"].map("{:.5f}".format)
    return rows, df


# Stage 5: Layered site-response model
def run_layered_site_case(config):
    """Run rock reference and layered-site simulations with identical source geometry.

    Dividing the layered surface spectrum by the reference surface spectrum gives
    the numerical site amplification function A(f).
    """
    from .fd_solver_1d import run_velocity_stress_solver

    receiver_depth_m = config.receiver_depth_m

    ref_material = homogeneous_rock_model(config)
    ref_result = run_velocity_stress_solver(
        config=config,
        material=ref_material,
        source_depth_m=config.source_depth_m,
        receiver_depth_m=receiver_depth_m,
    )

    layered_material = layered_soil_over_rock_model(config)
    layered_result = run_velocity_stress_solver(
        config=config,
        material=layered_material,
        source_depth_m=config.source_depth_m,
        receiver_depth_m=receiver_depth_m,
    )

    return ref_result, layered_result


# Stage 6: Site amplification analysis
@dataclass(frozen=True)
class AmplificationResult:
    """Spectral amplification result for a layered site."""
    frequencies_hz: np.ndarray
    ref_spectrum: np.ndarray
    layered_spectrum: np.ndarray
    amplification: np.ndarray
    resonance_frequencies_hz: np.ndarray
    f_min_hz: float
    f_max_hz: float

    def resonance_table(self):
        rows = [{"mode_n": n + 1, "resonance_frequency_hz": float(f)} for n, f in enumerate(self.resonance_frequencies_hz)]
        return pd.DataFrame(rows)


def surface_fourier_spectrum(time_s, surface_velocity):
    """Return one-sided Fourier amplitude spectrum of the surface velocity."""
    dt = float(time_s[1] - time_s[0])
    frequencies_hz = np.fft.rfftfreq(len(surface_velocity), d=dt)
    amplitude = np.abs(np.fft.rfft(surface_velocity))
    return frequencies_hz, amplitude


def compute_amplification(ref_result, layered_result, f_min_hz=None, f_max_hz=None, smoothing_octaves=0.1):
    """Compute A(f) = |V_layered(f)| / |V_ref(f)| within the source frequency band.

    The reference spectrum is floored at 1% of its in-band peak before dividing
    to avoid instability at frequencies with negligible source energy.
    """
    config = ref_result.config
    freq, spec_ref = surface_fourier_spectrum(ref_result.time_s, ref_result.surface_velocity_m_s)
    _, spec_lay = surface_fourier_spectrum(layered_result.time_s, layered_result.surface_velocity_m_s)

    if f_min_hz is None:
        f_min_hz = max(0.1, 0.5 * config.dominant_frequency_hz)
    if f_max_hz is None:
        f_max_hz = 2.0 * config.maximum_interpreted_frequency_hz

    if smoothing_octaves > 0.0:
        spec_ref = _octave_smooth(freq, spec_ref, smoothing_octaves)
        spec_lay = _octave_smooth(freq, spec_lay, smoothing_octaves)

    in_band = (freq >= f_min_hz) & (freq <= f_max_hz)
    ref_peak_in_band = spec_ref[in_band].max() if in_band.any() else spec_ref.max()
    floor = 0.01 * max(ref_peak_in_band, 1e-30)
    spec_ref_floored = np.maximum(spec_ref, floor)

    amplification = spec_lay / spec_ref_floored
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


def _octave_smooth(freq, spectrum, octaves):
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


def _quarter_wavelength_frequencies(config, n_modes=5):
    """Compute quarter-wavelength resonance frequencies for the soft layer."""
    H = config.soil_layer_thickness_m
    Vs_soil = config.soil_shear_velocity_m_s
    return np.array([(2 * n - 1) * Vs_soil / (4.0 * H) for n in range(1, n_modes + 1)])
