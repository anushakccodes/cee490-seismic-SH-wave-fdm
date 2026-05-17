"""Verification and analysis utilities for the 2D SH-wave model."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .fd_solver_2d import Simulation2DResult


@dataclass
class ArrivalCheck2D:
    """Arrival-time comparison for one receiver in a homogeneous medium."""
    receiver_name: str
    source_x_m: float
    source_z_m: float
    receiver_x_m: float
    receiver_z_m: float
    distance_m: float
    vs_m_s: float
    theoretical_peak_s: float
    numerical_peak_s: float
    absolute_error_s: float
    percent_error: float


def theoretical_2d_arrival_time(
    source_x: float,
    source_z: float,
    receiver_x: float,
    receiver_z: float,
    vs: float,
    peak_time_s: float,
) -> tuple[float, float]:
    """Return (distance_m, expected_peak_arrival_s) for a point source."""
    dist = float(np.sqrt((receiver_x - source_x) ** 2 + (receiver_z - source_z) ** 2))
    t_peak = peak_time_s + dist / vs
    return dist, t_peak


def _pick_peak_arrival(time_s: np.ndarray, signal: np.ndarray, expected_t: float, window_s: float = 0.4) -> float:
    """Return the time of the first significant local extremum within the search window.

    The Ricker wavelet produces a positive lobe followed by a larger negative trough.
    Taking the global |max| picks the trailing trough; instead we find the first local
    extremum that exceeds 30% of the windowed absolute maximum.
    """
    mask = np.abs(time_s - expected_t) <= window_s
    if not np.any(mask):
        return float(time_s[np.argmax(np.abs(signal))])

    t_win = time_s[mask]
    s_win = signal[mask]
    threshold = 0.30 * np.max(np.abs(s_win))

    # Find all local extrema in window (sign change in gradient)
    ds = np.diff(s_win)
    sign_change = (ds[:-1] * ds[1:]) < 0
    candidates = np.where(sign_change)[0] + 1  # index into s_win

    # Use a 60% threshold to skip small precursor lobes and find the main arrival
    threshold = 0.60 * np.max(np.abs(s_win))
    for ci in candidates:
        if np.abs(s_win[ci]) >= threshold:
            return float(t_win[ci])

    # Fallback: global |max| in window
    return float(t_win[np.argmax(np.abs(s_win))])


def arrival_time_verification_2d(
    result: Simulation2DResult,
    receiver_name: str,
) -> ArrivalCheck2D:
    """Compare numerical and theoretical arrival times for a named receiver.

    Assumes homogeneous medium — uses rock_shear_velocity_m_s from config.
    """
    config = result.config
    rx, rz = result.receiver_positions[receiver_name]
    sx, sz = result.source_x_m, result.source_z_m
    vs = config.rock_shear_velocity_m_s
    t0 = config.source_peak_time_s

    dist, t_theory = theoretical_2d_arrival_time(sx, sz, rx, rz, vs, t0)

    signal = result.receiver_records_m_s[receiver_name]
    t_num = _pick_peak_arrival(result.time_s, signal, t_theory)

    abs_err = abs(t_num - t_theory)
    pct_err = 100.0 * abs_err / t_theory if t_theory > 0 else float("nan")

    return ArrivalCheck2D(
        receiver_name=receiver_name,
        source_x_m=sx,
        source_z_m=sz,
        receiver_x_m=rx,
        receiver_z_m=rz,
        distance_m=dist,
        vs_m_s=vs,
        theoretical_peak_s=t_theory,
        numerical_peak_s=t_num,
        absolute_error_s=abs_err,
        percent_error=pct_err,
    )


def run_arrival_verification_table(result: Simulation2DResult) -> list[ArrivalCheck2D]:
    """Run arrival-time checks for all receivers that have a valid distance > 0."""
    checks = []
    for name, (rx, rz) in result.receiver_positions.items():
        dist = np.sqrt((rx - result.source_x_m) ** 2 + (rz - result.source_z_m) ** 2)
        if dist < 1.0:
            continue
        checks.append(arrival_time_verification_2d(result, name))
    return checks


def format_arrival_table(checks: list[ArrivalCheck2D]) -> str:
    """Return a formatted arrival-time verification table."""
    header = (
        f"{'Receiver':<14} {'dist(m)':>8} {'t_theory(s)':>12} "
        f"{'t_num(s)':>10} {'err(s)':>8} {'err(%)':>7}\n"
        + "-" * 65
    )
    rows = [header]
    for c in checks:
        rows.append(
            f"{c.receiver_name:<14} {c.distance_m:>8.1f} {c.theoretical_peak_s:>12.4f} "
            f"{c.numerical_peak_s:>10.4f} {c.absolute_error_s:>8.4f} {c.percent_error:>7.2f}"
        )
    return "\n".join(rows)


def compute_surface_peak_amplitudes(result: Simulation2DResult) -> dict[str, float]:
    """Return peak absolute velocity at each receiver."""
    return {
        name: float(np.max(np.abs(signal)))
        for name, signal in result.receiver_records_m_s.items()
    }


def surface_fourier_spectrum(
    result: Simulation2DResult,
    receiver_name: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (frequencies_hz, amplitude_spectrum) for one receiver.

    Amplitude is normalized so the peak equals 1.
    """
    signal = result.receiver_records_m_s[receiver_name]
    dt = float(result.time_s[1] - result.time_s[0])
    freqs = np.fft.rfftfreq(len(signal), d=dt)
    amp = np.abs(np.fft.rfft(signal))
    if amp.max() > 0:
        amp = amp / amp.max()
    return freqs, amp


def compute_2d_amplification(
    result_ref: Simulation2DResult,
    result_layered: Simulation2DResult,
    receiver_name: str,
    f_min: float = 0.5,
    f_max: float = 15.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (frequencies, amplification) = |layered| / |ref| in [f_min, f_max].

    Values outside the band are set to NaN.
    """
    freqs, amp_ref = surface_fourier_spectrum(result_ref, receiver_name)
    _, amp_lay = surface_fourier_spectrum(result_layered, receiver_name)

    in_band = (freqs >= f_min) & (freqs <= f_max)
    floor = 0.01 * amp_ref[in_band].max() if amp_ref[in_band].max() > 0 else 1e-10
    denom = np.where(amp_ref > floor, amp_ref, np.nan)
    ratio = amp_lay / denom
    ratio[~in_band] = np.nan
    return freqs, ratio
