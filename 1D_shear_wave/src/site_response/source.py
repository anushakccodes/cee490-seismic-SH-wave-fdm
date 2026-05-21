"""Ricker source-time functions and source spectra."""
from __future__ import annotations

import numpy as np
from .parameters import SimulationConfig, make_time_axis


def ricker_wavelet(time_s, dominant_frequency_hz, peak_time_s):
    """Return s(t) = (1 - 2*pi^2*f0^2*(t-t0)^2) exp(-pi^2*f0^2*(t-t0)^2)."""
    shifted_time = time_s - peak_time_s
    argument = np.pi * dominant_frequency_hz * shifted_time
    return (1.0 - 2.0 * argument**2) * np.exp(-argument**2)


def build_source_time_function(config):
    time_s = make_time_axis(config)
    source_values = ricker_wavelet(time_s, config.dominant_frequency_hz, config.source_peak_time_s)
    return time_s, source_values


def normalized_amplitude_spectrum(time_s, signal):
    """Return a one-sided normalized Fourier amplitude spectrum."""
    if len(time_s) < 2:
        raise ValueError("At least two time samples are required.")
    dt = float(time_s[1] - time_s[0])
    frequencies_hz = np.fft.rfftfreq(len(signal), d=dt)
    amplitude = np.abs(np.fft.rfft(signal))
    if amplitude.max() > 0.0:
        amplitude = amplitude / amplitude.max()
    return frequencies_hz, amplitude