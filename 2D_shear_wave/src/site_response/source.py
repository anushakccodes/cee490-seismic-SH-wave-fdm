"""Ricker source-time functions and source spectra."""
from __future__ import annotations

import numpy as np


def ricker_wavelet(time_s, dominant_frequency_hz, peak_time_s):
    """Return unit-normalized Ricker wavelet; peak amplitude = 1."""
    shifted_time = time_s - peak_time_s
    argument = np.pi * dominant_frequency_hz * shifted_time
    wav = (1.0 - 2.0 * argument**2) * np.exp(-argument**2)
    peak = np.max(np.abs(wav))
    return wav / peak if peak > 0.0 else wav


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