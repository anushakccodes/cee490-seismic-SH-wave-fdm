"""Consistent plotting utilities for the site-response project."""
from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib import font_manager
import numpy as np

from .analysis import ArrivalTimeCheck
from .fd_solver_1d import SimulationResult
from .materials import MaterialProfile

PROJECT_COLORS = {
    "deep_blue": "#0B3954",
    "teal": "#087E8B",
    "orange": "#E36414",
    "crimson": "#C52233",
    "green": "#2E933C",
    "purple": "#4F345A",
    "gray": "#5C677D",
}


def _preferred_font_family() -> str:
    """Use Arial when available; otherwise fall back to DejaVu Sans."""
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    return "Arial" if "Arial" in available_fonts else "DejaVu Sans"


def apply_project_style() -> None:
    """Apply the shared Matplotlib style used by all project figures."""
    plt.rcParams.update({
        "font.family": _preferred_font_family(),
        "font.size": 13,
        "axes.titlesize": 17,
        "axes.labelsize": 15,
        "legend.fontsize": 12,
        "xtick.labelsize": 12,
        "ytick.labelsize": 12,
        "figure.titlesize": 19,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "lines.linewidth": 2.3,
        "figure.dpi": 130,
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
    })


def save_figure(fig: plt.Figure, output_path: Path | str | None) -> None:
    """Save a figure when an output path is supplied."""
    if output_path is None:
        return
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)


def plot_velocity_profile(material: MaterialProfile, output_path: Path | str | None = None) -> plt.Figure:
    """Plot shear-wave velocity versus depth."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(6.5, 7.0))
    ax.plot(material.shear_velocity_m_s, material.depth_m, color=PROJECT_COLORS["deep_blue"], label="Vs(z)")
    if material.interface_depth_m is not None:
        ax.axhline(material.interface_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", label="soil-rock interface")
    ax.set_xlabel("Shear-wave velocity, Vs (m/s)")
    ax.set_ylabel("Depth, z (m)")
    ax.set_title("Shear-Wave Velocity Profile")
    ax.invert_yaxis()
    ax.legend(loc="best")
    save_figure(fig, output_path)
    return fig


def plot_impedance_profile(material: MaterialProfile, output_path: Path | str | None = None) -> plt.Figure:
    """Plot shear impedance Z = rho * Vs versus depth."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(6.5, 7.0))
    ax.plot(material.impedance_kg_m2_s, material.depth_m, color=PROJECT_COLORS["teal"], label="Z(z) = rho Vs")
    if material.interface_depth_m is not None:
        ax.axhline(material.interface_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", label="soil-rock interface")
    ax.set_xlabel("Shear impedance, Z (kg/m^2/s)")
    ax.set_ylabel("Depth, z (m)")
    ax.set_title("Shear Impedance Profile")
    ax.invert_yaxis()
    ax.legend(loc="best")
    save_figure(fig, output_path)
    return fig


def plot_source_time_history(time_s: np.ndarray, source_values: np.ndarray, output_path: Path | str | None = None) -> plt.Figure:
    """Plot the Ricker source wavelet in the time domain."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.plot(time_s, source_values, color=PROJECT_COLORS["purple"], label="Ricker wavelet")
    ax.axhline(0.0, color=PROJECT_COLORS["gray"], linewidth=1.2)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Normalized source amplitude")
    ax.set_title("Input Source Time Function")
    ax.legend(loc="best")
    save_figure(fig, output_path)
    return fig


def plot_source_spectrum(frequencies_hz: np.ndarray, normalized_amplitude: np.ndarray, dominant_frequency_hz: float, output_path: Path | str | None = None) -> plt.Figure:
    """Plot the normalized Fourier amplitude spectrum of the source."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(9.0, 4.8))
    ax.plot(frequencies_hz, normalized_amplitude, color=PROJECT_COLORS["green"], label="|S(f)|")
    ax.axvline(dominant_frequency_hz, color=PROJECT_COLORS["orange"], linestyle="--", label=f"f0 = {dominant_frequency_hz:g} Hz")
    ax.set_xlim(0.0, 4.0 * dominant_frequency_hz)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Normalized amplitude")
    ax.set_title("Source Frequency Spectrum")
    ax.legend(loc="best")
    save_figure(fig, output_path)
    return fig


def plot_homogeneous_seismogram(result: SimulationResult, arrival_check: ArrivalTimeCheck, output_path: Path | str | None = None) -> plt.Figure:
    """Plot receiver motion with numerical and theoretical arrival markers."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(10.5, 5.0))
    ax.plot(result.time_s, result.receiver_velocity_m_s, color=PROJECT_COLORS["deep_blue"], label=f"receiver at {result.receiver_depth_m:.0f} m")
    ax.axvline(arrival_check.theoretical_peak_arrival_s, color=PROJECT_COLORS["orange"], linestyle="--", label="theoretical peak arrival")
    ax.axvline(arrival_check.numerical_peak_arrival_s, color=PROJECT_COLORS["crimson"], linestyle=":", label="numerical peak arrival")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Particle velocity (model units)")
    ax.set_title("Homogeneous-Medium Arrival-Time Verification")
    ax.legend(loc="best")
    save_figure(fig, output_path)
    return fig


def plot_wavefield_image(result: SimulationResult, output_path: Path | str | None = None) -> plt.Figure:
    """Plot a depth-time image of particle velocity."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(10.5, 6.2))
    maximum_abs = np.max(np.abs(result.velocity_wavefield_m_s)) or 1.0
    image = ax.imshow(
        result.velocity_wavefield_m_s.T,
        extent=[result.stored_time_s[0], result.stored_time_s[-1], result.depth_m[-1], result.depth_m[0]],
        aspect="auto",
        cmap="RdBu_r",
        vmin=-maximum_abs,
        vmax=maximum_abs,
    )
    ax.axhline(result.source_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", linewidth=1.8, label="source depth")
    ax.axhline(result.receiver_depth_m, color=PROJECT_COLORS["green"], linestyle=":", linewidth=1.8, label="receiver depth")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Depth, z (m)")
    ax.set_title("Homogeneous Velocity Wavefield")
    ax.legend(loc="best")
    fig.colorbar(image, ax=ax, label="Particle velocity (model units)")
    save_figure(fig, output_path)
    return fig
