"""Consistent plotting utilities for the site-response project."""
from __future__ import annotations

from pathlib import Path
import matplotlib.pyplot as plt
import matplotlib.animation as animation
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
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    return "Arial" if "Arial" in available_fonts else "DejaVu Sans"


def apply_project_style() -> None:
    plt.rcParams.update({
        "font.family": _preferred_font_family(),
        "font.size": 11,
        "axes.titlesize": 13,
        "axes.labelsize": 11,
        "legend.fontsize": 10,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "figure.titlesize": 14,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "lines.linewidth": 1.8,
        "figure.dpi": 120,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })


def save_figure(fig: plt.Figure, output_path: Path | str | None) -> None:
    if output_path is None:
        return
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)


def plot_velocity_profile(material: MaterialProfile, output_path: Path | str | None = None) -> plt.Figure:
    apply_project_style()
    fig, ax = plt.subplots(figsize=(4.5, 5.5))
    ax.plot(material.shear_velocity_m_s, material.depth_m, color=PROJECT_COLORS["deep_blue"], label="Vs(z)")
    if material.interface_depth_m is not None:
        ax.axhline(material.interface_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", label="soil-rock interface")
    ax.set_xlabel("Shear-wave velocity (m/s)")
    ax.set_ylabel("Depth (m)")
    ax.set_title("Vs Profile")
    ax.invert_yaxis()
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_impedance_profile(material: MaterialProfile, output_path: Path | str | None = None) -> plt.Figure:
    apply_project_style()
    fig, ax = plt.subplots(figsize=(4.5, 5.5))
    ax.plot(material.impedance_kg_m2_s, material.depth_m, color=PROJECT_COLORS["teal"], label="Z(z) = ρVs")
    if material.interface_depth_m is not None:
        ax.axhline(material.interface_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", label="soil-rock interface")
    ax.set_xlabel("Shear impedance (kg/m²/s)")
    ax.set_ylabel("Depth (m)")
    ax.set_title("Impedance Profile")
    ax.invert_yaxis()
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_source_time_history(time_s: np.ndarray, source_values: np.ndarray, output_path: Path | str | None = None) -> plt.Figure:
    apply_project_style()
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    ax.plot(time_s, source_values, color=PROJECT_COLORS["purple"], label="Ricker wavelet")
    ax.axhline(0.0, color=PROJECT_COLORS["gray"], linewidth=0.8)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Normalized amplitude")
    ax.set_title("Source Time Function")
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_source_spectrum(frequencies_hz: np.ndarray, normalized_amplitude: np.ndarray, dominant_frequency_hz: float, output_path: Path | str | None = None) -> plt.Figure:
    apply_project_style()
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    ax.plot(frequencies_hz, normalized_amplitude, color=PROJECT_COLORS["green"], label="|S(f)|")
    ax.axvline(dominant_frequency_hz, color=PROJECT_COLORS["orange"], linestyle="--", label=f"f₀ = {dominant_frequency_hz:g} Hz")
    ax.set_xlim(0.0, 4.0 * dominant_frequency_hz)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Normalized amplitude")
    ax.set_title("Source Spectrum")
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_homogeneous_seismogram(result: SimulationResult, arrival_check: ArrivalTimeCheck, output_path: Path | str | None = None) -> plt.Figure:
    apply_project_style()
    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    ax.plot(result.time_s, result.receiver_velocity_m_s, color=PROJECT_COLORS["deep_blue"], label=f"receiver at {result.receiver_depth_m:.0f} m")
    ax.axvline(arrival_check.theoretical_peak_arrival_s, color=PROJECT_COLORS["orange"], linestyle="--", label="theoretical arrival")
    ax.axvline(arrival_check.numerical_peak_arrival_s, color=PROJECT_COLORS["crimson"], linestyle=":", label="numerical arrival")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Particle velocity (m/s)")
    ax.set_title("Homogeneous-Medium Arrival-Time Verification")
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_wavefield_image(result: SimulationResult, output_path: Path | str | None = None, title: str = "Velocity Wavefield") -> plt.Figure:
    apply_project_style()
    fig, ax = plt.subplots(figsize=(8.0, 4.8))
    maximum_abs = np.max(np.abs(result.velocity_wavefield_m_s)) or 1.0
    image = ax.imshow(
        result.velocity_wavefield_m_s.T,
        extent=[result.stored_time_s[0], result.stored_time_s[-1], result.depth_m[-1], result.depth_m[0]],
        aspect="auto",
        cmap="RdBu_r",
        vmin=-maximum_abs,
        vmax=maximum_abs,
    )
    ax.axhline(result.source_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", linewidth=1.4, label="source depth")
    ax.axhline(result.receiver_depth_m, color=PROJECT_COLORS["green"], linestyle=":", linewidth=1.4, label="receiver depth")
    if result.material.interface_depth_m is not None:
        ax.axhline(result.material.interface_depth_m, color=PROJECT_COLORS["crimson"], linestyle="-.", linewidth=1.4, label="soil-rock interface")
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Depth (m)")
    ax.set_title(title)
    ax.legend(loc="upper right", fontsize=9)
    fig.colorbar(image, ax=ax, label="Particle velocity (m/s)", pad=0.02)
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_convergence(dz_values: np.ndarray, errors: np.ndarray, output_path: Path | str | None = None) -> plt.Figure:
    """Log-log convergence plot: arrival-time error vs grid spacing.

    The reference line has slope 1 because peak-picking quantises the error to
    ±Δt/2, so with constant Courant number the measured error scales as O(Δz).
    """
    apply_project_style()
    fig, ax = plt.subplots(figsize=(5.5, 4.0))
    positive = errors > 0
    if positive.any():
        ax.loglog(dz_values[positive], errors[positive], "o-",
                  color=PROJECT_COLORS["deep_blue"], label="numerical error")
        # Reference slope-1 line anchored at coarsest valid point.
        ref_x = np.array([dz_values[positive].min(), dz_values[positive].max()])
        coarse_idx = np.argmax(dz_values[positive])
        ref_y = errors[positive][coarse_idx] * (ref_x / dz_values[positive][coarse_idx]) ** 1
        ax.loglog(ref_x, ref_y, "--", color=PROJECT_COLORS["gray"], label="slope 1 reference")
    ax.set_xlabel("Grid spacing Δz (m)")
    ax.set_ylabel("Arrival-time error (%)")
    ax.set_title("Grid-Refinement Convergence")
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_surface_comparison(
    ref_result: SimulationResult,
    layered_result: SimulationResult,
    output_path: Path | str | None = None,
) -> plt.Figure:
    """Overlay surface velocity time histories for rock reference vs layered site."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    ax.plot(ref_result.time_s, ref_result.surface_velocity_m_s, color=PROJECT_COLORS["gray"], label="rock reference", linewidth=1.6)
    ax.plot(layered_result.time_s, layered_result.surface_velocity_m_s, color=PROJECT_COLORS["crimson"], label="layered site", linewidth=1.6)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Surface velocity (m/s)")
    ax.set_title("Surface Motion: Rock Reference vs Layered Site")
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_spectra_comparison(
    frequencies_hz: np.ndarray,
    ref_spectrum: np.ndarray,
    layered_spectrum: np.ndarray,
    output_path: Path | str | None = None,
) -> plt.Figure:
    """Plot Fourier amplitude spectra of rock reference and layered site surface motion."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(7.0, 3.8))
    ax.semilogy(frequencies_hz, ref_spectrum + 1e-20, color=PROJECT_COLORS["gray"], label="|V_ref(f)|")
    ax.semilogy(frequencies_hz, layered_spectrum + 1e-20, color=PROJECT_COLORS["crimson"], label="|V_layered(f)|")
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Spectral amplitude (m/s)")
    ax.set_title("Surface Motion Spectra")
    ax.set_xlim(0.0, 30.0)
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_amplification(
    frequencies_hz: np.ndarray,
    amplification: np.ndarray,
    resonance_frequencies_hz: np.ndarray,
    output_path: Path | str | None = None,
    f_min_hz: float | None = None,
    f_max_hz: float | None = None,
) -> plt.Figure:
    """Plot site amplification factor A(f) with resonance frequency markers.

    Regions outside the reliable bandwidth (NaN in amplification) are not
    plotted.  The y-axis is capped at 3× the in-band peak to prevent
    near-zero reference values from dominating the scale.
    """
    apply_project_style()
    fig, ax = plt.subplots(figsize=(7.0, 3.8))

    # Only plot where amplification is finite (in-band mask applied upstream).
    valid = np.isfinite(amplification)
    if valid.any():
        ax.plot(frequencies_hz[valid], amplification[valid],
                color=PROJECT_COLORS["deep_blue"], label="A(f)")

    ax.axhline(1.0, color=PROJECT_COLORS["gray"], linewidth=0.8, linestyle="--", label="A = 1")

    # Shade the out-of-band regions to communicate where the ratio is unreliable.
    if f_min_hz is not None:
        ax.axvspan(0.0, f_min_hz, color=PROJECT_COLORS["gray"], alpha=0.12)
    if f_max_hz is not None:
        ax.axvspan(f_max_hz, frequencies_hz[-1], color=PROJECT_COLORS["gray"], alpha=0.12)

    for n, fn in enumerate(resonance_frequencies_hz, start=1):
        label = f"f{n} = {fn:.2f} Hz" if n <= 3 else None
        ax.axvline(fn, color=PROJECT_COLORS["orange"], linestyle=":", linewidth=1.2, label=label)

    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Amplification A(f)")
    ax.set_title("Site Amplification Factor")

    x_max = f_max_hz * 1.2 if f_max_hz is not None else 30.0
    ax.set_xlim(0.0, min(x_max, frequencies_hz[-1]))

    # Cap y-axis at 3× in-band peak (or 10 as a hard ceiling).
    if valid.any():
        in_band_peak = np.nanmax(amplification[valid])
        ax.set_ylim(bottom=0.0, top=min(3.0 * in_band_peak, 10.0))
    else:
        ax.set_ylim(bottom=0.0)

    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def _first_active_frame(wavefield: np.ndarray, threshold_fraction: float = 0.01) -> int:
    """Return the stored-frame index when the wavefield first exceeds a noise threshold."""
    global_max = np.max(np.abs(wavefield))
    if global_max == 0.0:
        return 0
    threshold = threshold_fraction * global_max
    for i, frame in enumerate(wavefield):
        if np.max(np.abs(frame)) >= threshold:
            return i
    return 0


def make_wave_animation(
    result: SimulationResult,
    output_path: Path | str | None = None,
    fps: int = 20,
    max_frames: int = 200,
) -> animation.FuncAnimation:
    """Animate the 1D velocity profile evolving through time.

    Frames before the wave arrives (amplitude < 1 % of peak) are skipped so
    the animation does not start with a blank screen.
    """
    apply_project_style()
    n_stored = result.velocity_wavefield_m_s.shape[0]
    first = _first_active_frame(result.velocity_wavefield_m_s)
    usable = n_stored - first
    step = max(1, usable // max_frames)
    frame_indices = list(range(first, n_stored, step))

    depth_m = result.depth_m
    wavefield = result.velocity_wavefield_m_s
    vmax = np.max(np.abs(wavefield)) or 1.0

    fig, ax = plt.subplots(figsize=(4.5, 5.5))
    line, = ax.plot([], [], color=PROJECT_COLORS["deep_blue"], linewidth=1.6)
    ax.set_xlim(-vmax, vmax)
    ax.set_ylim(depth_m[-1], depth_m[0])
    ax.set_xlabel("Particle velocity (m/s)")
    ax.set_ylabel("Depth (m)")
    ax.set_title("Wave Propagation")
    ax.axhline(result.source_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", linewidth=1.0, label="source")
    if result.material.interface_depth_m is not None:
        ax.axhline(result.material.interface_depth_m, color=PROJECT_COLORS["crimson"], linestyle="-.", linewidth=1.0, label="interface")
    ax.legend(loc="upper right", fontsize=8)
    time_text = ax.text(0.02, 0.97, "", transform=ax.transAxes, fontsize=9, va="top")

    def init():
        line.set_data([], [])
        time_text.set_text("")
        return line, time_text

    def update(frame_idx):
        v = wavefield[frame_idx]
        line.set_data(v, depth_m)
        time_text.set_text(f"t = {result.stored_time_s[frame_idx]:.3f} s")
        return line, time_text

    anim = animation.FuncAnimation(fig, update, frames=frame_indices, init_func=init, blit=True, interval=1000 // fps)
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        anim.save(str(output_path), writer="pillow", fps=fps)
    plt.close(fig)
    return anim


def make_layered_reflection_animation(
    result: SimulationResult,
    output_path: Path | str | None = None,
    fps: int = 20,
    max_frames: int = 200,
) -> animation.FuncAnimation:
    """Animate 1D wave reflection/reverberation in the layered model.

    The depth window is focused on the soft layer and immediate surroundings
    to make reflection/reverberation visible.  Blank early frames are skipped.
    """
    apply_project_style()
    n_stored = result.velocity_wavefield_m_s.shape[0]
    first = _first_active_frame(result.velocity_wavefield_m_s)
    usable = n_stored - first
    step = max(1, usable // max_frames)
    frame_indices = list(range(first, n_stored, step))

    depth_m = result.depth_m
    wavefield = result.velocity_wavefield_m_s
    vmax = np.max(np.abs(wavefield)) or 1.0

    interface = result.material.interface_depth_m
    if interface is not None:
        zoom_bottom = min(depth_m[-1], interface * 3.0)
    else:
        zoom_bottom = depth_m[-1]

    fig, ax = plt.subplots(figsize=(4.5, 5.5))
    line, = ax.plot([], [], color=PROJECT_COLORS["deep_blue"], linewidth=1.6)
    ax.set_xlim(-vmax, vmax)
    ax.set_ylim(zoom_bottom, 0.0)
    ax.set_xlabel("Particle velocity (m/s)")
    ax.set_ylabel("Depth (m)")
    ax.set_title("Layered-Site Reflection/Reverberation")
    ax.axhline(0.0, color=PROJECT_COLORS["gray"], linewidth=0.8, linestyle="-", label="free surface")
    if interface is not None:
        ax.axhline(interface, color=PROJECT_COLORS["crimson"], linestyle="-.", linewidth=1.2, label=f"interface ({interface:.0f} m)")
    ax.legend(loc="upper right", fontsize=8)
    time_text = ax.text(0.02, 0.97, "", transform=ax.transAxes, fontsize=9, va="top")

    def init():
        line.set_data([], [])
        time_text.set_text("")
        return line, time_text

    def update(frame_idx):
        v = wavefield[frame_idx]
        line.set_data(v, depth_m)
        time_text.set_text(f"t = {result.stored_time_s[frame_idx]:.3f} s")
        return line, time_text

    anim = animation.FuncAnimation(fig, update, frames=frame_indices, init_func=init, blit=True, interval=1000 // fps)
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        anim.save(str(output_path), writer="pillow", fps=fps)
    plt.close(fig)
    return anim
