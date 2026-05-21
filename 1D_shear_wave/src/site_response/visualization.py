"""Consistent plotting utilities."""
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


def _preferred_font_family():
    available_fonts = {font.name for font in font_manager.fontManager.ttflist}
    return "Arial" if "Arial" in available_fonts else "DejaVu Sans"


def apply_project_style():
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
        "axes.facecolor": "white",
        "figure.facecolor": "white",
        "grid.alpha": 0.25,
        "lines.linewidth": 1.8,
        "figure.dpi": 120,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })


def save_figure(fig, output_path):
    if output_path is None:
        return
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path)
    plt.close(fig)


def plot_velocity_profile(material, output_path=None):
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


def plot_impedance_profile(material, output_path=None):
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


def plot_profiles(material, velocity_output_path=None, impedance_output_path=None, combined_output_path=None):
    apply_project_style()
    fig, axes = plt.subplots(2, 2, figsize=(9.0, 8.25), sharey=True)
    (ax1, ax2), (ax3, ax4) = axes

    y_min = material.depth_m.min()
    y_max = material.depth_m.max()
    for ax in (ax1, ax2, ax3, ax4):
        if material.interface_depth_m is not None:
            h = material.interface_depth_m
            ax.axhspan(y_min, h, color="#e8d5b0", zorder=0, label="Soil")
            ax.axhspan(h, y_max, color="#dcdcdc", zorder=0, label="Rock")

    ax1.plot(material.shear_velocity_m_s, material.depth_m, color=PROJECT_COLORS["deep_blue"], label="Vs(z)", zorder=2)
    if material.interface_depth_m is not None:
        ax1.axhline(material.interface_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", label="soil-rock interface", zorder=3)
    ax1.set_xlabel("Shear-wave velocity (m/s)")
    ax1.set_ylabel("Depth (m)")
    ax1.set_title("Shear-Wave Velocity Profile")
    ax1.invert_yaxis()
    ax1.legend(loc="best")

    ax2.plot(material.density_kg_m3, material.depth_m, color=PROJECT_COLORS["crimson"], label="ρ(z)", zorder=2)
    if material.interface_depth_m is not None:
        ax2.axhline(material.interface_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", label="soil-rock interface", zorder=3)
    ax2.set_xlabel("Density (kg/m³)")
    ax2.set_title("Density Profile")
    ax2.legend(loc="best")

    ax3.plot(material.impedance_kg_m2_s, material.depth_m, color=PROJECT_COLORS["teal"], label="Z(z) = ρVs", zorder=2)
    if material.interface_depth_m is not None:
        ax3.axhline(material.interface_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", label="soil-rock interface", zorder=3)
    ax3.set_xlabel("Shear impedance (kg/m²/s)")
    ax3.set_ylabel("Depth (m)")
    ax3.set_title("Impedance Profile")
    ax3.legend(loc="best")

    ax4.plot(material.shear_modulus_pa, material.depth_m, color=PROJECT_COLORS["purple"], label="μ(z) = ρVs²", zorder=2)
    if material.interface_depth_m is not None:
        ax4.axhline(material.interface_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", label="soil-rock interface", zorder=3)
    ax4.set_xlabel("Shear modulus (Pa)")
    ax4.set_title("Material Stiffness Profile")
    ax4.legend(loc="best")

    fig.tight_layout()
    save_figure(fig, combined_output_path)
    return fig


def plot_source_time_history(time_s, source_values, output_path=None):
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


def plot_source_spectrum(frequencies_hz, normalized_amplitude, dominant_frequency_hz, output_path=None):
    apply_project_style()
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    ax.plot(frequencies_hz, normalized_amplitude, color=PROJECT_COLORS["green"], label="|S(f)|")
    ax.axvline(dominant_frequency_hz, color=PROJECT_COLORS["orange"], linestyle="--", label=f"$f_0$ = {dominant_frequency_hz:g} Hz")
    ax.set_xlim(0.0, 4.0 * dominant_frequency_hz)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Normalized amplitude")
    ax.set_title("Source Spectrum")
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_homogeneous_seismogram(result, arrival_check, output_path=None):
    apply_project_style()
    fig, (ax, ax_zoom) = plt.subplots(1, 2, figsize=(13.0, 3.8), gridspec_kw={"width_ratios": [2, 1]})

    t_th = arrival_check.theoretical_peak_arrival_s
    t_num = arrival_check.numerical_peak_arrival_s

    for a in (ax, ax_zoom):
        a.plot(result.time_s, result.receiver_velocity_m_s, color=PROJECT_COLORS["deep_blue"],
               label=f"receiver at {result.receiver_depth_m:.0f} m")
        a.axvline(t_th,  color="black", linestyle="--", linewidth=1.5, label=f"theoretical  ({t_th:.4f} s)")
        a.axvline(t_num, color="red",   linestyle=":",  linewidth=1.5, label=f"numerical  ({t_num:.4f} s)")

    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Particle velocity (m/s)")
    ax.set_title("Homogeneous-Medium Arrival-Time Verification (z=0m)")
    ax.legend(loc="best")

    # zoom window: ±5 dt around the theoretical peak
    half_win = max(5 * result.config.time_step_s, abs(t_num - t_th) * 3, 0.005)
    ax_zoom.set_xlim(t_th - half_win, t_th + half_win)
    ax_zoom.set_xlabel("Time (s)")
    ax_zoom.set_title("Zoomed View — Arrival Difference")
    ax_zoom.legend(loc="best", fontsize=8)

    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_wavefield_image(result, output_path=None, title="Velocity Wavefield"):
    apply_project_style()
    from matplotlib.ticker import MultipleLocator

    # Stress is defined on the staggered grid (N-1 points); interpolate to velocity grid for display.
    stress_depth_m = 0.5 * (result.depth_m[:-1] + result.depth_m[1:])
    stress_on_vel_grid = np.apply_along_axis(
        lambda row: np.interp(result.depth_m, stress_depth_m, row),
        axis=1, arr=result.stress_wavefield_pa,
    )

    extent = [result.stored_time_s[0], result.stored_time_s[-1], result.depth_m[-1], result.depth_m[0]]

    fig, (ax, ax_s) = plt.subplots(1, 2, figsize=(15.0, 4.8), sharey=True)

    absorb_start = result.depth_m[-1] - result.config.absorbing_layer_thickness_m

    def _decorate(a):
        a.axhspan(absorb_start, result.depth_m[-1], color="lightgray", alpha=0.35, zorder=1, label="absorbing layer")
        a.axhline(result.source_depth_m, color=PROJECT_COLORS["orange"], linestyle="--", linewidth=1.4, label="source depth", zorder=3)
        a.axhline(result.receiver_depth_m, color=PROJECT_COLORS["green"], linestyle=":", linewidth=1.4, label="receiver depth", zorder=3)
        if result.material.interface_depth_m is not None:
            a.axhline(result.material.interface_depth_m, color=PROJECT_COLORS["crimson"], linestyle="-.", linewidth=1.4, label="soil-rock interface", zorder=3)
        a.set_xlabel("Time (s)")
        a.xaxis.set_minor_locator(MultipleLocator(0.1))
        a.yaxis.set_minor_locator(MultipleLocator(10))
        a.grid(which="major", color="gray", linewidth=0.6, alpha=0.5)
        a.grid(which="minor", color="gray", linewidth=0.3, alpha=0.25)
        a.legend(loc="upper right", fontsize=9)

    for a in (ax, ax_s):
        a.set_facecolor("white")

    vmax = np.max(np.abs(result.velocity_wavefield_m_s)) or 1.0
    im_v = ax.imshow(result.velocity_wavefield_m_s.T, extent=extent, aspect="auto",
                     cmap="seismic", vmin=-vmax, vmax=vmax)
    ax.set_ylabel("Depth (m)")
    ax.set_title(title)
    _decorate(ax)
    fig.colorbar(im_v, ax=ax, label="Particle velocity (m/s)", pad=0.02)

    smax = np.max(np.abs(stress_on_vel_grid)) or 1.0
    im_s = ax_s.imshow(stress_on_vel_grid.T, extent=extent, aspect="auto",
                       cmap="bwr", vmin=-smax, vmax=smax)
    ax_s.set_title("Shear-Stress Wavefield")
    _decorate(ax_s)
    fig.colorbar(im_s, ax=ax_s, label="Shear stress (Pa)", pad=0.02)

    # Extend y-axis 20 m above the surface as visual padding so the 0 m receiver
    # line is not clipped; suppress ticks in that padding region.
    ax.set_ylim(result.depth_m[-1], -20)
    ax.yaxis.set_major_locator(MultipleLocator(50))
    ax.yaxis.set_minor_locator(MultipleLocator(10))

    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_convergence(dz_values, errors, l2_errors=None, output_path=None):
    """Two-panel log-log convergence plot.

    Left panel  — peak-picking arrival-time error (%): slope ≈ 1 because the
                  measurement is quantised to ±Δt/2 and Δt ∝ Δz.
    Right panel — Richardson L2 receiver error (relative): slope ≈ 2 because
                  comparing against the finest-grid reference at the same receiver
                  bypasses the quantisation and exposes the true leapfrog truncation
                  error order.
    """
    apply_project_style()
    dz = np.asarray(dz_values, dtype=float)

    has_l2 = l2_errors is not None and np.any(np.isfinite(l2_errors))
    ncols = 2 if has_l2 else 1
    fig, axes = plt.subplots(1, ncols, figsize=(6.0 * ncols, 4.5))
    if ncols == 1:
        axes = [axes]
    ax_pp, *rest = axes
    ax_l2 = rest[0] if rest else None

    # --- Left: peak-picking error ---
    pos = errors > 0
    if pos.any():
        ax_pp.loglog(dz[pos], errors[pos], "o-",
                     color=PROJECT_COLORS["deep_blue"], alpha=0.9,
                     label="peak-picking error")
        ref_x = np.array([dz[pos].min(), dz[pos].max()])
        ref_y = errors[pos][np.argmax(dz[pos])] * (ref_x / dz[pos].max()) ** 1
        ax_pp.loglog(ref_x, ref_y, "--",
                     color=PROJECT_COLORS["crimson"], alpha=0.7,
                     label="slope 1 reference")
    ax_pp.set_xlabel("Grid spacing Δz (m)")
    ax_pp.set_ylabel("Arrival-time error (%)")
    ax_pp.set_title("Peak-Picking Error\n(measurement-limited, slope ≈ 1)")
    ax_pp.legend(loc="upper left")

    # --- Right: L2 receiver error ---
    if ax_l2 is not None:
        l2 = np.asarray(l2_errors, dtype=float)
        pos2 = np.isfinite(l2) & (l2 > 0)
        if pos2.any():
            ax_l2.loglog(dz[pos2], l2[pos2], "s-",
                         color=PROJECT_COLORS["orange"], alpha=0.9,
                         label="L2 receiver error")
            ref_x2 = np.array([dz[pos2].min(), dz[pos2].max()])
            ref_y2 = l2[pos2][np.argmax(dz[pos2])] * (ref_x2 / dz[pos2].max()) ** 2
            ax_l2.loglog(ref_x2, ref_y2, "--",
                         color=PROJECT_COLORS["crimson"], alpha=0.7,
                         label="slope 2 reference")
        ax_l2.set_xlabel("Grid spacing Δz (m)")
        ax_l2.set_ylabel("Relative L2 error (receiver time series)")
        ax_l2.set_title("L2 Field Error\n(Richardson comparison, slope ≈ 2)")
        ax_l2.legend(loc="upper left")

    fig.suptitle("Grid-Refinement Convergence Study", fontsize=13)
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_surface_comparison(ref_result, layered_result, output_path=None):
    """Overlay surface velocity time histories for rock reference vs layered site."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(8.0, 3.8))
    ax.plot(ref_result.time_s, ref_result.surface_velocity_m_s, color=PROJECT_COLORS["gray"], label="Homogeneous rock reference", linewidth=1.6)
    ax.plot(layered_result.time_s, layered_result.surface_velocity_m_s, color=PROJECT_COLORS["crimson"], label="Layered site", linewidth=1.6)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Surface velocity (m/s)")
    ax.set_title("Surface Motion Comparison")
    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_spectra_comparison(frequencies_hz, ref_spectrum, layered_spectrum, output_path=None):
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


def plot_amplification(frequencies_hz, amplification, resonance_frequencies_hz, output_path=None, f_min_hz=None, f_max_hz=None):
    """Plot site amplification factor A(f) with resonance frequency markers.
    Regions outside the reliable bandwidth (NaN in amplification) are not
    plotted.  The y-axis is capped at 3× the in-band peak to prevent
    near-zero reference values from dominating the scale.
    """
    apply_project_style()
    fig, ax = plt.subplots(figsize=(7.0, 3.8))

    valid = np.isfinite(amplification)
    if valid.any():
        ax.plot(frequencies_hz[valid], amplification[valid],
                color=PROJECT_COLORS["deep_blue"], label="A(f)")

    ax.axhline(1.0, color=PROJECT_COLORS["gray"], linewidth=0.8, linestyle="--", label="A = 1")

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

    if valid.any():
        in_band_peak = np.nanmax(amplification[valid])
        ax.set_ylim(bottom=0.0, top=min(3.0 * in_band_peak, 10.0))
    else:
        ax.set_ylim(bottom=0.0)

    ax.legend(loc="best")
    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def plot_spectra_and_amplification(frequencies_hz, ref_spectrum, layered_spectrum,
                                   amplification, resonance_frequencies_hz,
                                   output_path=None, f_max_hz=None):
    """Side-by-side: left = surface motion spectra, right = site amplification factor."""
    from matplotlib.ticker import MultipleLocator
    apply_project_style()
    fig, (ax_spec, ax_amp) = plt.subplots(1, 2, figsize=(14.0, 4.2))

    # Left: surface motion spectra
    ax_spec.semilogy(frequencies_hz, ref_spectrum + 1e-20,
                     color=PROJECT_COLORS["gray"], label="|V_ref(f)|")
    ax_spec.semilogy(frequencies_hz, layered_spectrum + 1e-20,
                     color=PROJECT_COLORS["crimson"], label="|V_layered(f)|")
    ax_spec.set_xlabel("Frequency (Hz)")
    ax_spec.set_ylabel("Spectral amplitude (m/s)")
    ax_spec.set_title("Surface Motion Spectra")
    ax_spec.set_xlim(0.0, 18.0)
    ax_spec.xaxis.set_major_locator(MultipleLocator(2))
    ax_spec.legend(loc="best")

    # Right: site amplification factor
    valid = np.isfinite(amplification)
    if valid.any():
        ax_amp.plot(frequencies_hz[valid], amplification[valid],
                    color=PROJECT_COLORS["deep_blue"], label="A(f)")
    ax_amp.axhline(1.0, color=PROJECT_COLORS["gray"], linewidth=0.8, linestyle="--", label="A = 1")
    for n, fn in enumerate(resonance_frequencies_hz, start=1):
        ax_amp.axvline(fn, color=PROJECT_COLORS["orange"], linestyle=":", linewidth=1.2, label=f"f{n} = {fn:.2f} Hz")
    ax_amp.set_xlabel("Frequency (Hz)")
    ax_amp.set_ylabel("Amplification A(f)")
    ax_amp.set_title("Site Amplification Factor")
    ax_amp.set_xlim(0.0, 18.0)
    ax_amp.xaxis.set_major_locator(MultipleLocator(2))
    if valid.any():
        in_band_peak = np.nanmax(amplification[valid])
        ax_amp.set_ylim(bottom=0.0, top=min(3.0 * in_band_peak, 10.0))
    else:
        ax_amp.set_ylim(bottom=0.0)
    ax_amp.legend(loc="best")

    fig.tight_layout()
    save_figure(fig, output_path)
    return fig


def make_comparison_animation(ref_result, layered_result, output_path=None, fps=20, max_frames=200):
    """Three-panel animation: homogeneous (left), layered (centre), overlay (right).

    All panels share the same depth axis and velocity scale so amplitude
    differences and wave trapping are directly comparable.
    """
    apply_project_style()

    ref_wf = ref_result.velocity_wavefield_m_s
    lay_wf = layered_result.velocity_wavefield_m_s
    n_stored = min(ref_wf.shape[0], lay_wf.shape[0])

    first = min(_first_active_frame(ref_wf), _first_active_frame(lay_wf))
    usable = n_stored - first
    step = max(1, usable // max_frames)
    frame_indices = list(range(first, n_stored, step))

    depth_m = ref_result.depth_m
    vmax = max(np.max(np.abs(ref_wf)), np.max(np.abs(lay_wf))) or 1.0

    fig, (ax_ref, ax_lay, ax_both) = plt.subplots(1, 3, figsize=(13.5, 4.4), sharey=True)

    line_ref,  = ax_ref.plot( [], [], color=PROJECT_COLORS["gray"],    linewidth=1.6, zorder=3)
    line_lay,  = ax_lay.plot( [], [], color=PROJECT_COLORS["crimson"], linewidth=1.6, zorder=3)
    line_ref2, = ax_both.plot([], [], color=PROJECT_COLORS["gray"],    linewidth=1.6, zorder=3, label="Homogeneous")
    line_lay2, = ax_both.plot([], [], color=PROJECT_COLORS["crimson"], linewidth=1.6, zorder=3, label="Layered")

    interface_m = layered_result.material.interface_depth_m
    y_top = depth_m[0]
    y_bot = depth_m[-1]

    for ax, title, show_soil in [
        (ax_ref,  "Homogeneous Rock", False),
        (ax_lay,  "Layered Site",     True),
        (ax_both, "Overlay",          True),
    ]:
        # Material background colouring
        if show_soil and interface_m is not None:
            ax.axhspan(y_top, interface_m, color="#e8d5b0", zorder=0, label="Soil")
            ax.axhspan(interface_m, y_bot, color="#dcdcdc", zorder=0, label="Rock")
        else:
            ax.axhspan(y_top, y_bot, color="#dcdcdc", zorder=0, label="Rock")

        ax.set_xlim(-vmax, vmax)
        ax.set_ylim(y_bot, y_top)
        ax.set_xlabel("Particle velocity (m/s)")
        ax.set_title(title)
        ax.axhline(ref_result.source_depth_m, color=PROJECT_COLORS["orange"],
                   linestyle="--", linewidth=1.0, zorder=2,
                   label=f"source ({ref_result.source_depth_m:.0f} m)")
        if interface_m is not None:
            ax.axhline(interface_m, color=PROJECT_COLORS["orange"],
                       linestyle="-.", linewidth=1.0, zorder=2,
                       label=f"interface ({interface_m:.0f} m)")
        ax.legend(loc="lower right", fontsize=8)

    ax_ref.set_ylabel("Depth (m)")
    time_text = fig.text(0.5, 0.01, "", ha="center", fontsize=10)
    fig.suptitle("SH-Wave Propagation Comparison", fontsize=13)
    fig.tight_layout(rect=[0, 0.04, 1, 1])

    def init():
        for ln in (line_ref, line_lay, line_ref2, line_lay2):
            ln.set_data([], [])
        time_text.set_text("")
        return line_ref, line_lay, line_ref2, line_lay2, time_text

    def update(frame_idx):
        line_ref.set_data( ref_wf[frame_idx], depth_m)
        line_lay.set_data( lay_wf[frame_idx], depth_m)
        line_ref2.set_data(ref_wf[frame_idx], depth_m)
        line_lay2.set_data(lay_wf[frame_idx], depth_m)
        time_text.set_text(f"t = {ref_result.stored_time_s[frame_idx]:.3f} s")
        return line_ref, line_lay, line_ref2, line_lay2, time_text

    anim = animation.FuncAnimation(
        fig, update, frames=frame_indices, init_func=init, blit=True, interval=1000 // fps
    )
    if output_path is not None:
        output_path = Path(output_path)
        output_path.parent.mkdir(parents=True, exist_ok=True)
        anim.save(str(output_path), writer="pillow", fps=fps)
    plt.close(fig)
    return anim


def _first_active_frame(wavefield, threshold_fraction=0.01):
    """Return the stored-frame index when the wavefield first exceeds a noise threshold."""
    global_max = np.max(np.abs(wavefield))
    if global_max == 0.0:
        return 0
    threshold = threshold_fraction * global_max
    for i, frame in enumerate(wavefield):
        if np.max(np.abs(frame)) >= threshold:
            return i
    return 0


def make_wave_animation(result, output_path=None, fps=20, max_frames=200):
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


def make_layered_reflection_animation(result, output_path=None, fps=20, max_frames=200):
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
