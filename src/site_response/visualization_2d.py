"""Plotting and animation utilities for the 2D SH-wave model."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

from .materials_2d import Material2D, Simulation2DConfig
from .fd_solver_2d import Simulation2DResult

# -------------------------------------------------------------------------
# Style
# -------------------------------------------------------------------------
_COLORS = {
    "deep_blue":  "#1B4F72",
    "teal":       "#148F77",
    "amber":      "#B7770D",
    "crimson":    "#922B21",
    "green":      "#1E8449",
    "purple":     "#6C3483",
    "gray":       "#5D6D7E",
    "white":      "#FFFFFF",
    "near_black": "#1A1A2E",
}

# Colormaps
_CMAP_WAVEFIELD = "coolwarm"   # professional diverging colormap for velocity
_CMAP_MATERIAL  = "cividis"    # perceptually uniform, colorblind-safe for Vs/Z maps


def _style() -> None:
    plt.rcParams.update({
        "figure.dpi": 150,
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.grid": True,
        "grid.alpha": 0.25,
        "grid.linewidth": 0.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.8,
    })


def _save(fig: plt.Figure, output_path: Path | str | None) -> None:
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")


# -------------------------------------------------------------------------
# Shared receiver plotting helper
# -------------------------------------------------------------------------

def _plot_receivers(
    ax: plt.Axes,
    receiver_positions: dict[str, tuple[float, float]],
    surface_z_m: float = 0.0,
    label_receivers: bool = True,
    fontsize: int = 7,
) -> list[Line2D]:
    """Plot all receivers and return legend handles for surface/internal types.

    Surface receivers (z ≈ 0) are shown as inverted triangles (v).
    Internal receivers are shown as circles (o).
    Each marker is annotated with a short label.
    """
    surf_color    = "#F0A500"   # warm gold — visible on dark/light backgrounds
    intern_color  = "#00B4D8"   # sky blue

    added_surf_legend    = False
    added_intern_legend  = False
    handles: list[Line2D] = []

    for name, (rx, rz) in receiver_positions.items():
        is_surface = abs(rz - surface_z_m) < 2.0
        marker = "v" if is_surface else "o"
        color  = surf_color if is_surface else intern_color
        ms     = 7 if is_surface else 6

        h = ax.plot(rx, rz, marker, color=color, ms=ms,
                    markeredgecolor="white", markeredgewidth=0.6,
                    zorder=8, label=None)[0]

        if label_receivers:
            # Short display label: strip common prefix "surf_" if present
            display = name.replace("surf_", "x=") if name.startswith("surf_") else name
            ax.annotate(
                display,
                xy=(rx, rz), xytext=(4, -10),
                textcoords="offset points",
                fontsize=fontsize, color=color,
                zorder=9,
            )

        if is_surface and not added_surf_legend:
            handles.append(Line2D([0], [0], marker="v", color="w",
                                  markerfacecolor=surf_color,
                                  markeredgecolor="white", markeredgewidth=0.6,
                                  markersize=7, label="Surface receiver"))
            added_surf_legend = True
        elif not is_surface and not added_intern_legend:
            handles.append(Line2D([0], [0], marker="o", color="w",
                                  markerfacecolor=intern_color,
                                  markeredgecolor="white", markeredgewidth=0.6,
                                  markersize=6, label="Internal receiver"))
            added_intern_legend = True

    return handles


# -------------------------------------------------------------------------
# Material maps
# -------------------------------------------------------------------------

def plot_vs_map(
    material: Material2D,
    output_path: Path | str | None = None,
    title: str = "Shear-wave velocity map",
) -> plt.Figure:
    """Plot the 2D Vs field with depth on the vertical axis (z positive down)."""
    _style()
    fig, ax = plt.subplots(figsize=(9, 5))
    xx, zz = np.meshgrid(material.x_m, material.z_m)
    im = ax.pcolormesh(xx, zz, material.shear_velocity_m_s,
                       cmap=_CMAP_MATERIAL, shading="auto")
    cbar = fig.colorbar(im, ax=ax, label="Vs  (m/s)", shrink=0.9)
    cbar.ax.tick_params(labelsize=9)
    if material.interface_depth_m is not None:
        ax.axhline(material.interface_depth_m, color=_COLORS["amber"], lw=1.5, ls="--",
                   label=f"Interface  {material.interface_depth_m:.0f} m")
        ax.legend(fontsize=9, framealpha=0.7)
    ax.set_xlabel("x  (m)")
    ax.set_ylabel("Depth  (m)")
    ax.invert_yaxis()
    ax.set_title(title, fontsize=12)
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_impedance_map(
    material: Material2D,
    output_path: Path | str | None = None,
) -> plt.Figure:
    """Plot the 2D impedance field Z = rho * Vs."""
    _style()
    fig, ax = plt.subplots(figsize=(9, 5))
    xx, zz = np.meshgrid(material.x_m, material.z_m)
    im = ax.pcolormesh(xx, zz, material.impedance_kg_m2_s,
                       cmap=_CMAP_MATERIAL, shading="auto")
    cbar = fig.colorbar(im, ax=ax, label="Z = ρVs  (kg/m²/s)", shrink=0.9)
    cbar.ax.tick_params(labelsize=9)
    if material.interface_depth_m is not None:
        ax.axhline(material.interface_depth_m, color=_COLORS["amber"], lw=1.5, ls="--",
                   label=f"Interface  {material.interface_depth_m:.0f} m")
        ax.legend(fontsize=9, framealpha=0.7)
    ax.set_xlabel("x  (m)")
    ax.set_ylabel("Depth  (m)")
    ax.invert_yaxis()
    ax.set_title("Impedance map  Z = ρVs", fontsize=12)
    fig.tight_layout()
    _save(fig, output_path)
    return fig


# -------------------------------------------------------------------------
# Source wavelet
# -------------------------------------------------------------------------

def plot_source_wavelet_2d(
    result: Simulation2DResult,
    output_path: Path | str | None = None,
) -> plt.Figure:
    """Plot the Ricker source-time function used in the 2D run."""
    _style()
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(result.time_s, result.source_time_function,
            color=_COLORS["deep_blue"], lw=1.5)
    ax.axvline(result.config.source_peak_time_s, color=_COLORS["amber"],
               ls="--", lw=1.2, label=f"Peak  t₀ = {result.config.source_peak_time_s:.3f} s")
    ax.set_xlabel("Time  (s)")
    ax.set_ylabel("Normalized amplitude")
    ax.set_title(f"Ricker source  f₀ = {result.config.dominant_frequency_hz:.1f} Hz",
                 fontsize=12)
    ax.legend(fontsize=9)
    fig.tight_layout()
    _save(fig, output_path)
    return fig


# -------------------------------------------------------------------------
# Velocity snapshot
# -------------------------------------------------------------------------

def plot_velocity_snapshot(
    result: Simulation2DResult,
    snap_index: int,
    output_path: Path | str | None = None,
    overlay_circle: bool = False,
) -> plt.Figure:
    """Plot one 2D velocity snapshot with labelled source and receivers."""
    _style()
    field = result.velocity_snapshots_m_s[snap_index]
    t     = result.stored_time_s[snap_index]
    vmax  = float(np.max(np.abs(field))) or 1.0

    fig, ax = plt.subplots(figsize=(10, 5))
    xx, zz  = np.meshgrid(result.x_m, result.z_m)
    norm    = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im      = ax.pcolormesh(xx, zz, field, cmap=_CMAP_WAVEFIELD, norm=norm, shading="auto")
    cbar    = fig.colorbar(im, ax=ax, label="vy  (m/s)", shrink=0.9)
    cbar.ax.tick_params(labelsize=9)

    # Source
    ax.plot(result.source_x_m, result.source_z_m, "*",
            color=_COLORS["white"], ms=12,
            markeredgecolor=_COLORS["near_black"], markeredgewidth=0.8,
            zorder=8, label="Source")

    # Receivers (labelled)
    rec_handles = _plot_receivers(ax, result.receiver_positions, label_receivers=True)

    # Interface
    interface_handles: list[Line2D] = []
    if result.material.interface_depth_m is not None:
        ax.axhline(result.material.interface_depth_m,
                   color=_COLORS["amber"], lw=1.4, ls="--", zorder=7)
        interface_handles = [
            Line2D([0], [0], color=_COLORS["amber"], lw=1.4, ls="--",
                   label=f"Interface  {result.material.interface_depth_m:.0f} m")
        ]

    # Theoretical wavefront
    circle_handles: list[Line2D] = []
    if overlay_circle:
        t0 = result.config.source_peak_time_s
        vs = result.config.rock_shear_velocity_m_s
        r  = vs * max(t - t0, 0.0)
        theta = np.linspace(0, 2 * np.pi, 360)
        cx = result.source_x_m + r * np.cos(theta)
        cz = result.source_z_m + r * np.sin(theta)
        ax.plot(cx, cz, color=_COLORS["green"], lw=1.3, ls="--", zorder=7)
        circle_handles = [
            Line2D([0], [0], color=_COLORS["green"], lw=1.3, ls="--",
                   label=f"Theory  r = {r:.0f} m")
        ]

    src_handle = Line2D([0], [0], marker="*", color="w",
                        markerfacecolor=_COLORS["white"],
                        markeredgecolor=_COLORS["near_black"],
                        markeredgewidth=0.8, markersize=11, label="Source")
    ax.legend(handles=[src_handle] + rec_handles + interface_handles + circle_handles,
              fontsize=8, loc="lower right", framealpha=0.8)

    ax.set_xlabel("x  (m)")
    ax.set_ylabel("Depth  (m)")
    ax.invert_yaxis()
    ax.set_title(f"Velocity snapshot  t = {t:.3f} s  —  {result.material.label}",
                 fontsize=12)
    fig.tight_layout()
    _save(fig, output_path)
    return fig


# -------------------------------------------------------------------------
# Receiver records
# -------------------------------------------------------------------------

def plot_receiver_records(
    result: Simulation2DResult,
    output_path: Path | str | None = None,
    max_receivers: int = 8,
) -> plt.Figure:
    """Plot time histories for all (or up to max_receivers) receivers."""
    _style()
    names = list(result.receiver_records_m_s.keys())[:max_receivers]
    n = len(names)
    fig, axes = plt.subplots(n, 1, figsize=(10, 2.0 * n), sharex=True)
    if n == 1:
        axes = [axes]

    palette = plt.cm.tab10(np.linspace(0, 1, max(n, 1)))
    for ax, name, color in zip(axes, names, palette):
        signal = result.receiver_records_m_s[name]
        rx, rz = result.receiver_positions[name]
        ax.plot(result.time_s, signal, color=color, lw=1.0)
        ax.set_ylabel("vy  (m/s)", fontsize=8)
        ax.set_title(f"{name}   x = {rx:.0f} m,  z = {rz:.0f} m", fontsize=9)
    axes[-1].set_xlabel("Time  (s)")
    fig.suptitle(f"Receiver records — {result.material.label}", fontsize=11, y=1.01)
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_surface_array(
    result: Simulation2DResult,
    output_path: Path | str | None = None,
    surface_z_m: float = 0.0,
) -> plt.Figure:
    """Plot all surface receivers as a wiggle-trace gather."""
    _style()
    surface_recs = {
        name: (rx, rz)
        for name, (rx, rz) in result.receiver_positions.items()
        if abs(rz - surface_z_m) < 1.0
    }
    if not surface_recs:
        fig, ax = plt.subplots()
        ax.set_title("No surface receivers found")
        return fig

    names_sorted = sorted(surface_recs, key=lambda n: surface_recs[n][0])
    fig, ax = plt.subplots(figsize=(10, 5))
    palette = plt.cm.viridis(np.linspace(0.1, 0.9, len(names_sorted)))

    for name, color in zip(names_sorted, palette):
        sig  = result.receiver_records_m_s[name]
        rx   = surface_recs[name][0]
        norm = np.max(np.abs(sig)) or 1.0
        ax.plot(result.time_s, sig / norm * 50 + rx, color=color, lw=0.9)

    ax.set_xlabel("Time  (s)")
    ax.set_ylabel("Receiver x  (m)  [traces normalized + offset]")
    ax.set_title(f"Surface receiver gather — {result.material.label}", fontsize=12)
    fig.tight_layout()
    _save(fig, output_path)
    return fig


# -------------------------------------------------------------------------
# Arrival-time check figure
# -------------------------------------------------------------------------

def plot_arrival_check(
    result: Simulation2DResult,
    receiver_name: str,
    theoretical_t: float,
    numerical_t: float,
    output_path: Path | str | None = None,
) -> plt.Figure:
    """Plot receiver seismogram with theoretical and numerical arrival marks."""
    _style()
    signal = result.receiver_records_m_s[receiver_name]
    rx, rz = result.receiver_positions[receiver_name]
    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(result.time_s, signal, color=_COLORS["deep_blue"], lw=1.2, label="Numerical")
    ax.axvline(theoretical_t, color=_COLORS["green"], lw=1.5, ls="--",
               label=f"Theory  t = {theoretical_t:.4f} s")
    ax.axvline(numerical_t, color=_COLORS["crimson"], lw=1.5, ls=":",
               label=f"Numerical peak  t = {numerical_t:.4f} s")
    ax.set_xlabel("Time  (s)")
    ax.set_ylabel("vy  (m/s)")
    ax.set_title(
        f"Arrival check — {receiver_name}   x = {rx:.0f} m,  z = {rz:.0f} m",
        fontsize=12,
    )
    ax.legend(fontsize=9)
    fig.tight_layout()
    _save(fig, output_path)
    return fig


# -------------------------------------------------------------------------
# 2D wavefield animation
# -------------------------------------------------------------------------

def make_2d_wave_animation(
    result: Simulation2DResult,
    output_path: Path | str | None = None,
    fps: int = 15,
    dpi: int = 120,
    skip_blank: bool = True,
) -> animation.FuncAnimation:
    """Create an animation of the 2D velocity snapshots with labelled markers.

    Surface receivers are shown as gold inverted triangles (v) with x-position
    labels.  Internal receivers are shown as sky-blue circles (o).
    The source is a white star.  An interface line is drawn when present.

    Parameters
    ----------
    skip_blank : bool
        Skip leading frames where the wavefield is below 1% of the global peak.
    """
    _style()
    snapshots = result.velocity_snapshots_m_s
    vmax = float(np.max(np.abs(snapshots))) or 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

    # Find first frame with visible signal
    start_frame = 0
    if skip_blank:
        threshold = 0.01 * vmax
        for i, snap in enumerate(snapshots):
            if np.max(np.abs(snap)) > threshold:
                start_frame = max(0, i - 2)
                break

    frames = list(range(start_frame, len(snapshots)))
    xx, zz = np.meshgrid(result.x_m, result.z_m)

    fig, ax = plt.subplots(figsize=(10, 5))
    im = ax.pcolormesh(xx, zz, snapshots[start_frame],
                       cmap=_CMAP_WAVEFIELD, norm=norm, shading="auto")
    cbar = fig.colorbar(im, ax=ax, label="vy  (m/s)", shrink=0.9)
    cbar.ax.tick_params(labelsize=9)

    # Source marker
    ax.plot(result.source_x_m, result.source_z_m, "*",
            color=_COLORS["white"], ms=12,
            markeredgecolor=_COLORS["near_black"], markeredgewidth=0.8,
            zorder=8)

    # Receivers — labelled, colour-coded by type
    rec_handles = _plot_receivers(ax, result.receiver_positions,
                                  label_receivers=True, fontsize=7)

    # Interface line
    interface_handles: list[Line2D] = []
    if result.material.interface_depth_m is not None:
        ax.axhline(result.material.interface_depth_m,
                   color=_COLORS["amber"], lw=1.4, ls="--", zorder=7)
        interface_handles = [
            Line2D([0], [0], color=_COLORS["amber"], lw=1.4, ls="--",
                   label=f"Interface  {result.material.interface_depth_m:.0f} m")
        ]

    # Build legend (static — drawn once)
    src_handle = Line2D([0], [0], marker="*", color="w",
                        markerfacecolor=_COLORS["white"],
                        markeredgecolor=_COLORS["near_black"],
                        markeredgewidth=0.8, markersize=11, label="Source")
    ax.legend(handles=[src_handle] + rec_handles + interface_handles,
              fontsize=8, loc="lower right", framealpha=0.8)

    ax.set_xlabel("x  (m)")
    ax.set_ylabel("Depth  (m)")
    ax.invert_yaxis()
    title = ax.set_title("")
    fig.tight_layout()

    def _update(fi: int) -> list:
        snap = snapshots[fi]
        im.set_array(snap.ravel())
        title.set_text(
            f"{result.material.label} — t = {result.stored_time_s[fi]:.3f} s"
        )
        return [im, title]

    ani = animation.FuncAnimation(
        fig, _update, frames=frames, interval=1000 // fps, blit=True,
    )

    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        ani.save(str(output_path), writer="pillow", fps=fps, dpi=dpi)

    return ani
