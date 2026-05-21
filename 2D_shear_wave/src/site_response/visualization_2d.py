"""Plotting and animation utilities for the 2D SH-wave model."""
from __future__ import annotations

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from matplotlib import font_manager
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D

from .materials_2d import Material2D, Simulation2DConfig
from .fd_solver_2d import Simulation2DResult

PROJECT_COLORS = {
    "deep_blue":  "#0B3954",
    "teal":       "#087E8B",
    "orange":     "#E36414",
    "crimson":    "#C52233",
    "green":      "#2E933C",
    "purple":     "#4F345A",
    "gray":       "#5C677D",
    "white":      "#FFFFFF",
    "near_black": "#1A1A2E",
}

_CMAP_WAVEFIELD = "seismic"
_CMAP_MATERIAL  = "cividis"

FONT_SIZES = {
    "base":        11,
    "axes_label":  12,
    "axes_title":  13,
    "fig_title":   14,
    "colorbar":    10,
    "annotation":   9,
    "small":        9,
}


def _preferred_font_family():
    available = {f.name for f in font_manager.fontManager.ttflist}
    return "Arial" if "Arial" in available else "DejaVu Sans"


def apply_project_style():
    """Apply consistent rcParams for all project figures."""
    plt.rcParams.update({
        "font.family":      _preferred_font_family(),
        "font.size":        FONT_SIZES["base"],
        "axes.titlesize":   FONT_SIZES["axes_title"],
        "axes.labelsize":   FONT_SIZES["axes_label"],
        "legend.fontsize":  FONT_SIZES["base"],
        "xtick.labelsize":  FONT_SIZES["base"],
        "ytick.labelsize":  FONT_SIZES["base"],
        "figure.titlesize": FONT_SIZES["fig_title"],
        "axes.grid":        True,
        "axes.facecolor":   "white",
        "figure.facecolor": "white",
        "grid.alpha":       0.25,
        "grid.linewidth":   0.5,
        "lines.linewidth":  1.8,
        "figure.dpi":       120,
        "savefig.dpi":      150,
        "savefig.bbox":     "tight",
    })


def _save(fig, output_path):
    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output_path, bbox_inches="tight")


def _draw_interface(ax, material, lw=1.5, ls="--"):
    """Draw the soil-rock boundary on ax. Returns a legend handle or None."""
    oc = PROJECT_COLORS["orange"]
    if material.interface_depth_profile_m is not None:
        ax.plot(material.x_m, material.interface_depth_profile_m,
                color=oc, lw=lw, ls=ls, zorder=7)
        from matplotlib.lines import Line2D as _L2D
        return _L2D([0], [0], color=oc, lw=lw, ls=ls, label="Basin boundary")
    if material.interface_depth_m is not None:
        ax.axhline(material.interface_depth_m, color=oc, lw=lw, ls=ls, zorder=7)
        from matplotlib.lines import Line2D as _L2D
        return _L2D([0], [0], color=oc, lw=lw, ls=ls,
                    label=f"Interface  {material.interface_depth_m:.0f} m")
    return None


def _plot_receivers(ax, receiver_positions, surface_z_m=0.0, label_receivers=True):
    surf_color   = PROJECT_COLORS["green"]
    intern_color = PROJECT_COLORS["teal"]

    added_surf_legend   = False
    added_intern_legend = False
    handles = []

    for name, (rx, rz) in receiver_positions.items():
        is_surface = abs(rz - surface_z_m) < 2.0
        marker = "v" if is_surface else "o"
        color  = surf_color if is_surface else intern_color
        ms     = 7 if is_surface else 6

        h = ax.plot(rx, rz, marker, color=color, ms=ms,
                    markeredgecolor="white", markeredgewidth=0.6,
                    zorder=8, label=None)[0]

        if label_receivers:
            display = name.replace("surf_", "x=") if name.startswith("surf_") else name
            ax.annotate(
                display,
                xy=(rx, rz), xytext=(4, -10),
                textcoords="offset points",
                fontsize=FONT_SIZES["annotation"], color=color,
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


def plot_material_overview(materials, config, titles=None):
    """3-row × 3-column material overview figure.

    Row 1 (maps):  Vs map for each model.
    Row 2 (maps):  Impedance map Z = ρVs for each model.
    Row 3 (maps):  Shear-modulus map μ = ρVs² for each model.

    Parameters
    ----------
    materials : sequence of Material2D, length 3
        [homogeneous, layered, basin] models.
    config : Simulation2DConfig
    titles : sequence of str, optional
        Column titles; defaults to material labels.
    """
    apply_project_style()
    if titles is None:
        titles = [m.label for m in materials]

    vs_min = config.soil_shear_velocity_m_s
    vs_max = config.rock_shear_velocity_m_s
    z_min  = min(m.impedance_kg_m2_s.min() for m in materials)
    z_max  = max(m.impedance_kg_m2_s.max() for m in materials)
    mu_min = min(m.shear_modulus_pa.min() for m in materials)
    mu_max = max(m.shear_modulus_pa.max() for m in materials)

    fig, axes = plt.subplots(3, 3, figsize=(18, 13), constrained_layout=True)

    for col, (mat, title) in enumerate(zip(materials, titles)):
        xx, zz = np.meshgrid(mat.x_m, mat.z_m)

        # ── Row 0: Vs map ──────────────────────────────────────────────────
        ax = axes[0, col]
        im = ax.pcolormesh(xx, zz, mat.shear_velocity_m_s,
                           cmap=_CMAP_MATERIAL, shading="auto",
                           vmin=vs_min, vmax=vs_max)
        cbar = fig.colorbar(im, ax=ax, label="$V_s$  (m/s)", shrink=0.9)
        cbar.ax.tick_params(labelsize=FONT_SIZES["colorbar"])
        h = _draw_interface(ax, mat)
        if h is not None:
            ax.legend(handles=[h], fontsize=FONT_SIZES["small"], framealpha=0.7)
        ax.invert_yaxis()
        ax.set_title(f"Shear Velocity Map ({title})")
        ax.set_xlabel("x  (m)")
        ax.set_ylabel("Depth  (m)")

        # ── Row 1: Impedance map ───────────────────────────────────────────
        ax = axes[1, col]
        im2 = ax.pcolormesh(xx, zz, mat.impedance_kg_m2_s,
                            cmap=_CMAP_MATERIAL, shading="auto",
                            vmin=z_min, vmax=z_max)
        cbar = fig.colorbar(im2, ax=ax, label="Z = ρ$V_s$  (kg/m²/s)", shrink=0.9)
        cbar.ax.tick_params(labelsize=FONT_SIZES["colorbar"])
        h = _draw_interface(ax, mat)
        if h is not None:
            ax.legend(handles=[h], fontsize=FONT_SIZES["small"], framealpha=0.7)
        ax.invert_yaxis()
        ax.set_xlabel("x  (m)")
        ax.set_ylabel("Depth  (m)")

        # ── Row 2: Shear-modulus map ───────────────────────────────────────
        ax = axes[2, col]
        im3 = ax.pcolormesh(xx, zz, mat.shear_modulus_pa,
                            cmap=_CMAP_MATERIAL, shading="auto",
                            vmin=mu_min, vmax=mu_max)
        cbar = fig.colorbar(im3, ax=ax, label="μ = ρ$V_s^2$  (Pa)", shrink=0.9)
        cbar.ax.tick_params(labelsize=FONT_SIZES["colorbar"])
        h = _draw_interface(ax, mat)
        if h is not None:
            ax.legend(handles=[h], fontsize=FONT_SIZES["small"], framealpha=0.7)
        ax.invert_yaxis()
        ax.set_xlabel("x  (m)")
        ax.set_ylabel("Depth  (m)")

    axes[0, 0].set_ylabel("Depth  (m)")
    axes[0, 0].set_xlabel("x  (m)")

    for col, title in enumerate(titles):
        axes[2, col].set_title(f"Shear modulus ({title})")
        axes[1, col].set_title(f"Impedance ({title})")

    return fig


def plot_vs_map(material, output_path=None, title="Shear-wave velocity map"):
    """Plot the 2D Vs field with depth on the vertical axis (z positive down)."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(9, 5))
    xx, zz = np.meshgrid(material.x_m, material.z_m)
    im = ax.pcolormesh(xx, zz, material.shear_velocity_m_s,
                       cmap=_CMAP_MATERIAL, shading="auto")
    cbar = fig.colorbar(im, ax=ax, label="Vs  (m/s)", shrink=0.9)
    cbar.ax.tick_params(labelsize=FONT_SIZES["colorbar"])
    h = _draw_interface(ax, material)
    if h is not None:
        ax.legend(handles=[h], fontsize=FONT_SIZES["small"], framealpha=0.7)
    ax.set_xlabel("x  (m)")
    ax.set_ylabel("Depth  (m)")
    ax.invert_yaxis()
    ax.set_title(title)
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_impedance_map(material, output_path=None):
    """Plot the 2D impedance field Z = rho * Vs."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(9, 5))
    xx, zz = np.meshgrid(material.x_m, material.z_m)
    im = ax.pcolormesh(xx, zz, material.impedance_kg_m2_s,
                       cmap=_CMAP_MATERIAL, shading="auto")
    cbar = fig.colorbar(im, ax=ax, label="Z = ρVs  (kg/m²/s)", shrink=0.9)
    cbar.ax.tick_params(labelsize=FONT_SIZES["colorbar"])
    h = _draw_interface(ax, material)
    if h is not None:
        ax.legend(handles=[h], fontsize=FONT_SIZES["small"], framealpha=0.7)
    ax.set_xlabel("x  (m)")
    ax.set_ylabel("Depth  (m)")
    ax.invert_yaxis()
    ax.set_title("Impedance map  Z = ρVs")
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_source_wavelet_2d(result, output_path=None):
    """Plot the Ricker source-time function used in the 2D run."""
    apply_project_style()
    fig, ax = plt.subplots(figsize=(8, 3))
    ax.plot(result.time_s, result.source_time_function,
            color=PROJECT_COLORS["deep_blue"], lw=1.5)
    ax.axvline(result.config.source_peak_time_s, color=PROJECT_COLORS["orange"],
               ls="--", lw=1.2, label=f"Peak  t0 = {result.config.source_peak_time_s:.3f} s")
    ax.set_xlabel("Time  (s)")
    ax.set_ylabel("Normalized amplitude")
    ax.set_title(f"Ricker source  f0 = {result.config.dominant_frequency_hz:.1f} Hz")
    ax.legend()
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_velocity_snapshot(result, snap_index, output_path=None, overlay_circle=False):
    """Plot one 2D velocity snapshot with labelled source and receivers."""
    apply_project_style()
    field = result.velocity_snapshots_m_s[snap_index]
    t     = result.stored_time_s[snap_index]
    vmax  = float(np.max(np.abs(field))) or 1.0

    fig, ax = plt.subplots(figsize=(10, 5))
    xx, zz  = np.meshgrid(result.x_m, result.z_m)
    norm    = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)
    im      = ax.pcolormesh(xx, zz, field, cmap=_CMAP_WAVEFIELD, norm=norm, shading="auto")
    cbar    = fig.colorbar(im, ax=ax, label="vy  (m/s)", shrink=0.9)
    cbar.ax.tick_params(labelsize=FONT_SIZES["colorbar"])

    ax.plot(result.source_x_m, result.source_z_m, "*",
            color=PROJECT_COLORS["white"], ms=12,
            markeredgecolor=PROJECT_COLORS["near_black"], markeredgewidth=0.8,
            zorder=8, label="Source")

    rec_handles = _plot_receivers(ax, result.receiver_positions, label_receivers=True)

    h = _draw_interface(ax, result.material, lw=1.4)
    interface_handles = [h] if h is not None else []

    circle_handles = []
    if overlay_circle:
        t0 = result.config.source_peak_time_s
        vs = result.config.rock_shear_velocity_m_s
        r  = vs * max(t - t0, 0.0)
        theta = np.linspace(0, 2 * np.pi, 360)
        cx = result.source_x_m + r * np.cos(theta)
        cz = result.source_z_m + r * np.sin(theta)
        ax.plot(cx, cz, color=PROJECT_COLORS["green"], lw=1.3, ls="--", zorder=7)
        circle_handles = [
            Line2D([0], [0], color=PROJECT_COLORS["green"], lw=1.3, ls="--",
                   label=f"Theory  r = {r:.0f} m")
        ]

    src_handle = Line2D([0], [0], marker="*", color="w",
                        markerfacecolor=PROJECT_COLORS["white"],
                        markeredgecolor=PROJECT_COLORS["near_black"],
                        markeredgewidth=0.8, markersize=11, label="Source")
    ax.legend(handles=[src_handle] + rec_handles + interface_handles + circle_handles,
              loc="lower right", framealpha=0.8)

    ax.set_xlabel("x  (m)")
    ax.set_ylabel("Depth  (m)")
    ax.invert_yaxis()
    ax.set_title(f"Velocity snapshot  t = {t:.3f} s  [{result.material.label}]")
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_receiver_records(result, output_path=None, max_receivers=8):
    """Plot time histories for all (or up to max_receivers) receivers."""
    apply_project_style()
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
        ax.set_ylabel("vy  (m/s)")
        ax.set_title(f"{name}   x = {rx:.0f} m,  z = {rz:.0f} m")
    axes[-1].set_xlabel("Time  (s)")
    fig.suptitle(f"Receiver records [{result.material.label}]", y=1.01)
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_surface_array(result, output_path=None, surface_z_m=0.0):
    """Plot all surface receivers as a wiggle-trace gather."""
    apply_project_style()
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
    ax.set_title(f"Surface receiver gather [{result.material.label}]")
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def plot_arrival_check(result, receiver_name, theoretical_t, numerical_t, output_path=None):
    """Plot receiver seismogram with theoretical and numerical arrival marks."""
    apply_project_style()
    signal = result.receiver_records_m_s[receiver_name]
    rx, rz = result.receiver_positions[receiver_name]
    fig, ax = plt.subplots(figsize=(9, 3.5))
    ax.plot(result.time_s, signal, color=PROJECT_COLORS["deep_blue"], lw=1.2, label="Numerical")
    ax.axvline(theoretical_t, color=PROJECT_COLORS["green"], lw=1.5, ls="--",
               label=f"Theory  t = {theoretical_t:.4f} s")
    ax.axvline(numerical_t, color=PROJECT_COLORS["crimson"], lw=1.5, ls=":",
               label=f"Numerical peak  t = {numerical_t:.4f} s")
    ax.set_xlabel("Time  (s)")
    ax.set_ylabel("vy  (m/s)")
    ax.set_title(f"Arrival check at x = {rx:.0f} m,  z = {rz:.0f} m")
    ax.legend()
    fig.tight_layout()
    _save(fig, output_path)
    return fig


def make_2d_wave_animation(result, output_path=None, fps=15, dpi=120, skip_blank=True):
    """Animate the 2D velocity snapshots; skip_blank trims leading zero frames."""
    apply_project_style()
    snapshots = result.velocity_snapshots_m_s
    vmax = float(np.max(np.abs(snapshots))) or 1.0
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0, vmax=vmax)

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
    cbar.ax.tick_params(labelsize=FONT_SIZES["colorbar"])

    ax.plot(result.source_x_m, result.source_z_m, "*",
            color=PROJECT_COLORS["white"], ms=12,
            markeredgecolor=PROJECT_COLORS["near_black"], markeredgewidth=0.8,
            zorder=8)

    rec_handles = _plot_receivers(ax, result.receiver_positions, label_receivers=True)

    h = _draw_interface(ax, result.material, lw=1.4)
    interface_handles = [h] if h is not None else []

    src_handle = Line2D([0], [0], marker="*", color="w",
                        markerfacecolor=PROJECT_COLORS["white"],
                        markeredgecolor=PROJECT_COLORS["near_black"],
                        markeredgewidth=0.8, markersize=11, label="Source")
    ax.legend(handles=[src_handle] + rec_handles + interface_handles,
              loc="lower right", framealpha=0.8)

    ax.set_xlabel("x  (m)")
    ax.set_ylabel("Depth  (m)")
    ax.invert_yaxis()
    title = ax.set_title("")
    fig.tight_layout()

    def _update(fi):
        snap = snapshots[fi]
        im.set_array(snap.ravel())
        title.set_text(
            f"{result.material.label} [t = {result.stored_time_s[fi]:.3f} s]"
        )
        return [im, title]

    ani = animation.FuncAnimation(
        fig, _update, frames=frames, interval=1000 // fps, blit=True,
    )

    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        ani.save(str(output_path), writer="pillow", fps=fps, dpi=dpi)

    return ani


def make_2d_wave_animation_comparison(
    results,
    output_path=None,
    fps=15,
    dpi=120,
    skip_blank=True,
):
    """Animate three 2D results side-by-side on the same row.

    Parameters
    ----------
    results : sequence of Simulation2DResult
        Exactly three results (e.g. homogeneous, layered, basin).
    """
    if len(results) != 3:
        raise ValueError("Exactly three results are required.")

    apply_project_style()

    # ---- shared colour scale across all three cases ----
    global_vmax = max(
        float(np.max(np.abs(r.velocity_snapshots_m_s))) or 1.0 for r in results
    )
    norm = TwoSlopeNorm(vmin=-global_vmax, vcenter=0, vmax=global_vmax)

    # ---- find first non-blank frame for each result ----
    start_frames = []
    for r in results:
        s = 0
        if skip_blank:
            threshold = 0.01 * float(np.max(np.abs(r.velocity_snapshots_m_s))) or 1e-9
            for i, snap in enumerate(r.velocity_snapshots_m_s):
                if np.max(np.abs(snap)) > threshold:
                    s = max(0, i - 2)
                    break
        start_frames.append(s)

    n_frames = min(
        len(r.velocity_snapshots_m_s) - s
        for r, s in zip(results, start_frames)
    )
    frame_indices = list(range(n_frames))

    _ANIM_FONT = FONT_SIZES["axes_title"] + 4   # larger than static plots

    # ---- build figure ----
    fig, axes = plt.subplots(1, 3, figsize=(24, 6), sharey=False)

    meshes, titles = [], []
    for ax, r, s in zip(axes, results, start_frames):
        xx, zz = np.meshgrid(r.x_m, r.z_m)
        im = ax.pcolormesh(
            xx, zz, r.velocity_snapshots_m_s[s],
            cmap=_CMAP_WAVEFIELD, norm=norm, shading="auto",
        )
        cbar = fig.colorbar(im, ax=ax, label="vy  (m/s)", shrink=0.85, pad=0.02)
        cbar.ax.tick_params(labelsize=_ANIM_FONT - 2)
        cbar.set_label("vy  (m/s)", fontsize=_ANIM_FONT - 2)

        ax.plot(
            r.source_x_m, r.source_z_m, "*",
            color=PROJECT_COLORS["white"], ms=10,
            markeredgecolor=PROJECT_COLORS["near_black"], markeredgewidth=0.7,
            zorder=8,
        )
        _plot_receivers(ax, r.receiver_positions, label_receivers=False)

        _draw_interface(ax, r.material, lw=1.2)

        ax.set_xlabel("x  (m)", fontsize=_ANIM_FONT)
        ax.set_ylabel("Depth  (m)", fontsize=_ANIM_FONT)
        ax.tick_params(labelsize=_ANIM_FONT - 1)
        ax.invert_yaxis()
        t = ax.set_title("", fontsize=_ANIM_FONT, pad=10)
        meshes.append(im)
        titles.append((t, r, s))

    fig.tight_layout()
    fig.subplots_adjust(top=0.88)

    def _update(fi):
        artists = []
        for im, (title, r, s) in zip(meshes, titles):
            snap = r.velocity_snapshots_m_s[s + fi]
            im.set_array(snap.ravel())
            title.set_text(
                f"{r.material.label}  [t = {r.stored_time_s[s + fi]:.3f} s]"
            )
            artists += [im, title]
        return artists

    ani = animation.FuncAnimation(
        fig, _update, frames=frame_indices, interval=1000 // fps, blit=True,
    )

    if output_path is not None:
        Path(output_path).parent.mkdir(parents=True, exist_ok=True)
        ani.save(str(output_path), writer="pillow", fps=fps, dpi=dpi)
        plt.close(fig)

    return ani
