from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from site_response.materials_2d import (
    Simulation2DConfig,
    layered_2d_model,
    basin_2d_model,
    ensure_2d_output_directories,
    format_2d_parameter_summary,
)
from site_response.fd_solver_2d import run_2d_solver
from site_response.visualization_2d import plot_material_overview
from site_response.utils import save_2d_result_npz


def main():
    cfg = Simulation2DConfig()
    ensure_2d_output_directories(cfg)
    print(format_2d_parameter_summary(cfg))

    src_x = cfg.model_width_m / 2
    src_z = 200.0

    mat_layered = layered_2d_model(cfg)
    mat_basin   = basin_2d_model(cfg, base_depth_m=20.0, extra_depth_m=60.0, basin_half_width_m=180.0)

    Z_soil = cfg.soil_density_kg_m3 * cfg.soil_shear_velocity_m_s
    Z_rock = cfg.rock_density_kg_m3 * cfg.rock_shear_velocity_m_s
    R = (Z_rock - Z_soil) / (Z_rock + Z_soil)
    print(f"Z_soil = {Z_soil:.0f} kg/(m²·s),  Z_rock = {Z_rock:.0f} kg/(m²·s)")
    print(f"Impedance ratio Z_rock/Z_soil = {Z_rock / Z_soil:.2f}")
    print(f"Reflection coefficient R = {R:.3f}")

    print("\nRunning layered 2D simulation...")
    result_layered = run_2d_solver(cfg, mat_layered, source_x_m=src_x, source_z_m=src_z,
                                   receivers={}, store_every=10, source_amplitude=1.0e-4)
    save_2d_result_npz(result_layered, cfg.results_2d_dir / "layered_2d_results.npz")

    print("Running basin 2D simulation...")
    result_basin = run_2d_solver(cfg, mat_basin, source_x_m=src_x, source_z_m=src_z,
                                 receivers={}, store_every=10, source_amplitude=1.0e-4)
    save_2d_result_npz(result_basin, cfg.results_2d_dir / "basin_2d_results.npz")

    from site_response.visualization_2d import _save
    fig = plot_material_overview([mat_layered, mat_basin], cfg,
                                 titles=["Layered soil over rock", "Gaussian basin"])
    _save(fig, cfg.figures_2d_dir / "vs_maps_layered_basin.png")

    print("\nOutputs saved to figures/ and results/")


if __name__ == "__main__":
    main()
