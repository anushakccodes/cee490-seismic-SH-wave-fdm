from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from site_response.materials_2d import (
    Simulation2DConfig,
    homogeneous_2d_model,
    layered_2d_model,
    basin_2d_model,
    ensure_2d_output_directories,
    format_2d_parameter_summary,
    cfl_check_2d,
)
from site_response.fd_solver_2d import run_2d_solver
from site_response.analysis_2d import run_arrival_verification_table, format_arrival_table
from site_response.visualization_2d import (
    _save,
    plot_material_overview,
    plot_arrival_check,
    make_2d_wave_animation_comparison,
)
from site_response.utils import save_2d_result_npz


def main():
    cfg = Simulation2DConfig()
    ensure_2d_output_directories(cfg)
    print(format_2d_parameter_summary(cfg))
    cfl_check_2d(cfg)

    src_x = cfg.model_width_m / 2
    src_z = 200.0

    mat_homo    = homogeneous_2d_model(cfg)
    mat_layered = layered_2d_model(cfg)
    mat_basin   = basin_2d_model(cfg, base_depth_m=20.0, extra_depth_m=60.0,
                                 basin_half_width_m=180.0)

    fig = plot_material_overview(
        [mat_homo, mat_layered, mat_basin], cfg,
        titles=["Homogeneous rock", "Layered soil over rock", "Gaussian basin"],
    )
    _save(fig, cfg.figures_2d_dir / "vs_maps_all.png")

    receivers = {
        'surf_100': (100.0, 0.0),
        'surf_200': (200.0, 0.0),
        'surf_300': (300.0, 0.0),
        'surf_400': (400.0, 0.0),
        'surf_500': (500.0, 0.0),
        'surf_600': (600.0, 0.0),
        'surf_700': (700.0, 0.0),
        'internal': (src_x + 100.0, src_z),
    }

    print("\nRunning homogeneous 2D verification...")
    result_homo = run_2d_solver(cfg, mat_homo, source_x_m=src_x, source_z_m=src_z,
                                receivers=receivers, store_every=10, source_amplitude=1.0e-4)
    save_2d_result_npz(result_homo, cfg.results_2d_dir / "homogeneous_2d_results.npz")

    checks = run_arrival_verification_table(result_homo)
    print(format_arrival_table(checks))
    c_int = next(c for c in checks if c.receiver_name == "internal")
    plot_arrival_check(result_homo, "internal", c_int.theoretical_peak_s, c_int.numerical_peak_s,
                       output_path=cfg.figures_2d_dir / "homogeneous_arrival_check.png")

    print("\nRunning layered 2D simulation...")
    result_layered = run_2d_solver(cfg, mat_layered, source_x_m=src_x, source_z_m=src_z,
                                   receivers={}, store_every=10, source_amplitude=1.0e-4)
    save_2d_result_npz(result_layered, cfg.results_2d_dir / "layered_2d_results.npz")

    print("Running basin 2D simulation...")
    result_basin = run_2d_solver(cfg, mat_basin, source_x_m=src_x, source_z_m=src_z,
                                 receivers={}, store_every=10, source_amplitude=1.0e-4)
    save_2d_result_npz(result_basin, cfg.results_2d_dir / "basin_2d_results.npz")

    print("\nGenerating wavefield animation...")
    make_2d_wave_animation_comparison(
        [result_homo, result_layered, result_basin],
        output_path=cfg.figures_2d_dir / "wavefield_2d_comparison.gif",
        fps=15, dpi=100,
    )

    print("\nAll figures written to:", cfg.figures_2d_dir)


if __name__ == "__main__":
    main()
