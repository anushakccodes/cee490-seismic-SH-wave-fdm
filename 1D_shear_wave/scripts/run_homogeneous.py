from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from site_response.analysis import run_homogeneous_verification_case
from site_response.materials import homogeneous_rock_model, material_summary
from site_response.parameters import default_config, format_parameter_summary
from site_response.utils import prepare_project_outputs, save_simulation_result_npz
from site_response.visualization import plot_homogeneous_seismogram, plot_velocity_profile, plot_wavefield_image


def main():
    config = default_config()
    prepare_project_outputs(config)
    print(format_parameter_summary(config))
    material = homogeneous_rock_model(config)
    print(material_summary(material))
    result, arrival_check = run_homogeneous_verification_case(config)
    print(arrival_check.as_dataframe().to_string(index=False))
    save_simulation_result_npz(result, config.results_dir / "homogeneous_results.npz")
    arrival_check.as_dataframe().to_csv(config.results_dir / "homogeneous_arrival_check.csv", index=False)
    plot_velocity_profile(material, output_path=config.figures_dir / "homogeneous_velocity_profile.png")
    plot_homogeneous_seismogram(result, arrival_check, output_path=config.figures_dir / "homogeneous_verification.png")
    plot_wavefield_image(result, output_path=config.figures_dir / "homogeneous_wavefield.png")


if __name__ == "__main__":
    main()
