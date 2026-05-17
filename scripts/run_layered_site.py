"""Run the layered site-response case (Stages 5–6) from the command line."""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

from site_response.analysis import compute_amplification, run_layered_site_case
from site_response.parameters import default_config, format_parameter_summary
from site_response.utils import prepare_project_outputs, save_simulation_result_npz
from site_response.visualization import (
    apply_project_style,
    plot_amplification,
    plot_spectra_comparison,
    plot_surface_comparison,
    plot_wavefield_image,
)


def main() -> None:
    config = default_config()
    prepare_project_outputs(config)
    apply_project_style()
    print(format_parameter_summary(config))

    print("Running rock reference and layered site simulations...")
    ref_result, layered_result = run_layered_site_case(config)

    # Impedance contrast summary.
    Z_soil = config.soil_density_kg_m3 * config.soil_shear_velocity_m_s
    Z_rock = config.rock_density_kg_m3 * config.rock_shear_velocity_m_s
    R = (Z_rock - Z_soil) / (Z_rock + Z_soil)
    print(f"Z_soil = {Z_soil:.0f} kg/(m^2 s),  Z_rock = {Z_rock:.0f} kg/(m^2 s)")
    print(f"Impedance ratio Z_rock/Z_soil = {Z_rock / Z_soil:.2f}")
    print(f"Reflection coefficient R = {R:.3f}")

    # Save simulation arrays.
    save_simulation_result_npz(ref_result,     config.results_dir / "ref_results.npz")
    save_simulation_result_npz(layered_result, config.results_dir / "layered_results.npz")

    # Stage 5 figures.
    plot_wavefield_image(layered_result,
                         output_path=config.figures_dir / "layered_wavefield.png",
                         title="Layered-Site Velocity Wavefield")
    plot_surface_comparison(ref_result, layered_result,
                             output_path=config.figures_dir / "surface_comparison.png")

    # Stage 6: amplification.
    amp_result = compute_amplification(ref_result, layered_result)
    print("\nQuarter-wavelength resonance frequencies:")
    print(amp_result.resonance_table().to_string(index=False))
    amp_result.resonance_table().to_csv(config.results_dir / "resonance_frequencies.csv", index=False)

    plot_spectra_comparison(amp_result.frequencies_hz,
                             amp_result.ref_spectrum,
                             amp_result.layered_spectrum,
                             output_path=config.figures_dir / "spectra_comparison.png")
    plot_amplification(amp_result.frequencies_hz,
                        amp_result.amplification,
                        amp_result.resonance_frequencies_hz,
                        output_path=config.figures_dir / "site_amplification.png",
                        f_min_hz=amp_result.f_min_hz,
                        f_max_hz=amp_result.f_max_hz)

    print("\nLayered site outputs saved to figures/ and results/")


if __name__ == "__main__":
    main()
