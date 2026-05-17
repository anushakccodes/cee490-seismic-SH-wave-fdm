"""Reproduce all project figures and animations from the command line.

Run from the project root:
    python scripts/make_all_figures.py

Execution order:
    1. Homogeneous verification (Stage 3)
    2. Grid-refinement convergence study (Stage 4)
    3. Layered site-response + amplification (Stages 5–6)
    4. Animations (Stage 7)
"""
from __future__ import annotations

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
SRC_DIR = PROJECT_ROOT / "src"
if str(SRC_DIR) not in sys.path:
    sys.path.insert(0, str(SRC_DIR))

import numpy as np

from site_response.analysis import (
    compute_amplification,
    run_convergence_study,
    run_homogeneous_verification_case,
    run_layered_site_case,
)
from site_response.materials import homogeneous_rock_model, layered_soil_over_rock_model, material_summary
from site_response.parameters import default_config, format_parameter_summary
from site_response.source import build_source_time_function, normalized_amplitude_spectrum
from site_response.utils import prepare_project_outputs, save_simulation_result_npz
from site_response.visualization import (
    apply_project_style,
    make_layered_reflection_animation,
    make_wave_animation,
    plot_amplification,
    plot_convergence,
    plot_homogeneous_seismogram,
    plot_impedance_profile,
    plot_source_spectrum,
    plot_source_time_history,
    plot_spectra_comparison,
    plot_surface_comparison,
    plot_velocity_profile,
    plot_wavefield_image,
)


def main() -> None:
    config = default_config()
    prepare_project_outputs(config)
    apply_project_style()
    print(format_parameter_summary(config))

    # ------------------------------------------------------------------
    # Stage 2: Material profiles and source
    # ------------------------------------------------------------------
    layered_material = layered_soil_over_rock_model(config)
    print(material_summary(layered_material))

    plot_velocity_profile(layered_material,  output_path=config.figures_dir / "velocity_profile.png")
    plot_impedance_profile(layered_material, output_path=config.figures_dir / "impedance_profile.png")

    time_s, source_values = build_source_time_function(config)
    plot_source_time_history(time_s, source_values, output_path=config.figures_dir / "source_wavelet.png")
    frequencies_hz, source_spectrum = normalized_amplitude_spectrum(time_s, source_values)
    plot_source_spectrum(frequencies_hz, source_spectrum, config.dominant_frequency_hz,
                          output_path=config.figures_dir / "source_spectrum.png")
    print("Stage 2 figures saved.")

    # ------------------------------------------------------------------
    # Stage 3: Homogeneous verification
    # ------------------------------------------------------------------
    print("\nRunning homogeneous verification...")
    homogeneous_result, arrival_check = run_homogeneous_verification_case(config)
    print(arrival_check.as_dataframe().to_string(index=False))
    save_simulation_result_npz(homogeneous_result, config.results_dir / "homogeneous_results.npz")
    arrival_check.as_dataframe().to_csv(config.results_dir / "homogeneous_arrival_check.csv", index=False)
    plot_homogeneous_seismogram(homogeneous_result, arrival_check,
                                 output_path=config.figures_dir / "homogeneous_verification.png")
    plot_wavefield_image(homogeneous_result,
                          output_path=config.figures_dir / "homogeneous_wavefield.png")
    print("Stage 3 figures saved.")

    # ------------------------------------------------------------------
    # Stage 4: Convergence study
    # ------------------------------------------------------------------
    print("\nRunning convergence study (4 grid sizes)...")
    _, convergence_df = run_convergence_study(config, dz_values_m=[8.0, 4.0, 2.0, 1.0])
    convergence_df.to_csv(config.results_dir / "convergence_results.csv", index=False)
    print(convergence_df.to_string(index=False))

    dz_arr = convergence_df["grid_spacing_m"].values
    err_arr = convergence_df["percent_error"].values
    if len(dz_arr) >= 2 and err_arr.min() > 0:
        log_slope = np.polyfit(np.log(dz_arr), np.log(err_arr), 1)[0]
        print(f"Observed convergence rate: {log_slope:.2f}  (expected ~2.0 for 2nd-order scheme)")
    plot_convergence(dz_arr, err_arr, output_path=config.figures_dir / "convergence_plot.png")
    print("Stage 4 figures saved.")

    # ------------------------------------------------------------------
    # Stages 5–6: Layered site response and amplification
    # ------------------------------------------------------------------
    print("\nRunning layered site-response simulations...")
    ref_result, layered_result = run_layered_site_case(config)
    save_simulation_result_npz(ref_result,     config.results_dir / "ref_results.npz")
    save_simulation_result_npz(layered_result, config.results_dir / "layered_results.npz")

    plot_wavefield_image(layered_result,
                          output_path=config.figures_dir / "layered_wavefield.png",
                          title="Layered-Site Velocity Wavefield")
    plot_surface_comparison(ref_result, layered_result,
                             output_path=config.figures_dir / "surface_comparison.png")

    amp_result = compute_amplification(ref_result, layered_result)
    amp_result.resonance_table().to_csv(config.results_dir / "resonance_frequencies.csv", index=False)
    plot_spectra_comparison(amp_result.frequencies_hz, amp_result.ref_spectrum,
                             amp_result.layered_spectrum,
                             output_path=config.figures_dir / "spectra_comparison.png")
    plot_amplification(amp_result.frequencies_hz, amp_result.amplification,
                        amp_result.resonance_frequencies_hz,
                        output_path=config.figures_dir / "site_amplification.png",
                        f_min_hz=amp_result.f_min_hz,
                        f_max_hz=amp_result.f_max_hz)
    print("Stages 5–6 figures saved.")

    # ------------------------------------------------------------------
    # Stage 7: Animations
    # ------------------------------------------------------------------
    print("\nGenerating animations (this may take a moment)...")
    make_wave_animation(layered_result,
                         output_path=config.figures_dir / "wave_animation.gif",
                         fps=20, max_frames=200)
    make_layered_reflection_animation(layered_result,
                                       output_path=config.figures_dir / "animation_layered_reflection.gif",
                                       fps=20, max_frames=200)
    print("Stage 7 animations saved.")
    print("\nAll figures written to:", config.figures_dir)


if __name__ == "__main__":
    main()
