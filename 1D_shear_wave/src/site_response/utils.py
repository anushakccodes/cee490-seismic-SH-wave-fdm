"""Small helper functions for scripts and notebooks."""
from __future__ import annotations

from pathlib import Path
import numpy as np

from .fd_solver_1d import SimulationResult
from .parameters import SimulationConfig, ensure_output_directories


def prepare_project_outputs(config):
    ensure_output_directories(config)


def save_simulation_result_npz(result, output_path):
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        time_s=result.time_s,
        depth_m=result.depth_m,
        stored_time_s=result.stored_time_s,
        velocity_wavefield_m_s=result.velocity_wavefield_m_s,
        surface_velocity_m_s=result.surface_velocity_m_s,
        receiver_velocity_m_s=result.receiver_velocity_m_s,
        source_time_function=result.source_time_function,
        source_depth_m=result.source_depth_m,
        receiver_depth_m=result.receiver_depth_m,
    )
