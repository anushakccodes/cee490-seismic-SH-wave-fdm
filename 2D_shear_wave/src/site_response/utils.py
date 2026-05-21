"""Small helper functions for scripts and notebooks."""
from __future__ import annotations

from pathlib import Path
import numpy as np

from .fd_solver_2d import Simulation2DResult


def save_2d_result_npz(result: Simulation2DResult, output_path) -> None:
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output_path,
        time_s=result.time_s,
        x_m=result.x_m,
        z_m=result.z_m,
        stored_time_s=result.stored_time_s,
        velocity_snapshots_m_s=result.velocity_snapshots_m_s,
        **{f'rec_{k}': v for k, v in result.receiver_records_m_s.items()},
    )
