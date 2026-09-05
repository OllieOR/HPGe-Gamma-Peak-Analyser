"""Read GammaVision-style .TKA spectrum files."""

from pathlib import Path

import numpy as np


# The project calibration numbers the first stored count as MCA channel 2.
FIRST_CALIBRATION_CHANNEL = 2


def read_tka(file_path):
    """Return channels, counts, live time and real time from a TKA file."""
    file_path = Path(file_path)
    if file_path.suffix.lower() != ".tka":
        raise ValueError("File must have the .TKA extension.")
    if not file_path.exists():
        raise FileNotFoundError(file_path)

    data = np.loadtxt(file_path, dtype=float)
    if data.ndim != 1 or len(data) < 3:
        raise ValueError("TKA file must contain two times followed by counts.")

    live_time = float(data[0])
    real_time = float(data[1])
    raw_counts = data[2:]

    if live_time <= 0 or real_time <= 0:
        raise ValueError("Acquisition times must be positive.")
    if live_time > real_time:
        raise ValueError("Live time cannot be greater than real time.")
    if not np.all(np.isfinite(raw_counts)):
        raise ValueError("Spectrum contains a non-finite count value.")
    if np.any(raw_counts < 0):
        raise ValueError("Spectrum contains negative raw counts.")
    if not np.allclose(raw_counts, np.rint(raw_counts)):
        raise ValueError("Raw TKA counts must be whole numbers.")

    counts = np.rint(raw_counts).astype(np.int64)
    channels = np.arange(
        FIRST_CALIBRATION_CHANNEL,
        FIRST_CALIBRATION_CHANNEL + len(counts),
        dtype=float,
    )
    return channels, counts, live_time, real_time
