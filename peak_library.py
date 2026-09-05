"""Reference gamma lines used by the analyser."""

from dataclasses import dataclass


@dataclass(frozen=True)
class GammaLine:
    isotope: str
    energy_kev: float
    emission_probability: float | None
    series: str


# Curated full-energy peaks used in the original environmental-sample analysis.
# Emission probabilities are stored as fractions per decay.
GAMMA_LINES = (
    GammaLine("Ac-228", 209.253, 0.0397, "Th-232 series"),
    GammaLine("Pb-212", 238.632, 0.436, "Th-232 series"),
    GammaLine("Pb-214", 241.995, 0.0726, "U-238 series"),
    GammaLine("Ac-228", 270.245, 0.0355, "Th-232 series"),
    GammaLine("Tl-208", 277.371, 0.0631, "Th-232 series"),
    GammaLine("Pb-214", 295.224, 0.1828, "U-238 series"),
    GammaLine("Pb-212", 300.087, 0.0328, "Th-232 series"),
    GammaLine("Ac-228", 328.000, 0.0295, "Th-232 series"),
    GammaLine("Ac-228", 338.320, 0.1127, "Th-232 series"),
    GammaLine("Pb-214", 351.932, 0.3534, "U-238 series"),
    GammaLine("Ac-228", 409.462, 0.0192, "Th-232 series"),
    GammaLine("Ac-228", 463.004, 0.0440, "Th-232 series"),
    GammaLine("Bi-214", 562.300, 0.0089, "U-238 series"),
    GammaLine("Tl-208", 583.187, 0.850, "Th-232 series"),
    GammaLine("Bi-214", 609.316, 0.4516, "U-238 series"),
    GammaLine("Cs-137", 661.657, 0.8499, "Anthropogenic"),
    GammaLine("Bi-214", 665.453, 0.01521, "U-238 series"),
    GammaLine("Bi-214", 703.110, 0.00472, "U-238 series"),
    GammaLine("Bi-212", 727.330, 0.0667, "Th-232 series"),
    GammaLine("Ac-228", 755.315, 0.0104, "Th-232 series"),
    GammaLine("Bi-214", 768.356, 0.0494, "U-238 series"),
    GammaLine("Ac-228", 772.290, 0.0145, "Th-232 series"),
    GammaLine("Pb-214", 785.960, 0.0106, "U-238 series"),
    GammaLine("Ac-228", 794.947, 0.0425, "Th-232 series"),
    GammaLine("Bi-214", 806.185, 0.01255, "U-238 series"),
    GammaLine("Ac-228", 835.710, 0.0161, "Th-232 series"),
    GammaLine("Pb-214", 839.070, 0.00583, "U-238 series"),
    GammaLine("Tl-208", 860.557, 0.1250, "Th-232 series"),
    GammaLine("Ac-228", 904.190, 0.0077, "Th-232 series"),
    GammaLine("Ac-228", 911.204, 0.258, "Th-232 series"),
    GammaLine("Bi-214", 934.061, 0.0303, "U-238 series"),
    GammaLine("Ac-228", 964.766, 0.0499, "Th-232 series"),
    GammaLine("Ac-228", 968.971, 0.1580, "Th-232 series"),
    GammaLine("Pa-234m", 1001.030, 0.00847, "U-238 series"),
    GammaLine("Bi-214", 1120.287, 0.1478, "U-238 series"),
    GammaLine("Bi-214", 1155.190, 0.01624, "U-238 series"),
    GammaLine("Bi-214", 1238.110, 0.0579, "U-238 series"),
    GammaLine("Bi-214", 1280.960, 0.01425, "U-238 series"),
    GammaLine("Bi-214", 1377.669, 0.0400, "U-238 series"),
    GammaLine("Bi-214", 1385.310, 0.00795, "U-238 series"),
    GammaLine("Bi-214", 1401.516, 0.01324, "U-238 series"),
    GammaLine("Bi-214", 1407.980, 0.0215, "U-238 series"),
    GammaLine("K-40", 1460.822, 0.1066, "Primordial"),
    GammaLine("Bi-214", 1509.228, 0.0211, "U-238 series"),
    GammaLine("Ac-228", 1588.200, 0.0322, "Th-232 series"),
    GammaLine("Bi-212", 1620.500, 0.01486, "Th-232 series"),
    GammaLine("Ac-228", 1630.630, 0.0151, "Th-232 series"),
    GammaLine("Bi-214", 1661.316, 0.01037, "U-238 series"),
    GammaLine("Bi-214", 1729.595, 0.0292, "U-238 series"),
    GammaLine("Bi-214", 1764.539, 0.1517, "U-238 series"),
    GammaLine("Bi-214", 1847.420, 0.0211, "U-238 series"),
    GammaLine("Bi-214", 2118.550, 0.0117, "U-238 series"),
    GammaLine("Bi-214", 2204.210, 0.0508, "U-238 series"),
    GammaLine("Bi-214", 2447.860, 0.0155, "U-238 series"),
    GammaLine("Tl-208", 2614.511, 0.9979, "Th-232 series"),
)


def nearest_gamma_line(energy_kev, tolerance_kev=1.5, lines=GAMMA_LINES):
    """Return the nearest reference line within the requested tolerance."""
    if tolerance_kev <= 0:
        raise ValueError("Identification tolerance must be positive.")
    if not lines:
        return None

    nearest = min(lines, key=lambda line: abs(energy_kev - line.energy_kev))
    if abs(energy_kev - nearest.energy_kev) <= tolerance_kev:
        return nearest
    return None
