"""Energy calibration, background scaling and activity corrections."""

import numpy as np


def validate_spectrum(sample_counts, background_counts):
    """Check that sample and background spectra can be compared directly."""
    sample = np.asarray(sample_counts)
    background = np.asarray(background_counts)

    if sample.ndim != 1 or background.ndim != 1:
        raise ValueError("Spectra must be one-dimensional arrays.")
    if len(sample) == 0 or len(background) == 0:
        raise ValueError("Spectra cannot be empty.")
    if len(sample) != len(background):
        raise ValueError("Sample and background spectra must have equal lengths.")
    if not np.all(np.isfinite(sample)) or not np.all(np.isfinite(background)):
        raise ValueError("Spectra cannot contain NaN or infinite values.")
    if np.any(sample < 0) or np.any(background < 0):
        raise ValueError("Raw spectra cannot contain negative counts.")


def calibrate_energies(channels, slope, offset):
    """Convert MCA channel numbers to energy in keV."""
    channels = np.asarray(channels, dtype=float)
    if slope <= 0:
        raise ValueError("The energy-calibration slope must be positive.")
    return slope * channels + offset


def subtract_background(
    sample_counts,
    background_counts,
    sample_live_time,
    background_live_time,
):
    """Subtract a live-time-scaled background spectrum in raw count space."""
    validate_spectrum(sample_counts, background_counts)
    if sample_live_time <= 0 or background_live_time <= 0:
        raise ValueError("Live times must be positive.")

    scale = sample_live_time / background_live_time
    return (
        np.asarray(sample_counts, dtype=float)
        - scale * np.asarray(background_counts, dtype=float)
    )


def calculate_sample_geometry_correction(
    detector,
    filled_height_cm,
    beaker_filled=True,
):
    """Calculate the sample-geometry correction used in the project."""
    if detector not in (1, 2, 3):
        raise ValueError("Detector must be 1, 2, or 3.")
    if filled_height_cm <= 0:
        raise ValueError("Filled height must be positive.")

    detector_depth_cm = 0.9
    total_height_cm = filled_height_cm + detector_depth_cm

    if detector == 3:
        if not beaker_filled:
            return 0.96885
        upper_correction = (
            (0.73205 * 12.1)
            / np.cos(np.arctan(total_height_cm / 12.1))
            / 5
        ) ** 2
        lower_correction = 0.96885
        upper_volume = 114.9901 * filled_height_cm
        lower_volume = 318.4043
    else:
        if not beaker_filled:
            return 1.38465
        upper_correction = (
            (0.73205 * 15.05)
            / np.cos(np.arctan(total_height_cm / 15.05))
            / 5
        ) ** 2
        lower_correction = 1.38465
        upper_volume = 711.57859 * filled_height_cm
        lower_volume = 620.480257

    return (
        upper_correction * upper_volume
        + lower_correction * lower_volume
    ) / (upper_volume + lower_volume)


def _small_beaker_upper_correction(height_cm):
    return (
        0.0214358881 * height_cm**2
        + 0.0385845986 * height_cm
        + 3.1557914466
    )


def _small_beaker_geometry_uncertainty(height_cm):
    """Uncertainty expression retained from the Detector 3 appendix."""
    term_1 = (
        (114.9901 * height_cm / (114.9901 * height_cm + 318.4043))
        * 0.00214358881
        * np.sqrt(146.41 + 2.0 * (height_cm + 0.9) ** 2)
    )
    term_2 = (
        318.4043
        * (_small_beaker_upper_correction(height_cm) - 0.96885)
        / (114.9901 * height_cm + 318.4043) ** 2
        * np.sqrt(0.9031304877 * height_cm**2 + 33.05683368)
    )
    term_3 = (
        318.4043
        / (114.9901 * height_cm + 318.4043)
        * 0.0139201
    )
    term_4 = (
        114.9901
        * height_cm
        * (0.96885 - _small_beaker_upper_correction(height_cm))
        / (114.9901 * height_cm + 318.4043) ** 2
        * 8.668472272
    )
    return float(np.sqrt(term_1**2 + term_2**2 + term_3**2 + term_4**2))


def calculate_sample_geometry_uncertainty(
    detector,
    filled_height_cm,
    beaker_filled=True,
):
    """Return the available geometry uncertainty for the selected detector."""
    if detector == 3 and not beaker_filled:
        return 0.0139201
    if detector == 3 and beaker_filled:
        return _small_beaker_geometry_uncertainty(filled_height_cm)

    # Only Detector 3 had a recorded geometry uncertainty in the project files.
    return None


def detector_photopeak_efficiency(energy_kev, detector_info):
    """Evaluate the detector photopeak efficiency at one energy."""
    if energy_kev <= 0:
        raise ValueError("Efficiency can only be evaluated at positive energy.")
    if detector_info.get("efficiency_unit") != "fraction":
        raise ValueError(
            "Detector-efficiency unit has not been verified as a fraction."
        )

    efficiency = (
        detector_info["eff_coefficient"]
        * energy_kev ** detector_info["eff_power"]
    )
    if not 0 < efficiency <= 1:
        raise ValueError("Calculated photopeak efficiency is outside (0, 1].")
    return float(efficiency)


def relative_detector_efficiency_uncertainty(energy_kev, detector_info):
    """Propagate the fitted efficiency-curve parameter uncertainties."""
    if energy_kev <= 0:
        raise ValueError("Energy must be positive.")

    coefficient = detector_info["eff_coefficient"]
    coefficient_uncertainty = detector_info["eff_coefficient_uncertainty"]
    power_uncertainty = detector_info["eff_power_uncertainty"]
    return float(
        np.sqrt(
            (coefficient_uncertainty / coefficient) ** 2
            + (np.log(energy_kev) * power_uncertainty) ** 2
        )
    )


def calculate_activity_estimate(
    corrected_area_counts,
    statistical_uncertainty_counts,
    sample_area_counts,
    sample_live_time_seconds,
    sample_mass_kg,
    energy_kev,
    emission_probability,
    geometry_uncertainty,
    detector_info,
):
    """Convert a corrected photopeak area to an activity estimate in Bq/kg."""
    if sample_live_time_seconds <= 0:
        raise ValueError("Sample live time must be positive.")
    if sample_mass_kg <= 0:
        raise ValueError("Sample mass must be positive.")
    if statistical_uncertainty_counts < 0:
        raise ValueError("Statistical uncertainty cannot be negative.")
    if not 0 < emission_probability <= 1:
        raise ValueError("Gamma emission probability must lie in (0, 1].")

    efficiency = detector_photopeak_efficiency(energy_kev, detector_info)
    scale = 1.0 / (
        sample_live_time_seconds
        * sample_mass_kg
        * efficiency
        * emission_probability
    )

    activity = corrected_area_counts * scale
    statistical_uncertainty = statistical_uncertainty_counts * scale
    efficiency_uncertainty = (
        abs(activity)
        * relative_detector_efficiency_uncertainty(energy_kev, detector_info)
    )

    # g multiplies the sample area, so this term uses S rather than the net area.
    if geometry_uncertainty is None:
        geometry_activity_uncertainty = None
        known_terms = (statistical_uncertainty, efficiency_uncertainty)
    else:
        geometry_activity_uncertainty = abs(
            sample_area_counts * geometry_uncertainty * scale
        )
        known_terms = (
            statistical_uncertainty,
            efficiency_uncertainty,
            geometry_activity_uncertainty,
        )

    return {
        "activity_bq_per_kg": float(activity),
        "known_standard_uncertainty_bq_per_kg": float(
            np.sqrt(np.sum(np.square(known_terms)))
        ),
        "statistical_uncertainty_bq_per_kg": float(statistical_uncertainty),
        "detector_efficiency_uncertainty_bq_per_kg": float(
            efficiency_uncertainty
        ),
        "geometry_uncertainty_bq_per_kg": (
            float(geometry_activity_uncertainty)
            if geometry_activity_uncertainty is not None
            else None
        ),
        "detector_photopeak_efficiency": efficiency,
        "uncertainty_complete": False,
    }
