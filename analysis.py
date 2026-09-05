"""High-level analysis of a sample spectrum and matched background."""

import numpy as np

from calibration import (
    calibrate_energies,
    calculate_activity_estimate,
    calculate_sample_geometry_correction,
    calculate_sample_geometry_uncertainty,
    subtract_background,
    validate_spectrum,
)
from peak_fitting import (
    discover_unlisted_candidates,
    fit_raw_peak,
    fit_resolution_model,
    resolution_fwhm_kev,
)
from peak_library import GAMMA_LINES, GammaLine
from tka_reader import FIRST_CALIBRATION_CHANNEL


KNOWN_INTERFERENCES_KEV = {
    238.632: 241.995,
    241.995: 238.632,
    295.224: 300.087,
    300.087: 295.224,
    661.657: 665.453,
    665.453: 661.657,
    835.710: 839.070,
    839.070: 835.710,
    964.766: 968.971,
    968.971: 964.766,
    # Ac-228 at 1588.2 keV overlaps the Tl-208 double-escape region.
    1588.200: 1592.513,
}
ANNIHILATION_ENERGY_KEV = 510.999


def _combine_peak_fits(
    sample_fit,
    background_fit,
    sample_live_time_seconds,
    background_live_time_seconds,
    sample_mass_kg,
    geometry_correction,
    geometry_uncertainty,
    detector_info,
    energy_kev,
    emission_probability,
    detection_sigma,
):
    """Combine the separately fitted sample and background peak areas."""
    usable = (
        sample_fit["success"]
        and background_fit["success"]
        and sample_fit["valid"]
        and background_fit["valid"]
    )
    if not usable:
        return {
            "quantifiable": False,
            "detected": False,
            "reason": "Sample or background peak fit was not valid.",
        }

    live_time_ratio = sample_live_time_seconds / background_live_time_seconds
    raw_net_area = (
        sample_fit["area_counts"]
        - live_time_ratio * background_fit["area_counts"]
    )
    raw_net_uncertainty = np.sqrt(
        sample_fit["area_uncertainty_counts"] ** 2
        + (
            live_time_ratio
            * background_fit["area_uncertainty_counts"]
        )
        ** 2
    )
    raw_net_significance = raw_net_area / raw_net_uncertainty

    corrected_area = (
        geometry_correction * sample_fit["area_counts"]
        - live_time_ratio * background_fit["area_counts"]
    )
    statistical_uncertainty = np.sqrt(
        (
            geometry_correction
            * sample_fit["area_uncertainty_counts"]
        )
        ** 2
        + (
            live_time_ratio
            * background_fit["area_uncertainty_counts"]
        )
        ** 2
    )
    significance = corrected_area / statistical_uncertainty
    sample_peak_significance = (
        sample_fit["area_counts"]
        / sample_fit["area_uncertainty_counts"]
    )

    live_time_minutes = sample_live_time_seconds / 60.0
    rate_factor = 1.0 / (live_time_minutes * sample_mass_kg)
    corrected_rate = corrected_area * rate_factor
    corrected_rate_uncertainty = statistical_uncertainty * rate_factor

    if emission_probability is None:
        activity = {
            "activity_bq_per_kg": None,
            "known_standard_uncertainty_bq_per_kg": None,
            "statistical_uncertainty_bq_per_kg": None,
            "detector_efficiency_uncertainty_bq_per_kg": None,
            "geometry_uncertainty_bq_per_kg": None,
            "detector_photopeak_efficiency": None,
            "uncertainty_complete": False,
        }
    else:
        activity = calculate_activity_estimate(
            corrected_area,
            statistical_uncertainty,
            sample_fit["area_counts"],
            sample_live_time_seconds,
            sample_mass_kg,
            energy_kev,
            emission_probability,
            geometry_uncertainty,
            detector_info,
        )

    return {
        "quantifiable": True,
        "detected": bool(
            significance >= detection_sigma
            and corrected_area > 0
        ),
        "sample_peak_detected": bool(
            sample_peak_significance >= detection_sigma
        ),
        "reason": "",
        "live_time_ratio": float(live_time_ratio),
        "sample_peak_significance_sigma": float(sample_peak_significance),
        "raw_net_area_counts": float(raw_net_area),
        "raw_net_area_uncertainty_counts": float(raw_net_uncertainty),
        "raw_net_significance_sigma": float(raw_net_significance),
        "net_area_counts": float(corrected_area),
        "net_area_uncertainty_counts": float(statistical_uncertainty),
        "significance_sigma": float(significance),
        "statistical_decision_threshold_counts": float(
            detection_sigma * statistical_uncertainty
        ),
        "geometry_corrected_rate_cpm_per_kg": float(corrected_rate),
        "geometry_corrected_rate_uncertainty_cpm_per_kg": float(
            corrected_rate_uncertainty
        ),
        **activity,
    }


def _near_energy(
    energy_kev,
    reference_energies_kev,
    resolution_model,
    width_multiplier=1.25,
):
    """Return whether an energy lies within a resolution-aware exclusion."""
    for reference in reference_energies_kev:
        tolerance = max(
            1.0,
            width_multiplier
            * resolution_fwhm_kev(reference, resolution_model),
        )
        if abs(energy_kev - reference) <= tolerance:
            return True
    return False


def _expected_response_feature_energies(line_results):
    """Return detector-response energies to exclude from candidate labels."""
    response_energies = [ANNIHILATION_ENERGY_KEV]
    for result in line_results:
        fit = result["sample_fit"]
        if not fit.get("valid"):
            continue
        area_uncertainty = fit.get("area_uncertainty_counts", np.nan)
        if not np.isfinite(area_uncertainty) or area_uncertainty <= 0:
            continue
        if fit["area_counts"] / area_uncertainty < 5.0:
            continue
        gamma_energy = fit["energy_kev"]
        if gamma_energy > 2.0 * ANNIHILATION_ENERGY_KEV:
            response_energies.append(
                gamma_energy - ANNIHILATION_ENERGY_KEV
            )
            response_energies.append(
                gamma_energy - 2.0 * ANNIHILATION_ENERGY_KEV
            )
    return response_energies


def analyse_spectrum(
    sample_counts,
    background_counts,
    sample_live_time_seconds,
    background_live_time_seconds,
    sample_mass_kg,
    detector,
    detector_info,
    filled_height_cm,
    beaker_filled=True,
    minimum_energy_kev=200.0,
    detection_sigma=3.0,
    gamma_lines=GAMMA_LINES,
    discover_candidates=True,
):
    """Analyse reference gamma lines and optionally search for extra peaks."""
    validate_spectrum(sample_counts, background_counts)
    if sample_live_time_seconds <= 0 or background_live_time_seconds <= 0:
        raise ValueError("Live times must be positive.")
    if sample_mass_kg <= 0:
        raise ValueError("Sample mass must be positive.")
    if detection_sigma <= 0:
        raise ValueError("Detection threshold must be positive.")
    if minimum_energy_kev < 0:
        raise ValueError("Minimum energy cannot be negative.")

    channels = np.arange(
        FIRST_CALIBRATION_CHANNEL,
        FIRST_CALIBRATION_CHANNEL + len(sample_counts),
        dtype=float,
    )
    energies = calibrate_energies(
        channels,
        detector_info["slope"],
        detector_info["offset"],
    )
    net_counts = subtract_background(
        sample_counts,
        background_counts,
        sample_live_time_seconds,
        background_live_time_seconds,
    )
    geometry_correction = calculate_sample_geometry_correction(
        detector,
        filled_height_cm,
        beaker_filled,
    )
    geometry_uncertainty = calculate_sample_geometry_uncertainty(
        detector,
        filled_height_cm,
        beaker_filled,
    )

    line_results = []
    maximum_energy_kev = float(energies[-1])
    for line in gamma_lines:
        if not (minimum_energy_kev <= line.energy_kev <= maximum_energy_kev):
            continue

        sample_fit = fit_raw_peak(
            sample_counts,
            channels,
            detector_info,
            line.energy_kev,
            companion_energy_kev=KNOWN_INTERFERENCES_KEV.get(line.energy_kev),
        )
        background_fit = fit_raw_peak(
            background_counts,
            channels,
            detector_info,
            line.energy_kev,
            companion_energy_kev=KNOWN_INTERFERENCES_KEV.get(line.energy_kev),
        )
        net = _combine_peak_fits(
            sample_fit,
            background_fit,
            sample_live_time_seconds,
            background_live_time_seconds,
            sample_mass_kg,
            geometry_correction,
            geometry_uncertainty,
            detector_info,
            line.energy_kev,
            line.emission_probability,
            detection_sigma,
        )
        line_results.append(
            {
                "line": line,
                "sample_fit": sample_fit,
                "background_fit": background_fit,
                "net": net,
                "candidate_scan_significance_sigma": None,
            }
        )

    resolution_model = fit_resolution_model(line_results)
    discovery = {
        "enabled": bool(discover_candidates),
        "resolution_model": resolution_model,
        "effective_trials": 0,
        "scan_threshold_sigma": None,
        "maximum_scan_score_sigma": None,
        "screened_candidate_count": 0,
        "fitted_candidate_count": 0,
        "response_feature_rejection_count": 0,
        "known_line_rejection_count": 0,
        "duplicate_rejection_count": 0,
        "detected_candidate_count": 0,
    }

    if discover_candidates:
        candidates, scan_diagnostics = discover_unlisted_candidates(
            sample_counts,
            energies,
            resolution_model,
            minimum_energy_kev,
            maximum_energy_kev,
            detection_sigma,
        )
        discovery.update(scan_diagnostics)
        discovery["screened_candidate_count"] = len(candidates)

        known_energies = [
            line.energy_kev
            for line in gamma_lines
            if minimum_energy_kev <= line.energy_kev <= maximum_energy_kev
        ]
        response_energies = _expected_response_feature_energies(line_results)
        accepted_candidates = []
        accepted_energies = []
        for candidate in sorted(
            candidates,
            key=lambda item: item["scan_significance_sigma"],
            reverse=True,
        ):
            energy = candidate["energy_kev"]
            if _near_energy(energy, known_energies, resolution_model):
                discovery["known_line_rejection_count"] += 1
                continue
            if _near_energy(
                energy,
                response_energies,
                resolution_model,
            ):
                discovery["response_feature_rejection_count"] += 1
                continue
            if _near_energy(
                energy,
                accepted_energies,
                resolution_model,
            ):
                discovery["duplicate_rejection_count"] += 1
                continue
            accepted_candidates.append(candidate)
            accepted_energies.append(energy)

        # Only run the slower full fit on the 50 strongest candidates.
        accepted_candidates = accepted_candidates[:50]
        discovery["fitted_candidate_count"] = len(accepted_candidates)

        for candidate in accepted_candidates:
            reference_energy = candidate["energy_kev"]
            sample_fit = fit_raw_peak(
                sample_counts,
                channels,
                detector_info,
                reference_energy,
            )
            background_fit = fit_raw_peak(
                background_counts,
                channels,
                detector_info,
                reference_energy,
            )
            line = GammaLine(
                "Unknown",
                reference_energy,
                None,
                "Unlisted full-energy candidate",
            )
            net = _combine_peak_fits(
                sample_fit,
                background_fit,
                sample_live_time_seconds,
                background_live_time_seconds,
                sample_mass_kg,
                geometry_correction,
                geometry_uncertainty,
                detector_info,
                reference_energy,
                None,
                detection_sigma,
            )
            line_results.append(
                {
                    "line": line,
                    "sample_fit": sample_fit,
                    "background_fit": background_fit,
                    "net": net,
                    "candidate_scan_significance_sigma": candidate[
                        "scan_significance_sigma"
                    ],
                }
            )

    corrected_counts_diagnostic = (
        geometry_correction * np.asarray(sample_counts, dtype=float)
        - (
            sample_live_time_seconds
            / background_live_time_seconds
        )
        * np.asarray(background_counts, dtype=float)
    )
    diagnostic_rate = (
        corrected_counts_diagnostic
        / (sample_live_time_seconds / 60.0)
        / sample_mass_kg
    )
    line_results.sort(key=lambda result: result["line"].energy_kev)
    detected = [result for result in line_results if result["net"]["detected"]]
    discovery["detected_candidate_count"] = sum(
        result["line"].emission_probability is None
        for result in detected
    )

    return {
        "channels": channels,
        "energies_kev": energies,
        "net_counts_diagnostic": net_counts,
        "corrected_counts_diagnostic": corrected_counts_diagnostic,
        "net_rate_diagnostic_cpm_per_kg": diagnostic_rate,
        "geometry_correction": float(geometry_correction),
        "geometry_uncertainty": geometry_uncertainty,
        "background_live_time_scale": float(
            sample_live_time_seconds / background_live_time_seconds
        ),
        "line_results": line_results,
        "detected_line_results": detected,
        "discovery": discovery,
        "detection_sigma": float(detection_sigma),
    }
