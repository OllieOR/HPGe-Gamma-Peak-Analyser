"""Photopeak fitting and optional search for unlisted peaks."""

import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from scipy.stats import norm

from calibration import calibrate_energies


FWHM_TO_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))
MINIMUM_RESOLUTION_FWHM_KEV = 0.5
MAXIMUM_RESOLUTION_FWHM_KEV = 5.0


def gaussian_with_background(
    shifted_channels,
    area,
    centre_shift,
    sigma,
    log_background_level,
    background_log_slope,
):
    """Gaussian peak plus a slowly varying positive background."""
    x = np.asarray(shifted_channels, dtype=float)
    x_scale = max(float(np.max(np.abs(x))), 1.0)

    gaussian = (
        area
        / (sigma * np.sqrt(2.0 * np.pi))
        * np.exp(-0.5 * ((x - centre_shift) / sigma) ** 2)
    )
    background = np.exp(log_background_level + background_log_slope * x / x_scale)
    return gaussian + background


def _doublet_with_background(
    shifted_channels,
    target_area,
    companion_area,
    target_centre_shift,
    companion_centre_shift,
    target_sigma,
    companion_sigma,
    log_background_level,
    background_log_slope,
    companion_separation_channels,
):
    """Two Gaussian peaks on a shared local background."""
    x = np.asarray(shifted_channels, dtype=float)
    x_scale = max(float(np.max(np.abs(x))), 1.0)
    target = (
        target_area
        / (target_sigma * np.sqrt(2.0 * np.pi))
        * np.exp(
            -0.5
            * ((x - target_centre_shift) / target_sigma) ** 2
        )
    )
    companion = (
        companion_area
        / (companion_sigma * np.sqrt(2.0 * np.pi))
        * np.exp(
            -0.5
            * (
                (
                    x
                    - companion_separation_channels
                    - companion_centre_shift
                )
                / companion_sigma
            )
            ** 2
        )
    )
    background = np.exp(log_background_level + background_log_slope * x / x_scale)
    return target + companion + background


def _poisson_deviance_residuals(
    parameters,
    x,
    observed_counts,
    companion_separation_channels=None,
):
    """Return signed square-root Poisson deviance residuals."""
    if companion_separation_channels is None:
        expected = gaussian_with_background(x, *parameters)
    else:
        expected = _doublet_with_background(
            x,
            *parameters,
            companion_separation_channels,
        )
    expected = np.clip(expected, 1e-12, None)
    observed = np.asarray(observed_counts, dtype=float)

    deviance_term = np.empty_like(observed)
    positive = observed > 0
    deviance_term[positive] = 2.0 * (
        expected[positive]
        - observed[positive]
        + observed[positive] * np.log(observed[positive] / expected[positive])
    )
    deviance_term[~positive] = 2.0 * expected[~positive]

    # Round-off can produce a tiny negative value at a perfect match.
    deviance_term = np.maximum(deviance_term, 0.0)
    return np.sign(observed - expected) * np.sqrt(deviance_term)


def _failed_fit(reason, reference_energy_kev):
    return {
        "success": False,
        "valid": False,
        "reason": reason,
        "reference_energy_kev": float(reference_energy_kev),
    }


def fit_raw_peak(
    counts,
    channels,
    detector_info,
    reference_energy_kev,
    roi_half_width_kev=3.0,
    centre_tolerance_kev=2.0,
    companion_energy_kev=None,
):
    """Fit one expected photopeak to a raw count spectrum."""
    counts = np.asarray(counts, dtype=float)
    channels = np.asarray(channels, dtype=float)
    if counts.ndim != 1 or channels.ndim != 1 or len(counts) != len(channels):
        raise ValueError("Counts and channels must be equal one-dimensional arrays.")
    if np.any(counts < 0):
        raise ValueError("Raw counts cannot be negative.")
    if roi_half_width_kev <= 0 or centre_tolerance_kev <= 0:
        raise ValueError("ROI widths and centre tolerance must be positive.")

    slope = detector_info["slope"]
    offset = detector_info["offset"]
    energies = calibrate_energies(channels, slope, offset)
    if companion_energy_kev is None:
        roi = np.abs(energies - reference_energy_kev) <= roi_half_width_kev
        companion_separation_channels = None
    else:
        if companion_energy_kev == reference_energy_kev:
            raise ValueError("Companion and target energies must differ.")
        lower_energy = min(reference_energy_kev, companion_energy_kev)
        upper_energy = max(reference_energy_kev, companion_energy_kev)
        roi = (
            (energies >= lower_energy - roi_half_width_kev)
            & (energies <= upper_energy + roi_half_width_kev)
        )
        companion_separation_channels = (
            companion_energy_kev - reference_energy_kev
        ) / slope

    if np.count_nonzero(roi) < 20:
        return _failed_fit(
            "Fewer than 20 channels fall inside the fit ROI.",
            reference_energy_kev,
        )

    roi_channels = channels[roi]
    roi_counts = counts[roi]
    expected_channel = (reference_energy_kev - offset) / slope
    x = roi_channels - expected_channel

    edge_count = max(4, len(roi_counts) // 6)
    edge_values = np.concatenate((roi_counts[:edge_count], roi_counts[-edge_count:]))
    background_guess = max(float(np.median(edge_values)), 0.05)
    excess_counts = np.maximum(roi_counts - background_guess, 0.0)
    target_neighbourhood = np.abs(x) <= 1.5 / slope
    target_excess = excess_counts[target_neighbourhood]
    target_x = x[target_neighbourhood]
    area_guess = max(float(np.sum(target_excess)), 1.0)

    if np.sum(target_excess) > 0:
        centre_guess = float(
            np.sum(target_x * target_excess) / np.sum(target_excess)
        )
    else:
        centre_guess = 0.0
    centre_limit_channels = centre_tolerance_kev / slope
    centre_guess = float(
        np.clip(centre_guess, -0.8 * centre_limit_channels, 0.8 * centre_limit_channels)
    )

    # The measured HPGe peaks are close to 2 keV wide.
    sigma_guess = 2.0 / FWHM_TO_SIGMA / slope
    if companion_separation_channels is None:
        initial = np.array(
            [area_guess, centre_guess, sigma_guess, np.log(background_guess), 0.0],
            dtype=float,
        )
        lower = np.array(
            [0.0, -centre_limit_channels, 0.5, np.log(1e-6), -2.0],
            dtype=float,
        )
        upper = np.array(
            [np.inf, centre_limit_channels, 12.0, np.log(1e9), 2.0],
            dtype=float,
        )
        area_parameter_index = 0
        centre_parameter_index = 1
        sigma_parameter_index = 2
    else:
        companion_neighbourhood = (
            np.abs(x - companion_separation_channels) <= 1.5 / slope
        )
        companion_area_guess = max(
            float(np.sum(excess_counts[companion_neighbourhood])),
            1.0,
        )
        initial = np.array(
            [
                area_guess,
                companion_area_guess,
                centre_guess,
                centre_guess,
                sigma_guess,
                sigma_guess,
                np.log(background_guess),
                0.0,
            ],
            dtype=float,
        )
        lower = np.array(
            [
                0.0,
                0.0,
                -centre_limit_channels,
                -centre_limit_channels,
                0.5,
                0.5,
                np.log(1e-6),
                -2.0,
            ],
            dtype=float,
        )
        upper = np.array(
            [
                np.inf,
                np.inf,
                centre_limit_channels,
                centre_limit_channels,
                12.0,
                12.0,
                np.log(1e9),
                2.0,
            ],
            dtype=float,
        )
        area_parameter_index = 0
        centre_parameter_index = 2
        sigma_parameter_index = 4

    try:
        fit = least_squares(
            _poisson_deviance_residuals,
            initial,
            bounds=(lower, upper),
            args=(x, roi_counts, companion_separation_channels),
            x_scale="jac",
            max_nfev=5000,
        )
    except (ValueError, RuntimeError, FloatingPointError) as error:
        return _failed_fit(str(error), reference_energy_kev)

    if not fit.success or not np.all(np.isfinite(fit.x)):
        return _failed_fit(fit.message, reference_energy_kev)

    area = fit.x[area_parameter_index]
    centre_shift = fit.x[centre_parameter_index]
    sigma = fit.x[sigma_parameter_index]
    if companion_separation_channels is None:
        companion_area = None
        log_background = fit.x[3]
        background_log_slope = fit.x[4]
    else:
        companion_area = float(fit.x[1])
        companion_energy = float(
            companion_energy_kev + fit.x[3] * slope
        )
        companion_fwhm = float(FWHM_TO_SIGMA * fit.x[5] * slope)
        log_background = fit.x[6]
        background_log_slope = fit.x[7]
    degrees_of_freedom = max(len(roi_counts) - len(fit.x), 1)
    poisson_deviance = float(np.sum(fit.fun**2))
    reduced_deviance = poisson_deviance / degrees_of_freedom

    try:
        covariance = np.linalg.pinv(fit.jac.T @ fit.jac)
        parameter_uncertainties = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    except np.linalg.LinAlgError:
        parameter_uncertainties = np.full(len(fit.x), np.nan)

    area_uncertainty = float(parameter_uncertainties[area_parameter_index])
    centre_uncertainty_channels = float(
        parameter_uncertainties[centre_parameter_index]
    )
    sigma_uncertainty_channels = float(
        parameter_uncertainties[sigma_parameter_index]
    )

    fitted_energy = reference_energy_kev + centre_shift * slope
    fitted_fwhm = FWHM_TO_SIGMA * sigma * slope
    fitted_fwhm_uncertainty = FWHM_TO_SIGMA * sigma_uncertainty_channels * slope

    at_centre_boundary = (
        abs(abs(centre_shift) - centre_limit_channels)
        <= max(1e-5, 1e-3 * centre_limit_channels)
    )
    quality_problems = []
    if not np.isfinite(area_uncertainty) or area_uncertainty <= 0:
        quality_problems.append("area uncertainty is not finite and positive")
    if not 0.3 <= fitted_fwhm <= 5.0:
        quality_problems.append("FWHM is outside 0.3--5.0 keV")
    if at_centre_boundary:
        quality_problems.append("centroid reached its allowed shift boundary")
    if reduced_deviance >= 5.0:
        quality_problems.append(
            f"reduced Poisson deviance is {reduced_deviance:.2f} (limit 5.00)"
        )
    valid = not quality_problems

    return {
        "success": True,
        "valid": bool(valid),
        "reason": "" if valid else "; ".join(quality_problems) + ".",
        "reference_energy_kev": float(reference_energy_kev),
        "energy_kev": float(fitted_energy),
        "energy_uncertainty_kev": float(centre_uncertainty_channels * slope),
        "fwhm_kev": float(fitted_fwhm),
        "fwhm_uncertainty_kev": float(fitted_fwhm_uncertainty),
        "area_counts": float(area),
        "area_uncertainty_counts": area_uncertainty,
        "background_counts_per_channel": float(np.exp(log_background)),
        "background_log_slope": float(background_log_slope),
        "companion_energy_kev": (
            float(companion_energy_kev)
            if companion_energy_kev is not None
            else None
        ),
        "companion_area_counts": companion_area,
        "companion_fitted_energy_kev": (
            companion_energy
            if companion_energy_kev is not None
            else None
        ),
        "companion_fwhm_kev": (
            companion_fwhm
            if companion_energy_kev is not None
            else None
        ),
        "poisson_deviance": poisson_deviance,
        "reduced_deviance": float(reduced_deviance),
        "roi_channel_count": int(len(roi_counts)),
    }


def fit_resolution_model(line_results):
    """Fit FWHM² = a + bE using strong, well-behaved fitted peaks."""
    energies = []
    squared_widths = []
    for result in line_results:
        fit = result["sample_fit"]
        if not fit.get("valid"):
            continue
        area_uncertainty = fit.get("area_uncertainty_counts", np.nan)
        if not np.isfinite(area_uncertainty) or area_uncertainty <= 0:
            continue
        gross_significance = fit["area_counts"] / area_uncertainty
        if gross_significance < 5.0:
            continue
        energies.append(fit["energy_kev"])
        squared_widths.append(fit["fwhm_kev"] ** 2)

    fallback = {
        "intercept_kev2": 2.50,
        "slope_kev": 0.0018,
        "point_count": len(energies),
        "source": "fallback",
    }
    if len(energies) < 3:
        return fallback

    energies = np.asarray(energies, dtype=float)
    squared_widths = np.asarray(squared_widths, dtype=float)
    keep = np.ones(len(energies), dtype=bool)
    parameters = np.array([2.50, 0.0018], dtype=float)

    for _ in range(4):
        design = np.column_stack(
            (np.ones(np.count_nonzero(keep)), energies[keep])
        )
        observed = squared_widths[keep]
        try:
            fitted = least_squares(
                lambda values: design @ values - observed,
                parameters,
                bounds=([0.25, 0.0], [25.0, 0.02]),
            )
        except (ValueError, RuntimeError, FloatingPointError):
            return fallback
        if not fitted.success:
            return fallback
        parameters = fitted.x

        residuals = (
            squared_widths
            - parameters[0]
            - parameters[1] * energies
        )
        centre = np.median(residuals[keep])
        mad = np.median(np.abs(residuals[keep] - centre))
        if mad <= 1e-12:
            break
        proposed = np.abs(residuals - centre) <= 3.5 * 1.4826 * mad
        if np.count_nonzero(proposed) < 3 or np.array_equal(proposed, keep):
            break
        keep = proposed

    return {
        "intercept_kev2": float(parameters[0]),
        "slope_kev": float(parameters[1]),
        "point_count": int(np.count_nonzero(keep)),
        "source": "fitted",
    }


def resolution_fwhm_kev(energy_kev, resolution_model):
    """Evaluate the fitted detector resolution at one or more energies."""
    energy = np.asarray(energy_kev, dtype=float)
    width_squared = (
        resolution_model["intercept_kev2"]
        + resolution_model["slope_kev"] * energy
    )
    widths = np.sqrt(np.maximum(width_squared, 0.0))
    widths = np.clip(
        widths,
        MINIMUM_RESOLUTION_FWHM_KEV,
        MAXIMUM_RESOLUTION_FWHM_KEV,
    )
    if widths.ndim == 0:
        return float(widths)
    return widths


def global_scan_threshold_sigma(local_sigma, effective_trials):
    """Convert a local sigma threshold to a Bonferroni scan threshold."""
    if local_sigma <= 0:
        raise ValueError("Local sigma threshold must be positive.")
    if effective_trials < 1:
        raise ValueError("Effective trial count must be at least one.")
    family_false_alarm_probability = norm.sf(local_sigma)
    return float(
        norm.isf(family_false_alarm_probability / effective_trials)
    )


def matched_filter_peak_scores(
    counts,
    energies_kev,
    resolution_model,
    minimum_energy_kev,
    maximum_energy_kev,
):
    """Calculate local resolution-matched peak scores across a spectrum."""
    counts = np.asarray(counts, dtype=float)
    energies = np.asarray(energies_kev, dtype=float)
    if counts.ndim != 1 or energies.ndim != 1 or len(counts) != len(energies):
        raise ValueError("Counts and energies must be equal 1-D arrays.")
    if len(energies) < 3 or np.any(np.diff(energies) <= 0):
        raise ValueError("Candidate search requires increasing energies.")
    if np.any(counts < 0):
        raise ValueError("Raw counts cannot be negative.")

    channel_width_kev = float(np.median(np.diff(energies)))
    scores = np.zeros(len(counts), dtype=float)
    eligible = np.flatnonzero(
        (energies >= minimum_energy_kev)
        & (energies <= maximum_energy_kev)
    )

    for centre_index in eligible:
        centre_energy = energies[centre_index]
        fwhm = resolution_fwhm_kev(centre_energy, resolution_model)
        sigma_kev = fwhm / FWHM_TO_SIGMA
        half_width_channels = max(
            int(np.ceil(5.0 * sigma_kev / channel_width_kev)),
            8,
        )
        lower = centre_index - half_width_channels
        upper = centre_index + half_width_channels + 1
        if lower < 0 or upper > len(counts):
            continue

        x = energies[lower:upper] - centre_energy
        observed = counts[lower:upper]
        left = x <= -3.0 * sigma_kev
        right = x >= 3.0 * sigma_kev
        if np.count_nonzero(left) < 3 or np.count_nonzero(right) < 3:
            continue

        left_level = max(float(np.median(observed[left])), 0.5)
        right_level = max(float(np.median(observed[right])), 0.5)
        left_x = float(np.mean(x[left]))
        right_x = float(np.mean(x[right]))
        log_slope = (
            np.log(right_level) - np.log(left_level)
        ) / (right_x - left_x)
        log_continuum = np.log(left_level) + log_slope * (x - left_x)
        expected_continuum = np.exp(log_continuum)

        template = np.exp(-0.5 * (x / sigma_kev) ** 2)
        template /= np.sum(template)
        x_scaled = x / max(float(np.max(np.abs(x))), channel_width_kev)
        nuisance = np.column_stack((np.ones(len(x)), x_scaled))
        weights = 1.0 / np.maximum(expected_continuum, 0.5)
        gram = nuisance.T @ (weights[:, None] * nuisance)
        cross = nuisance.T @ (weights * template)
        try:
            nuisance_component = nuisance @ np.linalg.solve(gram, cross)
        except np.linalg.LinAlgError:
            continue
        orthogonal_template = template - nuisance_component
        information = float(
            np.sum(weights * orthogonal_template**2)
        )
        if information <= 0 or not np.isfinite(information):
            continue
        numerator = float(
            np.sum(
                weights
                * orthogonal_template
                * (observed - expected_continuum)
            )
        )
        scores[centre_index] = max(numerator / np.sqrt(information), 0.0)

    return scores


def discover_unlisted_candidates(
    sample_counts,
    energies_kev,
    resolution_model,
    minimum_energy_kev,
    maximum_energy_kev,
    local_detection_sigma,
):
    """Return unlisted peak candidates that pass the global scan threshold."""
    representative_energy = 0.5 * (
        minimum_energy_kev + maximum_energy_kev
    )
    representative_fwhm = resolution_fwhm_kev(
        representative_energy,
        resolution_model,
    )
    # Use every scanned channel in the Bonferroni correction. Adjacent scores
    # are correlated, so this overcounts the independent trials.
    effective_trials = max(
        int(
            np.count_nonzero(
                (energies_kev >= minimum_energy_kev)
                & (energies_kev <= maximum_energy_kev)
            )
        ),
        1,
    )
    scan_threshold = global_scan_threshold_sigma(
        local_detection_sigma,
        effective_trials,
    )
    scores = matched_filter_peak_scores(
        sample_counts,
        energies_kev,
        resolution_model,
        minimum_energy_kev,
        maximum_energy_kev,
    )
    channel_width_kev = float(np.median(np.diff(energies_kev)))
    minimum_distance_channels = max(
        int(round(representative_fwhm / channel_width_kev)),
        1,
    )
    peak_indices, properties = find_peaks(
        scores,
        height=scan_threshold,
        distance=minimum_distance_channels,
    )
    candidates = [
        {
            "energy_kev": float(energies_kev[index]),
            "scan_significance_sigma": float(height),
        }
        for index, height in zip(
            peak_indices,
            properties["peak_heights"],
        )
    ]
    return candidates, {
        "effective_trials": effective_trials,
        "scan_threshold_sigma": scan_threshold,
        "maximum_scan_score_sigma": float(np.max(scores)),
    }
