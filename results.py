"""Console and CSV output for the analyser."""

import csv
from contextlib import redirect_stdout
from pathlib import Path



def _format_fit(fit):
    if not fit["success"]:
        return f"fit failed ({fit['reason']})"
    validity = "" if fit["valid"] else " [quality check failed]"
    return (
        f"{fit['energy_kev']:.3f} ± {fit['energy_uncertainty_kev']:.3f} keV, "
        f"FWHM {fit['fwhm_kev']:.3f} keV, "
        f"area {fit['area_counts']:.1f} ± "
        f"{fit['area_uncertainty_counts']:.1f} counts, "
        f"reduced deviance {fit['reduced_deviance']:.2f}{validity}"
    )


def _format_value_uncertainty(value, uncertainty, decimals=3):
    return f"{value:.{decimals}f}±{uncertainty:.{decimals}f}"


def _print_peak_table(detected):
    headings = (
        ("#", 2, ">"),
        ("Isotope", 7, "<"),
        ("Energy±u/keV", 14, ">"),
        ("FWHM", 5, ">"),
        ("Stat/σ", 6, ">"),
        ("Rate±u/cpmkg", 13, ">"),
        ("Act±u/Bqkg", 13, ">"),
    )

    header = " ".join(
        f"{label:{alignment}{width}}"
        for label, width, alignment in headings
    )
    print(header)
    print("-" * len(header))

    for index, result in enumerate(detected, start=1):
        line = result["line"]
        fit = result["sample_fit"]
        net = result["net"]
        fitted_energy = _format_value_uncertainty(
            fit["energy_kev"],
            fit["energy_uncertainty_kev"],
        )
        rate = _format_value_uncertainty(
            net["geometry_corrected_rate_cpm_per_kg"],
            net["geometry_corrected_rate_uncertainty_cpm_per_kg"],
        )
        if net["activity_bq_per_kg"] is None:
            activity = "n/a"
        else:
            activity = _format_value_uncertainty(
                net["activity_bq_per_kg"],
                net["known_standard_uncertainty_bq_per_kg"],
            )
        print(
            f"{index:>2} "
            f"{line.isotope:<7} "
            f"{fitted_energy:>14} "
            f"{fit['fwhm_kev']:>5.3f} "
            f"{net['significance_sigma']:>6.2f} "
            f"{rate:>13} "
            f"{activity:>13}"
        )


def _print_concise_nondetection(result):
    line = result["line"]
    fit = result["sample_fit"]
    net = result["net"]

    label = f"{line.isotope} ({line.energy_kev:.3f} keV)"
    if not net["quantifiable"]:
        print(f"  {label}: unavailable — {net['reason']}")
        return

    fitted_energy = (
        f", fitted {fit['energy_kev']:.3f} keV" if fit["success"] else ""
    )
    print(
        f"  {label}{fitted_energy}: "
        f"{net['net_area_counts']:.1f} ± "
        f"{net['net_area_uncertainty_counts']:.1f} counts "
        f"({net['significance_sigma']:.2f} sigma)"
    )


def _print_verbose_report(results, show_nondetections):
    detected = results["detected_line_results"]

    print("\nFIT DETAILS")
    print(
        "Sample geometry correction: "
        f"{results['geometry_correction']:.6f} (multiplicative)"
    )
    print(
        "Background live-time scale: "
        f"{results['background_live_time_scale']:.6f}"
    )
    print(
        "Detection rule: positive corrected fitted area "
        "(g*S - alpha*B) at or above "
        f"{results['detection_sigma']:.1f} statistical sigma"
    )
    discovery = results["discovery"]
    if discovery["enabled"]:
        resolution = discovery["resolution_model"]
        print(
            "Unlisted-candidate scan: "
            f"{discovery['effective_trials']} effective trials, "
            f"{discovery['scan_threshold_sigma']:.2f} sigma global threshold"
        )
        print(
            "Resolution model: FWHM^2 = "
            f"{resolution['intercept_kev2']:.4f} + "
            f"{resolution['slope_kev']:.6f}*E "
            f"({resolution['source']}, {resolution['point_count']} peaks)"
        )
        print(
            "Candidate scan counts "
            "(screened / fitted / corrected detections): "
            f"{discovery['screened_candidate_count']} / "
            f"{discovery['fitted_candidate_count']} / "
            f"{discovery['detected_candidate_count']}"
        )

    for index, result in enumerate(detected, start=1):
        line = result["line"]
        net = result["net"]
        print(
            f"\nPeak {index}: {line.isotope} "
            f"({line.energy_kev:.3f} keV reference)"
        )
        print(f"  Series/source: {line.series}")
        if line.emission_probability is None:
            print("  Identification: unlisted candidate; activity not assigned")
            print(
                "  Candidate-scan significance: "
                f"{result['candidate_scan_significance_sigma']:.2f} sigma"
            )
        else:
            print(
                "  Emission probability: "
                f"{100 * line.emission_probability:.2f}%"
            )
        print(f"  Sample fit: {_format_fit(result['sample_fit'])}")
        print(f"  Background fit: {_format_fit(result['background_fit'])}")
        print(
            "  Gross sample-peak significance: "
            f"{net['sample_peak_significance_sigma']:.2f} sigma"
        )
        print(
            "  Raw detector-space difference (S - alpha*B): "
            f"{net['raw_net_area_counts']:.1f} ± "
            f"{net['raw_net_area_uncertainty_counts']:.1f} counts "
            f"({net['raw_net_significance_sigma']:.2f} sigma)"
        )
        print(
            "  Corrected area (g*S - alpha*B): "
            f"{net['net_area_counts']:.1f} ± "
            f"{net['net_area_uncertainty_counts']:.1f} counts"
        )
        print(
            "  Corrected-area significance: "
            f"{net['significance_sigma']:.2f} sigma"
        )
        print(
            "  Local statistical decision threshold: "
            f"{net['statistical_decision_threshold_counts']:.1f} counts"
        )
        print(
            "  Corrected count rate: "
            f"{net['geometry_corrected_rate_cpm_per_kg']:.3f} ± "
            f"{net['geometry_corrected_rate_uncertainty_cpm_per_kg']:.3f} "
            "counts/min/kg"
        )
        if net["activity_bq_per_kg"] is None:
            print(
                "  Activity estimate: n/a until the radionuclide and "
                "emission probability are independently assigned"
            )
        else:
            print(
                "  Activity estimate: "
                f"{net['activity_bq_per_kg']:.3f} ± "
                f"{net['known_standard_uncertainty_bq_per_kg']:.3f} Bq/kg "
                "(known uncertainty terms only)"
            )
            geometry_component = net["geometry_uncertainty_bq_per_kg"]
            geometry_text = (
                f"{geometry_component:.3f}"
                if geometry_component is not None
                else "not available"
            )
            print(
                "  Known uncertainty components "
                "(fit+background / detector efficiency / geometry): "
                f"{net['statistical_uncertainty_bq_per_kg']:.3f} / "
                f"{net['detector_efficiency_uncertainty_bq_per_kg']:.3f} / "
                f"{geometry_text} Bq/kg"
            )
            print(
                "  Photopeak efficiency used: "
                f"{100 * net['detector_photopeak_efficiency']:.3f}%"
            )

    if show_nondetections:
        print("\nDETAILED NON-DETECTIONS")
        for result in results["line_results"]:
            if result["net"]["detected"]:
                continue
            line = result["line"]
            net = result["net"]
            print(f"\n{line.isotope} at {line.energy_kev:.3f} keV")
            print(f"  Sample fit: {_format_fit(result['sample_fit'])}")
            print(f"  Background fit: {_format_fit(result['background_fit'])}")
            if net["quantifiable"]:
                print(
                    "  Raw detector-space difference (S - alpha*B): "
                    f"{net['raw_net_area_counts']:.1f} ± "
                    f"{net['raw_net_area_uncertainty_counts']:.1f} counts "
                    f"({net['raw_net_significance_sigma']:.2f} sigma)"
                )
                print(
                    "  Corrected area (g*S - alpha*B): "
                    f"{net['net_area_counts']:.1f} ± "
                    f"{net['net_area_uncertainty_counts']:.1f} counts "
                    f"({net['significance_sigma']:.2f} sigma)"
                )
            else:
                print(f"  Net result unavailable: {net['reason']}")


def print_analysis_report(
    results,
    show_nondetections=False,
    verbose=False,
):
    """Print the compact table and optional fit details."""
    detected = results["detected_line_results"]

    if not detected:
        print(
            f"No peaks detected at >= "
            f"{results['detection_sigma']:.1f} sigma."
        )
    else:
        print(
            f"{len(detected)} peaks detected at >= "
            f"{results['detection_sigma']:.1f} sigma\n"
        )
        _print_peak_table(detected)

    if show_nondetections:
        print("\nNON-DETECTIONS")
        for result in results["line_results"]:
            if not result["net"]["detected"]:
                _print_concise_nondetection(result)

    if verbose:
        _print_verbose_report(results, show_nondetections)


def save_results_csv(results, output_path, include_nondetections=False):
    """Write fitted peak results to a CSV file and return its path."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "isotope",
        "series",
        "reference_energy_kev",
        "fitted_energy_kev",
        "energy_uncertainty_kev",
        "fwhm_kev",
        "significance_sigma",
        "corrected_rate_cpm_per_kg",
        "corrected_rate_uncertainty_cpm_per_kg",
        "activity_bq_per_kg",
        "activity_uncertainty_bq_per_kg",
        "detected",
    ]

    rows = (
        results["line_results"]
        if include_nondetections
        else results["detected_line_results"]
    )
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for result in rows:
            line = result["line"]
            fit = result["sample_fit"]
            net = result["net"]
            writer.writerow(
                {
                    "isotope": line.isotope,
                    "series": line.series,
                    "reference_energy_kev": line.energy_kev,
                    "fitted_energy_kev": fit.get("energy_kev"),
                    "energy_uncertainty_kev": fit.get("energy_uncertainty_kev"),
                    "fwhm_kev": fit.get("fwhm_kev"),
                    "significance_sigma": net.get("significance_sigma"),
                    "corrected_rate_cpm_per_kg": net.get(
                        "geometry_corrected_rate_cpm_per_kg"
                    ),
                    "corrected_rate_uncertainty_cpm_per_kg": net.get(
                        "geometry_corrected_rate_uncertainty_cpm_per_kg"
                    ),
                    "activity_bq_per_kg": net.get("activity_bq_per_kg"),
                    "activity_uncertainty_bq_per_kg": net.get(
                        "known_standard_uncertainty_bq_per_kg"
                    ),
                    "detected": net.get("detected", False),
                }
            )

    return output_path


def save_analysis_report(
    results,
    output_path,
    *,
    header_lines=(),
    show_nondetections=False,
):
    """Save a detailed text report and return its path."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open("w", encoding="utf-8") as handle:
        with redirect_stdout(handle):
            for line in header_lines:
                print(line)
            if header_lines:
                print()
            print_analysis_report(
                results,
                show_nondetections=show_nondetections,
                verbose=True,
            )

    return output_path
