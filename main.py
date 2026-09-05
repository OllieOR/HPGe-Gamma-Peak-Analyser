"""Command-line entry point for the Gamma Peak Analyser."""

import argparse
from pathlib import Path

from analysis import analyse_spectrum
from results import (
    print_analysis_report,
    save_analysis_report,
    save_results_csv,
)
from sample_data import BACKGROUND_INFO, DETECTOR_INFO, SAMPLE_INFO
from tka_reader import read_tka


PROJECT_FOLDER = Path(__file__).parent
DATA_FOLDER = PROJECT_FOLDER / "data"
OUTPUT_FOLDER = PROJECT_FOLDER / "outputs"
DEFAULT_SAMPLE = "BARROWLEEK.TKA"


def _parse_arguments():
    parser = argparse.ArgumentParser(
        description=(
            "Fit environmental gamma lines in an HPGe sample spectrum using "
            "the detector-matched background."
        )
    )
    parser.add_argument(
        "sample",
        nargs="?",
        default=DEFAULT_SAMPLE,
        help=f"registered TKA filename (default: {DEFAULT_SAMPLE})",
    )
    parser.add_argument(
        "--detection-sigma",
        type=float,
        default=3.0,
        help="minimum corrected-area significance (default: 3.0)",
    )
    parser.add_argument(
        "--minimum-energy",
        type=float,
        default=200.0,
        help="lowest reference energy to analyse in keV (default: 200)",
    )
    parser.add_argument(
        "--show-nondetections",
        action="store_true",
        help="also print lines that did not pass the detection threshold",
    )
    parser.add_argument(
        "--verbose",
        action="store_true",
        help="print the underlying sample/background fits and uncertainties",
    )
    parser.add_argument(
        "--no-discovery",
        action="store_true",
        help="skip the optional search for unlisted peaks",
    )
    parser.add_argument(
        "--save-figure",
        type=Path,
        metavar="PATH",
        help=(
            "override the automatic background-subtracted spectrum path "
            "(default: outputs/<sample>_spectrum.png)"
        ),
    )
    parser.add_argument(
        "--save-report",
        type=Path,
        metavar="PATH",
        help=(
            "override the automatic detailed-report path "
            "(default: outputs/<sample>_report.txt)"
        ),
    )
    parser.add_argument(
        "--save-csv",
        type=Path,
        metavar="PATH",
        help="save the fitted peak table as CSV",
    )
    parser.add_argument(
        "--no-output-files",
        action="store_true",
        help="print results only; do not create the automatic report or figure",
    )
    parser.add_argument(
        "--list-samples",
        action="store_true",
        help="list registered sample filenames and exit",
    )
    return parser.parse_args()



def _output_stem(sample_filename):
    stem = Path(sample_filename).stem.strip().lower()
    return "_".join(stem.replace("-", " ").split())


def _report_header(run):
    metadata = run["sample_metadata"]
    return (
        "GAMMA PEAK ANALYSER - DETAILED REPORT",
        f"Sample: {run['sample_filename']}",
        f"Detector: {run['detector']}",
        f"Sample mass: {metadata['sample_mass_g']:.1f} g",
        f"Sample live time: {run['sample_live_time']:.1f} s",
        f"Background: {run['background_filename']}",
        f"Background live time: {run['background_live_time']:.1f} s",
    )


def _read_and_check(file_path, expected_live_time):
    channels, counts, live_time, real_time = read_tka(file_path)
    if abs(live_time - expected_live_time) > 0.5:
        raise ValueError(
            f"{file_path.name}: live time {live_time:g} s does not match "
            f"the registered value {expected_live_time:g} s."
        )
    return channels, counts, live_time, real_time


def analyse_file(
    sample_filename,
    detection_sigma=3.0,
    minimum_energy_kev=200.0,
    discover_candidates=True,
):
    """Load one registered sample and run the full analysis."""
    if sample_filename not in SAMPLE_INFO:
        available = ", ".join(sorted(SAMPLE_INFO))
        raise ValueError(
            f"Unknown sample {sample_filename!r}. Registered samples: {available}"
        )

    sample_metadata = SAMPLE_INFO[sample_filename]
    detector = sample_metadata["detector"]
    background_metadata = BACKGROUND_INFO[detector]

    sample_file = DATA_FOLDER / sample_filename
    background_file = DATA_FOLDER / background_metadata["filename"]

    sample_channels, sample_counts, sample_live_time, sample_real_time = (
        _read_and_check(sample_file, sample_metadata["expected_live_time_s"])
    )
    (
        background_channels,
        background_counts,
        background_live_time,
        background_real_time,
    ) = _read_and_check(
        background_file,
        background_metadata["expected_live_time_s"],
    )

    if len(sample_channels) != len(background_channels):
        raise ValueError("Sample and background contain different channel counts.")

    analysis = analyse_spectrum(
        sample_counts=sample_counts,
        background_counts=background_counts,
        sample_live_time_seconds=sample_live_time,
        background_live_time_seconds=background_live_time,
        sample_mass_kg=sample_metadata["sample_mass_g"] / 1000.0,
        detector=detector,
        detector_info=DETECTOR_INFO[detector],
        filled_height_cm=sample_metadata["filled_height_cm"],
        beaker_filled=True,
        minimum_energy_kev=minimum_energy_kev,
        detection_sigma=detection_sigma,
        discover_candidates=discover_candidates,
    )

    run_information = {
        "sample_filename": sample_filename,
        "sample_metadata": sample_metadata,
        "sample_channels": sample_channels,
        "sample_counts": sample_counts,
        "sample_live_time": sample_live_time,
        "sample_real_time": sample_real_time,
        "background_filename": background_metadata["filename"],
        "background_channels": background_channels,
        "background_counts": background_counts,
        "background_live_time": background_live_time,
        "background_real_time": background_real_time,
        "detector": detector,
    }
    return run_information, analysis


def main():
    args = _parse_arguments()

    if args.list_samples:
        print("Registered samples:")
        for filename in sorted(SAMPLE_INFO):
            print(f"  {filename}")
        return

    try:
        run, analysis = analyse_file(
            args.sample,
            detection_sigma=args.detection_sigma,
            minimum_energy_kev=args.minimum_energy,
            discover_candidates=not args.no_discovery,
        )
    except (FileNotFoundError, ValueError) as error:
        raise SystemExit(str(error)) from error

    print("\nGAMMA PEAK ANALYSER")
    print(
        f"Sample: {run['sample_filename']} | Detector: {run['detector']} | "
        f"Mass: {run['sample_metadata']['sample_mass_g']:.1f} g"
    )
    print_analysis_report(
        analysis,
        show_nondetections=args.show_nondetections,
        verbose=args.verbose,
    )

    if not args.no_output_files:
        from plotting import save_spectrum_figure

        stem = _output_stem(run["sample_filename"])
        report_path = args.save_report or OUTPUT_FOLDER / f"{stem}_report.txt"
        figure_path = args.save_figure or OUTPUT_FOLDER / f"{stem}_spectrum.png"

        report_output = save_analysis_report(
            analysis,
            report_path,
            header_lines=_report_header(run),
            show_nondetections=args.show_nondetections,
        )
        figure_output = save_spectrum_figure(run, analysis, figure_path)

        print(f"\nSaved detailed report: {report_output}")
        print(f"Saved background-subtracted figure: {figure_output}")

    if args.save_csv:
        output = save_results_csv(
            analysis,
            args.save_csv,
            include_nondetections=args.show_nondetections,
        )
        print(f"Saved CSV: {output}")


if __name__ == "__main__":
    main()
