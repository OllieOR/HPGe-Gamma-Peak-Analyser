"""Regression tests for the analyser. Run with unittest discovery."""

import tempfile
import unittest
from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from main import DATA_FOLDER
from peak_library import GAMMA_LINES, GammaLine, nearest_gamma_line
from results import (
    print_analysis_report,
    save_analysis_report,
    save_results_csv,
)
from sample_data import BACKGROUND_INFO, DETECTOR_INFO, SAMPLE_INFO
from analysis import _combine_peak_fits
from calibration import (
    calibrate_energies,
    calculate_activity_estimate,
    calculate_sample_geometry_correction,
    subtract_background,
    validate_spectrum,
)
from peak_fitting import (
    _doublet_with_background,
    discover_unlisted_candidates,
    fit_raw_peak,
    fit_resolution_model,
    gaussian_with_background,
    global_scan_threshold_sigma,
    matched_filter_peak_scores,
    resolution_fwhm_kev,
)
from tka_reader import FIRST_CALIBRATION_CHANNEL, read_tka


SYNTHETIC_DETECTOR = {
    "slope": 0.25,
    "offset": 0.0,
    "eff_coefficient": 1.0,
    "eff_power": 0.0,
    "eff_coefficient_uncertainty": 0.0,
    "eff_power_uncertainty": 0.0,
    "efficiency_unit": "fraction",
}


def synthetic_peak(area, energy_shift_kev=0.0, background=8.0):
    channels = np.arange(FIRST_CALIBRATION_CHANNEL, 4098, dtype=float)
    reference_energy = 500.0
    expected_channel = reference_energy / SYNTHETIC_DETECTOR["slope"]
    shifted_channels = channels - expected_channel
    model = gaussian_with_background(
        shifted_channels,
        area,
        energy_shift_kev / SYNTHETIC_DETECTOR["slope"],
        3.2,
        np.log(background),
        0.0,
    )
    return channels, model, reference_energy


def _fit_pair(sample_counts, background_counts, channels, reference):
    sample_fit = fit_raw_peak(
        sample_counts,
        channels,
        SYNTHETIC_DETECTOR,
        reference,
    )
    background_fit = fit_raw_peak(
        background_counts,
        channels,
        SYNTHETIC_DETECTOR,
        reference,
    )
    return sample_fit, background_fit


def _net_peak(
    sample_fit,
    background_fit,
    reference,
    geometry=1.0,
    emission_probability=1.0,
):
    return _combine_peak_fits(
        sample_fit,
        background_fit,
        sample_live_time_seconds=1000.0,
        background_live_time_seconds=1000.0,
        sample_mass_kg=1.0,
        geometry_correction=geometry,
        geometry_uncertainty=0.0,
        detector_info=SYNTHETIC_DETECTOR,
        energy_kev=reference,
        emission_probability=emission_probability,
        detection_sigma=3.0,
    )


class InputTests(unittest.TestCase):
    def test_energy_calibration(self):
        channels = np.array([2.0, 3.0, 4.0])
        actual = calibrate_energies(channels, slope=2.0, offset=10.0)
        np.testing.assert_allclose(actual, [14.0, 16.0, 18.0])

    def test_tka_channel_offset(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "tiny.TKA"
            path.write_text("10\n12\n1\n2\n3\n", encoding="utf-8")
            channels, counts, live, real = read_tka(path)
        np.testing.assert_array_equal(channels, [2.0, 3.0, 4.0])
        np.testing.assert_array_equal(counts, [1, 2, 3])
        self.assertEqual(live, 10.0)
        self.assertEqual(real, 12.0)

    def test_background_scaling(self):
        actual = subtract_background(
            np.array([10.0, 10.0, 10.0]),
            np.array([20.0, 10.0, 0.0]),
            sample_live_time=10.0,
            background_live_time=20.0,
        )
        np.testing.assert_allclose(actual, [0.0, 5.0, 10.0])

    def test_invalid_spectrum(self):
        with self.assertRaises(ValueError):
            validate_spectrum([1.0, 2.0], [1.0])
        with self.assertRaises(ValueError):
            validate_spectrum([1.0, -2.0], [1.0, 2.0])


class PeakFitTests(unittest.TestCase):
    def test_poisson_peak_area(self):
        channels, counts, reference = synthetic_peak(area=2400.0)
        fit = fit_raw_peak(
            counts,
            channels,
            SYNTHETIC_DETECTOR,
            reference,
        )
        self.assertTrue(fit["valid"])
        self.assertAlmostEqual(fit["area_counts"], 2400.0, delta=2.0)
        self.assertAlmostEqual(fit["energy_kev"], reference, delta=0.002)

    def test_identical_spectra(self):
        channels, counts, reference = synthetic_peak(area=2400.0)
        sample_fit, background_fit = _fit_pair(
            counts,
            counts.copy(),
            channels,
            reference,
        )
        net = _net_peak(sample_fit, background_fit, reference)
        self.assertTrue(net["quantifiable"])
        self.assertAlmostEqual(net["net_area_counts"], 0.0, delta=1e-7)
        self.assertAlmostEqual(net["significance_sigma"], 0.0, delta=1e-7)
        self.assertFalse(net["detected"])

    def test_shifted_background_peak(self):
        channels, sample_counts, reference = synthetic_peak(
            area=3600.0,
            energy_shift_kev=-0.30,
        )
        _, background_counts, _ = synthetic_peak(
            area=2400.0,
            energy_shift_kev=0.35,
        )
        sample_fit, background_fit = _fit_pair(
            sample_counts,
            background_counts,
            channels,
            reference,
        )
        net = _net_peak(sample_fit, background_fit, reference)
        self.assertAlmostEqual(net["net_area_counts"], 1200.0, delta=4.0)
        self.assertTrue(net["detected"])

    def test_geometry_before_background(self):
        channels, sample_counts, reference = synthetic_peak(area=3000.0)
        _, background_counts, _ = synthetic_peak(area=2000.0)
        sample_fit, background_fit = _fit_pair(
            sample_counts,
            background_counts,
            channels,
            reference,
        )
        net = _net_peak(
            sample_fit,
            background_fit,
            reference,
            geometry=2.0,
        )
        self.assertAlmostEqual(
            net["raw_net_area_counts"],
            1000.0,
            delta=4.0,
        )
        self.assertAlmostEqual(
            net["net_area_counts"],
            4000.0,
            delta=6.0,
        )

    def test_238_242_doublet(self):
        channels = np.arange(FIRST_CALIBRATION_CHANNEL, 4098, dtype=float)
        target_energy = 238.632
        companion_energy = 241.995
        expected_channel = target_energy / SYNTHETIC_DETECTOR["slope"]
        x = channels - expected_channel
        counts = _doublet_with_background(
            x,
            3000.0,
            1200.0,
            0.4,
            -0.2,
            3.0,
            3.3,
            np.log(10.0),
            0.0,
            (companion_energy - target_energy) / SYNTHETIC_DETECTOR["slope"],
        )
        fit = fit_raw_peak(
            counts,
            channels,
            SYNTHETIC_DETECTOR,
            target_energy,
            companion_energy_kev=companion_energy,
        )
        self.assertTrue(fit["valid"])
        self.assertAlmostEqual(fit["area_counts"], 3000.0, delta=4.0)
        self.assertAlmostEqual(fit["companion_area_counts"], 1200.0, delta=4.0)

    def test_lower_energy_companion(self):
        channels = np.arange(FIRST_CALIBRATION_CHANNEL, 4098, dtype=float)
        lower_energy = 238.632
        target_energy = 241.995
        expected_channel = target_energy / SYNTHETIC_DETECTOR["slope"]
        x = channels - expected_channel
        counts = _doublet_with_background(
            x,
            1200.0,
            3000.0,
            -0.2,
            0.4,
            3.3,
            3.0,
            np.log(10.0),
            0.0,
            (lower_energy - target_energy)
            / SYNTHETIC_DETECTOR["slope"],
        )
        fit = fit_raw_peak(
            counts,
            channels,
            SYNTHETIC_DETECTOR,
            target_energy,
            companion_energy_kev=lower_energy,
        )
        self.assertTrue(fit["valid"])
        self.assertAlmostEqual(fit["area_counts"], 1200.0, delta=4.0)
        self.assertAlmostEqual(fit["companion_area_counts"], 3000.0, delta=4.0)

    def test_unknown_peak_has_no_activity(self):
        channels, sample_counts, reference = synthetic_peak(area=3200.0)
        _, background_counts, _ = synthetic_peak(area=1000.0)
        sample_fit, background_fit = _fit_pair(
            sample_counts,
            background_counts,
            channels,
            reference,
        )
        net = _net_peak(
            sample_fit,
            background_fit,
            reference,
            emission_probability=None,
        )
        self.assertTrue(net["detected"])
        self.assertGreater(net["geometry_corrected_rate_cpm_per_kg"], 0.0)
        self.assertIsNone(net["activity_bq_per_kg"])


class CandidateScanTests(unittest.TestCase):
    @staticmethod
    def _resolution_model():
        return {
            "intercept_kev2": 2.5,
            "slope_kev": 0.0018,
            "point_count": 0,
            "source": "test",
        }

    def test_resolution_fit(self):
        line_results = []
        for energy in (250.0, 500.0, 1000.0, 1500.0, 2500.0):
            fwhm = np.sqrt(2.8 + 0.0015 * energy)
            line_results.append(
                {
                    "line": GammaLine("Test", energy, 1.0, "Synthetic"),
                    "sample_fit": {
                        "valid": True,
                        "energy_kev": energy,
                        "fwhm_kev": fwhm,
                        "area_counts": 1000.0,
                        "area_uncertainty_counts": 20.0,
                    },
                }
            )
        model = fit_resolution_model(line_results)
        self.assertEqual(model["source"], "fitted")
        self.assertAlmostEqual(model["intercept_kev2"], 2.8, places=5)
        self.assertAlmostEqual(model["slope_kev"], 0.0015, places=7)
        self.assertAlmostEqual(
            resolution_fwhm_kev(1000.0, model),
            np.sqrt(4.3),
            places=5,
        )

    def test_scan_threshold(self):
        threshold = global_scan_threshold_sigma(3.0, effective_trials=10000)
        self.assertGreater(threshold, 5.0)

    def test_scan_finds_peak(self):
        energies = np.arange(200.0, 800.0, 0.25)
        model = self._resolution_model()
        fwhm = resolution_fwhm_kev(503.25, model)
        sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
        channel_width = energies[1] - energies[0]
        counts = (
            np.full(len(energies), 20.0)
            + 500.0
            * channel_width
            / (sigma * np.sqrt(2.0 * np.pi))
            * np.exp(-0.5 * ((energies - 503.25) / sigma) ** 2)
        )
        candidates, diagnostics = discover_unlisted_candidates(
            counts,
            energies,
            model,
            minimum_energy_kev=220.0,
            maximum_energy_kev=780.0,
            local_detection_sigma=3.0,
        )
        self.assertTrue(
            any(abs(item["energy_kev"] - 503.25) < 0.5 for item in candidates)
        )
        self.assertGreater(diagnostics["scan_threshold_sigma"], 4.0)

    def test_flat_continuum(self):
        energies = np.arange(200.0, 800.0, 0.25)
        counts = np.full(len(energies), 20.0)
        scores = matched_filter_peak_scores(
            counts,
            energies,
            self._resolution_model(),
            minimum_energy_kev=220.0,
            maximum_energy_kev=780.0,
        )
        self.assertAlmostEqual(float(np.max(scores)), 0.0, places=10)


class MetadataTests(unittest.TestCase):
    def test_nearest_line(self):
        line = nearest_gamma_line(609.40, tolerance_kev=1.5)
        self.assertIsNotNone(line)
        self.assertEqual(line.isotope, "Bi-214")
        self.assertIsNone(nearest_gamma_line(700.0, tolerance_kev=1.5))

    def test_gamma_library(self):
        energies = [line.energy_kev for line in GAMMA_LINES]
        self.assertEqual(energies, sorted(energies))
        self.assertEqual(len(energies), len(set(energies)))
        self.assertGreaterEqual(len(energies), 40)
        for line in GAMMA_LINES:
            self.assertGreater(line.emission_probability, 0.0)
            self.assertLessEqual(line.emission_probability, 1.0)

    def test_geometry_values(self):
        correction = calculate_sample_geometry_correction(
            detector=1,
            filled_height_cm=2.4,
            beaker_filled=True,
        )
        self.assertAlmostEqual(correction, 4.101590, places=6)

    def test_tka_metadata(self):
        for filename, metadata in SAMPLE_INFO.items():
            with self.subTest(filename=filename):
                channels, counts, live, real = read_tka(DATA_FOLDER / filename)
                self.assertEqual(len(channels), len(counts))
                self.assertAlmostEqual(
                    live,
                    metadata["expected_live_time_s"],
                    delta=0.5,
                )
                self.assertLessEqual(live, real)

    def test_detector_assignments(self):
        reference_energy = 1460.822
        expected_channels = {
            detector: (
                reference_energy - detector_info["offset"]
            )
            / detector_info["slope"]
            for detector, detector_info in DETECTOR_INFO.items()
        }

        def k40_excess(counts, expected_channel):
            centre = int(
                round(expected_channel - FIRST_CALIBRATION_CHANNEL)
            )
            peak_mean = np.mean(counts[centre - 8:centre + 9])
            side_mean = np.mean(
                np.r_[
                    counts[centre - 50:centre - 25],
                    counts[centre + 25:centre + 50],
                ]
            )
            return peak_mean - side_mean

        registered = [
            (filename, metadata["detector"])
            for filename, metadata in SAMPLE_INFO.items()
        ]
        registered.extend(
            (metadata["filename"], detector)
            for detector, metadata in BACKGROUND_INFO.items()
        )

        for filename, recorded_detector in registered:
            with self.subTest(filename=filename):
                _, counts, _, _ = read_tka(DATA_FOLDER / filename)
                scores = {
                    detector: k40_excess(counts, channel)
                    for detector, channel in expected_channels.items()
                }
                inferred_detector = max(scores, key=scores.get)
                self.assertEqual(inferred_detector, recorded_detector)

    def test_parsnip_mass(self):
        self.assertEqual(SAMPLE_INFO["Parsnip .TKA"]["sample_mass_g"], 173.8)

        for metadata in BACKGROUND_INFO.values():
            filename = metadata["filename"]
            with self.subTest(filename=filename):
                _, _, live, real = read_tka(DATA_FOLDER / filename)
                self.assertAlmostEqual(
                    live,
                    metadata["expected_live_time_s"],
                    delta=0.5,
                )
                self.assertLessEqual(live, real)

    def test_efficiency_units(self):
        for detector_info in DETECTOR_INFO.values():
            self.assertEqual(detector_info["efficiency_unit"], "fraction")

    def test_activity_calculation(self):
        result = calculate_activity_estimate(
            corrected_area_counts=2000.0,
            statistical_uncertainty_counts=20.0,
            sample_area_counts=1000.0,
            sample_live_time_seconds=100.0,
            sample_mass_kg=2.0,
            energy_kev=500.0,
            emission_probability=0.5,
            geometry_uncertainty=0.0,
            detector_info={
                "eff_coefficient": 0.25,
                "eff_power": 0.0,
                "eff_coefficient_uncertainty": 0.0,
                "eff_power_uncertainty": 0.0,
                "efficiency_unit": "fraction",
            },
        )
        self.assertAlmostEqual(result["activity_bq_per_kg"], 80.0)


class ReportTests(unittest.TestCase):
    @staticmethod
    def _default_analysis():
        from main import analyse_file

        _, analysis = analyse_file("CARROTSDET1.TKA")
        return analysis

    def test_default_report(self):
        stream = StringIO()
        with redirect_stdout(stream):
            print_analysis_report(self._default_analysis())
        output = stream.getvalue()

        self.assertIn("41 peaks detected at >= 3.0 sigma", output)
        self.assertIn("Energy±u/keV", output)
        self.assertIn("Stat/σ", output)
        self.assertIn("Act±u/Bqkg", output)
        self.assertIn("K-40", output)
        self.assertTrue(
            all(len(line) <= 66 for line in output.splitlines()),
            "Default report must fit a 66-column terminal without wrapping.",
        )
        self.assertNotIn("Sample fit:", output)
        self.assertNotIn("Background fit:", output)
        self.assertNotIn("Raw detector-space difference", output)
        self.assertNotIn("Known uncertainty components", output)
        self.assertNotIn("Uncertainty:", output)
        self.assertNotIn("Activities from different gamma lines", output)

    def test_verbose_report(self):
        stream = StringIO()
        with redirect_stdout(stream):
            print_analysis_report(self._default_analysis(), verbose=True)
        output = stream.getvalue()

        self.assertIn("FIT DETAILS", output)
        self.assertIn("Sample fit:", output)
        self.assertIn("Background fit:", output)
        self.assertIn("Raw detector-space difference", output)
        self.assertIn("Known uncertainty components", output)


class FileOutputTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from main import analyse_file

        cls.run_info, cls.analysis = analyse_file(
            "BARROWLEEK.TKA", discover_candidates=False
        )

    def test_csv_output(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "peaks.csv"
            save_results_csv(self.analysis, path)
            text = path.read_text(encoding="utf-8")

        self.assertIn("reference_energy_kev", text)
        self.assertIn("K-40", text)

    def test_text_report(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "report.txt"
            save_analysis_report(
                self.analysis,
                path,
                header_lines=("GAMMA PEAK ANALYSER - DETAILED REPORT",),
            )
            text = path.read_text(encoding="utf-8")

        self.assertIn("GAMMA PEAK ANALYSER - DETAILED REPORT", text)
        self.assertIn("FIT DETAILS", text)
        self.assertIn("K-40", text)

    def test_spectrum_figure(self):
        from plotting import make_spectrum_figure, save_spectrum_figure

        figure = make_spectrum_figure(self.run_info, self.analysis)
        plotted_rate = figure.axes[0].lines[0].get_ydata()
        energies = self.analysis["energies_kev"]
        mask = (energies >= 200.0) & (energies <= 2800.0)
        expected_rate = (
            self.analysis["corrected_counts_diagnostic"][mask]
            / (self.run_info["sample_live_time"] / 60.0)
        )
        np.testing.assert_allclose(plotted_rate, expected_rate)
        lower_limit, upper_limit = figure.axes[0].get_ylim()
        self.assertLess(lower_limit, np.min(expected_rate))
        self.assertGreater(upper_limit, np.max(expected_rate))
        plt.close(figure)

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "spectrum.png"
            save_spectrum_figure(self.run_info, self.analysis, path)
            self.assertTrue(path.exists())
            self.assertGreater(path.stat().st_size, 10_000)

    def test_background_figure(self):
        from plotting import save_sample_background_figure

        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "sample_vs_background.png"
            save_sample_background_figure(self.run_info, self.analysis, path)
            self.assertTrue(path.exists())
            self.assertGreater(path.stat().st_size, 10_000)


if __name__ == "__main__":
    unittest.main(verbosity=2)
