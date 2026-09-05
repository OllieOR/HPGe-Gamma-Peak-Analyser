"""Plotting helpers for the HPGe spectra."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_MINIMUM_ENERGY_KEV = 200.0
DEFAULT_MAXIMUM_ENERGY_KEV = 2800.0


def _count_rate(counts, live_time_seconds):
    return np.asarray(counts, dtype=float) / (live_time_seconds / 60.0)


def _plot_mask(energies, minimum_energy_kev, maximum_energy_kev=None):
    if maximum_energy_kev is None:
        maximum_energy_kev = min(DEFAULT_MAXIMUM_ENERGY_KEV, float(energies[-1]))
    mask = (energies >= minimum_energy_kev) & (energies <= maximum_energy_kev)
    return mask, maximum_energy_kev


def _select_peak_labels(results, max_labels=8, minimum_spacing_kev=70.0):
    candidates = [
        result
        for result in results["detected_line_results"]
        if result["line"].emission_probability is not None
        and result["sample_fit"].get("valid", False)
    ]
    candidates.sort(
        key=lambda result: result["net"]["sample_peak_significance_sigma"],
        reverse=True,
    )

    selected = []
    for result in candidates:
        energy = result["sample_fit"]["energy_kev"]
        if any(
            abs(energy - other["sample_fit"]["energy_kev"])
            < minimum_spacing_kev
            for other in selected
        ):
            continue
        selected.append(result)
        if len(selected) == max_labels:
            break

    return sorted(selected, key=lambda result: result["sample_fit"]["energy_kev"])


def _sample_name(run):
    return Path(run["sample_filename"]).stem.replace("_", " ").title()


def _add_peak_labels(ax, analysis, minimum_energy_kev, maximum_energy_kev, max_labels):
    labels = _select_peak_labels(analysis, max_labels=max_labels)
    label_heights = (0.96, 0.80)
    for index, result in enumerate(labels):
        fit = result["sample_fit"]
        energy = fit["energy_kev"]
        if not minimum_energy_kev <= energy <= maximum_energy_kev:
            continue
        ax.axvline(energy, linewidth=0.7, alpha=0.25)
        ax.text(
            energy,
            label_heights[index % len(label_heights)],
            f"{result['line'].isotope}\n{energy:.0f} keV",
            transform=ax.get_xaxis_transform(),
            rotation=90,
            ha="right",
            va="top",
            fontsize=8,
        )


def make_spectrum_figure(
    run,
    analysis,
    minimum_energy_kev=DEFAULT_MINIMUM_ENERGY_KEV,
    maximum_energy_kev=None,
    max_labels=8,
):
    """Plot the corrected sample spectrum with its matched background removed."""
    energies = np.asarray(analysis["energies_kev"], dtype=float)
    net_rate = (
        np.asarray(analysis["corrected_counts_diagnostic"], dtype=float)
        / (run["sample_live_time"] / 60.0)
    )

    mask, maximum_energy_kev = _plot_mask(
        energies,
        minimum_energy_kev,
        maximum_energy_kev,
    )

    masked_energies = energies[mask]
    masked_net_rate = net_rate[mask]

    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.step(
        masked_energies,
        masked_net_rate,
        where="mid",
        linewidth=0.9,
        label="Corrected sample - matched background",
    )
    ax.axhline(0.0, linewidth=0.8, alpha=0.45)

    _add_peak_labels(ax, analysis, minimum_energy_kev, maximum_energy_kev, max_labels)

    upper_limit = max(float(np.max(masked_net_rate)), 1.0)
    lower_limit = min(float(np.min(masked_net_rate)), 0.0)
    padding = 0.08 * max(abs(upper_limit), abs(lower_limit), 1e-6)

    ax.set_ylim(lower_limit - padding, upper_limit + padding)
    ax.set_xlim(minimum_energy_kev, maximum_energy_kev)
    ax.set_title(
        f"{_sample_name(run)} - HPGe detector {run['detector']}\n"
        "Background-subtracted spectrum"
    )
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Corrected net count rate (counts min$^{-1}$)")
    ax.grid(axis="y", alpha=0.2)
    ax.legend(frameon=False, loc="upper right")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def make_sample_background_figure(
    run,
    analysis,
    minimum_energy_kev=DEFAULT_MINIMUM_ENERGY_KEV,
    maximum_energy_kev=None,
    max_labels=8,
):
    """Create a diagnostic overlay of the raw sample and matched background."""
    energies = np.asarray(analysis["energies_kev"], dtype=float)
    sample_rate = _count_rate(run["sample_counts"], run["sample_live_time"])
    background_rate = _count_rate(
        run["background_counts"], run["background_live_time"]
    )

    mask, maximum_energy_kev = _plot_mask(
        energies,
        minimum_energy_kev,
        maximum_energy_kev,
    )

    fig, ax = plt.subplots(figsize=(11, 6.2))
    ax.step(
        energies[mask],
        np.ma.masked_less_equal(sample_rate[mask], 0),
        where="mid",
        linewidth=0.9,
        label="Sample",
    )
    ax.step(
        energies[mask],
        np.ma.masked_less_equal(background_rate[mask], 0),
        where="mid",
        linewidth=0.8,
        alpha=0.75,
        label="Matched background",
    )

    _add_peak_labels(ax, analysis, minimum_energy_kev, maximum_energy_kev, max_labels)

    ax.set_title(
        f"{_sample_name(run)} - HPGe detector {run['detector']}\n"
        "Sample and matched background"
    )
    ax.set_xlabel("Energy (keV)")
    ax.set_ylabel("Count rate (counts min$^{-1}$)")
    ax.set_yscale("log")
    ax.set_xlim(minimum_energy_kev, maximum_energy_kev)
    ax.grid(axis="y", which="both", alpha=0.2)
    ax.legend(frameon=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    fig.tight_layout()
    return fig


def save_spectrum_figure(run, analysis, output_path, **plot_options):
    """Save the main background-subtracted spectrum figure."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig = make_spectrum_figure(run, analysis, **plot_options)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path



def save_sample_background_figure(run, analysis, output_path, **plot_options):
    """Save the diagnostic sample-vs-background figure."""
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    fig = make_sample_background_figure(run, analysis, **plot_options)
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    return output_path
