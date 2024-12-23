"""
Main deconvolution function that uses the graph operations and navigation
methods to identify isotopic envelopes among the provided peaks.
"""

from typing import List, Tuple, Literal
from .graph_ops import construct_graph
from .navigation import navigate_left, navigate_right
from .data_structures import DeconvolutedPeak


def deconvolute(
        peaks: List[Tuple[float, float]],
        tolerance: float = 50,
        tolerance_type: Literal['ppm', 'da'] = 'ppm',
        charge_range: Tuple[int, int] = (1, 3),
        max_left_decrease: float = 0.6,
        max_right_decrease: float = 0.9,
) -> List[DeconvolutedPeak]:
    """
    Performs the main deconvolution procedure:
      1. Builds a graph of peaks connected by isotope spacing.
      2. Separates the graph by charge.
      3. Iterates over peaks in descending intensity,
         navigating left and right to form isotopic envelopes.

    Args:
        peaks: A list of (mz, intensity) tuples.
        tolerance: Numeric tolerance (in ppm or da).
        tolerance_type: Either 'ppm' or 'da'.
        charge_range: (min_charge, max_charge) to consider.
        max_left_decrease: Max fraction drop allowed to go left.
        max_right_decrease: Max fraction drop allowed to go right.
    Returns:
        A list of DeconvolutedPeak objects containing identified isotopic envelopes.
    """
    graph = construct_graph(peaks, tolerance, tolerance_type, charge_range)

    # For quick peak lookup
    indexed_peaks = {i: p for i, p in enumerate(peaks)}
    dpeaks = []

    # Sort peaks by intensity descending
    sorted_peaks = sorted(range(len(peaks)), key=lambda x: indexed_peaks[x][1], reverse=True)

    for peak_idx in sorted_peaks:
        # Skip if peak is already visited
        if graph[peak_idx].seen:
            continue

        # Attempt each charge from highest to lowest,
        # build a set of candidate peaks using navigate_left/right
        results = {}
        for charge in range(charge_range[1], charge_range[0] - 1, -1):
            left_peaks = navigate_left(graph, peak_idx, charge, max_left_decrease)
            right_peaks = navigate_right(graph, peak_idx, charge, max_right_decrease)
            peaks_in_profile = sorted(set(left_peaks + right_peaks))

            decon_peak = DeconvolutedPeak(
                peaks=[graph[p].value for p in peaks_in_profile],
                charge=charge,
            )
            results[charge] = decon_peak

        # Pick the best (highest total intensity) result among tested charges
        best_charge = max(results, key=lambda c: results[c].total_intensity)
        best_result = results[best_charge]

        # If only a single peak is found, we treat charge as unknown (None).
        if best_result.num_peaks == 1:
            best_result.charge = None

        # Mark all peaks in the best_result as seen
        for p_idx in best_result.peaks:
            graph[p_idx.index].seen = True

        dpeaks.append(best_result)

    return dpeaks
