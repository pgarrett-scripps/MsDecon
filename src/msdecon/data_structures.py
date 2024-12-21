"""
Data structure classes and constants used throughout the package.
"""

import dataclasses
from typing import List, Tuple, Optional

NEUTRON_MASS = 1.00866491578

@dataclasses.dataclass
class Peak:
    """
    Represents a single peak with its m/z (mass-to-charge)
    and intensity values.
    """
    mz: float
    intensity: float


@dataclasses.dataclass
class IsotopeGap:
    """
    Stores information about isotope spacing between two peaks.
    """
    offset: float
    charge: int
    mz_error: float
    ppm_error: float


@dataclasses.dataclass
class DeconvolutedPeak:
    """
    Represents the result of deconvolution on an isotopic envelope,
    capturing the monoisotopic peak, the largest peak, and overall info.
    """
    monoisotopic_peak: Optional[Peak]
    largest_peak: Peak
    total_intensity: float
    charge: Optional[int]
    mz_window: Tuple[float, float]
    peaks: List[int]  # Indices of peaks contributing to this envelope

    @property
    def neutral_mass(self) -> Optional[float]:
        """
        Returns the neutral (uncharged) mass if charge is known.
        """
        if self.charge is None or self.monoisotopic_peak is None:
            return None
        return self.monoisotopic_peak.mz * self.charge - self.charge * 1.007276466812


class GraphNode:
    """
    Node wrapper for rustworkx PyGraph nodes, tagging each node with
    the original peak data and an index.
    """
    def __init__(self, value: Peak):
        self.index = None
        self.value: Peak = value
        self.seen = False

    def __str__(self) -> str:
        return f"GraphNode: {self.value} @ index: {self.index}"


class GraphEdge:
    """
    Edge wrapper for rustworkx PyGraph edges, carrying isotope gap info.
    """
    def __init__(self, value: IsotopeGap):
        self.index = None
        self.value: IsotopeGap = value

    def __str__(self) -> str:
        return f"EdgeNode: {self.value} @ index: {self.index}"
