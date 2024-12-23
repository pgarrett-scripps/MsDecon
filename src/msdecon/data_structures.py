"""
Data structure classes and constants used throughout the package.
"""

import dataclasses
from typing import List, Tuple, Optional

NEUTRON_MASS = 1.00866491578
PROTON_MASS = 1.007276466812


@dataclasses.dataclass
class Peak:
    """
    Represents a single peak with its m/z (mass-to-charge)
    and intensity values.
    """
    mz: float
    intensity: float
    index: Optional[int] = None


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
    peaks: List[Peak]
    charge: Optional[int]

    @property
    def base_peak(self) -> Peak:
        """
        Returns the monoisotopic peak from the list of peaks.
        """
        return self.peaks[0]

    @property
    def largest_peak(self) -> Peak:
        """
        Returns the largest peak from the list of peaks.
        """
        return max(self.peaks, key=lambda x: x.intensity)

    @property
    def mz_window(self) -> Tuple[float, float]:
        """
        Returns the m/z window of the isotopic envelope.
        """
        return self.peaks[0].mz, self.peaks[-1].mz

    @property
    def intensity_window(self) -> Tuple[float, float]:
        """
        Returns the intensity window of the isotopic envelope.
        """
        return self.peaks[0].intensity, self.peaks[-1].intensity

    @property
    def total_intensity(self) -> float:
        """
        Returns the total intensity of the isotopic envelope.
        """
        return sum(p.intensity for p in self.peaks)

    @property
    def num_peaks(self) -> int:
        """
        Returns the number of peaks in the isotopic envelope.
        """
        return len(self.peaks)

    def base_peak_neutral_mass(self, charge_carrier: float = PROTON_MASS) -> Optional[float]:
        """
        Returns the neutral (uncharged) mass if charge is known.
        """
        if self.charge is None or self.base_peak.mz is None:
            return None
        return self.base_peak.mz * self.charge - self.charge * charge_carrier

    def largest_peak_neutral_mass(self, charge_carrier: float = PROTON_MASS) -> Optional[float]:
        """
        Returns the neutral (uncharged) mass if charge is known.
        """
        if self.charge is None or self.largest_peak.mz is None:
            return None
        return self.largest_peak.mz * self.charge - self.charge * charge_carrier

    @property
    def isotope_gaps_mz(self) -> List[float]:
        """
        Returns a list of IsotopeGap objects between peaks.
        """
        gaps = []
        for i in range(1, len(self.peaks)):
            gap = self.peaks[i].mz - self.peaks[i - 1].mz
            gaps.append(gap)

        return gaps

    @property
    def isotope_gaps_neutral_mass(self) -> List[float]:
        """
        Returns a list of IsotopeGap objects between peaks.
        """
        isotope_gaps = self.isotope_gaps_mz
        gaps = []
        for gap in isotope_gaps:
            gaps.append(gap * self.charge)

        return gaps

    @property
    def isotope_gaps_ppm_error(self) -> List[float]:
        """
        Returns a list of ppm errors associated with isotope gaps.
        """
        gap_errors = []
        for i in range(1, len(self.peaks)):
            gap_errors.append((self.peaks[i - 1].mz + (NEUTRON_MASS / self.charge) - self.peaks[i].mz) /
                              self.peaks[i].mz * 1e6)

        return gap_errors


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
