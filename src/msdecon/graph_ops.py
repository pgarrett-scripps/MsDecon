"""
Graph construction and subgraph separation utilities.
"""

import rustworkx as rx
from typing import List, Tuple, Literal, Dict
from .data_structures import Peak, GraphNode, GraphEdge, NEUTRON_MASS


def get_tolerance(mz: float, tolerance: float, tolerance_type: Literal['ppm', 'da']) -> float:
    """
    Returns the absolute tolerance in Da, given a base m/z and either
    a parts-per-million (ppm) or Dalton (da) setting.
    """
    if tolerance_type == 'ppm':
        return mz * tolerance / 1e6
    elif tolerance_type == 'da':
        return tolerance


def construct_graph(
        peaks: List[Tuple[float, float]],
        tolerance: float,
        tolerance_type: Literal['ppm', 'da'],
        charge_range: Tuple[int, int]
) -> rx.PyGraph:
    """
    Constructs a rustworkx PyGraph from a list of (mz, intensity) peaks,
    adding edges between peaks if they fall within the expected isotope spacing.
    """
    graph = rx.PyGraph()
    graph.add_nodes_from([GraphNode(Peak(mz=mz, intensity=intensity, index=i)) for i, (mz, intensity)
                          in enumerate(peaks)])

    for index in graph.node_indices():
        # Assign node indices for easy reference
        graph[index].index = index

    # Determine the smallest and largest offset for the given charge range
    min_isotope_offset = NEUTRON_MASS / charge_range[1]
    max_isotope_offset = NEUTRON_MASS / charge_range[0]

    # Precompute valid offsets for each charge
    valid_isotope_offsets = []
    for charge in range(charge_range[0], charge_range[1] + 1):
        valid_isotope_offsets.append((charge, NEUTRON_MASS / charge))

    # Loop over peaks and connect them if they match an isotope offset
    for i in range(len(peaks)):
        for j in range(i + 1, len(peaks)):
            mz_i = peaks[i][0]
            mz_j = peaks[j][0]
            mz_diff = abs(mz_j - mz_i)

            tol = get_tolerance(mz_i, tolerance, tolerance_type)

            # If difference is too small or too large, skip early
            if mz_diff < min_isotope_offset - tol:
                continue
            if mz_diff > max_isotope_offset + tol:
                break

            # Check potential matches for each possible charge
            for charge, offset in valid_isotope_offsets:
                if abs(mz_diff - offset) <= tol:
                    from .data_structures import IsotopeGap  # local import is fine if needed
                    ppm_error = (mz_diff - offset) / mz_i * 1e6
                    peak_edge = GraphEdge(
                        IsotopeGap(
                            offset=offset,
                            charge=charge,
                            mz_error=mz_diff - offset,
                            ppm_error=ppm_error
                        )
                    )
                    graph.add_edge(i, j, peak_edge)

    # Assign edge indices for easy reference
    for index, data in graph.edge_index_map().items():
        data[2].index = index

    return graph
