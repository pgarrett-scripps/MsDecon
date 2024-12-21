# msdecon

**msdecon** is a Python package for performing simple isotopic envelope deconvolution in mass spectrometry data. It uses a graph-based approach to identify and group peaks belonging to the same isotopic distribution. By mapping peaks and their relationships into a graph, and then navigating through that graph, **msdecon** helps you automatically find monoisotopic peaks, estimate charge states, and retrieve isotopic clusters.

## Features

- **Graph Construction**: Builds a graph of peaks connected by isotope spacing.
- **Charge State Filtering**: Supports user-defined charge ranges to focus on relevant isotopes.
- **Automated Traversal**: Navigates through the graph to find isotopic envelopes.
- **Flexible Tolerances**: Choose between ppm or Da tolerances for peak matching.
- **Returns Detailed Results**: Provides monoisotopic peaks, charge states, total envelope intensity, and more.

## Installation

You can install **msdecon** from source or integrate it into your existing Python environment. For example:

```bash
pip install .
```

(Ensure you are in the directory with `setup.py` or `pyproject.toml`.)

## Usage

```python
from msdecon.deconvolution import deconvolute

# Example input: list of (mz, intensity) pairs
peaks = [
    (689.6649, 52.0),
    (689.93, 9.0),
    (689.9836, 71.0),
    # ...
]

results = deconvolute(
    peaks,
    tolerance=0.01,
    tolerance_type='da',
    charge_range=(1, 3)
)

for r in results:
    print(f"Monoisotopic Peak: {r.monoisotopic_peak}, Charge: {r.charge}, Total Intensity: {r.total_intensity}")
```

## Documentation

- **`deconvolute(peaks, tolerance, tolerance_type, charge_range, ...)`**  
  The primary function for isotope envelope deconvolution. Returns a list of `DeconvolutedPeak` objects.
  
- **`DeconvolutedPeak`**  
  A dataclass representing the identified isotopic envelope, including monoisotopic peak, charge state, and total intensity.

- **`Peak` / `IsotopeGap`**  
  Internal data structures representing basic MS entities.

- **`GraphNode` / `GraphEdge`**  
  Wrappers around the peaks and isotope gap data for use in the underlying graph structure.

## Dependencies

- [rustworkx](https://github.com/Qiskit/rustworkx) for efficient graph operations
- Python 3.8+ (recommended)

## Contributing

Pull requests and issues are welcome! Before contributing, please ensure you follow best practices in code formatting and provide tests or examples as needed.

## License

This project is licensed under the MIT License. See [LICENSE](LICENSE) for details.