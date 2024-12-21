# msdecon

**msdecon** is a Python package for performing simple isotopic deconvolution. It uses a graph-based approach 
to identify and group peaks belonging to the same isotopic distribution. It doesn't use any isotope distribution
scoring, rather it uses very simple logic when assigning isotopic distribution.

## Algorithm

- Start with the most abunant peak in the spectra.
- look left (lower m/z) and right (higher m/z) for the next peak in the isotopic distribution which is +/- a valid 
isotope offset (Neutron / charge), ~1 for +1, and ~0.5 for +2...
- Isotope distributions Must be sequential (no skipping peaks) and best be monotonically decreasing in intensity.
- Additionally to monotonically decreasing intensity, the intensity of the next peak should also be less than a certain 
scale factor of the previous peak (0.6 for left and 0.8 for right).

## Installation

```bash
pip install .
```

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

## Streamlit app

### Hosted on Streamlit Cloud (Try Me): [msdecon](https://msdecon.streamlit.app/)

![Deconvolution Example](example.gif)

This repo also contains the source code for a streamlit app which can be used to visualize the 
deconvolution results. To run the app, use the following command:

```bash
pip install requirements.txt
streamlit run streamlit_app.py
```