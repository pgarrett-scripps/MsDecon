"""
This Streamlit application allows you to deconvolute mass spectra using the msdecon package.
It takes input spectra in (mz intensity) format, filters by intensity, and then uses the
deconvolute function from msdecon to identify monoisotopic peaks and assign charge states.
"""

import numpy as np
import streamlit as st
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import plotly.graph_objects as go
import pandas as pd

from src.msdecon.deconvolution import deconvolute

st.set_page_config(
    page_title="MsDecon",
    page_icon=":bar_chart:",
    menu_items={
        'Get Help': 'https://github.com/pgarrett-scripps/MsDecon',
        'Report a bug': "https://github.com/pgarrett-scripps/MsDecon/issues",
        'About': "# This is a Streamlit app for deconvoluting mass spectra."
    }
)

st.title('Deconvolute Mass Sepctra :bar_chart:')

st.caption("""
This app deconvolutes mass spectra into monoisotopic peaks by combining intensities from isotopic clusters. 
When possible, charge states are assigned to peaks. The input spectra must be centroided.
""")

c1, c2 = st.columns(2)

# read default spectra 'default_spectra.txt'
default_spectra = open('default_spectra.txt', 'r').read()

# mz intensity\n
spectra = c1.text_area('Paste spectra here', default_spectra, height=125,
                       help='Paste the spectra in MS2 format: mz intensity\n')

deconvolute_tolerance = c2.number_input('Tolerance', value=50.0, step=0.001,
                                        help='The tolerance to apply when assigning isotopic peaks')
tolerance_type = c2.selectbox('Tolerance type', ['ppm', 'da'], index=0,
                              help='The type of tolerance to use for deconvolution')
min_charge, max_charge = st.slider('Charge range', 1, 10, (1, 3),
                                   help='The min and max charge to consider for deconvolution.')

with st.expander('Advanced Options'):
    c1, c2 = st.columns(2)
    min_left_intensity_decrease = c1.number_input('Min left intensity decrease', value=0.65,
                                                  help='Minimum sequential intensity decrease (%) when navigating left (lower m/z)')
    min_right_intensity_decrease = c2.number_input('Min right intensity decrease', value=0.95,
                                                   help='Minimum sequential intensity decrease (%) when navigating right (higher m/z)')
    min_intensity = st.number_input('Min intensity', value=0)


if tolerance_type == 'ppm':
    # sure it is less than 100 ppm
    if deconvolute_tolerance > 100:
        st.warning('Tolerance is greater than 100 ppm')
        st.stop()

if tolerance_type == 'da':
    # sure it is less than 0.1
    if deconvolute_tolerance > 0.1:
        st.warning('Tolerance is greater than 0.1 Da')
        st.stop()
#deconvolute_min_intensity = st.number_input('Deconvolution min intensity', value=0.1)

# Validate and parse input
lines = [line.strip() for line in spectra.split('\n') if line.strip()]
try:
    mz_array = [float(line.split()[0]) for line in lines]
    intensity_array = [float(line.split()[1]) for line in lines]
except (ValueError, IndexError):
    st.error("Invalid input format. Please ensure each line contains 'mz intensity' separated by a space.")
    st.stop()

peaks = [(mz, intensity) for mz, intensity in zip(mz_array, intensity_array)]

# Filter by minimum intensity and sort
peaks = list(filter(lambda x: x[1] >= min_intensity, peaks))
if len(peaks) == 0:
    st.stop()

if len(peaks) > 10_000:
    st.warning('Too many peaks. Please reduce the number of peaks')
    st.stop()

# sort peaks by mz
peaks = sorted(peaks, key=lambda x: x[0])


@st.cache_data(show_spinner=False)
def run_deconvolution(peaks, charge_range, tolerance, tolerance_type, min_left_decrease, min_right_decrease):
    return deconvolute(
        peaks,
        charge_range=charge_range,
        tolerance=tolerance,
        tolerance_type=tolerance_type,
        min_left_intensity_decrease=min_left_decrease,
        min_right_intensity_decrease=min_right_decrease
    )


dpeaks = run_deconvolution(
    peaks,
    (min_charge, max_charge),
    deconvolute_tolerance,
    tolerance_type,
    min_left_intensity_decrease,
    min_right_intensity_decrease
)

peaks_df = pd.DataFrame()
peaks_df['monoisotopic_mz'] = [dpeak.monoisotopic_peak.mz for dpeak in dpeaks]
peaks_df['monoisotopic_intensity'] = [dpeak.monoisotopic_peak.intensity for dpeak in dpeaks]
peaks_df['largest_mz'] = [dpeak.largest_peak.mz for dpeak in dpeaks]
peaks_df['largest_intensity'] = [dpeak.largest_peak.intensity for dpeak in dpeaks]
peaks_df['charge'] = [peak.charge for peak in dpeaks]
peaks_df['envelope_start_mz'] = [peak.mz_window[0] for peak in dpeaks]
peaks_df['envelope_end_mz'] = [peak.mz_window[1] for peak in dpeaks]
peaks_df['peaks'] = [len(peak.peaks) for peak in dpeaks]
peaks_df['neutral_mass'] = [peak.neutral_mass for peak in dpeaks]

# Handle missing charges by filling with 0 or a placeholder value
peaks_df['charge'] = peaks_df['charge'].fillna(0)  # Use 0 or a designated value for None

fig = go.Figure()

fig.add_trace(go.Scatter(
    x=sum([[mz, mz, None] for mz, intensity in peaks], []),
    y=sum([[0, intensity, None] for mz, intensity in peaks], []),
    mode='lines',
    name='Raw Spectrum',
    line=dict(color='gray'),
    opacity=0.3,
))

# Normalize charge states to map to a colormap
unique_charges = np.unique(peaks_df['charge'])
norm = mcolors.Normalize(vmin=min(unique_charges), vmax=max(unique_charges))
colormap = plt.get_cmap('viridis')  # Choose a colormap ('viridis', 'plasma', etc.)

# Assign a distinct color for None or 0 (unassigned charges)
special_color = '#B0B0B0'  # Light gray for None/NaN charges

# Generate colors for each charge state
charge_colors = {charge: mcolors.rgb2hex(colormap(norm(charge))) for charge in unique_charges if charge != 0}
charge_colors[0] = special_color  # Explicitly assign color for charge = 0 (was None)

# Map colors to charge states
colors = [charge_colors.get(charge, special_color) for charge in peaks_df['charge']]

# Plot each charge state as a separate trace
for charge in unique_charges:
    mask = peaks_df['charge'] == charge
    mz_values = peaks_df['monoisotopic_mz'][mask]
    intensity_values = peaks_df['monoisotopic_intensity'][mask]

    fig.add_trace(go.Scatter(
        x=sum([[mz, mz, None] for mz in mz_values], []),
        y=sum([[0, intensity, None] for intensity in intensity_values], []),
        mode='lines',
        name=f'+{int(charge)} Peaks' if charge > 0 else '+? Peaks',
        line=dict(color=charge_colors[charge], width=2),
        hovertext=sum([
            [
                f'Charge: {charge}<br>m/z: {mz}<br>Intensity: {intensity}<br>Window Start: {start}<br>Window End: {end}',
                f'Charge: {charge}<br>m/z: {mz}<br>Intensity: {intensity}<br>Window Start: {start}<br>Window End: {end}',
                None
            ]
            for mz, intensity, start, end in zip(
                mz_values,
                intensity_values,
                peaks_df['envelope_start_mz'][mask],
                peaks_df['envelope_end_mz'][mask]
            )
        ], [])
    ))

fig.update_layout(
    title='Deconvoluted Mass Spectrum',
    width=800,
    height=500,
    legend=dict(
        x=0.02,  # Horizontal position (0 = left, 1 = right)
        y=0.98,  # Vertical position (0 = bottom, 1 = top)
        xanchor='left',
        yanchor='top',
        bgcolor='rgba(255,255,255,0.5)',  # Semi-transparent background
    )

)

# show peaks in original spectrum (count) and then in the deconvoluted spectrum

c1, c2, c3 = st.columns(3)
c1.metric('Original Peaks', len(peaks),
          help='Number of peaks in the original spectrum (After filtering by min intensity)')
c2.metric('Deconvoluted Peaks', len(peaks_df), help='Number of peaks in the deconvoluted spectrum')

# peaks with valid charge
valid_peaks = peaks_df[peaks_df['charge'] != 0]
c3.metric('Valid Peaks', len(valid_peaks), help='Number of peaks with a valid charge state')

@st.fragment
def display():
    min_peak_mz = min([peak[0] for peak in peaks])
    max_peak_mz = max([peak[0] for peak in peaks])

    filter_min_mz, filter_max_mz = st.slider('Select a range of m/z', min_peak_mz, max_peak_mz,
                                             (min_peak_mz, max_peak_mz))

    peaks_df_filtered = peaks_df[
        (peaks_df['monoisotopic_mz'] >= filter_min_mz) & (peaks_df['monoisotopic_mz'] <= filter_max_mz)]
    max_intensity = max(peaks_df_filtered['monoisotopic_intensity'])
    # zoom into plotly plot and filter df

    fig.update_layout(xaxis=dict(range=[filter_min_mz, filter_max_mz]), yaxis=dict(range=[0, max_intensity]))

    st.plotly_chart(fig, use_container_width=True)

    st.caption('Peaks Table')

    # set 0 chareg state to None:

    peaks_df_filtered['charge'] = peaks_df_filtered['charge'].replace(0, None)

    st.dataframe(peaks_df_filtered, use_container_width=True, hide_index=True)


display()
