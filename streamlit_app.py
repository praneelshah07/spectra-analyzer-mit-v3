import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.style import use
use('fast')
import json
import random
import zipfile
import io
from scipy.signal import find_peaks
from scipy.spatial.distance import pdist, squareform
from scipy.cluster.hierarchy import linkage
import requests

# Preload ZIP file from GitHub and extract CSV inside it
ZIP_URL = 'https://raw.githubusercontent.com/praneelshah07/MIT-Project/main/ASM_Vapor_Spectra.csv.zip'

def load_data_from_zip(zip_url):
    try:
        # Download the zip file
        response = requests.get(zip_url)
        if response.status_code != 200:
            st.error("Error downloading the zip file from the server.")
            return None

        # Open the downloaded zip file in memory
        zip_file = zipfile.ZipFile(io.BytesIO(response.content))

        # Find the first CSV file inside the zip
        csv_file = None
        for file in zip_file.namelist():
            if file.endswith('.csv'):
                csv_file = file
                break

        if csv_file:
            # Read the CSV file from the zip file
            with zip_file.open(csv_file) as f:
                df = pd.read_csv(f)
            return df
        else:
            st.error("No CSV file found inside the ZIP from the server.")
            return None
    except Exception as e:
        st.error(f"Error extracting CSV from ZIP: {e}")
        return None

# Set up the Streamlit app
st.title("Spectra Visualization App")

# Try to load data from the zip file on GitHub
data = load_data_from_zip(ZIP_URL)
if data is not None:
    st.write("Using preloaded data from GitHub zip file.")

# File uploader for manual CSV or ZIP file
uploaded_file = st.file_uploader("If you would like to enter another dataset, insert it here", type=["csv", "zip"])

# If a file is uploaded, it will override the preloaded data
if uploaded_file is not None:
    # Handle ZIP file upload
    if uploaded_file.name.endswith('.zip'):
        # Open the uploaded zip file
        with zipfile.ZipFile(uploaded_file, 'r') as z:
            # Get a list of all files in the zip
            file_list = z.namelist()

            # Find the first CSV file in the zip
            csv_file = None
            for file in file_list:
                if file.endswith('.csv'):
                    csv_file = file
                    break

            # If a CSV file is found, read it
            if csv_file:
                with z.open(csv_file) as f:
                    data = pd.read_csv(f)
                st.write(f"Extracted and reading: {csv_file}")
            else:
                st.error("No CSV file found inside the uploaded ZIP.")
    
    # Handle CSV file upload
    elif uploaded_file.name.endswith('.csv'):
        data = pd.read_csv(uploaded_file)
        st.write("Using uploaded CSV file data.")

# Check if data exists (either preloaded or uploaded)
if data is not None:
    # Convert JSON string to lists and normalize the spectra
    data['Raw_Spectra_Intensity'] = data['Raw_Spectra_Intensity'].apply(json.loads)
    data['Raw_Spectra_Intensity'] = data['Raw_Spectra_Intensity'].apply(np.array)
    data['Normalized_Spectra_Intensity'] = data['Raw_Spectra_Intensity'].apply(lambda x: x / max(x))

    # Preview the dataframe to ensure data is loaded correctly
    st.write(data.head())






    

    # Select SMILES for molecules you want to highlight
    unique_smiles = data['SMILES'].unique()
    selected_smiles = st.multiselect('Select molecules by SMILES to highlight:', unique_smiles)

    # Add a checkbox to enable or disable peak finding
    peak_finding_enabled = st.checkbox('Enable Peak Finding and Labeling', value=False)

    # Initialize plot with adjusted DPI for better resolution
    fig, ax = plt.subplots(figsize=(16, 6.5), dpi=100)

    # Calculate wavelength
    wavenumber = np.arange(4000, 500, -1)
    wavelength = 10000 / wavenumber  # in microns

    # Color palette for highlighted spectra
    color_palette = ['r', 'g', 'b', 'c', 'm', 'y']  # Add more colors if needed
    random.shuffle(color_palette)  # Shuffle colors to randomize highlights

    # Plot the spectra
    target_spectra = {}  # Store selected spectra for highlighting
    for smiles, spectra in data[['SMILES', 'Normalized_Spectra_Intensity']].values:
        if smiles in selected_smiles:
            target_spectra[smiles] = spectra  # Store for highlighting later
        else:
            ax.fill_between(wavelength, 0, spectra, color="k", alpha=0.01)  # Plot all other spectra

    # Highlight the selected spectra with different colors and annotate peaks if enabled
    for i, smiles in enumerate(target_spectra):
        spectra = target_spectra[smiles]
        ax.fill_between(wavelength, 0, spectra, color=color_palette[i % len(color_palette)], 
                        alpha=0.5, label=f"{smiles}")
        
        # If peak finding is enabled, find and annotate peaks
        if peak_finding_enabled:
            peaks, _ = find_peaks(spectra, height=0.05)  # Adjust height parameter for sensitivity
            for peak in peaks:
                peak_wavelength = wavelength[peak]
                peak_intensity = spectra[peak]
                ax.text(peak_wavelength, peak_intensity + 0.05, f'{round(peak_wavelength, 1)}', 
                        fontsize=10, ha='center', color=color_palette[i % len(color_palette)])

    # Customize plot axes and ticks
    ax.set_xscale('log')
    ax.set_xlim([2.5, 20])

    major_ticks = [3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 20]
    ax.set_xticks(major_ticks, minor=False)

    ax.tick_params(axis="x", labelsize=16)
    ax.tick_params(axis="y", labelsize=16)

    ax.set_xticks([3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 15, 20])
    ax.set_xticklabels(["3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "15", "20"])

    ax.tick_params(direction="in",
                   labelbottom=True, labeltop=False, labelleft=True, labelright=False,
                   bottom=True, top=True, left=True, right=True)

    ax.set_xlabel("Wavelength ($\mu$m)", fontsize=22)
    ax.set_ylabel("Absorbance (Normalized to 1)", fontsize=22)

    # Show legend
    if selected_smiles:
        ax.legend()

    # Display the plot in Streamlit
    st.pyplot(fig)
