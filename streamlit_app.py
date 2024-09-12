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
import requests

# preloaded zip
ZIP_URL = 'https://raw.githubusercontent.com/praneelshah07/MIT-Project/main/ASM_Vapor_Spectra.csv.zip'

def load_data_from_zip(zip_url):
    try:
        response = requests.get(zip_url)
        if response.status_code != 200:
            st.error("Error downloading the zip file from the server.")
            return None

        zip_file = zipfile.ZipFile(io.BytesIO(response.content))
        csv_file = None
        for file in zip_file.namelist():
            if file.endswith('.csv'):
                csv_file = file
                break

        if csv_file:
            with zip_file.open(csv_file) as f:
                df = pd.read_csv(f)
            return df
        else:
            st.error("No CSV file found inside the ZIP from the server.")
            return None
    except Exception as e:
        st.error(f"Error extracting CSV from ZIP: {e}")
        return None

# set up app
st.title("Spectra Visualization App")

# load data from zip
data = load_data_from_zip(ZIP_URL)
if data is not None:
    st.write("Using preloaded data from GitHub zip file.")

# file uploader
uploaded_file = st.file_uploader("If you would like to enter another dataset, insert it here", type=["csv", "zip"])

if uploaded_file is not None:
    if uploaded_file.name.endswith('.zip'):
        with zipfile.ZipFile(uploaded_file, 'r') as z:
            file_list = z.namelist()
            csv_file = None
            for file in file_list:
                if file.endswith('.csv'):
                    csv_file = file
                    break
            if csv_file:
                with z.open(csv_file) as f:
                    data = pd.read_csv(f)
                st.write(f"Extracted and reading: {csv_file}")
            else:
                st.error("No CSV file found inside the uploaded ZIP.")
    
    elif uploaded_file.name.endswith('.csv'):
        data = pd.read_csv(uploaded_file)
        st.write("Using uploaded CSV file data.")

if data is not None:
    data['Raw_Spectra_Intensity'] = data['Raw_Spectra_Intensity'].apply(json.loads)
    data['Raw_Spectra_Intensity'] = data['Raw_Spectra_Intensity'].apply(np.array)
    data['Normalized_Spectra_Intensity'] = data['Raw_Spectra_Intensity'].apply(lambda x: x / max(x))

    columns_to_display = ["Formula", "IUPAC chemical name", "SMILES", "Molecular Weight", "Boiling Point (oC)"]
    headerdata = data[columns_to_display]
    st.write(headerdata)

    unique_smiles = data['SMILES'].unique()
    selected_smiles = st.multiselect('Select molecules by SMILES to highlight:', unique_smiles)

    peak_finding_enabled = st.checkbox('Enable Peak Finding and Labeling', value=False)

    confirm_button = st.button('Confirm Selection and Start Plotting')

    if confirm_button:
        st.write("The code will take some time to run, please wait...")

        fig, ax = plt.subplots(figsize=(16, 6.5), dpi=100)
        wavenumber = np.arange(4000, 500, -1)
        wavelength = 10000 / wavenumber

        color_palette = ['r', 'g', 'b', 'c', 'm', 'y']
        random.shuffle(color_palette)

        target_spectra = {}
        for smiles, spectra in data[['SMILES', 'Normalized_Spectra_Intensity']].values:
            if smiles in selected_smiles:
                target_spectra[smiles] = spectra
            else:
                ax.fill_between(wavelength, 0, spectra, color="k", alpha=0.01)

        for i, smiles in enumerate(target_spectra):
            spectra = target_spectra[smiles]
            ax.fill_between(wavelength, 0, spectra, color=color_palette[i % len(color_palette)], 
                            alpha=0.5, label=f"{smiles}")
            
            if peak_finding_enabled:
                peaks, _ = find_peaks(spectra, height=0.05)
                for peak in peaks:
                    peak_wavelength = wavelength[peak]
                    peak_intensity = spectra[peak]
                    ax.text(peak_wavelength, peak_intensity + 0.05, f'{round(peak_wavelength, 1)}', 
                            fontsize=10, ha='center', color=color_palette[i % len(color_palette)])

        # customize plot
        ax.set_xscale('log')
        ax.set_xlim([2.5, 20])

        major_ticks = [3, 4, 5, 6, 7, 8, 9, 11, 12, 15, 20]
        ax.set_xticks(major_ticks)

        # number of label matches number of ticks
        ax.set_xticklabels([str(tick) for tick in major_ticks])

        ax.tick_params(direction="in",
            labelbottom=True, labeltop=False, labelleft=True, labelright=False,
            bottom=True, top=True, left=True, right=True)

        ax.set_xlabel("Wavelength ($\mu$m)", fontsize=22)
        ax.set_ylabel("Absorbance (Normalized to 1)", fontsize=22)

        if selected_smiles:
            ax.legend()

        st.pyplot(fig)

        # download button
        buf = io.BytesIO()
        fig.savefig(buf, format='png')
        buf.seek(0)
        st.download_button(label="Download Plot as PNG", data=buf, file_name="spectra_plot.png", mime="image/png")
