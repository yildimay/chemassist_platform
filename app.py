import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import streamlit as st
from gaussian.gaussian_input_ui import gaussian_input_ui
from gaussian.gaussian_fixer_ui import gaussian_fixer_ui
from smiles_to_mol import smiles_ui

# Set Streamlit page configuration
st.set_page_config(page_title="ChemAssist Platform", layout="wide")

# Sidebar module selection
st.sidebar.title("ChemAssist Modules")
selected_module = st.sidebar.selectbox("Choose a tool:", ["Gaussian", "Molecule Builder", "VASP (coming soon)", "GROMACS (coming soon)", "Environmental Hazard Panel"])

if selected_module == "Gaussian":
    gaussian_tool = st.sidebar.radio("Gaussian Tool", ["Error Fixer", "Input Creator"])
    if gaussian_tool == "Error Fixer":
        gaussian_fixer_ui()
    elif gaussian_tool == "Input Creator":
        gaussian_input_ui()

elif selected_module == "Molecule Builder":
    smiles_ui()

elif selected_module == "VASP (coming soon)":
    st.markdown("### VASP module is under development.")

elif selected_module == "GROMACS (coming soon)":
    st.markdown("### GROMACS module is under development.")
elif choice == "Environmental Hazard Panel":
    env_hazard_panel.main()
