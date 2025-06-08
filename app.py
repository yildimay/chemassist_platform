import streamlit as st

# Page config
st.set_page_config(
    page_title="Chem Assist",
    layout="wide",
    initial_sidebar_state="expanded"
)

# Imports
from gaussian.gaussian_fixer_ui import gaussian_fixer_ui
from gaussian.gaussian_input_ui import gaussian_input_ui
from smiles_to_mol import smiles_ui
from env_hazard_panel import main as hazard_ui
from modules.gromacs_mdp_generator import (
    generate_em_mdp,
    generate_nvt_mdp,
    generate_npt_mdp,
    generate_md_mdp
)
#from modules.gromacs_error_fixer_ai import gromacs_error_fixer_ai_ui

# Sidebar module selection
selected_module = st.sidebar.radio(
    "Select a module:",
    [
        "Gaussian",
        "GROMACS",
        "Molecule Builder",
        "Environmental Hazard Panel",
        "VASP (coming soon)"
    ]
)

# ========== MODULE LOGIC ==========

# Gaussian module
if selected_module == "Gaussian":
    gaussian_tool = st.sidebar.radio("Select Tool", ["Error Fixer", "Input Creator"])
    if gaussian_tool == "Error Fixer":
        gaussian_fixer_ui()
    elif gaussian_tool == "Input Creator":
        gaussian_input_ui()

# Molecule builder
elif selected_module == "Molecule Builder":
    smiles_ui()

# Environmental hazard panel
elif selected_module == "Environmental Hazard Panel":
    hazard_ui()

# VASP placeholder
elif selected_module == "VASP (coming soon)":
    st.markdown("### VASP module is under development.")

# GROMACS module
elif selected_module == "GROMACS":
    st.header("🧬 GROMACS Tools")

    gromacs_tool = st.sidebar.radio("Select Tool", ["MDP Generator", "Error Fixer", "Simulation Builder"])

    if gromacs_tool == "MDP Generator":
        gsub_tool = st.selectbox("Select Stage", ["Energy Minimization", "NVT", "NPT", "Production"])

        if gsub_tool == "Energy Minimization":
            generate_em_mdp()
        elif gsub_tool == "NVT":
            generate_nvt_mdp()
        elif gsub_tool == "NPT":
            generate_npt_mdp()
        elif gsub_tool == "Production":
            generate_md_mdp()

    elif gromacs_tool == "Error Fixer":
            gromacs_error_fixer_ai_ui()

    elif gromacs_tool == "Simulation Builder":
        st.markdown("### Simulation builder coming soon.")
