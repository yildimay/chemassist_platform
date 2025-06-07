import streamlit as st
st.set_page_config(
    page_title="Chem Assist",
    layout="wide",
    initial_sidebar_state="expanded"
)
from gaussian.gaussian_fixer_ui import gaussian_fixer_ui
from gaussian.gaussian_input_ui import gaussian_input_ui
from smiles_to_mol import smiles_ui
from env_hazard_panel import main as hazard_ui
from modules.gromacs_mdp_generator import generate_em_mdp

# … any other imports …

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

if selected_module == "Gaussian":
    gaussian_tool = st.sidebar.radio("Select Tool", ["Error Fixer", "Input Creator"])
    if gaussian_tool == "Error Fixer":
        gaussian_fixer_ui()
    elif gaussian_tool == "Input Creator":
        gaussian_input_ui()

elif selected_module == "Molecule Builder":
    smiles_ui()

elif selected_module == "Environmental Hazard Panel":
    hazard_ui()

elif selected_module == "VASP (coming soon)":
    st.markdown("### VASP module is under development.")

elif selected_module == "GROMACS":
      gromacs_tool = st.sidebar.radio("Select Tool", ["MDP Generator", "Error Fixer", "Simulation Builder"])
    if gromacs_tool == "MDP Generator":
        gsub_tool = st.selectbox("Select Stage", ["Energy Minimization", "NVT", "NPT", "Production"])

            if gsub_tool == "Energy Minimization":
                generate_em_mdp()
            if gsub_tool == "NVT":
                st.markdown("### this tool is under development.")
            if gsub_tool == "NPT":
                st.markdown("### this tool is under development.")
            if gsub_tool == "Production":
                st.markdown("### this tool is under development.")



if selected_module == "GROMACS":
    st.header("GROMACS Tools")

    selected_tool = st.selectbox("Select Tool", ["MDP Generator", "Error Fixer", "Simulation Builder"])

   
        # Next up: add similar generators for NVT, NPT, etc.
