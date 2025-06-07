import streamlit as st
st.set_page_config(
    page_title="Chem Assist",
    layout="wide",
    initial_sidebar_state="expanded"
)
from gaussian.gaussian_fixer_ui import gaussian_fixer_ui
from gaussian.gaussian_input_ui import gaussian_input_ui
from smiles_to_mol import smiles_ui
from env_hazard_panel import main as hazard_ui   # import your panel

# … any other imports …

selected_module = st.sidebar.radio(
    "Select a module:",
    [
        "Gaussian",
        "Molecule Builder",
        "Environmental Hazard Panel",
        "VASP (coming soon)",
        "GROMACS (coming soon)"
    ]
)

if selected_module == "Gaussian":
    gaussian_tool = st.sidebar.radio("Gaussian Tool", ["Error Fixer", "Input Creator"])
    if gaussian_tool == "Error Fixer":
        gaussian_fixer_ui()
    elif gaussian_tool == "Input Creator":
        gaussian_input_ui()

elif selected_module == "Molecule Builder":
    smiles_ui()

elif selected_module == "Environmental Hazard Panel":
    # Directly call your panel's main() function here:
    hazard_ui()

elif selected_module == "VASP (coming soon)":
    st.markdown("### VASP module is under development.")

elif selected_module == "GROMACS (coming soon)":
    st.markdown("### GROMACS module is under development.")


from modules.gromacs_mdp_generator import generate_em_mdp

if selected_software == "GROMACS":
    st.header("GROMACS Tools")

    selected_tool = st.selectbox("Select Tool", ["MDP Generator", "Error Fixer", "Simulation Builder"])

    if selected_tool == "MDP Generator":
        sim_stage = st.selectbox("Select Stage", ["Energy Minimization", "NVT", "NPT", "Production"])

        if sim_stage == "Energy Minimization":
            generate_em_mdp()
        # Next up: add similar generators for NVT, NPT, etc.
