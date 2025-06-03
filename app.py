import streamlit as st
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
