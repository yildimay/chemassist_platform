import streamlit as st
from gaussian_fixer import gaussian_ui
from smiles_to_mol import smiles_ui
from rdkit.Chem import Draw

st.set_page_config(page_title="ChemAssist Platform", layout="wide")

# Sidebar setup
st.sidebar.title("🧪 ChemAssist Modules")
selected_module = st.sidebar.selectbox(
    "Choose a tool:", 
    ["Gaussian Error Fixer", "Molecule Builder"]
)

# Module routing
if selected_module == "Gaussian Error Fixer":
    gaussian_ui()

elif selected_module == "Molecule Builder":
    smiles_ui()

# Footer
st.sidebar.markdown("---")
st.sidebar.markdown("© 2023 ChemAssist Platform")
