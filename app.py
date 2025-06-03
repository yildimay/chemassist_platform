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
selected_module = st.sidebar.selectbox("Choose a tool:", ["Gaussian", "Molecule Builder", "VASP (coming soon)", "GROMACS (coming soon)"])

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
    
# env_hazard_panel.py

import streamlit as st
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

def draw_molecule(mol, width=300, height=200):
    """
    Draw an RDKit molecule as SVG and render it in Streamlit.
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')
    st.write(f'<div>{svg}</div>', unsafe_allow_html=True)

def calculate_physicochemical_properties(smiles):
    """
    Stub function for later: compute molecular weight, logP, TPSA, etc.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {}
    from rdkit.Chem import Descriptors
    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "logP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        # Add more properties later...
    }

def main():
    st.set_page_config(page_title="Chem Assist", layout="wide")
    st.sidebar.title("Tools")
    choice = st.sidebar.radio(
        "Select a module:",
        ["Gaussian Fixer", "DFT Helper", "Environmental Hazard Panel", "Other Modules…"]
    )

    if choice == "Environmental Hazard Panel":
        st.header("Environmental Hazard Panel")
        st.write(
            "Upload your molecule below (SMILES or .mol/.sdf). "
            "We’ll visualize it and, in the next step, compute key properties and risk indicators."
        )

        # --- 1. Molecule Input ---
        smiles_input = st.text_input("Enter molecule SMILES string:", "")
        file_input = st.file_uploader("Or upload a .mol / .sdf file:", type=["mol", "sdf"])

        mol_obj = None
        if smiles_input:
            mol_obj = Chem.MolFromSmiles(smiles_input)
            if mol_obj is None:
                st.error("Invalid SMILES string. Please check your input.")
        elif file_input is not None:
            try:
                block = file_input.read().decode("utf-8")
                mol_obj = Chem.MolFromMolBlock(block)
                if mol_obj is None:
                    st.error("Could not parse the uploaded file. Make sure it’s a valid .mol or .sdf.")
            except Exception as e:
                st.error(f"Error reading file: {e}")

        # --- 2. Display Molecule Structure ---
        if mol_obj is not None:
            st.subheader("Molecule Structure")
            draw_molecule(mol_obj)

            # --- 3. Placeholder: Basic Properties ---
            st.subheader("Physicochemical Properties (preliminary)")
            props = calculate_physicochemical_properties(smiles_input if smiles_input else "")
            if props:
                st.table(props)
            else:
                st.write("Properties will appear here once a valid molecule is provided.")

            # --- 4. Placeholder: Environmental Risk Calculations ---
            st.subheader("Environmental Risk Indicators")
            st.write(
                "⚙️ Calculation modules (ECOSAR, BCF, LC₅₀, etc.) will go here. "
                "For now, this section is under construction."
            )

            # Example placeholder graphics area—fill in later
            st.info("Graphical risk profiles will be plotted here in a future iteration.")

            # --- 5. Placeholder: Risk Score ---
            st.subheader("Overall Risk Score")
            st.write("Risk score computation coming soon…")

            # (Later: implement “Download Report” button)
            st.download_button(
                label="Download Full Report (PDF/JSON)",
                data=None,
                file_name="env_risk_report.pdf",
                mime="application/pdf",
                disabled=True
            )

    else:
        # Handle other modules or show a placeholder
        st.write(f"You selected: **{choice}**")
        st.write("This module is not implemented yet in this prototype.")

if __name__ == "__main__":
    main()
