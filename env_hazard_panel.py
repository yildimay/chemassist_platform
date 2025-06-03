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
        # …add more later…
    }

def main():
    st.header("Environmental Hazard Panel")
    st.write(
        "Upload a molecule (SMILES or .mol/.sdf), see its 2D structure, "
        "and (soon) compute environmental risk metrics."
    )

    # 1) Molecule input
    smiles_input = st.text_input("Enter SMILES string:", "")
    file_input = st.file_uploader("Or upload a .mol/.sdf file:", type=["mol", "sdf"])

    mol_obj = None
    if smiles_input:
        mol_obj = Chem.MolFromSmiles(smiles_input)
        if mol_obj is None:
            st.error("Invalid SMILES. Please correct it.")
    elif file_input is not None:
        try:
            block = file_input.read().decode("utf-8")
            mol_obj = Chem.MolFromMolBlock(block)
            if mol_obj is None:
                st.error("Can't parse this file. Make sure it's a valid .mol/.sdf.")
        except Exception as e:
            st.error(f"Error reading file: {e}")

    # 2) Show structure
    if mol_obj:
        st.subheader("Molecule Structure")
        draw_molecule(mol_obj)

        # 3) Show basic RDKit properties
        st.subheader("Physicochemical Properties (demo)")
        props = calculate_physicochemical_properties(smiles_input if smiles_input else "")
        if props:
            st.table(props)
        else:
            st.write("Will show properties here once molecule is valid.")

        # 4) Placeholder: Environmental risk calculations
        st.subheader("Environmental Risk Indicators")
        st.write("⚙️ Advanced QSAR/EPI Suite/ECOSAR modules will appear here (coming soon).")
        st.info("Graphs & scores will be rendered here in a future update.")

        # 5) Placeholder: Overall risk score
        st.subheader("Overall Risk Score")
        st.write("Computation of a final risk score is under construction.")

        # Download‐button stub
        st.download_button(
            label="Download Full Report (PDF/JSON)",
            data=None,
            file_name="env_risk_report.pdf",
            mime="application/pdf",
            disabled=True
        )
    else:
        st.write("Awaiting a valid molecule to analyze…")
