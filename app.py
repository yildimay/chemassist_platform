import streamlit as st
import requests
import difflib
import base64
from PIL import Image
import os
import py3Dmol
from rdkit import Chem
from rdkit.Chem import Draw, AllChem

st.set_page_config(page_title="ChemAssist Platform", layout="centered")

# API setup
GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
EXPLAIN_MODEL = "llama3-8b-8192"
FIX_MODEL = "llama3-70b-8192"
GROQ_API_KEY = os.environ.get("GROQ_API_KEY", "")

if "fix_prompt" not in st.session_state:
    st.session_state.fix_prompt = ""

st.sidebar.title("🧪 ChemAssist Tools")
selected_software = st.sidebar.selectbox("Which software are you using?", ["Select...", "Gaussian", "Modeling"])

if selected_software == "Modeling":
    st.title("🧪 Molecule Builder (SMILES → MOL)")
    st.markdown("Enter a SMILES string to generate and download a .mol file.")

    smiles_input = st.text_input("💬 Enter SMILES", value="CCO")

    if st.button("🛠 Convert to MOL"):
        try:
            mol = Chem.MolFromSmiles(smiles_input)
            mol_filename = "molecule.mol"
            Chem.MolToMolFile(mol, mol_filename)

            st.success("MOL file generated!")
            col1, col2 = st.columns([1, 1])

            with col1:
                st.download_button("💾 Download .mol file", open(mol_filename, "rb").read(), file_name=mol_filename)
                img = Draw.MolToImage(mol, size=(300, 300))
                st.image(img, caption="Molecule Preview")

            with col2:
                st.subheader("🧬 3D Molecule Viewer")
                try:
                    mol3d = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol3d, AllChem.ETKDG())
                    xyz_block = Chem.MolToXYZBlock(mol3d)
                    viewer = py3Dmol.view(width=400, height=300)
                    viewer.addModel(xyz_block, "xyz")
                    viewer.setStyle({"stick": {}})
                    viewer.zoomTo()
                    viewer.show()
                    st.components.v1.html(viewer._make_html(), height=300)
                except Exception as e:
                    st.warning("3D structure could not be generated.")
