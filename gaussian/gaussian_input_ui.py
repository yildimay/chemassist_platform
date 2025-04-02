import streamlit as st
import os
import requests
import json
from rdkit import Chem
from rdkit.Chem import AllChem
import fitz  # PyMuPDF


def truncate_structure_input(structure: str, max_chars: int = 3000) -> str:
    if len(structure) <= max_chars:
        return structure
    lines = structure.splitlines()
    truncated_lines = []
    total_len = 0
    for line in lines:
        if total_len + len(line) > max_chars:
            break
        truncated_lines.append(line)
        total_len += len(line) + 1
    return "\n".join(truncated_lines) + "\n... [structure truncated for safety]"


def smiles_to_xyz(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)
    conf = mol.GetConformer()
    xyz = f"{mol.GetNumAtoms()}\nGenerated from SMILES\n"
    for atom in mol.GetAtoms():
        pos = conf.GetAtomPosition(atom.GetIdx())
        xyz += f"{atom.GetSymbol()} {pos.x:.4f} {pos.y:.4f} {pos.z:.4f}\n"
    return xyz


def extract_text_from_pdf(uploaded_file) -> str:
    doc = fitz.open(stream=uploaded_file.read(), filetype="pdf")
    text = "\n".join([page.get_text() for page in doc])
    return text


def call_groq_for_gjf_from_paper(smiles_xyz, paper_text):
    GROQ_API_KEY = os.environ.get("GROQ_API_KEY")
    if not GROQ_API_KEY:
        raise ValueError("GROQ_API_KEY not set in environment variables.")

    prompt = f"""
Given the following molecule in XYZ format:
{truncate_structure_input(smiles_xyz)}

And this paper excerpt:
{paper_text[:4000]}

Extract the method, basis set, charge, multiplicity, and generate a Gaussian input (.gjf) file. Output only the .gjf content.
"""

    headers = {
        "Authorization": f"Bearer {GROQ_API_KEY}",
        "Content-Type": "application/json"
    }

    payload = {
        "model": "llama3-70b-8192",
        "messages": [{"role": "user", "content": prompt}],
        "temperature": 0.5
    }

    response = requests.post(
        "https://api.groq.com/openai/v1/chat/completions",
        headers=headers,
        data=json.dumps(payload)
    )
    response.raise_for_status()
    return response.json()["choices"][0]["message"]["content"]


def gaussian_input_ui():
    st.title("Gaussian Input from SMILES + Paper")

    st.markdown("### Upload SMILES and Reference Paper")
    smiles = st.text_input("Enter SMILES string")
    uploaded_pdf = st.file_uploader("Upload reference paper (PDF)", type="pdf")

    if st.button("Validate & Generate .gjf from Paper"):
        if not smiles or not uploaded_pdf:
            st.warning("Please provide both a SMILES and a paper PDF.")
            return

        try:
            with st.spinner("Processing paper and molecule..."):
                xyz_data = smiles_to_xyz(smiles)
                paper_text = extract_text_from_pdf(uploaded_pdf)
                gjf_output = call_groq_for_gjf_from_paper(xyz_data, paper_text)
                st.success(".gjf file generated from SMILES + paper!")
                st.code(gjf_output, language="gjf")
                st.download_button("Download .gjf File", data=gjf_output, file_name="generated_input.gjf")
        except Exception as e:
            st.error(f"Error: {e}")
