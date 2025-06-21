# -*- coding: utf-8 -*-
"""smiles_to_mol.py
Streamlit app that:
  • Visualises a SMILES string in 2D, 3D, or XYZ form.
  • Extracts chemical names from an uploaded research‑article PDF and lets the
    user fetch the selected molecule's structure from PubChem.

Add to requirements.txt (or pyproject.toml):
    streamlit>=1.34.0
    rdkit-pypi>=2023.9.4
    pymupdf==1.24.2             # for PDF parsing (fitz)
    chemdataextractor==1.5.0    # optional, better name‑NER
    pubchempy>=1.0.4

Usage:
    streamlit run smiles_to_mol.py
"""

###############################################################################
# Imports
###############################################################################
import io
import re
import requests

import streamlit as st
import streamlit.components.v1 as components

from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdmolfiles import MolToMolBlock, MolToXYZBlock

# Optional dependencies -------------------------------------------------------
try:
    import fitz  # PyMuPDF
except ImportError:  # pragma: no cover
    fitz = None

try:
    from chemdataextractor import Document as CDE_Document
except ImportError:  # pragma: no cover
    CDE_Document = None

###############################################################################
# Helper functions
###############################################################################

def render_2d_molecule(mol, size=(350, 350)) -> None:
    """Render RDKit Mol as PNG in Streamlit."""
    img = Draw.MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    st.image(buf.getvalue(), use_column_width=False)


def render_3d_molecule(mol) -> None:
    """Render RDKit Mol in an interactive 3Dmol.js viewer."""
    mol_block = MolToMolBlock(mol).replace("\n", "\\n")

    html = f"""
    <div id=\"viewer\" style=\"width:100%;height:400px\"></div>
    <script src=\"https://3Dmol.org/build/3Dmol.js\"></script>
    <script>
      const viewer = $3Dmol.createViewer("viewer", {{backgroundColor:"white"}});
      viewer.addModel(`{mol_block}`, "mol");
      viewer.setStyle({{}}, {{stick:{{}}}});
      viewer.zoomTo();
      viewer.render();
    </script>
    """
    components.html(html, height=420, scrolling=False)

###############################################################################
# PDF → molecule‑name extraction utilities
###############################################################################

def extract_molecule_names_from_text(text: str):
    """Return a deduplicated list of chemical names detected in raw text."""
    names = set()

    # 1) Named‑entity recognition via ChemDataExtractor, if available
    if CDE_Document is not None:
        try:
            doc = CDE_Document(text)
            names.update(c.text for c in doc.cems)
        except Exception:
            pass  # fallback to regex below

    # 2) Regex heuristic as fallback / complement
    regex_hits = re.findall(r"\b([A-Z][a-z][A-Za-z0-9\-]{2,})\b", text)
    names.update(regex_hits)

    # 3) Filter obvious non‑chemical noise
    noise = re.compile(
        r"^(Table|Figure|Received|Revised|Accepted|Science|Journal|Available|"
        r"Vertex|Explorer|Changsha|Lianyungang|Shelx)$",
        re.I,
    )
    return sorted({n for n in names if not noise.match(n)}, key=str.lower)


def extract_molecule_names_from_pdf(file_like):
    """Read a PDF (file‑like object) and harvest candidate molecule names."""
    if fitz is None:
        st.error("PyMuPDF (pymupdf) is not installed; cannot read PDF files.")
        return []

    text = ""
    try:
        with fitz.open(stream=file_like.read(), filetype="pdf") as doc:
            for page in doc:
                text += page.get_text()
    except Exception as err:
        st.error(f"PDF reading error: {err}")
        return []

    return extract_molecule_names_from_text(text)


def fetch_smiles_from_pubchem(name: str):
    """Return the first PubChem Isomeric SMILES for *name*, or None."""
    try:
        url = (
            "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/"
            f"{name}/property/IsomericSMILES/TXT"
        )
        resp = requests.get(url, timeout=10)
        if resp.status_code == 200 and resp.text.strip():
            return resp.text.strip()
    except Exception:
        pass
    return None

###############################################################################
# Streamlit UI
###############################################################################

def main():
    st.title("SMILES → 3D / 2D Viewer & PDF Molecule Extractor")

    # -------------------- SMILES viewer sidebar --------------------
    st.sidebar.header("SMILES Input")
    smiles_in = st.sidebar.text_input("SMILES:", "C1=CC=CC=C1")
    disp_mode = st.sidebar.radio("Display", ["2D", "3D", "XYZ"])

    if st.sidebar.button("Render"):
        if not smiles_in:
            st.error("Please enter a SMILES string.")
        else:
            mol = Chem.MolFromSmiles(smiles_in)
            if mol is None:
                st.error("❌ Invalid SMILES.")
            else:
                if disp_mode != "2D":
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    try:
                        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
                    except Exception:
                        pass  # ignore minor optimisation failures

                if disp_mode == "2D":
                    render_2d_molecule(mol)
                elif disp_mode == "3D":
                    render_3d_molecule(mol)
                else:  # XYZ coordinates
                    st.code(MolToXYZBlock(mol), language="xyz")

    # ----------------------- PDF extraction -----------------------
    st.markdown("---")
    st.header("Extract Molecules from PDF")

    pdf = st.file_uploader("Upload a research‑article PDF", type=["pdf"])
    if pdf and st.button("Scan PDF"):
        names = extract_molecule_names_from_pdf(pdf)
        if names:
            st.success(f"Found {len(names)} candidate names.")
            chosen = st.selectbox("Select a molecule", names)

            if st.button("Fetch & Render from PubChem"):
                smi = fetch_smiles_from_pubchem(chosen)
                if smi:
                    st.info(f"SMILES from PubChem: {smi}")
                    mol = Chem.MolFromSmiles(smi)
                    mol = Chem.AddHs(mol)
                    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                    AllChem.UFFOptimizeMolecule(mol)
                    render_3d_molecule(mol)
                    st.download_button(
                        "Download MOL",
                        MolToMolBlock(mol),
                        file_name=f"{chosen}.mol",
                    )
                else:
                    st.error("No PubChem entry found for this name.")
        else:
            st.warning("No molecule‑like names detected in the PDF.")


if __name__ == "__main__":
    main()
