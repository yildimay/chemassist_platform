# -*- coding: utf-8 -*-
"""
smiles_to_mol.py
Streamlit uygulaması
  • SMILES metnini 2 D / 3 D / XYZ olarak gösterir
  • PDF yüklenirse gerçek kimyasal isimleri çıkarıp listeler
     (saf-Python RegEx + PubChem doğrulaması – Python 3.11’de derleme yok)

requirements.txt  ➜  yalnızca bunlar:
    streamlit>=1.34.0
    rdkit-pypi>=2023.9.4
    pymupdf==1.24.2          # fitz
    pubchempy==1.0.4
"""
import re

RE_CANDIDATE = re.compile(r'\b([A-Z][a-z]+[a-z0-9() \-]{2,})\b')
import streamlit as st
import io
import re
from typing import Iterable, List

import requests
from rdkit import Chem
from rdkit.Chem import AllChem, MolToMolBlock, Draw
from rdkit.Chem.rdmolfiles import MolToXYZBlock
import streamlit.components.v1 as components

import fitz                 # PyMuPDF
import pubchempy as pcp      # PubChem API wrapper
# ----------------------------------------------------------------------
# Helper • PDF → kimyasal isimler (saf-Python, 3.11-uyumlu)
STOP_WORDS = {
    "table", "figure", "online", "revised", "received",
    "accepted", "science", "explorer", "vertex"
}

def smiles_ui():
    """Backward-compatibility alias — optional."""
    main()


def _extract_text(pdf_file: io.BufferedReader | io.BytesIO) -> str:
    "Return raw text concatenated from all pages"
    with fitz.open(stream=pdf_file.read(), filetype="pdf") as doc:
        return "\n".join(page.get_text("text") for page in doc)


def _is_real_molecule(name: str) -> bool:
    "True if PubChem has at least one record for this common name"
    try:
        return bool(pcp.get_compounds(name, "name", record_type="3d"))
    except Exception:
        # ağ hatası vs. durumunda false deme – mümkünse tut
        return False


def extract_molecules_from_pdf(pdf_file: io.BufferedReader | io.BytesIO) -> List[str]:
    """Return list of validated chemical names from the given PDF file-like object."""
    text = _extract_text(pdf_file)
    candidates = {m.group(1) for m in RE_CANDIDATE.finditer(text)}
    # kaba kara-liste & uzunluk filtresi
    filtered = {
        c for c in candidates
        if c.lower() not in STOP_WORDS and len(c) <= 30
    }
    # PubChem onayı
    molecules = [n for n in filtered if _is_real_molecule(n)]
    return sorted(set(molecules), key=str.lower)
# ----------------------------------------------------------------------
# RDKit görselleştirme yardımcıları
def render_2d_molecule(mol):
    img = Draw.MolToImage(mol, size=(300, 300))
    st.image(img, use_column_width=False)


def render_3d_molecule(mol):
    mol_block = MolToMolBlock(mol)
    html = f"""
        <div id="viewer" style="width:100%;height:400px;"></div>
        <script src="https://3Dmol.org/build/3Dmol.js"></script>
        <script>
          const viewer = $3Dmol.createViewer(
              "viewer", {{backgroundColor: "white"}}
          );
          viewer.addModel(`{mol_block}`, "mol");
          viewer.setStyle({{}}, {{stick:{{}}}});
          viewer.zoomTo();
          viewer.render();
        </script>
    """
    components.html(html, height=420, scrolling=False)
# ----------------------------------------------------------------------
# Streamlit UI
def main():
    st.title("SMILES ↔ Molecule & PDF Molecule Extractor")

    # ==== PDF Bölümü ====
    st.subheader("📄  PDF’den molekül isimleri çek")
    pdf_bytes = st.file_uploader("Araştırma makalesi PDF’i seç (isteğe bağlı)", type=["pdf"])
    if st.button("Molekülleri bul") and pdf_bytes:
        with st.spinner("PDF taranıyor…"):
            names = extract_molecules_from_pdf(pdf_bytes)
        if names:
            st.success(f"Bulunan moleküller ({len(names)}):")
            st.write(", ".join(names))
        else:
            st.warning("Geçerli kimyasal isim bulunamadı.")
    st.markdown("---")

    # ==== SMILES Görselleştirme ====
    st.subheader("🧪  SMILES ➜ MOL görselleştir")
    smi = st.text_input("SMILES:", "C1=CC=CC=C1")  # varsayılan benzen
    mode = st.radio("Görüntü Modu", ("2D", "3D", "XYZ koordinatları"))

    if st.button("Göster") and smi:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                st.error("❌ Geçersiz SMILES")
                return
            if mode != "2D":
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.UFFOptimizeMolecule(mol)

            if mode == "2D":
                render_2d_molecule(mol)
            elif mode == "3D":
                render_3d_molecule(mol)
            else:
                st.code(MolToXYZBlock(mol), language="xyz")
        except Exception as exc:
            st.error(f"⚠️ {mode} gösterimi başarısız: {exc}")

if __name__ == "__main__":
    main()



# 📄 PDF'den Molekül Seç ve Modelle - entegre modül
with st.expander("📄 PDF'den Molekül Seç ve Modelle"):
    example_molecules = ["Paracetamol", "Imatinib", "Quantum", "The", "Water", "Diclofenac"]
    stopwords = {"the", "water", "quantum", "access", "open", "clean", "all"}

    valid_mols = [m for m in example_molecules if m.lower() not in stopwords]
    st.write("**PDF'den çekilen moleküller:**")
    for mol in valid_mols:
        col1, col2 = st.columns([4, 1])
        col1.markdown(f"- **{mol}**")
        if col2.button("Modelle", key=f"btn_{mol}"):
            results = pcp.get_compounds(mol, 'name')
            if results:
                compound = results[0]
                st.session_state["selected_smiles"] = compound.isomeric_smiles
                st.session_state["selected_name"] = mol

    if "selected_smiles" in st.session_state:
        st.subheader(f"🧪 {st.session_state['selected_name']} Molekülü")
        mol = Chem.MolFromSmiles(st.session_state["selected_smiles"])
        if mol:
            AllChem.Compute2DCoords(mol)
            st.image(Draw.MolToImage(mol, size=(300, 300)))
            st.code(st.session_state["selected_smiles"], language="smiles")
        else:
            st.error("Molekül görselleştirilemedi.")
