# --- import section ---
import streamlit as st
import io                                     # <── eklendi
from rdkit import Chem
import requests 
from rdkit.Chem import AllChem, MolToMolBlock, Draw
from rdkit.Chem.rdmolfiles import MolToXYZBlock   # <── eklendi
import streamlit.components.v1 as components

# ---------- YARDIMCI GÖRSELLEŞTİRME FONKSİYONLARI ----------
def render_2d_molecule(mol):
    """RDKit 2D çizim çıktısını Streamlit'e gönderir."""
    img = Draw.MolToImage(mol, size=(300, 300))
    st.image(img, use_column_width=False)

def render_3d_molecule(mol):
    """py3Dmol + 3Dmol.js ile 3D gösterim yapar (HTML/JS gömülü)."""
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

# ---------- ANA UI ----------
def smiles_ui():
    st.title("SMILES ➜ MOL Dönüştürücü")
    smi = st.text_input("Enter SMILES string:", "C1=CC=CC=C1")

    display_mode = st.radio(
        "Display Mode", ["2D Viewer", "3D Viewer", "XYZ Coordinates"]
    )

    if st.button("Convert to MOL") and smi:
        try:
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                st.error("❌ Geçersiz SMILES girdisi")
                return

            # 3D gereken durumlarda koordinat ekle + opsiyonel UFF optimizasyon
            if display_mode != "2D Viewer":
                AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                AllChem.UFFOptimizeMolecule(mol)

            # Görselleştir
            if display_mode == "2D Viewer":
                render_2d_molecule(mol)
            elif display_mode == "3D Viewer":
                render_3d_molecule(mol)
            else:  # XYZ Coordinates
                xyz_block = Chem.MolToXYZBlock(mol)
                st.code(xyz_block, language="xyz")

        except Exception as e:
            st.error(f"⚠️ {display_mode} failed: {e}")

if __name__ == "__main__":
    smiles_ui()
    
def smiles_ui():
    st.title("🧬 SMILES to MOL Converter")

    st.markdown("### 📩 Input")
    smiles = st.text_input("Enter SMILES string:", "")
    display_mode = st.radio("Display Mode", ["2D Viewer", "3D Viewer", "XYZ Coordinates"])

    if st.button("Convert to MOL") and smiles:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                st.error("❌ Invalid SMILES string.")
                return

            mol = Chem.AddHs(mol)
            embedded = False
            for method in [AllChem.ETKDG(), AllChem.ETKDGv2(), AllChem.ETKDGv3()]:
                if AllChem.EmbedMolecule(mol, method) == 0:
                    embedded = True
                    break

            if not embedded:
                raise ValueError("3D embedding failed.")

            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=100)
            except:
                try:
                    props = AllChem.MMFFGetMoleculeProperties(mol)
                    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
                    ff.Minimize(maxIts=150)
                except:
                    raise ValueError("3D optimization with UFF and MMFF failed.")

            if display_mode == "XYZ Coordinates":
                xyz_block = MolToXYZBlock(mol)
                st.code(xyz_block, language="xyz")
            elif display_mode == "2D Viewer":
                mol_2d = Chem.MolFromSmiles(smiles)
                img = Draw.MolToImage(mol_2d, size=(400, 400))
                buf = io.BytesIO()
                img.save(buf, format="PNG")
                st.image(buf.getvalue(), caption="2D Structure")
            elif display_mode == "3D Viewer":
                render_3d_molecule(mol)

        except Exception as e:
            st.error(f"3D optimization failed: {e}")

    # --- PDF Extraction Section ---
    st.markdown("---")
    st.header("📄 Extract Molecule from Research Paper (PDF)")
    pdf_file = st.file_uploader("Upload a PDF research paper:", type=["pdf"])

    if pdf_file is not None:
        import fitz  # PyMuPDF
        import re
        text = ""
        try:
            with fitz.open(stream=pdf_file.read(), filetype="pdf") as doc:
                for page in doc:
                    text += page.get_text()
        except Exception as e:
            st.error(f"PDF reading error: {e}")

        def find_molecule_names(text):
            candidates = re.findall(r"\b([A-Z][a-z]+[a-z0-9\-() ]{2,})\b", text)
            filtered = [name for name in candidates if len(name.split()) < 4 and any(c.isdigit() for c in name)]
            return list(set(filtered))

        molecule_names = find_molecule_names(text)

        if molecule_names:
            st.success(f"Detected molecule names: {', '.join(molecule_names)}")
            selected_name = st.selectbox("Select a molecule:", molecule_names)

            if st.button("Fetch Structure from PubChem"):
                try:
                    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{selected_name}/property/IsomericSMILES/TXT"
                    response = requests.get(url)
                    if response.status_code == 200:
                        extracted_smiles = response.text.strip()
                        mol = Chem.MolFromSmiles(extracted_smiles)
                        mol = Chem.AddHs(mol)
                        AllChem.EmbedMolecule(mol, AllChem.ETKDG())
                        AllChem.UFFOptimizeMolecule(mol)
                        mol_block = MolToMolBlock(mol)
                        st.download_button("Download .mol file", mol_block, file_name=f"{selected_name}.mol")
                        st.info(f"SMILES used: {extracted_smiles}")
                        if display_mode == "3D Viewer":
                            render_3d_molecule(mol)
                    else:
                        st.error("Could not fetch SMILES from PubChem.")
                except Exception as e:
                    st.error(f"Error fetching structure: {e}")
        else:
            st.warning("No molecule names detected in the PDF.")
