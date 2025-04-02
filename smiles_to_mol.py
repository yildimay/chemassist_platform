import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
from rdkit.Chem.rdmolfiles import MolToMolBlock, MolToXYZBlock
import streamlit.components.v1 as components
from PIL import Image
import io


def render_3d_molecule(mol):
    try:
        mol_block = MolToMolBlock(mol)
        mol_block = mol_block.replace("\n", "\\n")  # preserve formatting for JS
        viewer_html = f"""
        <div style='height: 500px; width: 100%; position: relative;'>
        <script src='https://3Dmol.csb.pitt.edu/build/3Dmol-min.js'></script>
        <div id='viewer' style='height: 100%; width: 100%;'></div>
        <script>
            let viewer = $3Dmol.createViewer('viewer', {{backgroundColor: "white"}});
            viewer.addModel("{mol_block}", "mol");
            viewer.setStyle({{}}, {{stick: {{}}}});
            viewer.zoomTo();
            viewer.render();
        </script>
        </div>
        """
        components.html(viewer_html, height=500)
    except Exception as e:
        st.error(f"Failed to render 3D model: {e}")


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

            # Try multiple embedding methods
            embedded = False
            for method in [AllChem.ETKDG(), AllChem.ETKDGv2(), AllChem.ETKDGv3()]:
                if AllChem.EmbedMolecule(mol, method) == 0:
                    embedded = True
                    break

            if not embedded:
                raise ValueError("3D embedding failed.")

            # Soft optimization
            try:
                AllChem.UFFOptimizeMolecule(mol, maxIters=100)
            except:
                try:
                    props = AllChem.MMFFGetMoleculeProperties(mol)
                    ff = AllChem.MMFFGetMoleculeForceField(mol, props)
                    ff.Minimize(maxIts=150)
                except:
                    raise ValueError("3D optimization with UFF and MMFF failed.")

            # Show outputs
            if display_mode == "XYZ Coordinates":
                xyz_data = MolToXYZBlock(mol)
                st.text_area("XYZ Coordinates", value=xyz_data, height=300)

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
