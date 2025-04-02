import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import base64
import os
import tempfile
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D

def smiles_ui():
    st.title("🧬 SMILES to MOL Converter")

    st.header("📩 Input")
    smiles = st.text_input("Enter SMILES string:", "")

    view_mode = st.radio("Display Mode", ["3D Viewer", "XYZ Coordinates"])

    if st.button("Convert to MOL") and smiles:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            mol = Chem.AddHs(mol)
            embed_result = AllChem.EmbedMolecule(mol)
            if embed_result != 0:
                st.warning("3D embedding failed.")
                return

            opt_result = AllChem.UFFOptimizeMolecule(mol)
            if opt_result != 0:
                st.warning("3D optimization failed.")
                return

            mol_block = Chem.MolToMolBlock(mol)
            xyz_block = Chem.MolToXYZBlock(mol)

            st.success("MOL file generated!")
            b64 = base64.b64encode(mol_block.encode()).decode()
            href = f'<a href="data:file/mol;base64,{b64}" download="molecule.mol">Download .mol file</a>'
            st.markdown(href, unsafe_allow_html=True)

            col1, col2 = st.columns([1, 1])

            with col1:
                rdDepictor.Compute2DCoords(mol)
                drawer = rdMolDraw2D.MolDraw2DCairo(300, 300)
                drawer.drawOptions().addStereoAnnotation = True
                drawer.drawOptions().explicitMethyl = True
                drawer.drawOptions().addAtomIndices = False
                drawer.DrawMolecule(mol)
                drawer.FinishDrawing()
                png = drawer.GetDrawingText()
                st.image(png, caption="Molecule Preview")

            with col2:
                if view_mode == "XYZ Coordinates":
                    st.subheader("📄 XYZ Coordinates")
                    st.code(xyz_block, language="xyz")
                elif view_mode == "3D Viewer":
                    st.subheader("🧬 3D Molecule Viewer")
                    try:
                        import py3Dmol
                        viewer = py3Dmol.view(width=400, height=300)
                        viewer.addModel(xyz_block, "xyz")
                        viewer.setStyle({"stick": {}})
                        viewer.zoomTo()
                        st.components.v1.html(viewer._make_html(), height=300)
                    except Exception as e:
                        st.warning(f"3D structure could not be generated: {e}")
        else:
            st.error("Invalid SMILES string.")

def gaussian_input_ui():
    pass
