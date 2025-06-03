# env_hazard_panel.py

import streamlit as st
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors
import requests

def draw_molecule(mol, width=300, height=200):
    """
    Draw an RDKit molecule as SVG and render it in Streamlit.
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')
    st.write(f'<div>{svg}</div>', unsafe_allow_html=True)

def calculate_physicochemical_properties(mol):
    """
    Compute RDKit‐based physicochemical properties for a molecule object.
    Returns a dict with Molecular Weight, logP, TPSA, etc.
    """
    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "logP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
        # You can append more RDKit descriptors here if needed
    }

def get_pubchem_cid(smiles: str):
    """
    Given a SMILES string, query PubChem PUG-REST to get the first CID.
    Returns an integer CID or None if not found / error.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{smiles}/cids/JSON"
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            return cids[0]
        else:
            return None
    except Exception as e:
        st.error(f"Error fetching CID from PubChem: {e}")
        return None

def get_ecosar_data(cid: int):
    """
    Given a PubChem CID, fetch ECOSAR predictions via PUG-REST.
    Returns a dict with keys like 'Fish LC50', 'Daphnia EC50', 'Algae EC50', 
    or None on failure.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/ECOSAR/JSON"
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        # The JSON structure usually looks like:
        # {
        #   "InformationList": {
        #     "Information": [
        #       {
        #         "CID": 1234,
        #         "ECOSAR": [
        #           { "Category": "Neutral Organics", "Fish_LC50_mg_L": "1.2", "Daphnia_EC50_mg_L": "2.3", "Algae_EC50_mg_L": "0.5", ... },
        #           { ... if multiple classes ... }
        #         ]
        #       }
        #     ]
        #   }
        # }
        info_list = data.get("InformationList", {}).get("Information", [])
        if not info_list:
            return None

        ecosar_entries = info_list[0].get("ECOSAR", [])
        # Pick the first entry (often the “best” category). You could iterate if you want to combine them.
        if not ecosar_entries:
            return None

        first = ecosar_entries[0]
        # Extract numeric values if present, else fallback to “N/A”
        fish_lc50 = first.get("Fish_LC50_mg_L", "N/A")
        daphnia_ec50 = first.get("Daphnia_EC50_mg_L", "N/A")
        algae_ec50 = first.get("Algae_EC50_mg_L", "N/A")
        category = first.get("Category", "Unknown")

        # Convert to float where possible
        def to_float(val):
            try:
                return float(val)
            except:
                return None

        return {
            "ECOSAR Category": category,
            "Fish LC₅₀ (mg/L)": to_float(fish_lc50),
            "Daphnia EC₅₀ (mg/L)": to_float(daphnia_ec50),
            "Algae EC₅₀ (mg/L)": to_float(algae_ec50),
        }
    except Exception as e:
        st.error(f"Error fetching ECOSAR data from PubChem: {e}")
        return None

def calculate_risk_score(props: dict, ecosar: dict):
    """
    Simple risk‐scoring function:
      - +30 points if logP ≥ 4.0
      - +30 points if Fish LC₅₀ ≤ 1 mg/L (high acute fish toxicity)
      - +20 points if Daphnia EC₅₀ ≤ 1 mg/L
      - +20 points if Algae EC₅₀ ≤ 1 mg/L
    Caps at 100.
    Returns an integer score between 0 and 100.
    """
    score = 0

    # 1) Bioaccumulation proxy (logP)
    logp = props.get("logP")
    if logp is not None:
        if logp >= 4.0:
            score += 30

    # 2) Fish acute toxicity
    fish_lc50 = ecosar.get("Fish LC₅₀ (mg/L)") if ecosar else None
    if fish_lc50 is not None:
        if fish_lc50 <= 1.0:
            score += 30

    # 3) Daphnia acute toxicity
    daph_ec50 = ecosar.get("Daphnia EC₅₀ (mg/L)") if ecosar else None
    if daph_ec50 is not None:
        if daph_ec50 <= 1.0:
            score += 20

    # 4) Algae acute toxicity
    algae_ec50 = ecosar.get("Algae EC₅₀ (mg/L)") if ecosar else None
    if algae_ec50 is not None:
        if algae_ec50 <= 1.0:
            score += 20

    return min(score, 100)

def main():
    st.header("Environmental Hazard Panel")
    st.write(
        "Upload a molecule (SMILES or .mol/.sdf), visualize its structure, "
        "and obtain preliminary environmental risk indicators."
    )

    # -----------------------------
    # 1) Molecule input & RDKit draw
    # -----------------------------
    smiles_input = st.text_input("Enter SMILES string:", "")
    file_input = st.file_uploader("Or upload a .mol / .sdf file:", type=["mol", "sdf"])

    mol_obj = None
    actual_smiles = None

    if smiles_input:
        mol_obj = Chem.MolFromSmiles(smiles_input)
        actual_smiles = smiles_input
        if mol_obj is None:
            st.error("Invalid SMILES. Please check your input.")
    elif file_input is not None:
        try:
            block = file_input.read().decode("utf-8")
            mol_obj = Chem.MolFromMolBlock(block)
            if mol_obj:
                # Get SMILES back from the Mol object for later API calls
                actual_smiles = Chem.MolToSmiles(mol_obj)
            else:
                st.error("Could not parse this file. Make sure it’s a valid .mol/.sdf.")
        except Exception as e:
            st.error(f"Error reading file: {e}")

    if mol_obj is None:
        st.info("Awaiting a valid molecule to analyze…")
        return

    # ---------------
    # 2) Show structure
    # ---------------
    st.subheader("Molecule Structure")
    draw_molecule(mol_obj)

    # -----------------------------------
    # 3) Compute & display basic RDKit props
    # -----------------------------------
    st.subheader("Physicochemical Properties")
    props = calculate_physicochemical_properties(mol_obj)
    st.table(props)

    # -----------------------------------
    # 4) Fetch PubChem CID & ECOSAR data
    # -----------------------------------
    st.subheader("ECOSAR Predictions (via PubChem)")

    with st.spinner("Querying PubChem for ECOSAR data…"):
        cid = get_pubchem_cid(actual_smiles)
        if cid is None:
            st.error("No PubChem CID found for this SMILES. Skipping ECOSAR.")
            ecosar = None
        else:
            ecosar = get_ecosar_data(cid)

    if ecosar is None:
        st.write("ECOSAR data not available.")
    else:
        # Display a small table with the key ECOSAR values
        ecosar_display = {
            "Category": ecosar.get("ECOSAR Category", "N/A"),
            "Fish LC₅₀ (mg/L)": ecosar.get("Fish LC₅₀ (mg/L)", "N/A"),
            "Daphnia EC₅₀ (mg/L)": ecosar.get("Daphnia EC₅₀ (mg/L)", "N/A"),
            "Algae EC₅₀ (mg/L)": ecosar.get("Algae EC₅₀ (mg/L)", "N/A"),
        }
        st.table(ecosar_display)

    # ---------------------------
    # 5) Calculate & show risk score
    # ---------------------------
    st.subheader("Overall Risk Score")

    if ecosar:
        risk_score = calculate_risk_score(props, ecosar)
        delta_color = "normal"
        if risk_score >= 75:
            delta_color = "inverse"  # red-ish
        elif risk_score >= 40:
            delta_color = "off"      # yellow-ish
        else:
            delta_color = "normal"   # green-ish

        st.metric(label="Risk Score (0–100)", value=f"{risk_score}", delta_color=delta_color)

        # Add a brief interpretation
        if risk_score >= 75:
            st.error("⚠️ High environmental hazard predicted!")
        elif risk_score >= 40:
            st.warning("⚠️ Moderate environmental hazard predicted.")
        else:
            st.success("✅ Low environmental hazard predicted.")
    else:
        st.write("Cannot compute overall risk score without ECOSAR data.")

    # ---------------------------
    # 6) Download report (stub)
    # ---------------------------
    st.subheader("Download Full Report")
    st.download_button(
        label="Download Report (JSON)",
        data=None,  # later: `json.dumps({...all results...})`
        file_name="env_risk_report.json",
        mime="application/json",
        disabled=True
    )
