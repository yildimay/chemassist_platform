# env_hazard_panel.py

import os
import json
import requests
import streamlit as st
import pandas as pd

from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from urllib.parse import quote
from dotenv import load_dotenv

# ─── Load .env (for local dev) ────────────────────────────────────────────────
load_dotenv()  # ensures that if you have a local .env file with no keys in Git, they get loaded

# -------------------------------------------------------------------------------
#  Helpers to fetch public, database-driven toxicity and environmental data
# -------------------------------------------------------------------------------

def draw_molecule(mol: Chem.Mol, width: int = 300, height: int = 200):
    """
    Render an RDKit Mol as SVG in Streamlit.
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")
    st.write(f"<div>{svg}</div>", unsafe_allow_html=True)


def calculate_physicochemical_properties(mol: Chem.Mol) -> dict:
    """
    Compute basic RDKit descriptors:
      - Molecular Weight
      - logP
      - TPSA (Topological Polar Surface Area)
    Returns a dict.
    """
    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "logP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
    }


def get_pubchem_cid(smiles: str) -> int | None:
    """
    Given a SMILES string, URL-encode it and query PubChem PUG-REST to get the first CID.
    Returns the integer CID or None if not found / on error.
    """
    encoded = quote(smiles, safe="")
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded}/cids/JSON"
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        return cids[0] if cids else None
    except Exception as e:
        st.error(f"Error fetching CID from PubChem: {e}")
        return None


def get_ecosar_data(cid: int) -> dict | None:
    """
    Given a PubChem CID, fetch ECOSAR predictions (Fish, Daphnia, Algae LC₅₀/EC₅₀) via PUG-REST.
    Returns a dict with keys:
       - ECOSAR Category
       - Fish LC₅₀ (mg/L)
       - Daphnia EC₅₀ (mg/L)
       - Algae EC₅₀ (mg/L)
    or None if no data is available.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/ECOSAR/JSON"
    try:
        resp = requests.get(url, timeout=10)
        # If PubChem replies 400 BadRequest, treat as “no ECOSAR data”
        if resp.status_code == 400:
            return None
        resp.raise_for_status()
        data = resp.json()
        info_list = data.get("InformationList", {}).get("Information", [])
        if not info_list:
            return None
        ecosar_entries = info_list[0].get("ECOSAR", [])
        if not ecosar_entries:
            return None
        first = ecosar_entries[0]
        def to_float(x):
            try:
                return float(x)
            except:
                return None
        return {
            "ECOSAR Category": first.get("Category", "Unknown"),
            "Fish LC₅₀ (mg/L)": to_float(first.get("Fish_LC50_mg_L")),
            "Daphnia EC₅₀ (mg/L)": to_float(first.get("Daphnia_EC50_mg_L")),
            "Algae EC₅₀ (mg/L)": to_float(first.get("Algae_EC50_mg_L")),
        }
    except Exception as e:
        st.error(f"Error fetching ECOSAR data: {e}")
        return None


def get_predicted_bcf(cid: int) -> float | None:
    """
    Fetch the Predicted_BCF (Bioconcentration Factor) from PubChem Property endpoint.
    Returns the float BCF (L/kg) or None if missing.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Predicted_BCF/JSON"
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if not props:
            return None
        val = props[0].get("Predicted_BCF")
        try:
            return float(val)
        except:
            return None
    except Exception:
        return None  # don’t spam errors if BCF is simply not available


def get_predicted_water_solubility(cid: int) -> float | None:
    """
    Fetch the Predicted_Water_Solubility (mg/L) from PubChem Property endpoint.
    Returns the float solubility or None if missing.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Predicted_Water_Solubility/JSON"
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        props = data.get("PropertyTable", {}).get("Properties", [])
        if not props:
            return None
        val = props[0].get("Predicted_Water_Solubility")
        try:
            return float(val)
        except:
            return None
    except Exception:
        return None


def get_assay_toxicity(cid: int) -> dict:
    """
    Query PubChem’s assaysummary endpoint for this CID and filter out any assays
    with “Fish”, “Daphnia”, or “Algae” in their description. Returns a dict:
      {
        "Fish Assay Count": int,
        "Daphnia Assay Count": int,
        "Algae Assay Count": int
      }
    (We simply count how many active assays exist for each category,
     as a rough proxy for “active” / “toxic” flags.)
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    result = {"Fish Assay Count": 0, "Daphnia Assay Count": 0, "Algae Assay Count": 0}
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        summaries = data.get("AssaySummaries", {}).get("AssaySummary", [])
        for entry in summaries:
            desc = entry.get("Description", "").lower()
            # We count only assays flagged “Active” (i.e. likely to show toxicity)
            outcome = entry.get("Outcome", "").lower()
            if outcome != "active":
                continue
            if "fish" in desc:
                result["Fish Assay Count"] += 1
            if "daphnia" in desc:
                result["Daphnia Assay Count"] += 1
            if "algae" in desc:
                result["Algae Assay Count"] += 1
        return result
    except Exception:
        return result  # if the call fails, just return zeros


def calculate_risk_score_from_db(props: dict,
                                 ecosar: dict | None,
                                 bcf: float | None,
                                 assay_counts: dict) -> int:
    """
    Heuristic risk score (0–100) based on database‐derived endpoints:
      - +20 if logP ≥ 4.0 (possible bioaccumulation)
      - +20 if ECOSAR Fish LC₅₀ ≤ 1 mg/L
      - +15 if Predicted_BCF ≥ 2000 L/kg
      - +10 if Predicted Water Solubility < 1 mg/L (very insoluble → more persistent)
      - +10 if Fish Assay Count ≥ 1 (active fish‐toxicity assays)
      - +10 if Daphnia Assay Count ≥ 1
      - +10 if Algae Assay Count ≥ 1
    Caps at 100.
    """
    score = 0

    logp = props.get("logP")
    if logp is not None and logp >= 4.0:
        score += 20

    if ecosar:
        fish_lc50 = ecosar.get("Fish LC₅₀ (mg/L)")
        if fish_lc50 is not None and fish_lc50 <= 1.0:
            score += 20

    if bcf is not None and bcf >= 2000:
        score += 15

    # if it’s extremely insoluble, that can correlate with persistence
    if assay_counts.get("Fish Assay Count", 0) >= 1:
        score += 10
    if assay_counts.get("Daphnia Assay Count", 0) >= 1:
        score += 10
    if assay_counts.get("Algae Assay Count", 0) >= 1:
        score += 10

    # penalize very low water solubility (makes it more persistent in sediments)
    if assay_counts.get("Algae Assay Count", 0) == 0 and bcf is None:
        # If we have no assay count, we can check water sol:
        ws = assay_counts.get("Predicted Water Solubility (mg/L)", None)
        if ws is not None and ws < 1.0:
            score += 10

    return min(score, 100)


# -------------------------------------------------------------------------------
#  Streamlit UI
# -------------------------------------------------------------------------------

def main():
    # **No** st.set_page_config() here; that belongs in your root script exactly once

    st.header("Environmental Hazard Panel (Database-Driven)")
    st.write("""
        Upload a molecule (SMILES or .mol/.sdf), see its 2D structure and basic properties, 
        then fetch publicly available toxicity/eco-hazard data from PubChem.  
        Finally, compute a simple risk score based on those endpoints.
    """)

    # ─── Molecule Input ──────────────────────────────────────────────────────────
    smiles_input = st.text_input("Enter SMILES string:", "")
    file_input   = st.file_uploader("Or upload a .mol / .sdf file:", type=["mol", "sdf"])

    mol_obj       = None
    actual_smiles = None

    if smiles_input:
        mol_obj       = Chem.MolFromSmiles(smiles_input)
        actual_smiles = smiles_input
        if not mol_obj:
            st.error("Invalid SMILES. Please check your input.")
    elif file_input is not None:
        try:
            block   = file_input.read().decode("utf-8")
            mol_obj = Chem.MolFromMolBlock(block)
            if mol_obj:
                # Convert uploaded Mol back to SMILES for API queries
                actual_smiles = Chem.MolToSmiles(mol_obj)
            else:
                st.error("Could not parse the file. Make sure it’s a valid .mol/.sdf.")
        except Exception as e:
            st.error(f"Error reading file: {e}")

    if mol_obj is None:
        st.info("Awaiting a valid molecule (SMILES or .mol/.sdf)…")
        return

    # ─── Draw the Molecule ────────────────────────────────────────────────────────
    st.subheader("Molecule Structure")
    draw_molecule(mol_obj)

    # ─── Basic RDKit Properties ─────────────────────────────────────────────────
    st.subheader("Physicochemical Properties")
    props = calculate_physicochemical_properties(mol_obj)
    st.table(props)

    # ─── Get PubChem CID ─────────────────────────────────────────────────────────
    st.subheader("PubChem Identifier & Data Fetch")
    cid = get_pubchem_cid(actual_smiles)
    if cid is None:
        st.error("No PubChem CID found for this SMILES. Unable to fetch database data.")
        return
    st.write(f"**PubChem CID:** {cid}")

    # ─── 1) ECOSAR Predictions ───────────────────────────────────────────────────
    ecosar = get_ecosar_data(cid)
    if ecosar is None:
        st.write("**ECOSAR data**: Not available for this compound.")
    else:
        st.write("**ECOSAR Predictions:**")
        st.table({
            "Category": [ecosar.get("ECOSAR Category", "N/A")],
            "Fish LC₅₀ (mg/L)": [ecosar.get("Fish LC₅₀ (mg/L)", "N/A")],
            "Daphnia EC₅₀ (mg/L)": [ecosar.get("Daphnia EC₅₀ (mg/L)", "N/A")],
            "Algae EC₅₀ (mg/L)": [ecosar.get("Algae EC₅₀ (mg/L)", "N/A")],
        })

    # ─── 2) Predicted Bioconcentration Factor (BCF) ──────────────────────────────
    bcf = get_predicted_bcf(cid)
    if bcf is None:
        st.write("**Predicted BCF**: Not available.")
    else:
        st.write(f"**Predicted BCF (L/kg)**: {bcf:.1f}")

    # ─── 3) Predicted Water Solubility ────────────────────────────────────────────
    wc = get_predicted_water_solubility(cid)
    if wc is None:
        st.write("**Predicted Water Solubility (mg/L)**: Not available.")
    else:
        st.write(f"**Predicted Water Solubility (mg/L)**: {wc:.2f}")

    # ─── 4) Assay Summary: Fish, Daphnia, Algae Counts ────────────────────────────
    assay_counts = get_assay_toxicity(cid)
    st.write("**PubChem BioAssay ‘Active’ Counts:**")
    df_counts = pd.DataFrame([
        { "Endpoint": k, "Active Assay Count": v } 
        for k, v in assay_counts.items()
        if k.endswith("Assay Count")
    ])
    st.table(df_counts)

    # ─── Compute & Display Risk Score ────────────────────────────────────────────
    st.subheader("Overall Risk Score")
    # Note: we’ll feed the ‘Predicted Water Solubility’ into our scoring logic as well
    assay_counts["Predicted Water Solubility (mg/L)"] = wc
    risk_score = calculate_risk_score_from_db(props, ecosar, bcf, assay_counts)

    if risk_score >= 75:
        st.metric(label="Risk Score (0–100)", value=risk_score, delta_color="inverse")
        st.error("⚠️ High environmental hazard predicted!")
    elif risk_score >= 40:
        st.metric(label="Risk Score (0–100)", value=risk_score, delta_color="off")
        st.warning("⚠️ Moderate environmental hazard predicted.")
    else:
        st.metric(label="Risk Score (0–100)", value=risk_score, delta_color="normal")
        st.success("✅ Low environmental hazard predicted.")

    # ─── Download All Data as JSON Report ────────────────────────────────────────
    st.subheader("Download Full Report")
    report_data = {
        "SMILES": actual_smiles,
        "PubChem CID": cid,
        "Physicochemical Properties": props,
        "ECOSAR Predictions": ecosar or {"Category": None, "Fish LC₅₀ (mg/L)": None, "Daphnia EC₅₀ (mg/L)": None, "Algae EC₅₀ (mg/L)": None},
        "Predicted BCF (L/kg)": bcf,
        "Predicted Water Solubility (mg/L)": wc,
        "Assay Counts": assay_counts,
        "Risk Score": risk_score
    }
    st.download_button(
        label="Download JSON Report",
        data=json.dumps(report_data, indent=2),
        file_name="env_hazard_report.json",
        mime="application/json"
    )
