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

# ─── Load .env (for local development) ────────────────────────────────────────
load_dotenv()

# -------------------------------------------------------------------------------
#  Helpers to fetch public, database-driven toxicity data from PubChem
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
      - TPSA
    Returns a dictionary.
    """
    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "logP":           round(Descriptors.MolLogP(mol), 2),
        "TPSA":           round(Descriptors.TPSA(mol), 2),
    }


def count_halogens(mol: Chem.Mol) -> int:
    """
    Count the number of halogen atoms (Cl, Br, F, I) in the molecule.
    """
    return sum(1 for atom in mol.GetAtoms() 
               if atom.GetSymbol() in ("F", "Cl", "Br", "I"))


def get_pubchem_cid(smiles: str) -> int | None:
    """
    URL-encode SMILES, query PubChem to get the first CID.
    Returns integer CID, or None on error.
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
    Fetch ECOSAR data (Fish LC₅₀, Daphnia EC₅₀, Algae EC₅₀) from PubChem.
    Returns a dict or None if no data available.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/xrefs/ECOSAR/JSON"
    try:
        resp = requests.get(url, timeout=10)
        # 400 = “no ECOSAR data” → return None
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
    Fetch Predicted_BCF from PubChem. Returns float (L/kg) or None.
    """
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/Predicted_BCF/JSON"
    )
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
        return None


def get_predicted_water_solubility(cid: int) -> float | None:
    """
    Fetch Predicted_Water_Solubility (mg/L) from PubChem. Returns float or None.
    """
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/"
        f"{cid}/property/Predicted_Water_Solubility/JSON"
    )
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
    Query PubChem’s assaysummary endpoint. Count the number of “Active” assays
    for Fish, Daphnia, and Algae. Returns a dict with those counts.
    """
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/assaysummary/JSON"
    result = {
        "Fish Assay Count": 0,
        "Daphnia Assay Count": 0,
        "Algae Assay Count": 0
    }
    try:
        resp = requests.get(url, timeout=10)
        resp.raise_for_status()
        data = resp.json()
        summaries = data.get("AssaySummaries", {}).get("AssaySummary", [])
        for entry in summaries:
            desc = entry.get("Description", "").lower()
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
        return result


def calculate_risk_score_from_db(
    mol: Chem.Mol,
    props: dict,
    ecosar: dict | None,
    bcf: float | None,
    water_sol: float | None,
    assay_counts: dict
) -> int:
    """
    Refined risk‐scoring logic:
      1) logP-based (heavier penalty if very high logP):
         - +30 if logP ≥ 5.0
         - +15 if 4.0 ≤ logP < 5.0
      2) If ECOSAR Fish LC₅₀ ≤ 1 mg/L: +20
      3) BCF:
         - if (bcf is not None and bcf ≥ 2000): +15
         - else if (bcf is None AND logP available): estimate bcf = 10^(logP – 1)
             and if est_bcf ≥ 2000: +15
      4) Halogen penalty if data is missing:
         - if ecosar is None and count_halogens(mol) ≥ 4: +15
      5) Active assay counts:
         - +10 if Fish Assay Count ≥ 1
         - +10 if Daphnia Assay Count ≥ 1
         - +10 if Algae Assay Count ≥ 1
      6) Very low solubility → persistence:
         - if (water_sol is not None and water_sol < 1.0) AND (no assays flagged): +10
    Caps at 100.
    """
    score = 0
    logp = props.get("logP")

    # 1) logP penalty
    if logp is not None:
        if logp >= 5.0:
            score += 30
        elif logp >= 4.0:
            score += 15

    # 2) ECOSAR fish acute toxicity
    if ecosar:
        fish_lc50 = ecosar.get("Fish LC₅₀ (mg/L)")
        if fish_lc50 is not None and fish_lc50 <= 1.0:
            score += 20
    else:
        # 4) Halogen penalty if no ECOSAR data
        if count_halogens(mol) >= 4:
            score += 15

    # 3) BCF penalty
    if bcf is not None:
        if bcf >= 2000:
            score += 15
    else:
        # fallback estimate from logP
        if logp is not None:
            est_bcf = 10 ** (logp - 1)  # rough free‐energy estimate
            if est_bcf >= 2000:
                score += 15

    # 5) Active assay counts
    if assay_counts.get("Fish Assay Count", 0) >= 1:
        score += 10
    if assay_counts.get("Daphnia Assay Count", 0) >= 1:
        score += 10
    if assay_counts.get("Algae Assay Count", 0) >= 1:
        score += 10

    # 6) Very low water solubility → persistence
    if assay_counts.get("Fish Assay Count", 0) == 0 \
       and assay_counts.get("Daphnia Assay Count", 0) == 0 \
       and assay_counts.get("Algae Assay Count", 0) == 0:
        if water_sol is not None and water_sol < 1.0:
            score += 10

    return min(score, 100)


# -------------------------------------------------------------------------------
#  Streamlit UI
# -------------------------------------------------------------------------------

def main():
    # (No set_page_config here—your dispatcher script does that once.)

    st.header("Environmental Hazard Panel (Database-Driven, Refined Scoring)")
    st.write(
        "Upload a molecule (SMILES or .mol/.sdf), see its structure & basic properties, "
        "and fetch publicly available toxicity & environmental data from PubChem.  "
        "A refined risk score flags highly chlorinated/persistent compounds even if data is missing."
    )

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
                actual_smiles = Chem.MolToSmiles(mol_obj)
            else:
                st.error("Could not parse the file. Make sure it’s valid .mol/.sdf.")
        except Exception as e:
            st.error(f"Error reading file: {e}")

    if mol_obj is None:
        st.info("Awaiting a valid molecule (SMILES or .mol/.sdf)…")
        return

    # ─── 1) Draw the molecule ────────────────────────────────────────────────────
    st.subheader("Molecule Structure")
    draw_molecule(mol_obj)

    # ─── 2) Basic RDKit properties ───────────────────────────────────────────────
    st.subheader("Physicochemical Properties")
    props = calculate_physicochemical_properties(mol_obj)
    st.table(props)

    # ─── 3) Get PubChem CID ─────────────────────────────────────────────────────
    st.subheader("PubChem Identifier & Data Fetch")
    cid = get_pubchem_cid(actual_smiles)
    if cid is None:
        st.error("No PubChem CID found for this SMILES. Cannot fetch downstream data.")
        return
    st.write(f"**PubChem CID:** {cid}")

    # ─── 4) ECOSAR Predictions ───────────────────────────────────────────────────
    ecosar = get_ecosar_data(cid)
    if ecosar is None:
        st.write("**ECOSAR data**: Not available on PubChem.")
    else:
        st.write("**ECOSAR Predictions**")
        st.table({
            "Category":           [ecosar.get("ECOSAR Category", "N/A")],
            "Fish LC₅₀ (mg/L)":   [ecosar.get("Fish LC₅₀ (mg/L)", "N/A")],
            "Daphnia EC₅₀ (mg/L)": [ecosar.get("Daphnia EC₅₀ (mg/L)", "N/A")],
            "Algae EC₅₀ (mg/L)":   [ecosar.get("Algae EC₅₀ (mg/L)", "N/A")],
        })

    # ─── 5) Predicted BCF ────────────────────────────────────────────────────────
    bcf = get_predicted_bcf(cid)
    if bcf is None:
        st.write("**Predicted BCF**: Not available.")
    else:
        st.write(f"**Predicted BCF (L/kg)**: {bcf:.1f}")

    # ─── 6) Predicted Water Solubility ───────────────────────────────────────────
    water_sol = get_predicted_water_solubility(cid)
    if water_sol is None:
        st.write("**Predicted Water Solubility (mg/L)**: Not available.")
    else:
        st.write(f"**Predicted Water Solubility (mg/L)**: {water_sol:.2f}")

    # ─── 7) PubChem BioAssay “Active” Counts ────────────────────────────────────
    assay_counts = get_assay_toxicity(cid)
    st.write("**PubChem BioAssay ‘Active’ Counts** (Fish / Daphnia / Algae)")
    df_counts = pd.DataFrame([
        {"Endpoint": k, "Active Assay Count": v}
        for k, v in assay_counts.items()
        if k.endswith("Assay Count")
    ])
    st.table(df_counts)

    # ─── 8) Compute & display risk score ─────────────────────────────────────────
    st.subheader("Overall Risk Score")
    risk_score = calculate_risk_score_from_db(
        mol=mol_obj,
        props=props,
        ecosar=ecosar,
        bcf=bcf,
        water_sol=water_sol,
        assay_counts=assay_counts
    )

    if risk_score >= 75:
        st.metric(label="Risk Score (0–100)", value=risk_score, delta_color="inverse")
        st.error("⚠️ High environmental hazard predicted!")
    elif risk_score >= 40:
        st.metric(label="Risk Score (0–100)", value=risk_score, delta_color="off")
        st.warning("⚠️ Moderate environmental hazard predicted.")
    else:
        st.metric(label="Risk Score (0–100)", value=risk_score, delta_color="normal")
        st.success("✅ Low environmental hazard predicted.")

    # ─── 9) Download full JSON report ───────────────────────────────────────────
    st.subheader("Download Full Report")
    report_data = {
        "SMILES": actual_smiles,
        "PubChem CID": cid,
        "Physicochemical": props,
        "ECOSAR Predictions": ecosar or {
            "ECOSAR Category": None,
            "Fish LC₅₀ (mg/L)": None,
            "Daphnia EC₅₀ (mg/L)": None,
            "Algae EC₅₀ (mg/L)": None
        },
        "Predicted BCF (L/kg)": bcf,
        "Predicted Water Solubility (mg/L)": water_sol,
        "Assay Counts": assay_counts,
        "Risk Score": risk_score
    }
    st.download_button(
        label="Download JSON Report",
        data=json.dumps(report_data, indent=2),
        file_name="env_hazard_report.json",
        mime="application/json"
    )
