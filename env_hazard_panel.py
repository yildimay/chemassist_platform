# env_hazard_panel.py

import os
import json
import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D

def draw_molecule(mol, width=300, height=200):
    """
    Draw an RDKit molecule object as SVG and render it in Streamlit.
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace('svg:', '')
    st.write(f'<div>{svg}</div>', unsafe_allow_html=True)

def calculate_physicochemical_properties(mol):
    """
    Compute basic RDKit-based physicochemical properties:
    - Molecular Weight
    - logP (Octanol-water partition coefficient)
    - TPSA (Topological Polar Surface Area)
    Returns a dictionary.
    """
    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "logP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
    }

def query_groq_toxicity(smiles: str) -> dict | None:
    """
    Send a prompt to a Groq-hosted LLM (e.g., a fine-tuned LLaMA3-70B) that predicts
    Fish LC50, Daphnia EC50, and Algae EC50 in mg/L. Returns a dict or None on error.
    """
    groq_api_url = "https://api.groq.example.com/v1/completions"  # <-- replace with your actual Groq endpoint
    api_key = os.getenv("GROQ_API_KEY")
    if not api_key:
        st.error("GROQ_API_KEY environment variable not set.")
        return None

    prompt = (
        f"SMILES: {smiles}\n"
        "Predict the following endpoints in mg/L (numbers only) and return valid JSON:\n"
        "{\n"
        '  "fish_lc50_mg_L": ,\n'
        '  "daphnia_ec50_mg_L": ,\n'
        '  "algae_ec50_mg_L": \n'
        "}\n"
    )

    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }
    payload = {
        "model": "llama3-70b-fine-tuned-toxicity",  # <-- your fine-tuned model name
        "prompt": prompt,
        "max_tokens": 64,
        "temperature": 0.0,
        "stop": ["}"]  # ensure it stops after closing brace
    }

    try:
        resp = requests.post(groq_api_url, json=payload, headers=headers, timeout=20)
        resp.raise_for_status()
        result = resp.json()
        text = result["choices"][0]["text"].strip()
        # Ensure a closing brace if the stop token cut it off
        if not text.endswith("}"):
            text = text + "}"
        preds = json.loads(text)
        return preds
    except Exception as e:
        st.error(f"Error querying Groq toxicity model: {e}")
        return None

def calculate_risk_score(props: dict, preds: dict) -> int:
    """
    Simple heuristic-based risk scoring:
      - +30 if logP ≥ 4.0
      - +30 if Fish LC50 ≤ 1 mg/L
      - +20 if Daphnia EC50 ≤ 1 mg/L
      - +20 if Algae EC50 ≤ 1 mg/L
    Caps at 100.
    """
    score = 0

    # 1) logP bioaccumulation proxy
    logp = props.get("logP")
    if logp is not None and logp >= 4.0:
        score += 30

    # 2) Fish acute toxicity
    fish_val = preds.get("fish_lc50_mg_L")
    if fish_val is not None and fish_val <= 1.0:
        score += 30

    # 3) Daphnia acute toxicity
    daph_val = preds.get("daphnia_ec50_mg_L")
    if daph_val is not None and daph_val <= 1.0:
        score += 20

    # 4) Algae acute toxicity
    algae_val = preds.get("algae_ec50_mg_L")
    if algae_val is not None and algae_val <= 1.0:
        score += 20

    return min(score, 100)

def main():
    # **NO** st.set_page_config(...) here!
    st.header("Environmental Hazard Panel (Groq-Powered QSAR)")
    st.write(
        "Upload a molecule (SMILES or .mol/.sdf) to visualize its structure, "
        "compute basic physicochemical properties, and obtain ecotoxicity predictions "
        "from a Groq-hosted LLM."
    )

    # --- 1) Molecule Input ---
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
                actual_smiles = Chem.MolToSmiles(mol_obj)
            else:
                st.error("Could not parse this file. Make sure it’s a valid .mol or .sdf.")
        except Exception as e:
            st.error(f"Error reading file: {e}")

    if mol_obj is None:
        st.info("Awaiting a valid molecule to analyze…")
        return

    # --- 2) Draw the Molecule ---
    st.subheader("Molecule Structure")
    draw_molecule(mol_obj)

    # --- 3) Basic RDKit Properties ---
    st.subheader("Physicochemical Properties")
    props = calculate_physicochemical_properties(mol_obj)
    st.table(props)

    # --- 4) Query Groq for Toxicity Predictions ---
    st.subheader("Groq LLM Toxicity Predictions")
    with st.spinner("Sending SMILES to Groq model…"):
        preds = query_groq_toxicity(actual_smiles)

    if preds is None:
        st.warning("No toxicity predictions returned. Please check your Groq setup or input.")
        return

    # Display the three predicted endpoints
    fish_val = preds.get("fish_lc50_mg_L", "N/A")
    daph_val = preds.get("daphnia_ec50_mg_L", "N/A")
    algae_val = preds.get("algae_ec50_mg_L", "N/A")

    tox_table = {
        "Endpoint": ["Fish LC₅₀ (mg/L)", "Daphnia EC₅₀ (mg/L)", "Algae EC₅₀ (mg/L)"],
        "Predicted Value": [fish_val, daph_val, algae_val]
    }
    st.table(tox_table)

    # --- 5) Compute & Show Overall Risk Score ---
    st.subheader("Overall Risk Score")
    score = calculate_risk_score(props, preds)
    if score >= 75:
        delta_color = "inverse"  # red
        interpretation = "⚠️ High environmental hazard predicted!"
    elif score >= 40:
        delta_color = "off"      # yellow
        interpretation = "⚠️ Moderate environmental hazard predicted."
    else:
        delta_color = "normal"   # green
        interpretation = "✅ Low environmental hazard predicted."

    st.metric(label="Risk Score (0–100)", value=f"{score}", delta_color=delta_color)
    st.write(interpretation)

    # --- 6) Download Full Report ---
    st.subheader("Download Full Report")
    report_data = {
        "SMILES": actual_smiles,
        "Physicochemical": props,
        "Toxicity Predictions": {
            "fish_lc50_mg_L": fish_val,
            "daphnia_ec50_mg_L": daph_val,
            "algae_ec50_mg_L": algae_val
        },
        "Risk Score": score
    }
    st.download_button(
        label="Download Report (JSON)",
        data=json.dumps(report_data, indent=2),
        file_name="env_toxicity_report.json",
        mime="application/json"
    )
