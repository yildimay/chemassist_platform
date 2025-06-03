# env_hazard_panel.py

import os
import json
import streamlit as st
import openai                       # ensure openai>=1.0.0 is installed
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Draw import rdMolDraw2D
from dotenv import load_dotenv

# ─── Load .env for local development (ignored by Git via .gitignore) ─────────
load_dotenv()

# ─── Configure openai to point at Groq’s OpenAI-compatible endpoint ──────────
openai.api_key  = os.getenv("GROQ_API_KEY")
openai.api_base = os.getenv("OPENAI_API_BASE", "https://api.groq.com/openai")

def draw_molecule(mol, width: int = 300, height: int = 200):
    """
    Render an RDKit Mol as SVG in Streamlit.
    """
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")
    st.write(f"<div>{svg}</div>", unsafe_allow_html=True)

def calculate_physicochemical_properties(mol) -> dict:
    """
    Compute:
      - Molecular Weight
      - logP
      - TPSA
    via RDKit Descriptors.
    """
    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "logP": round(Descriptors.MolLogP(mol), 2),
        "TPSA": round(Descriptors.TPSA(mol), 2),
    }

def query_groq_toxicity(smiles: str) -> dict | None:
    """
    Query the Groq-hosted LLM (via OpenAI-compatible API) for Fish LC₅₀,
    Daphnia EC₅₀, and Algae EC₅₀. Returns a dict or None on error.
    """
    if not openai.api_key:
        st.error("OPENAI_API_KEY not set in environment. Cannot query Groq model.")
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

    try:
        # ─── NEW STYLE (openai>=1.0.0) ──────────────────────────────────────────
        resp = openai.chat.completions.create(
            model="llama3-70b-fine-tuned-toxicity",  # your actual fine-tuned model name
            messages=[{"role": "user", "content": prompt}],
            temperature=0.0,
            max_tokens=64,
            stop=["}"],  # ensures JSON closes properly
        )
        text = resp.choices[0].message.content.strip()
        if not text.endswith("}"):
            text += "}"
        return json.loads(text)

    except Exception as e:
        st.error(f"Error querying Groq toxicity model: {e}")
        return None

def calculate_risk_score(props: dict, preds: dict) -> int:
    """
    Heuristic risk score (0–100):
      +30 if logP ≥ 4.0 (bioaccumulation concern)
      +30 if Fish LC₅₀ ≤ 1 mg/L (high fish toxicity)
      +20 if Daphnia EC₅₀ ≤ 1 mg/L
      +20 if Algae EC₅₀ ≤ 1 mg/L
    """
    score = 0

    logp = props.get("logP")
    if logp is not None and logp >= 4.0:
        score += 30

    fish_val  = preds.get("fish_lc50_mg_L")
    daph_val  = preds.get("daphnia_ec50_mg_L")
    algae_val = preds.get("algae_ec50_mg_L")

    if fish_val  is not None and fish_val  <= 1.0: score += 30
    if daph_val  is not None and daph_val  <= 1.0: score += 20
    if algae_val is not None and algae_val <= 1.0: score += 20

    return min(score, 100)

def main():
    # ─── STREAMLIT UI ──────────────────────────────────────────────────────────
    st.header("Environmental Hazard Panel (Groq-Powered QSAR)")
    st.write(
        "Upload a molecule (SMILES or .mol/.sdf), visualize it, compute basic RDKit properties, "
        "and retrieve ecotoxicity predictions from a Groq-hosted LLM via the OpenAI-compatible API."
    )

    smiles_input = st.text_input("Enter SMILES string:", "")
    file_input   = st.file_uploader("Or upload a .mol / .sdf file:", type=["mol", "sdf"])

    mol_obj       = None
    actual_smiles = None

    if smiles_input:
        mol_obj       = Chem.MolFromSmiles(smiles_input)
        actual_smiles = smiles_input
        if not mol_obj:
            st.error("Invalid SMILES. Please double-check your input.")
    elif file_input is not None:
        try:
            block   = file_input.read().decode("utf-8")
            mol_obj = Chem.MolFromMolBlock(block)
            if mol_obj:
                actual_smiles = Chem.MolToSmiles(mol_obj)
            else:
                st.error("Could not parse the uploaded file. Ensure it’s a valid .mol/.sdf.")
        except Exception as e:
            st.error(f"Error reading file: {e}")

    if mol_obj is None:
        st.info("Awaiting a valid molecule (SMILES or .mol/.sdf)…")
        return

    # ─── 1) Draw the molecule ───────────────────────────────────────────────────
    st.subheader("Molecule Structure")
    draw_molecule(mol_obj)

    # ─── 2) Show basic RDKit properties ───────────────────────────────────────
    st.subheader("Physicochemical Properties")
    props = calculate_physicochemical_properties(mol_obj)
    st.table(props)

    # ─── 3) Query Groq-hosted LLM for toxicity ────────────────────────────────
    st.subheader("Groq LLM Toxicity Predictions")
    with st.spinner("Sending SMILES to Groq model…"):
        preds = query_groq_toxicity(actual_smiles)

    if preds is None:
        st.warning("No toxicity predictions returned. Check your API key or model configuration.")
        return

    fish_val  = preds.get("fish_lc50_mg_L", "N/A")
    daph_val  = preds.get("daphnia_ec50_mg_L", "N/A")
    algae_val = preds.get("algae_ec50_mg_L", "N/A")

    st.table({
        "Endpoint": ["Fish LC₅₀ (mg/L)", "Daphnia EC₅₀ (mg/L)", "Algae EC₅₀ (mg/L)"],
        "Predicted Value": [fish_val, daph_val, algae_val]
    })

    # ─── 4) Compute & show overall risk score ─────────────────────────────────
    st.subheader("Overall Risk Score")
    score = calculate_risk_score(props, preds)
    if score >= 75:
        st.metric(label="Risk Score (0–100)", value=score, delta_color="inverse")
        st.error("⚠️ High environmental hazard predicted!")
    elif score >= 40:
        st.metric(label="Risk Score (0–100)", value=score, delta_color="off")
        st.warning("⚠️ Moderate environmental hazard predicted.")
    else:
        st.metric(label="Risk Score (0–100)", value=score, delta_color="normal")
        st.success("✅ Low environmental hazard predicted.")

    # ─── 5) Download a JSON report ─────────────────────────────────────────────
    st.subheader("Download Full Report")
    report_data = {
        "SMILES": actual_smiles,
        "Physicochemical": props,
        "Toxicity Predictions": {
            "fish_lc50_mg_L":   fish_val,
            "daphnia_ec50_mg_L":daph_val,
            "algae_ec50_mg_L":  algae_val
        },
        "Risk Score": score
    }
    st.download_button(
        label="Download Report (JSON)",
        data=json.dumps(report_data, indent=2),
        file_name="env_hazard_report.json",
        mime="application/json"
    )
