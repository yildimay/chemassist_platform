# env_hazard_panel.py

import os
import json
import streamlit as st
import openai
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors
from dotenv import load_dotenv

# ─── Load local .env if present ───────────────────────────────────────────────
# (This will do nothing on Render, because .env isn't checked in. Locally, it loads your .env.)
load_dotenv()

# ─── Configure OpenAI client (point to Groq) ─────────────────────────────────
openai.api_key  = os.getenv("GROQ_API_KEY")
openai.api_base = os.getenv("OPENAI_API_BASE", "https://api.groq.com/openai")

def draw_molecule(mol, width=300, height=200):
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText().replace("svg:", "")
    st.write(f"<div>{svg}</div>", unsafe_allow_html=True)

def calculate_physicochemical_properties(mol):
    return {
        "Molecular Weight": round(Descriptors.MolWt(mol), 2),
        "logP":           round(Descriptors.MolLogP(mol), 2),
        "TPSA":           round(Descriptors.TPSA(mol), 2),
    }

def query_groq_toxicity(smiles: str) -> dict | None:
    """  
    Ask the Groq‐hosted LLM (via openai.ChatCompletion) for Fish/Daphnia/Algae endpoints.
    """
    if not openai.api_key:
        st.error("OPENAI_API_KEY is not set. Cannot query Groq.")
        return None

    prompt = (
        f"SMILES: {smiles}\n"
        "Predict the following ecotoxicity endpoints in mg/L (return exactly valid JSON):\n"
        "{\n"
        '  "fish_lc50_mg_L": ,\n'
        '  "daphnia_ec50_mg_L": ,\n'
        '  "algae_ec50_mg_L": \n'
        "}\n"
    )

    try:
        completion = openai.ChatCompletion.create(
            model="llama3-70b-fine-tuned-toxicity",  # your fine-tuned model name on Groq
            messages=[{"role": "user", "content": prompt}],
            temperature=0.0,
            max_tokens=64,
            stop=["}"],   # ensure the JSON ends cleanly
        )
        text = completion.choices[0].message.content.strip()
        if not text.endswith("}"):
            text += "}"
        return json.loads(text)
    except Exception as e:
        st.error(f"Error querying Groq toxicity model: {e}")
        return None

def calculate_risk_score(props: dict, preds: dict) -> int:
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
    st.header("Environmental Hazard Panel (Groq-Powered)")

    # 1) Input
    smiles_input = st.text_input("Enter SMILES string:", "")
    file_input   = st.file_uploader("Upload a .mol/.sdf file:", type=["mol","sdf"])

    mol_obj      = None
    actual_smiles= None

    if smiles_input:
        mol_obj       = Chem.MolFromSmiles(smiles_input)
        actual_smiles = smiles_input
        if not mol_obj:
            st.error("Invalid SMILES. Please fix it.")
    elif file_input:
        try:
            molblock = file_input.read().decode("utf-8")
            mol_obj  = Chem.MolFromMolBlock(molblock)
            if mol_obj:
                actual_smiles = Chem.MolToSmiles(mol_obj)
            else:
                st.error("Unable to parse the uploaded file.")
        except Exception as e:
            st.error(f"Error reading file: {e}")

    if not mol_obj:
        st.info("Waiting for a valid molecule…")
        return

    # 2) Draw
    st.subheader("Molecule Structure")
    draw_molecule(mol_obj)

    # 3) PhysChem Props
    st.subheader("Physicochemical Properties")
    props = calculate_physicochemical_properties(mol_obj)
    st.table(props)

    # 4) Query Groq LLM
    st.subheader("Groq LLM Toxicity Predictions")
    with st.spinner("Querying Groq model..."):
        preds = query_groq_toxicity(actual_smiles)

    if preds is None:
        st.warning("No toxicity predictions (check your API key or model).")
        return

    fish_val  = preds.get("fish_lc50_mg_L", "N/A")
    daph_val  = preds.get("daphnia_ec50_mg_L", "N/A")
    algae_val = preds.get("algae_ec50_mg_L", "N/A")

    st.table({
        "Endpoint":        ["Fish LC₅₀ (mg/L)", "Daphnia EC₅₀ (mg/L)", "Algae EC₅₀ (mg/L)"],
        "Predicted Value": [fish_val,             daph_val,            algae_val]
    })

    # 5) Risk Score
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

    # 6) Download report
    st.subheader("Download Full Report")
    report = {
        "SMILES": actual_smiles,
        "Physicochemical": props,
        "Predictions": {
            "fish_lc50_mg_L":   fish_val,
            "daphnia_ec50_mg_L":daph_val,
            "algae_ec50_mg_L":  algae_val
        },
        "Risk Score": score
    }
    st.download_button(
        label="Download Report (JSON)",
        data=json.dumps(report, indent=2),
        file_name="env_hazard_report.json",
        mime="application/json"
    )

# Note: do NOT put set_page_config() here.
