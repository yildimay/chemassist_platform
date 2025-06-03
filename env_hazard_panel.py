# env_hazard_panel.py (excerpt)

import os
import streamlit as st
import requests
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem import Descriptors
from urllib.parse import quote

# … other helper functions as before (draw_molecule, calculate_physicochemical_properties) …

def query_groq_toxicity(smiles: str) -> dict | None:
    """
    Send a prompt to your Groq-hosted LLaMA3-70B (fine-tuned) to get toxicity predictions.
    Returns a dict like:
      {
        "Fish_LC50_mg_L": float,
        "Daphnia_EC50_mg_L": float,
        "Algae_EC50_mg_L": float
      }
    or None on error.
    """
    groq_api_url = "https://api.groq.example.com/v1/completions"  # <— your actual endpoint
    api_key = os.getenv("GROQ_API_KEY")  # assume you set this in your environment

    # Construct a clear, minimal prompt:
    prompt = (
        f"SMILES: {smiles}\n"
        "Predict the following endpoints in mg/L (just give numbers):\n"
        "Fish LC50: \nDaphnia EC50: \nAlgae EC50: \n\n"
        "Return JSON with keys fish_lc50_mg_L, daphnia_ec50_mg_L, algae_ec50_mg_L."
    )

    headers = {
        "Content-Type": "application/json",
        "Authorization": f"Bearer {api_key}"
    }
    payload = {
        "model": "llama3-70b-fine-tuned-toxicity",
        "prompt": prompt,
        "max_tokens": 64,
        "temperature": 0.0,
        "stop": ["\n\n"]  # or whatever stops your model
    }

    try:
        resp = requests.post(groq_api_url, json=payload, headers=headers, timeout=20)
        resp.raise_for_status()
        result = resp.json()
        # Extract the model’s “generated text”:
        text = result["choices"][0]["text"].strip()
        # Example text might be:
        #   {
        #     "fish_lc50_mg_L": 2.3,
        #     "daphnia_ec50_mg_L": 5.6,
        #     "algae_ec50_mg_L": 1.2
        #   }
        # So we can parse JSON directly:
        pred = None
        try:
            pred = json.loads(text)
        except ValueError:
            # If the model returned plain numbers line by line, you might have to parse manually
            # e.g. "Fish LC50: 2.3 mg/L\nDaphnia EC50: 5.6 mg/L\nAlgae EC50: 1.2 mg/L"
            # You can write quick regex or splitlines to grab floats.
            pass

        return pred
    except Exception as e:
        st.error(f"Error querying Groq toxicity model: {e}")
        return None


def main():
    st.header("Environmental Hazard Panel (Groq-Powered QSAR)")
    st.write(
        "Upload a molecule (SMILES or .mol/.sdf), visualize its structure, "
        "and obtain toxicity predictions from a Groq-hosted LLM model."
    )

    # 1) Molecule input & RDKit draw
    smiles_input = st.text_input("Enter SMILES string:", "")
    file_input = st.file_uploader("Or upload a .mol /.sdf file:", type=["mol", "sdf"])

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
                st.error("Could not parse this file. Make sure it’s a valid .mol/.sdf.")
        except Exception as e:
            st.error(f"Error reading file: {e}")

    if mol_obj is None:
        st.info("Awaiting a valid molecule to analyze…")
        return

    # 2) Show structure
    st.subheader("Molecule Structure")
    draw_molecule(mol_obj)

    # 3) Basic RDKit props
    st.subheader("Physicochemical Properties")
    props = calculate_physicochemical_properties(mol_obj)
    st.table(props)

    # 4) Query Groq-based toxicity
    st.subheader("Groq LLM Toxicity Predictions")
    with st.spinner("Sending SMILES to Groq model…"):
        preds = query_groq_toxicity(actual_smiles)

    if preds is None:
        st.write("No toxicity predictions returned.")
        return

    # 5) Display the predicted endpoints
    fish_val = preds.get("fish_lc50_mg_L")
    daph_val = preds.get("daphnia_ec50_mg_L")
    algae_val = preds.get("algae_ec50_mg_L")

    tox_table = {
        "Endpoint": ["Fish LC₅₀ (mg/L)", "Daphnia EC₅₀ (mg/L)", "Algae EC₅₀ (mg/L)"],
        "Predicted Value": [fish_val, daph_val, algae_val]
    }
    st.table(tox_table)

    # 6) Compute a simple risk score from Groq outputs
    st.subheader("Overall Risk Score")
    score = 0
    if fish_val is not None and fish_val <= 1.0:
        score += 30
    if daph_val is not None and daph_val <= 1.0:
        score += 20
    if algae_val is not None and algae_val <= 1.0:
        score += 20
    logp = props.get("logP")
    if logp is not None and logp >= 4.0:
        score += 30
    score = min(score, 100)

    if score >= 75:
        delta_color = "inverse"
        interpretation = "⚠️ High environmental hazard predicted!"
    elif score >= 40:
        delta_color = "off"
        interpretation = "⚠️ Moderate environmental hazard predicted."
    else:
        delta_color = "normal"
        interpretation = "✅ Low environmental hazard predicted."

    st.metric(label="Risk Score (0–100)", value=f"{score}", delta_color=delta_color)
    st.write(interpretation)

    # 7) Download JSON report
    st.subheader("Download Report")
    report_data = {
        "SMILES": actual_smiles,
        "Physicochemical": props,
        "Predictions": {
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
