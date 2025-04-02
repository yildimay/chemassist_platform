import streamlit as st
import tempfile
import base64
import os
import requests

def call_ai(prompt):
    GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
    GROQ_API_KEY = os.environ.get("GROQ_API_KEY", "")
    AI_MODEL = "llama3-70b-8192"

    headers = {
        "Authorization": f"Bearer {GROQ_API_KEY}",
        "Content-Type": "application/json"
    }
    data = {
        "model": AI_MODEL,
        "messages": [{"role": "user", "content": prompt}],
        "max_tokens": 1500,
        "temperature": 0.3
    }
    try:
        response = requests.post(GROQ_API_URL, headers=headers, json=data)
        response.raise_for_status()
        return response.json()["choices"][0]["message"]["content"]
    except Exception as e:
        return f"[ERROR] AI processing failed: {e}"

def gaussian_input_ui():
    st.title("Gaussian Input Creator")
    st.markdown("This tool generates a Gaussian input (.gjf) file using your molecule's SMILES and a reference paper.")

    st.subheader("1. Upload Reference Paper")
    uploaded_pdf = st.file_uploader("Upload a reference PDF file", type=["pdf"])

    st.subheader("2. Enter Molecule SMILES")
    smiles_input = st.text_input("Enter SMILES string of the molecule")

    st.subheader("3. Select Gaussian Version")
    selected_version = st.selectbox("Choose Gaussian version", ["G16", "G09", "Other"])

    st.subheader("4. Generate Input File")
    if st.button("Create .gjf File"):
        if not uploaded_pdf or not smiles_input:
            st.warning("Please provide both a PDF and a SMILES string.")
        else:
            with tempfile.NamedTemporaryFile(delete=False, suffix=".pdf") as tmp:
                tmp.write(uploaded_pdf.read())
                tmp_path = tmp.name
            try:
                with open(tmp_path, "rb") as f:
                    encoded_pdf = base64.b64encode(f.read()).decode()
                ai_prompt = f"""
You are a computational chemistry assistant. A user uploaded a reference paper (PDF, base64 encoded) and provided the SMILES of a molecule. Based on the content of the reference and the molecular structure, generate a Gaussian input (.gjf) file with the appropriate route section, method, basis set, and coordinates.

SMILES: {smiles_input}
Gaussian Version: {selected_version}
PDF (base64): {encoded_pdf}

Respond with only the final .gjf file content, no commentary.
"""
                result = call_ai(ai_prompt)
                st.code(result, language="text")
                st.download_button("Download .gjf File", result, file_name="generated.gjf", mime="text/plain")
            except Exception as e:
                st.error(f"Error reading or processing file: {e}")
