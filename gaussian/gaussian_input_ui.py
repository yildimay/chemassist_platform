import os
import requests
import json
import streamlit as st


def truncate_structure_input(structure: str, max_chars: int = 3000) -> str:
    """
    Trims long molecule structure data for API safety.
    Adds a note to the prompt so the AI knows it's partial.
    """
    if len(structure) <= max_chars:
        return structure
    else:
        lines = structure.splitlines()
        truncated_lines = []
        total_len = 0

        for line in lines:
            if total_len + len(line) > max_chars:
                break
            truncated_lines.append(line)
            total_len += len(line) + 1  # for \n

        return "\n".join(truncated_lines) + "\n... [structure truncated for safety]"


def call_groq_for_gjf(job_name, method_basis, charge, multiplicity, structure):
    GROQ_API_KEY = os.environ.get("GROQ_API_KEY")  # Make sure it's set in Render

    if not GROQ_API_KEY:
        raise ValueError("GROQ_API_KEY not found in environment variables.")

    safe_structure = truncate_structure_input(structure)

    prompt = f"""
Write a Gaussian input (.gjf) file using the following data:
- Job name: {job_name}
- Method/Basis: {method_basis}
- Charge: {charge}
- Multiplicity: {multiplicity}
- Molecular structure (may be truncated):\n{safe_structure}
Only output the .gjf content.
"""

    headers = {
        "Authorization": f"Bearer {GROQ_API_KEY}",
        "Content-Type": "application/json"
    }

    payload = {
        "model": "llama3-70b-8192",
        "messages": [{"role": "user", "content": prompt}],
        "temperature": 0.5
    }

    response = requests.post("https://api.groq.com/openai/v1/chat/completions",
                             headers=headers, data=json.dumps(payload))
    response.raise_for_status()
    return response.json()["choices"][0]["message"]["content"]


# 🔁 Use this inside your Streamlit UI:
if st.button("Create .gjf File"):
    with st.spinner("Generating input file using AI..."):
        try:
            gjf_result = call_groq_for_gjf(job_name, method_basis, charge, multiplicity, structure)
            st.code(gjf_result.strip(), language='gjf')
        except Exception as e:
            st.error(f"[ERROR] AI processing failed: {e}")
