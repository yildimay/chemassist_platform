import streamlit as st
import requests
import difflib
import base64
from PIL import Image
import os

# === AI CONFIG ===
GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
EXPLAIN_MODEL = "llama3-8b-8192"
FIX_MODEL = "llama3-70b-8192"
GROQ_API_KEY = os.environ.get("GROQ_API_KEY", "")

# === HELPER FUNCTIONS ===
def call_model(prompt, model="groq-explain", openai_api_key=None):
    if openai_api_key:
        headers = {
            "Authorization": f"Bearer {openai_api_key}",
            "Content-Type": "application/json"
        }
        payload = {
            "model": "gpt-4",
            "messages": [{"role": "user", "content": prompt}],
            "max_tokens": 1000,
            "temperature": 0.3
        }
        res = requests.post("https://api.openai.com/v1/chat/completions", headers=headers, json=payload)
        try:
            return res.json()["choices"][0]["message"]["content"]
        except:
            return "[❌ Error from OpenAI GPT-4]"
    else:
        llama_model = EXPLAIN_MODEL if model == "groq-explain" else FIX_MODEL
        headers = {
            "Authorization": f"Bearer {GROQ_API_KEY}",
            "Content-Type": "application/json"
        }
        payload = {
            "model": llama_model,
            "messages": [{"role": "user", "content": prompt}],
            "max_tokens": 1000,
            "temperature": 0.3
        }
        try:
            res = requests.post(GROQ_API_URL, headers=headers, json=payload)
            res.raise_for_status()
            return res.json()["choices"][0]["message"]["content"]
        except Exception as e:
            st.error(f"❌ Error: {e}")
            return None

def read_uploaded_file(file):
    return file.read().decode("utf-8", errors="ignore")

def generate_diff(original, fixed):
    return "\n".join(difflib.unified_diff(
        original.strip().splitlines(),
        fixed.strip().splitlines(),
        fromfile="original.gjf", tofile="fixed.gjf", lineterm=""
    ))

def gaussian_ui():
    st.title("Gaussian Error Assistant")

    st.subheader("📂 File Upload")
    use_test = st.checkbox("Use built-in sample files")
    gjf_content = log_content = ""

    if use_test:
        gjf_content = "%chk=broken.chk\n\nGaussian Broken Job\n\n0 1\nC       0.000000    0.000000    0.000000\nH       0.000000    0.000000    1.089000\nH       1.026719    0.000000   -0.363000\nH      -0.513360   -0.889165   -0.363000\nH      -0.513360    0.889165   -0.363000"
        log_content = " Optimization failed to converge\n Displacement too small, but forces too large\n Error termination via Lnk1e in /g16/l101.exe"
    else:
        gjf_file = st.file_uploader("Upload broken .gjf file", type=["gjf", "com"])
        log_file = st.file_uploader("Upload related .log or .out file", type=["log", "out"])
        if gjf_file: gjf_content = read_uploaded_file(gjf_file)
        if log_file: log_content = read_uploaded_file(log_file)

    if st.button("🔍 Analyze & Fix"):
        if not gjf_content or not log_content:
            st.warning("Missing input files.")
            st.stop()

        with st.spinner("Analyzing..."):
            explain_prompt = f"""You're a Gaussian expert.\n\nAnalyze this input and log file.\n\n### Problem\n### Likely Cause\n### Recommended Solution\n\n-- .gjf file --\n{gjf_content}\n\n-- .log file --\n{log_content}"""
            explanation = call_model(explain_prompt, model="groq-explain")
            st.subheader("📘 Explanation")
            st.markdown(explanation)

            fix_prompt = f"""Fix the broken .gjf using the .log file.\nReturn only the corrected .gjf file.\n\n-- .gjf --\n{gjf_content}\n\n-- .log --\n{log_content}"""
            fixed_gjf = call_model(fix_prompt, model="groq-fix")
            if fixed_gjf:
                st.subheader("✅ Fixed .gjf")
                st.code(fixed_gjf, language="text")
                st.download_button("Download Fixed .gjf", fixed_gjf, file_name="fixed.gjf")
                st.subheader("🔍 Difference")
                st.code(generate_diff(gjf_content, fixed_gjf), language="diff")
