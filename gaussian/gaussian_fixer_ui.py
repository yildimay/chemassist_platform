import streamlit as st
import requests
import difflib
import os

def gaussian_fixer_ui():
    st.title("Gaussian Error Fixer")

    st.subheader("Upload Files")
    gjf_file = st.file_uploader("Upload broken .gjf file", type=["gjf", "com"])
    log_file = st.file_uploader("Upload related .log or .out file", type=["log", "out"])

    if not gjf_file or not log_file:
        st.info("Please upload both a .gjf and a .log/.out file to proceed.")
        return

    gjf_content = gjf_file.read().decode("utf-8", errors="ignore")
    log_content = log_file.read().decode("utf-8", errors="ignore")

GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
GROQ_API_KEY = os.getenv("GROQ_API_KEY", "")  # Now it works with Render env vars
MODEL = "llama3-70b-8192"

    def call_groq(prompt):
        headers = {
            "Authorization": f"Bearer {GROQ_API_KEY}",
            "Content-Type": "application/json"
        }
        data = {
            "model": MODEL,
            "messages": [{"role": "user", "content": prompt}],
            "max_tokens": 1000,
            "temperature": 0.3
        }
        try:
            res = requests.post(GROQ_API_URL, headers=headers, json=data)
            res.raise_for_status()
            return res.json()["choices"][0]["message"]["content"]
        except Exception as e:
            st.error(f"AI request failed: {e}")
            return None

    if st.button("Analyze and Fix"):
        with st.spinner("Analyzing error and generating fix..."):
            prompt = f"""
You're a Gaussian expert.
Given the following input file (.gjf) and log file (.log), explain what went wrong and generate a corrected .gjf file.

Input (.gjf):
{gjf_content}

Log (.log):
{log_content}

Output ONLY the corrected .gjf file.
"""
            response = call_groq(prompt)
            if response:
                st.subheader("Fixed .gjf File")
                st.code(response, language="text")
                st.download_button("Download Fixed .gjf", response, file_name="fixed_input.gjf", mime="text/plain")
            else:
                st.error("Could not generate a fix.")
