import streamlit as st
import requests
import difflib

st.set_page_config(page_title="Gaussian Fixer", layout="centered")

# === Custom CSS ===
st.markdown("""
<style>
.stApp {
    background-color: #0f0f0f;
    padding: 2rem;
}
h1 {
    color: #5B8DEE;
    font-size: 2.4rem;
    font-weight: 800;
    margin-bottom: 1.2rem;
}
h2, h3 {
    color: #ffffff;
    margin-top: 1.5rem;
    font-weight: 600;
}
div[data-testid="stMarkdownContainer"] p {
    font-size: 1.05rem;
    color: #cccccc;
}
.stButton > button {
    background-color: #5B8DEE;
    color: white;
    padding: 0.6rem 1.2rem;
    border-radius: 8px;
    font-weight: bold;
    border: none;
    transition: 0.3s ease;
}
.stButton > button:hover {
    background-color: #3A68D8;
    transform: scale(1.02);
}
.stTextInput input, textarea {
    background-color: #1f1f1f;
    color: #eee;
    border-radius: 8px;
}
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# === Title ===
st.markdown("""<h1>&#128640; Gaussian Error Fixer <span style='color:#aaa'>— Free. Anonymous. Clean AF.</span></h1>""", unsafe_allow_html=True)

# === Config ===
GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
EXPLAIN_MODEL = "llama3-8b-8192"
FIX_MODEL = "llama3-70b-8192"
GROQ_API_KEY = st.secrets["GROQ_API_KEY"]

def read_uploaded_file(file):
    return file.read().decode("utf-8", errors="ignore")

def call_groq(prompt, model):
    headers = {
        "Authorization": f"Bearer {GROQ_API_KEY}",
        "Content-Type": "application/json"
    }
    data = {
        "model": model,
        "messages": [{"role": "user", "content": prompt}],
        "max_tokens": 900,
        "temperature": 0.3
    }
    try:
        res = requests.post(GROQ_API_URL, headers=headers, json=data)
        res.raise_for_status()
        return res.json()["choices"][0]["message"]["content"]
    except Exception as e:
        st.error(f"❌ Error: {e}")
        return None

def generate_diff(original, fixed):
    original_lines = original.strip().splitlines()
    fixed_lines = fixed.strip().splitlines()
    return "\n".join(difflib.unified_diff(original_lines, fixed_lines, fromfile="original.gjf", tofile="fixed.gjf", lineterm=""))

# === Test Data ===
test_gjf = "%chk=broken.chk\n\nGaussian Broken Job\n\n0 1\nC       0.000000    0.000000    0.000000\nH       0.000000    0.000000    1.089000\nH       1.026719    0.000000   -0.363000\nH      -0.513360   -0.889165   -0.363000\nH      -0.513360    0.889165   -0.363000"

test_log = """ Optimization failed to converge
 Displacement too small, but forces too large
 Error termination via Lnk1e in /g16/l101.exe
"""

# === Upload UI ===
st.subheader("🧪 Upload Files")
use_test = st.checkbox("Use built-in test files")
gjf_content, log_content = "", ""

if use_test:
    gjf_content = test_gjf
    log_content = test_log
else:
    gjf_file = st.file_uploader("Upload .gjf file", type=["gjf", "com"])
    log_file = st.file_uploader("Upload .log or .out file", type=["log", "out"])
    if gjf_file: gjf_content = read_uploaded_file(gjf_file)
    if log_file: log_content = read_uploaded_file(log_file)

# === Analyze & Fix ===
if st.button("🔍 Analyze & Fix"):
    if not gjf_content or not log_content:
        st.warning("Missing files. Scroll down to enter your problem manually.")
        st.stop()

    with st.spinner("Analyzing..."):
        explain_prompt = f"""You're a Gaussian expert.\n\nAnalyze the input and log file. Format your answer as:\n\n### 🔍 Problem\n### ❓ Why It Happens\n### 🛠 How to Fix\n\n-- .gjf --\n{gjf_content}\n\n-- .log --\n{log_content}"""
        explanation = call_groq(explain_prompt, EXPLAIN_MODEL)
        st.subheader("📘 Explanation")
        st.markdown(explanation)
        st.download_button("📄 Download Explanation", explanation, file_name="explanation.txt", mime="text/plain")

        fix_prompt = f"""Fix the broken Gaussian .gjf using the .log info.\nOnly return a valid, fixed .gjf file.\n\n-- .gjf --\n{gjf_content}\n\n-- .log --\n{log_content}"""
        fixed_gjf = call_groq(fix_prompt, FIX_MODEL)
        if fixed_gjf:
            st.subheader("✅ Fixed .gjf File")
            st.code(fixed_gjf, language="text")
            st.download_button("💾 Download .gjf", fixed_gjf, file_name="fixed.gjf", mime="text/plain")
            st.subheader("🔍 Difference")
            st.code(generate_diff(gjf_content, fixed_gjf), language="diff")

            # === Retry / Refine Options ===
            st.subheader("🔁 Refine or Retry Fix")
            col1, col2, col3 = st.columns(3)

            with col1:
                if st.button("🔬 Tighten Optimization"):
                    refine_prompt = fix_prompt + "\n\nPlease tighten geometry optimization settings if possible (e.g., tighter convergence or opt=tight)."
                    refined = call_groq(refine_prompt, FIX_MODEL)
                    if refined:
                        st.subheader("🔧 Refined .gjf (Tighter Opt)")
                        st.code(refined, language="text")
                        st.download_button("💾 Download Refined .gjf", refined, file_name="refined_opt.gjf", mime="text/plain")

            with col2:
                if st.button("🚀 Increase Memory/CPUs"):
                    memory_prompt = fix_prompt + "\n\nPlease increase memory and processor count if current job settings are low."
                    refined = call_groq(memory_prompt, FIX_MODEL)
                    if refined:
                        st.subheader("🔧 More Resources .gjf")
                        st.code(refined, language="text")
                        st.download_button("💾 Download More Resources .gjf", refined, file_name="more_resources.gjf", mime="text/plain")

            with col3:
                if st.button("🧪 Alternate Method"):
                    retry_prompt = fix_prompt + "\n\nTry using a different method or functional to improve stability (e.g., switch from B3LYP to M06)."
                    refined = call_groq(retry_prompt, FIX_MODEL)
                    if refined:
                        st.subheader("🔧 Alternate Method .gjf")
                        st.code(refined, language="text")
                        st.download_button("💾 Download Alternate .gjf", refined, file_name="alternate_method.gjf", mime="text/plain")

# === Manual Fallback ===
st.divider()
st.subheader("💬 No files? Just describe your problem:")
manual = st.text_area("Describe issue or paste partial input")

if st.button("🧠 Manual Explain"):
    if manual.strip():
        with st.spinner("Thinking..."):
            manual_prompt = f"""You're a Gaussian troubleshooting assistant.\n\nAnalyze the following text. Format response as:\n\n### 🔍 Problem\n### ❓ Why It Happens\n### 🛠 How to Fix\n\nIf it's not related to Gaussian or computational chemistry, respond:\n\"This tool is for Gaussian input/log troubleshooting. Please provide a chemistry-related issue.\"\n\n-- User message --\n{manual}"""
            result = call_groq(manual_prompt, EXPLAIN_MODEL)
            st.markdown(result)
    else:
        st.warning("Type something first, bruh.")
