import streamlit as st
import requests
import difflib

# -------------------🔥 STYLES 🔥-------------------
st.set_page_config(page_title="Gaussian Fixer", layout="centered")
st.markdown("""
<style>
/* Whole app background */
.stApp {
    background-color: #0f0f0f;
    padding: 2rem;
}

/* App title */
h1 {
    color: #5B8DEE;
    font-size: 2.4rem;
    font-weight: 800;
    margin-bottom: 1.2rem;
}

/* Subheaders */
h2, h3 {
    color: #ffffff;
    margin-top: 1.5rem;
    font-weight: 600;
}

/* Markdown */
div[data-testid="stMarkdownContainer"] p {
    font-size: 1.05rem;
    color: #cccccc;
}

/* Buttons */
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

/* File upload */
.css-1n76uvr {
    color: #dddddd;
}

/* Inputs */
.stTextInput input, textarea {
    background-color: #1f1f1f;
    color: #eee;
    border-radius: 8px;
}

/* Hide hamburger and footer */
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# -------------------🧠 LOGIC -------------------

st.markdown("<h1>🚀 Gaussian Error Fixer <span style='color:#aaa'>— Free. Anonymous. Clean AF.</span></h1>", unsafe_allow_html=True)

# Config
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

# Test data
test_gjf = """%chk=broken.chk

Gaussian Broken Job

0 1
C       0.000000    0.000000    0.000000
H       0.000000    0.000000    1.089000
H       1.026719    0.000000   -0.363000
H      -0.513360   -0.889165   -0.363000
H      -0.513360    0.889165   -0.363000
"""

test_log = """ Optimization failed to converge
 Displacement too small, but forces too large
 Error termination via Lnk1e in /g16/l101.exe
"""

# Upload UI
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

if st.button("🔍 Analyze & Fix"):
    if not gjf_content or not log_content:
        st.warning("Missing files. Scroll down to enter your problem manually.")
        st.stop()

    with st.spinner("Analyzing..."):
        explain_prompt = f"""You're a Gaussian expert.

Analyze the input and log file. Format your answer as:

### 🔍 Problem
### ❓ Why It Happens
### 🛠 How to Fix

-- .gjf --
{gjf_content}

-- .log --
{log_content}
"""
        explanation = call_groq(explain_prompt, EXPLAIN_MODEL)
        st.subheader("📘 Explanation")
        st.markdown(explanation)
        st.download_button("📄 Download Explanation", explanation, file_name="explanation.txt", mime="text/plain")

        fix_prompt = f"""Fix the broken Gaussian .gjf using the .log info.
Only return a valid, fixed .gjf file.

-- .gjf --
{gjf_content}

-- .log --
{log_content}
"""
        fixed_gjf = call_groq(fix_prompt, FIX_MODEL)
        if fixed_gjf:
            st.subheader("✅ Fixed .gjf File")
            st.code(fixed_gjf, language="text")
            st.download_button("💾 Download .gjf", fixed_gjf, file_name="fixed.gjf", mime="text/plain")
            st.subheader("🔍 Diff")
            st.code(generate_diff(gjf_content, fixed_gjf), language="diff")

# Manual fallback
st.divider()
st.subheader("💬 No files? Just describe your problem:")
manual = st.text_area("Describe issue or paste partial input")

if st.button("🧠 Manual Explain"):
    if manual.strip():
        with st.spinner("Thinking..."):
            manual_prompt = f"""You're a Gaussian troubleshooting assistant.

Analyze the following text. Format response as:

### 🔍 Problem
### ❓ Why It Happens
### 🛠 How to Fix

If it's not related to Gaussian or computational chemistry, respond:
"This tool is for Gaussian input/log troubleshooting. Please provide a chemistry-related issue."

-- User message --
{manual}
"""
            result = call_groq(manual_prompt, EXPLAIN_MODEL)
            st.markdown(result)
    else:
        st.warning("Type something first, bruh.")
