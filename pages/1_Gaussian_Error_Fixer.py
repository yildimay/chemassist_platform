import streamlit as st
import requests
import difflib

st.title("⚛️ Gaussian Error Fixer — Free & Anonymous")

# GROQ setup
GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
EXPLAIN_MODEL = "llama3-8b-8192"
FIX_MODEL = "llama3-70b-8192"
GROQ_API_KEY = st.secrets.get("GROQ_API_KEY", "")

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
        st.error(f"❌ Error contacting model: {e}")
        return None

def generate_diff(original, fixed):
    original_lines = original.strip().splitlines()
    fixed_lines = fixed.strip().splitlines()
    diff = difflib.unified_diff(original_lines, fixed_lines, lineterm="", fromfile="original.gjf", tofile="fixed.gjf")
    return "\n".join(diff)

# Built-in sample
test_gjf = """%chk=broken.chk

Gaussian Broken Job

0 1
C       0.000000    0.000000    0.000000
H       0.000000    0.000000    1.089000
H       1.026719    0.000000   -0.363000
H      -0.513360   -0.889165   -0.363000
H      -0.513360    0.889165   -0.363000
"""

test_log = """ Initial command:
 %chk=broken.chk
 #p B3LYP/6-31G(d) opt freq

 Gaussian Test Job

 0 1
 C 0.000000 0.000000 0.000000
 H 0.000000 0.000000 1.089000
 H 1.026719 0.000000 -0.363000
 H -0.513360 -0.889165 -0.363000
 H -0.513360 0.889165 -0.363000

 Optimization failed to converge
 Displacement too small, but forces too large
 Error termination via Lnk1e in /g16/l101.exe
"""

# UI
st.subheader("🧪 Upload Files")
use_test_mode = st.checkbox("Use built-in test files")
gjf_content = ""
log_content = ""

if use_test_mode:
    gjf_content = test_gjf
    log_content = test_log
else:
    gjf_file = st.file_uploader("Upload broken .gjf file", type=["gjf", "com"])
    log_file = st.file_uploader("Upload .log or .out file", type=["log", "out"])
    if gjf_file and log_file:
        gjf_content = read_uploaded_file(gjf_file)
        log_content = read_uploaded_file(log_file)

# Main button
if st.button("🔍 Analyze & Fix"):
    if not gjf_content or not log_content:
        st.warning("Missing files. Scroll down to enter your problem manually.")
        st.stop()

    with st.spinner("Analyzing the problem..."):
        explain_prompt = f"""You're a Gaussian expert.

Analyze the following input and log file.

Separate your response into:

### 🔍 Problem
### ❓ Why It Happens
### 🛠 How to Fix

-- .gjf file --
{gjf_content}

-- .log file --
{log_content}
"""
        explanation = call_groq(explain_prompt, EXPLAIN_MODEL)
        if explanation:
            st.subheader("📘 Explanation & Fix Suggestion")
            st.markdown(explanation)
            st.download_button("📄 Download Explanation", explanation, file_name="explanation.txt", mime="text/plain")

        fix_prompt = f"""You are an expert in Gaussian input files.

Fix the broken .gjf file using info from the log. Output only the corrected .gjf.

-- .gjf file --
{gjf_content}

-- .log file --
{log_content}
"""
        fixed_gjf = call_groq(fix_prompt, FIX_MODEL)
        if fixed_gjf:
            st.subheader("✅ Fixed .gjf File")
            st.code(fixed_gjf, language="text")
            st.download_button("💾 Download Fixed .gjf", fixed_gjf, file_name="fixed_input.gjf", mime="text/plain")

            diff = generate_diff(gjf_content, fixed_gjf)
            st.subheader("🔍 Difference")
            st.code(diff, language="diff")
        else:
            st.info("Couldn't generate fixed .gjf.")

st.divider()
st.subheader("💬 No files? Describe your issue here:")
manual_prompt = st.text_area("Type your problem or paste part of a Gaussian input file")

if st.button("🧠 Explain Manually"):
    if not manual_prompt.strip():
        st.warning("Please enter something first.")
    else:
        with st.spinner("Thinking..."):
            manual_wrap_prompt = f"""You're a professional Gaussian and quantum chemistry troubleshooting assistant.

A user is describing a problem they're having with a Gaussian input file, a quantum chemistry job, or an error message.

Analyze their message as if it's part of a Gaussian calculation issue.

Respond with:

### 🔍 Problem  
(Explain what might be wrong)  

### ❓ Why It Happens  
(Possible cause)  

### 🛠 How to Fix  
(How to solve or debug it)

-- User message --
{manual_prompt}
"""

result = call_groq(manual_wrap_prompt, EXPLAIN_MODEL)manual_wrap_prompt = f"""You're a professional Gaussian and quantum chemistry troubleshooting assistant.

A user is describing a problem they're having with a Gaussian input file, a quantum chemistry job, or an error message.

Analyze their message as if it's part of a Gaussian calculation issue.

Respond with:

### 🔍 Problem  
(Explain what might be wrong)  

### ❓ Why It Happens  
(Possible cause)  

### 🛠 How to Fix  
(How to solve or debug it)

-- User message --
{manual_prompt}
"""

manual_wrap_prompt = f"""You're a professional Gaussian and quantum chemistry troubleshooting assistant.

A user is describing a problem they're having with a Gaussian input file, a quantum chemistry job, or an error message.

Analyze their message as if it's part of a Gaussian calculation issue.

Respond with:

### 🔍 Problem  
(Explain what might be wrong)  

### ❓ Why It Happens  
(Possible cause)  

### 🛠 How to Fix  
(How to solve or debug it)

-- User message --
{manual_prompt}
"""

result = call_groq(manual_wrap_prompt, EXPLAIN_MODEL)
            if result:
                st.subheader("📘 Explanation")
                st.markdown(result)
                st.download_button("📄 Download Explanation", result, file_name="manual_explanation.txt", mime="text/plain")
            else:
                st.error("Couldn’t analyze the input.")
