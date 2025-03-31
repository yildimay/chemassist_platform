import streamlit as st
import requests

st.title("⚛️ Gaussian Error Fixer + GJF Generator")

GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
EXPLAIN_MODEL = "llama3-8b-8192"
FIX_MODEL = "llama3-70b-8192"
GROQ_API_KEY = st.secrets.get("GROQ_API_KEY", "")

# Login and tier selection
st.subheader("🔐 Login")
user_email = st.text_input("Enter your email:")
is_paid = st.checkbox("I'm a paid user", value=False)

if not user_email:
    st.stop()

st.divider()
st.subheader("🧪 File Input")

use_test_mode = st.checkbox("Use built-in sample files (for testing)", value=False)

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
    except requests.exceptions.HTTPError as e:
        if res.status_code == 401:
            st.error("🔐 Unauthorized: Check your GROQ_API_KEY.")
        elif res.status_code == 429:
            st.error("🚫 Quota exceeded: You've hit the Groq usage limit. Try again later.")
        else:
            st.error(f"❌ HTTP error {res.status_code}: {res.text}")
        return None
    except Exception as e:
        st.error(f"⚠️ Unexpected error: {e}")
        return None

# Built-in test content
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

# File input
gjf_content = ""
log_content = ""

if use_test_mode:
    gjf_content = test_gjf
    log_content = test_log
else:
    gjf_file = st.file_uploader("Upload broken .gjf file", type=["gjf", "com"])
    log_file = st.file_uploader("Upload related .log or .out file", type=["log", "out"])
    if gjf_file and log_file:
        gjf_content = read_uploaded_file(gjf_file)
        log_content = read_uploaded_file(log_file)

if st.button("Analyze / Fix"):
    if not gjf_content or not log_content:
        st.warning("Missing required input files.")
        st.stop()

    with st.spinner("Processing..."):
        explain_prompt = f"""You're a Gaussian error expert.

A user submitted this Gaussian input file and log file.

Separate your response into 3 sections:

### 🔍 Problem:
(Explain what's wrong)

### ❓ Why It Happens:
(Explain the likely cause)

### 🛠 How to Fix:
(Explain how to resolve it manually)

-- .gjf file --
{gjf_content}

-- .log file --
{log_content}
"""
        explanation = call_groq(explain_prompt, EXPLAIN_MODEL)
        if explanation:
            st.subheader("📘 Explanation & Suggested Fix")
            st.markdown(explanation)
            st.download_button("📄 Download Explanation", explanation, file_name="explanation.txt", mime="text/plain")
        else:
            st.info("Something went wrong — please try again shortly.")

        if is_paid:
            fix_prompt = f"""You are an expert in Gaussian input files.

Fix the following broken Gaussian input file (.gjf) using information from the log file.
Output only the corrected .gjf file (no explanation), with proper route section, charge/multiplicity, and fixed atom coordinates.

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
            else:
                st.info("GJF generation failed — please try again later.")
        else:
            st.info("Upgrade to a paid plan to unlock .gjf file generation.")
