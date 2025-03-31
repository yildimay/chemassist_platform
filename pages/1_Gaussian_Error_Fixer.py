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

gjf_file = st.file_uploader("Upload broken .gjf file", type=["gjf", "com"])
log_file = st.file_uploader("Upload related .log or .out file", type=["log", "out"])

if st.button("Analyze / Fix"):
    if not gjf_file or not log_file:
        st.warning("Please upload both a .gjf and a .log/.out file.")
        st.stop()

    with st.spinner("Processing..."):
        gjf_content = read_uploaded_file(gjf_file)
        log_content = read_uploaded_file(log_file)

        # Always show explanation
        explain_prompt = f"""You're a Gaussian error expert.

A user submitted this Gaussian input file and log file.

Explain:
1. What is wrong
2. Why it likely happened
3. How they can fix it manually

-- .gjf file --
{gjf_content}

-- .log file --
{log_content}
"""
        explanation = call_groq(explain_prompt, EXPLAIN_MODEL)
        st.subheader("📘 Explanation & Suggested Fix")
        st.markdown(explanation)

        # Only paid users get fixed .gjf
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
            st.subheader("✅ Fixed .gjf File")
            st.code(fixed_gjf, language="text")
            st.download_button("💾 Download Fixed .gjf", fixed_gjf, file_name="fixed_input.gjf", mime="text/plain")
        else:
            st.info("Upgrade to a paid plan to unlock .gjf file generation.")
