import streamlit as st
import requests
import difflib
import base64
from PIL import Image
import os

# 🧠 AI model toggle (sidebar)
st.sidebar.subheader("🧠 AI Model Settings")
use_gpt4 = st.sidebar.checkbox("Use GPT-4 (OpenAI)", value=False)
openai_api_key = ""
if use_gpt4:
    openai_api_key = st.sidebar.text_input("Enter your OpenAI API Key", type="password")

# 🌐 API settings
GROQ_API_URL = "https://api.groq.com/openai/v1/chat/completions"
EXPLAIN_MODEL = "llama3-8b-8192"
FIX_MODEL = "llama3-70b-8192"
GROQ_API_KEY = os.environ.get("GROQ_API_KEY", "")

st.set_page_config(page_title="Gaussian Error Assistant", layout="centered")

# === Init session ===
if "fix_prompt" not in st.session_state:
    st.session_state.fix_prompt = ""

# === API Call Function ===
def call_model(prompt, model="groq-explain"):
    if use_gpt4 and openai_api_key:
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
        response = requests.post("https://api.openai.com/v1/chat/completions", headers=headers, json=payload)
        try:
            return response.json()["choices"][0]["message"]["content"]
        except:
            return "[❌ Error from OpenAI GPT-4]"
    else:
        llama_model = EXPLAIN_MODEL if model == "groq-explain" else FIX_MODEL
        headers = {
            "Authorization": f"Bearer {GROQ_API_KEY}",
            "Content-Type": "application/json"
        }
        data = {
            "model": llama_model,
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

def read_uploaded_file(file):
    return file.read().decode("utf-8", errors="ignore")

def generate_diff(original, fixed):
    return "\n".join(difflib.unified_diff(
        original.strip().splitlines(),
        fixed.strip().splitlines(),
        fromfile="original.gjf", tofile="fixed.gjf", lineterm=""
    ))

# === OCR + Cleaner ===
def extract_text_ocr_space(image_file):
    api_key = "helloworld"
    base64_img = base64.b64encode(image_file.read()).decode("utf-8")
    payload = {
        "base64Image": f"data:image/png;base64,{base64_img}",
        "language": "eng",
        "isOverlayRequired": False,
    }
    headers = {"apikey": api_key}
    res = requests.post("https://api.ocr.space/parse/image", data=payload, headers=headers)
    try:
        return res.json()["ParsedResults"][0]["ParsedText"]
    except:
        return "[OCR failed to extract text]"

def clean_ocr_text(text):
    junk = ["JOBID", "SON", "@", "Terminal", "Tabs", "Help", "File", "Edit", "View"]
    lines = text.splitlines()
    return "\n".join([l for l in lines if not any(j in l for j in junk)]).strip()

# === Interface ===
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
        explain_prompt = f"""You're a Gaussian expert.

Analyze this input and log file.

### Problem
### Likely Cause
### Recommended Solution

-- .gjf file --
{gjf_content}

-- .log file --
{log_content}
"""
        explanation = call_model(explain_prompt, model="groq-explain")
        st.subheader("📘 Explanation")
        st.markdown(explanation)

        fix_prompt = f"""Fix the broken .gjf using the .log file.
Return only the corrected .gjf file.

-- .gjf --
{gjf_content}

-- .log --
{log_content}
"""
        st.session_state.fix_prompt = fix_prompt
        fixed_gjf = call_model(fix_prompt, model="groq-fix")
        if fixed_gjf:
            st.subheader("✅ Fixed .gjf")
            st.code(fixed_gjf, language="text")
            st.download_button("Download Fixed .gjf", fixed_gjf, file_name="fixed.gjf")
            st.subheader("🔍 Difference")
            st.code(generate_diff(gjf_content, fixed_gjf), language="diff")

# === Retry Buttons ===
if st.session_state.fix_prompt:
    st.subheader("🔁 Refine or Retry")
    col1, col2, col3 = st.columns(3)
    with col1:
        if st.button("Tighten Optimization"):
            refined = call_model(st.session_state.fix_prompt + "\nUse tighter convergence criteria.", model="groq-fix")
            if refined: st.code(refined)
    with col2:
        if st.button("Increase Resources"):
            refined = call_model(st.session_state.fix_prompt + "\nIncrease memory and processor count.", model="groq-fix")
            if refined: st.code(refined)
    with col3:
        if st.button("Alternate Method"):
            refined = call_model(st.session_state.fix_prompt + "\nTry a different method or functional.", model="groq-fix")
            if refined: st.code(refined)

# === Manual Entry ===
st.subheader("✍️ Manual Error Entry")
manual_input = st.text_area("Paste Gaussian error text here:")
if st.button("Analyze Manual Text"):
    if manual_input.strip():
        with st.spinner("Analyzing..."):
            result = call_model(f"Analyze the following Gaussian output:\n{manual_input}", model="groq-explain")
            st.markdown(result)

# === Image OCR Analysis ===
st.subheader("🖼️ Image-Based Error Analysis")
img = st.file_uploader("Upload Screenshot (.png/.jpg)", type=["png", "jpg", "jpeg"])
if img and st.button("Analyze Image Error"):
    with st.spinner("Extracting and analyzing..."):
        raw = extract_text_ocr_space(img)
        cleaned = clean_ocr_text(raw)
        st.image(img, caption="Uploaded Screenshot", use_container_width=True)
        st.code(cleaned)
        analysis = call_model(f"Analyze this extracted Gaussian error:\n{cleaned}", model="groq-explain")
        st.markdown(analysis)
