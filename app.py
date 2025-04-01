import streamlit as st
import requests
import difflib
from PIL import Image
import base64

st.set_page_config(page_title="Gaussian Error Assistant", layout="centered")

if "fix_prompt" not in st.session_state:
    st.session_state.fix_prompt = ""

# === GROQ API CONFIG ===
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

def extract_text_ocr_space(image_file):
    api_key = "helloworld"  # Free OCR.space API key
    image_bytes = image_file.read()
    base64_image = base64.b64encode(image_bytes).decode("utf-8")
    payload = {
        "base64Image": f"data:image/png;base64,{base64_image}",
        "language": "eng",
        "isOverlayRequired": False,
    }
    headers = {"apikey": api_key}
    response = requests.post("https://api.ocr.space/parse/image", data=payload, headers=headers)
    try:
        return response.json()["ParsedResults"][0]["ParsedText"]
    except:
        return "[OCR failed to extract text]"

# === UI ===
st.title("Gaussian Error Assistant")

st.subheader("File Upload")
use_test = st.checkbox("Use built-in sample files")
gjf_content, log_content = "", ""

if use_test:
    gjf_content = "%chk=broken.chk\n\nGaussian Broken Job\n\n0 1\nC       0.000000    0.000000    0.000000\nH       0.000000    0.000000    1.089000\nH       1.026719    0.000000   -0.363000\nH      -0.513360   -0.889165   -0.363000\nH      -0.513360    0.889165   -0.363000"
    log_content = """ Optimization failed to converge\n Displacement too small, but forces too large\n Error termination via Lnk1e in /g16/l101.exe\n"""
else:
    gjf_file = st.file_uploader("Upload .gjf file", type=["gjf", "com"])
    log_file = st.file_uploader("Upload .log or .out file", type=["log", "out"])
    if gjf_file: gjf_content = read_uploaded_file(gjf_file)
    if log_file: log_content = read_uploaded_file(log_file)

if st.button("Analyze & Fix"):
    if not gjf_content or not log_content:
        st.warning("Missing input files.")
        st.stop()

    with st.spinner("Analyzing..."):
        explain_prompt = f"""You're a Gaussian expert.\n\nAnalyze the input and log file. Respond with:\n\n### Problem\n### Likely Cause\n### Recommended Solution\n\n-- .gjf --\n{gjf_content}\n\n-- .log --\n{log_content}"""
        explanation = call_groq(explain_prompt, EXPLAIN_MODEL)
        st.subheader("Explanation")
        st.markdown(explanation)
        st.download_button("Download Explanation", explanation, file_name="explanation.txt")

        fix_prompt = f"""Fix this broken Gaussian .gjf using the .log file.\nReturn only the fixed .gjf file.\n\n-- .gjf --\n{gjf_content}\n\n-- .log --\n{log_content}"""
        st.session_state.fix_prompt = fix_prompt
        fixed_gjf = call_groq(fix_prompt, FIX_MODEL)
        if fixed_gjf:
            st.subheader("Fixed .gjf File")
            st.code(fixed_gjf, language="text")
            st.download_button("Download .gjf", fixed_gjf, file_name="fixed.gjf")
            st.subheader("Diff vs Original")
            st.code(generate_diff(gjf_content, fixed_gjf), language="diff")

# === Retry / Refine ===
if st.session_state.fix_prompt:
    st.subheader("Refine or Retry")
    col1, col2, col3 = st.columns(3)
    with col1:
        if st.button("Tighten Optimization"):
            prompt = st.session_state.fix_prompt + "\n\nUse tighter convergence criteria."
            refined = call_groq(prompt, FIX_MODEL)
            if refined:
                st.code(refined)
                st.download_button("Download Refined .gjf", refined, file_name="tight_opt.gjf")
    with col2:
        if st.button("Increase Resources"):
            prompt = st.session_state.fix_prompt + "\n\nIncrease memory and processor count."
            refined = call_groq(prompt, FIX_MODEL)
            if refined:
                st.code(refined)
                st.download_button("Download .gjf", refined, file_name="more_mem.gjf")
    with col3:
        if st.button("Alternate Method"):
            prompt = st.session_state.fix_prompt + "\n\nTry a different method or functional."
            refined = call_groq(prompt, FIX_MODEL)
            if refined:
                st.code(refined)
                st.download_button("Download .gjf", refined, file_name="alternate_method.gjf")

# === Manual ===
st.subheader("Manual Text Input")
man_text = st.text_area("Describe the problem or paste error text")
if st.button("Analyze Manual Text"):
    if man_text.strip():
        with st.spinner("Thinking..."):
            manual_prompt = f"You're a Gaussian assistant. Analyze the following:\n\n{man_text}"
            result = call_groq(manual_prompt, EXPLAIN_MODEL)
            st.markdown(result)

# === Image OCR ===
st.subheader("Screenshot Error Analysis")
image_file = st.file_uploader("Upload screenshot (.png/.jpg)", type=["png", "jpg", "jpeg"])
if image_file and st.button("Analyze Image Error"):
    with st.spinner("Extracting text & analyzing..."):
        ocr_text = extract_text_ocr_space(image_file)
        st.image(image_file, caption="Uploaded Screenshot", use_column_width=True)
        st.code(ocr_text, language="text")
        prompt = f"Analyze this Gaussian error text extracted from a screenshot:\n\n{ocr_text}"
        result = call_groq(prompt, EXPLAIN_MODEL)
        st.markdown(result)
