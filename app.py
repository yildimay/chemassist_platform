import streamlit as st
import requests
import difflib
import pytesseract
from PIL import Image
import pytesseract

st.set_page_config(page_title="Gaussian Error Assistant", layout="centered")

# === Init session state for prompts ===
if "fix_prompt" not in st.session_state:
    st.session_state.fix_prompt = ""

# === Custom CSS ===
st.markdown("""
<style>
.stApp {
    background-color: #f8f9fa;
    padding: 2rem;
    font-family: 'Segoe UI', sans-serif;
}
h1 {
    color: #2c3e50;
    font-size: 2.2rem;
    font-weight: 700;
    margin-bottom: 1rem;
}
h2, h3 {
    color: #2c3e50;
    margin-top: 1.5rem;
    font-weight: 600;
}
div[data-testid="stMarkdownContainer"] p {
    font-size: 1rem;
    color: #2c3e50;
}
.stButton > button {
    background-color: #0069d9;
    color: white;
    padding: 0.5rem 1rem;
    border-radius: 5px;
    font-weight: 500;
    border: none;
    transition: 0.3s ease;
}
.stButton > button:hover {
    background-color: #0056b3;
    transform: scale(1.01);
}
.stTextInput input, textarea {
    background-color: #ffffff;
    color: #2c3e50;
    border-radius: 6px;
    border: 1px solid #ced4da;
}
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
header {visibility: hidden;}
</style>
""", unsafe_allow_html=True)

# === Title ===
st.markdown("""<h1>Gaussian Error Assistant <span style='color:#6c757d'>– AI-powered diagnostic and correction tool</span></h1>""", unsafe_allow_html=True)

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

def extract_text_from_image(uploaded_file):
    image = Image.open(uploaded_file)
    return pytesseract.image_to_string(image)

# === Upload UI ===
st.subheader("File Upload")
use_test = st.checkbox("Use demonstration files for testing purposes")
gjf_content, log_content = "", ""

if use_test:
    gjf_content = "%chk=broken.chk\n\nGaussian Broken Job\n\n0 1\nC       0.000000    0.000000    0.000000\nH       0.000000    0.000000    1.089000\nH       1.026719    0.000000   -0.363000\nH      -0.513360   -0.889165   -0.363000\nH      -0.513360    0.889165   -0.363000"
    log_content = """ Optimization failed to converge\n Displacement too small, but forces too large\n Error termination via Lnk1e in /g16/l101.exe\n"""
else:
    gjf_file = st.file_uploader("Upload a .gjf file", type=["gjf", "com"])
    log_file = st.file_uploader("Upload the corresponding .log or .out file", type=["log", "out"])
    if gjf_file: gjf_content = read_uploaded_file(gjf_file)
    if log_file: log_content = read_uploaded_file(log_file)

# === Analyze & Fix ===
if st.button("Analyze & Fix"):
    if not gjf_content or not log_content:
        st.warning("Required files are missing. Please ensure both input and log files are provided.")
        st.stop()

    with st.spinner("Analyzing input and log files..."):
        explain_prompt = f"""You're an expert in Gaussian computational chemistry.\n\nAnalyze the provided input and log files. Respond with:\n\n### Problem\n### Likely Cause\n### Recommended Solution\n\n-- Input File --\n{gjf_content}\n\n-- Log File --\n{log_content}"""
        explanation = call_groq(explain_prompt, EXPLAIN_MODEL)
        st.subheader("Diagnostic Explanation")
        st.markdown(explanation)
        st.download_button("Download Explanation", explanation, file_name="explanation.txt", mime="text/plain")

        fix_prompt = f"""Based on the input and log files, correct the Gaussian input file.\nOnly return the corrected .gjf content.\n\n-- Input File --\n{gjf_content}\n\n-- Log File --\n{log_content}"""
        st.session_state.fix_prompt = fix_prompt
        fixed_gjf = call_groq(fix_prompt, FIX_MODEL)
        if fixed_gjf:
            st.subheader("Corrected Input File (.gjf)")
            st.code(fixed_gjf, language="text")
            st.download_button("Download Fixed File", fixed_gjf, file_name="fixed.gjf", mime="text/plain")
            st.subheader("Comparison: Original vs. Corrected")
            st.code(generate_diff(gjf_content, fixed_gjf), language="diff")

# === Retry / Refine Options ===
if st.session_state.fix_prompt:
    st.subheader("Optional Refinements")
    col1, col2, col3 = st.columns(3)

    with col1:
        if st.button("Tighten Optimization"):
            refine_prompt = st.session_state.fix_prompt + "\n\nPlease apply tighter convergence criteria for optimization."
            refined = call_groq(refine_prompt, FIX_MODEL)
            if refined:
                st.subheader("Tighter Optimization")
                st.code(refined, language="text")
                st.download_button("Download Refined File", refined, file_name="refined_opt.gjf", mime="text/plain")

    with col2:
        if st.button("Allocate More Resources"):
            memory_prompt = st.session_state.fix_prompt + "\n\nPlease increase memory and processor usage for the job."
            refined = call_groq(memory_prompt, FIX_MODEL)
            if refined:
                st.subheader("Enhanced Resources")
                st.code(refined, language="text")
                st.download_button("Download Modified File", refined, file_name="more_resources.gjf", mime="text/plain")

    with col3:
        if st.button("Alternate Method"):
            retry_prompt = st.session_state.fix_prompt + "\n\nConsider using a different method or functional to improve convergence."
            refined = call_groq(retry_prompt, FIX_MODEL)
            if refined:
                st.subheader("Alternate Computational Method")
                st.code(refined, language="text")
                st.download_button("Download Alternate File", refined, file_name="alternate_method.gjf", mime="text/plain")

# === Manual Fallback ===
st.divider()
st.subheader("Manual Input")
manual = st.text_area("Describe your issue or paste Gaussian input/output text")

if st.button("Analyze Manual Entry"):
    if manual.strip():
        with st.spinner("Analyzing manual entry..."):
            manual_prompt = f"""You are a Gaussian error diagnostic assistant.\n\nAnalyze the following user-provided content. Respond with:\n\n### Problem\n### Likely Cause\n### Recommended Solution\n\nIf irrelevant, kindly respond: 'This tool is designed to assist with Gaussian-related issues.'\n\n-- User Message --\n{manual}"""
            result = call_groq(manual_prompt, EXPLAIN_MODEL)
            st.markdown(result)
    else:
        st.warning("Please enter a valid message for analysis.")

# === Image OCR Upload ===
st.divider()
st.subheader("Image-Based Error Analysis")
image_file = st.file_uploader("Upload a screenshot of the Gaussian error message (PNG or JPG)", type=["png", "jpg", "jpeg"])

if image_file and st.button("Analyze Image Error"):
    with st.spinner("Extracting text from image and analyzing..."):
        ocr_text = extract_text_from_image(image_file)
        st.image(image_file, caption="Uploaded Image", use_column_width=True)
        st.code(ocr_text, language="text")

        image_prompt = f"""You are a Gaussian troubleshooting expert.\nAnalyze the extracted text from a user-provided screenshot.\n\n{ocr_text}\n\nRespond with:\n### Problem\n### Likely Cause\n### Recommended Solution\n"""
        result = call_groq(image_prompt, EXPLAIN_MODEL)
        st.markdown(result)
