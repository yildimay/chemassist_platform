import streamlit as st
from PIL import Image
import pytesseract
from utils.mdp_editor import apply_mdp_fixes
from modules.groq_client import ask_groq  # your Groq-backed AI call

def extract_text_from_image(image_file):
    image = Image.open(image_file)
    return pytesseract.image_to_string(image)

def gromacs_error_fixer_ai_ui():
    st.subheader("🧠 GROMACS AI Error Fixer (with Screenshot & Auto-Fix)")

    uploaded_file = st.file_uploader("📄 Upload `.log`, `.mdp`, or `.txt` file", type=["log", "mdp", "txt"])
    uploaded_image = st.file_uploader("🕼 Or upload screenshot of error", type=["png", "jpg", "jpeg"])

    file_content = ""
    is_mdp = False

    if uploaded_file:
        file_content = uploaded_file.read().decode("utf-8", errors="ignore")
        is_mdp = uploaded_file.name.endswith(".mdp")
        st.text_area("📃 File Content", file_content, height=250)

    elif uploaded_image:
        file_content = extract_text_from_image(uploaded_image)
        st.text_area("�퓾 Extracted Text from Image", file_content, height=250)

    if file_content:
        # Guard against non-GROMACS input
        if not any(keyword in file_content.lower() for keyword in ["gromacs", "lincs", "mdp", "mdrun", "tpr", "gro"]):
            st.warning("This doesn’t appear to be a GROMACS-related file or error. Please upload a relevant log or .mdp file.")
            return

        if st.button("🤖 Analyze with AI"):
            with st.spinner("Asking Groq for help..."):
                response = ask_groq(file_content)
                st.markdown("### 💡 AI Diagnosis & Fixes")
                st.markdown(response)

                if is_mdp and "change" in response.lower():
                    if st.button("�� Apply Fix to .mdp"):
                        fixed_mdp = apply_mdp_fixes(file_content, response)
                        st.text_area("✅ Fixed .mdp", fixed_mdp, height=300)
                        st.download_button("📅 Download Fixed .mdp", fixed_mdp, file_name="fixed.mdp", mime="text/plain")
import openai  # Replace with Groq client if different

def ask_groq(log_text):
    system_prompt = """
You are a strict and specialized assistant trained only in **GROMACS troubleshooting**.
You are not allowed to answer general questions or unrelated topics.

Your job is to:
1. Read error logs, `.mdp` settings, or simulation outputs from GROMACS.
2. Identify common issues: LINCS warnings, bad timestep, unstable barostat, corrupted .tpr, etc.
3. Suggest specific parameter or setup changes based on best GROMACS practices.

🚡 Reject anything that is not related to GROMACS with this exact message:
"I'm only trained to help with GROMACS troubleshooting. Please upload a related file or error."
"""

    prompt = f"""
{system_prompt}

--- FILE CONTENT START ---
{log_text[:4000]}
--- FILE CONTENT END ---
"""

    response = openai.ChatCompletion.create(
        model="gpt-4o",  # Or your Groq LLaMA3 endpoint
        messages=[
            {"role": "user", "content": prompt}
        ],
        temperature=0.1
    )

    return response.choices[0].message.content
