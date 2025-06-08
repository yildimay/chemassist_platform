# 📁 modules/gromacs_fixer_ui.py
import streamlit as st
from PIL import Image
import pytesseract
import requests
import os

def gromacs_fixer_ui():
    st.header("🧠 GROMACS Error Fixer Chat")

    st.markdown("""
    Welcome to the **GROMACS Error Fixer Chat**.
    - Describe your GROMACS issue in the chatbox below.
    - Optionally, upload a screenshot of the error.
    - Optionally, upload relevant GROMACS files (`.mdp`, `.gro`, `.log`, `.tpr`).

    **Note**: Text input is required. Other inputs are optional.
    """)

    if "messages" not in st.session_state:
        st.session_state.messages = []

    for msg in st.session_state.messages:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])

    with st.form("gromacs_chat_form"):
        user_input = st.text_area("Describe your GROMACS error:", height=100, placeholder="e.g. LINCS WARNING, unstable pressure...")
        uploaded_image = st.file_uploader("Optional: Upload error screenshot", type=["png", "jpg", "jpeg"])
        uploaded_mdp = st.file_uploader("Optional: Upload MDP file (.mdp)", type=["mdp"])
        uploaded_log = st.file_uploader("Optional: Upload log file (.log)", type=["log"])
        uploaded_gro = st.file_uploader("Optional: Upload GRO file (.gro)", type=["gro"])
        submitted = st.form_submit_button("Send")

    if submitted:
        if user_input.strip() == "":
            st.warning("Text input is required.")
        else:
            st.session_state.messages.append({"role": "user", "content": user_input})
            with st.chat_message("user"):
                st.markdown(user_input)

            image_text = ""
            if uploaded_image:
                image = Image.open(uploaded_image)
                image_text = pytesseract.image_to_string(image)
                st.markdown("*Extracted from image:*\n" + image_text)

            mdp_text = uploaded_mdp.read().decode("utf-8", errors="ignore") if uploaded_mdp else ""
            log_text = uploaded_log.read().decode("utf-8", errors="ignore") if uploaded_log else ""
            gro_text = uploaded_gro.read().decode("utf-8", errors="ignore") if uploaded_gro else ""

            full_input = f"USER TEXT:\n{user_input}\n\nIMAGE TEXT:\n{image_text}\n\nMDP FILE:\n{mdp_text}\n\nLOG FILE:\n{log_text}\n\nGRO FILE:\n{gro_text}"

            try:
                headers = {
                    "Authorization": f"Bearer {os.getenv('GROQ_API_KEY')}",
                    "Content-Type": "application/json"
                }

                data = {
                    "model": "llama3-70b-8192",
                    "messages": [
                        {
                            "role": "system",
                            "content": """
You are a GROMACS error fixer. Only answer GROMACS-related questions.
If the message is unrelated, respond: 'This tool is for GROMACS troubleshooting only. Please provide a GROMACS-related error or file.'
"""
                        },
                        {"role": "user", "content": full_input}
                    ]
                }

                response = requests.post("https://api.groq.com/openai/v1/chat/completions", headers=headers, json=data)
                response.raise_for_status()
                ai_response = response.json()['choices'][0]['message']['content']
            except Exception as e:
                ai_response = f"\u274c Error while contacting Groq API: {e}"

            st.session_state.messages.append({"role": "assistant", "content": ai_response})
            with st.chat_message("assistant"):
                st.markdown(ai_response)
