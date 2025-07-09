# gromacs_fixer_ui.py

def gromacs_fixer_ai_ui():
    import streamlit as st
    from PIL import Image
    import pytesseract
    import requests
    import os

    st.header("🧠 GROMACS Error Fixer Chat")

    st.markdown("""
    Welcome to the **GROMACS Error Fixer Chat**.
    - Describe your GROMACS issue in the chatbox below.
    - Optionally, upload a screenshot of the error.
    - Optionally, upload GROMACS input (`.gjf`) and/or output (`.log` / `.out`) files.

    **Note**: Text input is required. Other inputs are optional.
    """)

    # Initialize chat history
    if "messages" not in st.session_state:
        st.session_state.messages = []

    # Display chat history
    for msg in st.session_state.messages:
        with st.chat_message(msg["role"]):
            st.markdown(msg["content"])

    # User input and uploads
    with st.form("chat_form"):
        user_input = st.text_area("Describe your GROMACS error:", height=100, placeholder="e.g. Link 9999 or SCF not converging...")
        uploaded_image = st.file_uploader("Optional: Upload error screenshot", type=["png", "jpg", "jpeg"])
        uploaded_input_file = st.file_uploader("Optional: Upload GROMACS Input File (.gjf)", type=["gjf"])
        uploaded_output_file = st.file_uploader("Optional: Upload GROMACS Output File (.log, .out)", type=["log", "out"])
        submitted = st.form_submit_button("Send")

    if submitted:
        if user_input.strip() == "":
            st.warning("Text input is required.")
        else:
            # Save user message
            st.session_state.messages.append({"role": "user", "content": user_input})
            with st.chat_message("user"):
                st.markdown(user_input)

            # Process image if available
            image_text = ""
            if uploaded_image is not None:
                image = Image.open(uploaded_image)
                image_text = pytesseract.image_to_string(image)
                st.markdown("*Extracted from image:* \n" + image_text)

            # Read file contents
            input_file_text = ""
            output_file_text = ""
            if uploaded_input_file is not None:
                input_file_text = uploaded_input_file.read().decode("utf-8", errors="ignore")
            if uploaded_output_file is not None:
                output_file_text = uploaded_output_file.read().decode("utf-8", errors="ignore")

            # Combine and process all inputs
            final_input = f"USER MESSAGE:\n{user_input}\n\nIMAGE TEXT:\n{image_text}\n\nINPUT FILE:\n{input_file_text}\n\nOUTPUT FILE:\n{output_file_text}"

            # Call Groq API
            try:
                headers = {
                    "Authorization": f"Bearer {os.getenv('GROQ_API_KEY')}",
                    "Content-Type": "application/json"
                }

                data = {
                    "model": "llama3-70b-8192",
                    "messages": [
                        {"role": "system", "content": "You are a helpful assistant for GROMACS error fixing. If the question is not related to GROMACS software or GROMACS errors, reply with: 'This tool is specifically for GROMACS error fixing. Please ask a GROMACS-related question.'"},
                        {"role": "user", "content": final_input}
                    ]
                }

                response = requests.post("https://api.groq.com/openai/v1/chat/completions", headers=headers, json=data)
                response.raise_for_status()
                ai_response = response.json()['choices'][0]['message']['content']
            except Exception as e:
                ai_response = f"❌ Error while contacting Groq API: {e}"

            st.session_state.messages.append({"role": "assistant", "content": ai_response})
            with st.chat_message("assistant"):
                st.markdown(ai_response)
