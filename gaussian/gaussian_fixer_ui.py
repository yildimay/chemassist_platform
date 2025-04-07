import streamlit as st
from PIL import Image
import pytesseract
import requests
import os

st.set_page_config(page_title="Gaussian Error Fixer Chat", layout="centered")
st.title("🧠 Gaussian Error Fixer Chat")

st.markdown("""
Welcome to the **Gaussian Error Fixer Chat**.
- Describe your Gaussian issue in the chatbox below.
- Optionally, upload a screenshot of the error.

**Note**: Text input is required. Image upload is optional.
""")

# Initialize chat history
if "messages" not in st.session_state:
    st.session_state.messages = []

# Display chat history
for msg in st.session_state.messages:
    with st.chat_message(msg["role"]):
        st.markdown(msg["content"])

# User input and image upload
with st.form("chat_form"):
    user_input = st.text_area("Describe your Gaussian error:", height=100, placeholder="e.g. Link 9999 or SCF not converging...")
    uploaded_image = st.file_uploader("Optional: Upload error screenshot", type=["png", "jpg", "jpeg"])
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

        # Combine and process both inputs
        final_input = user_input + "\n" + image_text

        # Call Groq API
        try:
            headers = {
                "Authorization": f"Bearer {os.getenv('GROQ_API_KEY')}",
                "Content-Type": "application/json"
            }

            data = {
                "model": "llama3-70b-8192",
                "messages": [
                    {"role": "system", "content": "You are a helpful assistant for fixing Gaussian errors."},
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
