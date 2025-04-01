# Add post-fix button options
st.subheader("🔁 Refine or Retry Fix")
col1, col2, col3 = st.columns(3)

with col1:
    if st.button("🔬 Tighten Optimization"):
        refine_prompt = fix_prompt + "\n\nAlso tighten the geometry optimization settings if possible."
        refined = call_groq(refine_prompt, FIX_MODEL)
        if refined:
            st.subheader("🔧 Refined .gjf")
            st.code(refined, language="text")
            st.download_button("💾 Download Refined .gjf", refined, file_name="refined.gjf")

with col2:
    if st.button("🚀 Increase Resources"):
        memory_prompt = fix_prompt + "\n\nAlso increase memory and processor usage if it's too low."
        refined = call_groq(memory_prompt, FIX_MODEL)
        if refined:
            st.subheader("💾 Fixed with More Resources")
            st.code(refined, language="text")
            st.download_button("📄 Download Fixed .gjf", refined, file_name="more_resources.gjf")

with col3:
    if st.button("🔁 Retry with Different Functional"):
        retry_prompt = fix_prompt + "\n\nTry using a different functional or method if it looks unstable."
        refined = call_groq(retry_prompt, FIX_MODEL)
        if refined:
            st.subheader("🧪 Alternate Fix (.gjf)")
            st.code(refined, language="text")
            st.download_button("📄 Download Alternate .gjf", refined, file_name="alternate_method.gjf")
