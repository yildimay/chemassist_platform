Gaussian Refine Retry
92
93
94
95
96
97
98
99
100
101
102
103
104
105
106
107
108
109
110
111
112
113
114
115
116
117
118
119
120
121
122
123
124
125
126
127
128
129
130
131
132
133
134
135
136
137
138
139
140
141
142
143
144
145
146
147
148
149
150
151
152
153
154
155
156
157
158
159
160
161
162
163
164
165
166
167
168
169
170
171
172
173
174
175
176
177
test_log = """ Optimization failed to converge
 Displacement too small, but forces too large
 Error termination via Lnk1e in /g16/l101.exe
"""

# === Upload UI ===
st.subheader("\ud83e\uddea Upload Files")
use_test = st.checkbox("Use built-in test files")
gjf_content, log_content = "", ""

if use_test:
    gjf_content = test_gjf
    log_content = test_log
else:
    gjf_file = st.file_uploader("Upload .gjf file", type=["gjf", "com"])
    log_file = st.file_uploader("Upload .log or .out file", type=["log", "out"])
    if gjf_file: gjf_content = read_uploaded_file(gjf_file)
    if log_file: log_content = read_uploaded_file(log_file)

# === Analyze & Fix ===
if st.button("\ud83d\udd0d Analyze & Fix"):
    if not gjf_content or not log_content:
        st.warning("Missing files. Scroll down to enter your problem manually.")
        st.stop()

    with st.spinner("Analyzing..."):
        explain_prompt = f"""You're a Gaussian expert.\n\nAnalyze the input and log file. Format your answer as:\n\n### \ud83d\udd0d Problem\n### \u2753 Why It Happens\n### \ud83d\udee0 How to Fix\n\n-- .gjf --\n{gjf_content}\n\n-- .log --\n{log_content}"""
        explanation = call_groq(explain_prompt, EXPLAIN_MODEL)
        st.subheader("\ud83d\udcd8 Explanation")
        st.markdown(explanation)
        st.download_button("\ud83d\udcc4 Download Explanation", explanation, file_name="explanation.txt", mime="text/plain")

        fix_prompt = f"""Fix the broken Gaussian .gjf using the .log info.\nOnly return a valid, fixed .gjf file.\n\n-- .gjf --\n{gjf_content}\n\n-- .log --\n{log_content}"""
        fixed_gjf = call_groq(fix_prompt, FIX_MODEL)
        if fixed_gjf:
            st.subheader("\u2705 Fixed .gjf File")
            st.code(fixed_gjf, language="text")
            st.download_button("\ud83d\udcc4 Download .gjf", fixed_gjf, file_name="fixed.gjf", mime="text/plain")
            st.subheader("\ud83d\udd0d Difference")
            st.code(generate_diff(gjf_content, fixed_gjf), language="diff")

            # === Retry / Refine Options ===
            st.subheader("\ud83d\udd01 Refine or Retry Fix")
            col1, col2, col3 = st.columns(3)

            with col1:
                if st.button("\ud83d\udd2c Tighten Optimization"):
                    refine_prompt = fix_prompt + "\n\nPlease tighten geometry optimization settings if possible (e.g., tighter convergence or opt=tight)."
                    refined = call_groq(refine_prompt, FIX_MODEL)
                    if refined:
                        st.subheader("\ud83d\udd27 Refined .gjf (Tighter Opt)")
                        st.code(refined, language="text")
                        st.download_button("\ud83d\udcc4 Download Refined .gjf", refined, file_name="refined_opt.gjf", mime="text/plain")

            with col2:
                if st.button("\ud83d\ude80 Increase Memory/CPUs"):
                    memory_prompt = fix_prompt + "\n\nPlease increase memory and processor count if current job settings are low."
                    refined = call_groq(memory_prompt, FIX_MODEL)
                    if refined:
                        st.subheader("\ud83d\udd27 More Resources .gjf")
                        st.code(refined, language="text")
                        st.download_button("\ud83d\udcc4 Download More Resources .gjf", refined, file_name="more_resources.gjf", mime="text/plain")

            with col3:
                if st.button("\ud83e\uddea Alternate Method"):
                    retry_prompt = fix_prompt + "\n\nTry using a different method or functional to improve stability (e.g., switch from B3LYP to M06)."
                    refined = call_groq(retry_prompt, FIX_MODEL)
                    if refined:
                        st.subheader("\ud83d\udd27 Alternate Method .gjf")
                        st.code(refined, language="text")
                        st.download_button("\ud83d\udcc4 Download Alternate .gjf", refined, file_name="alternate_method.gjf", mime="text/plain")

# === Manual Fallback ===
st.divider()
st.subheader("\ud83d\udcac No files? Just describe your problem:")
manual = st.text_area("Describe issue or paste partial input")

if st.button("\ud83e\udde0 Manual Explain"):
    if manual.strip():
        with st.spinner("Thinking..."):
            manual_prompt = f"""You're a Gaussian troubleshooting assistant.\n\nAnalyze the following text. Format response as:\n\n### \ud83d\udd0d Problem\n### \u2753 Why It Happens\n### \ud83d\udee0 How to Fix\n\nIf it's not related to Gaussian or computational chemistry, respond:\n\"This tool is for Gaussian input/log troubleshooting. Please provide a chemistry-related issue.\"\n\n-- User message --\n{manual}"""
            result = call_groq(manual_prompt, EXPLAIN_MODEL)
            st.markdown(result)
    else:
        st.warning("Type something first, bruh.")


