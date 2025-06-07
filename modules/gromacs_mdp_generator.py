import streamlit as st

def generate_em_mdp():
    st.subheader("Energy Minimization (.mdp)")
    
    emtol = st.number_input("EM Tolerance (emtol)", value=1000.0)
    emstep = st.number_input("EM Step Size (emstep)", value=0.01)
    nsteps = st.number_input("Number of Steps (nsteps)", value=50000, step=1000)
    cutoff_scheme = st.selectbox("Cutoff Scheme", ["Verlet", "group"], index=0)
    coulombtype = st.selectbox("Coulomb Type", ["PME", "Cut-off", "Reaction-Field"], index=0)
    rcoulomb = st.number_input("Coulomb Radius (rcoulomb)", value=1.0)
    rvdw = st.number_input("VDW Radius (rvdw)", value=1.0)

    mdp_content = f"""integrator    = steep
emtol         = {emtol}
emstep        = {emstep}
nsteps        = {nsteps}
cutoff-scheme = {cutoff_scheme}
coulombtype   = {coulombtype}
rcoulomb      = {rcoulomb}
rvdw          = {rvdw}
pbc           = xyz
"""

    st.text_area("Preview MDP", mdp_content, height=250)

    st.download_button(
        label="📥 Download em.mdp",
        data=mdp_content,
        file_name="em.mdp",
        mime="text/plain"
    )
