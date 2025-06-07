import streamlit as st

# ───────────────────────────────────────────────
# ENERGY MINIMIZATION
# ───────────────────────────────────────────────

def generate_em_mdp():
    st.subheader("Energy Minimization (.mdp)")

    emtol = st.number_input("EM Tolerance (`emtol`)", value=1000.0)
    emstep = st.number_input("EM Step Size (`emstep`)", value=0.01)
    nsteps = st.number_input("Max Number of Steps (`nsteps`)", value=50000, step=1000)
    cutoff_scheme = st.selectbox("Cutoff Scheme", ["Verlet", "group"], index=0)
    coulombtype = st.selectbox("Coulomb Type", ["PME", "Cut-off", "Reaction-Field"], index=0)
    rcoulomb = st.number_input("Coulomb Radius (`rcoulomb`)", value=1.0)
    rvdw = st.number_input("VDW Radius (`rvdw`)", value=1.0)

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

    st.text_area("Preview MDP", mdp_content, height=280)
    st.download_button("📥 Download em.mdp", data=mdp_content, file_name="em.mdp", mime="text/plain")


# ───────────────────────────────────────────────
# NVT EQUILIBRATION
# ───────────────────────────────────────────────

def generate_nvt_mdp():
    st.subheader("NVT Equilibration (.mdp)")

    nsteps = st.number_input("Number of Steps", value=50000, step=1000)
    dt = st.number_input("Timestep (ps)", value=0.002)
    tcoupl = st.selectbox("Thermostat", ["V-rescale", "Berendsen", "Nose-Hoover"], index=0)
    ref_t = st.number_input("Temperature (K)", value=300.0)
    tau_t = st.number_input("Temperature coupling constant", value=0.1)

    mdp_content = f"""integrator    = md
nsteps        = {nsteps}
dt            = {dt}
cutoff-scheme = Verlet
nstxout       = 1000
nstvout       = 1000
nstenergy     = 1000
nstlog        = 1000
continuation  = no
tcoupl        = {tcoupl}
tc-grps       = System
tau_t         = {tau_t}
ref_t         = {ref_t}
pcoupl        = no
constraints   = h-bonds
pbc           = xyz
"""

    st.text_area("Preview MDP", mdp_content, height=280)
    st.download_button("📥 Download nvt.mdp", data=mdp_content, file_name="nvt.mdp", mime="text/plain")


# ───────────────────────────────────────────────
# NPT EQUILIBRATION
# ───────────────────────────────────────────────

def generate_npt_mdp():
    st.subheader("NPT Equilibration (.mdp)")

    nsteps = st.number_input("Number of Steps", value=50000, step=1000)
    dt = st.number_input("Timestep (ps)", value=0.002)
    tcoupl = st.selectbox("Thermostat", ["V-rescale", "Berendsen", "Nose-Hoover"], index=0)
    ref_t = st.number_input("Temperature (K)", value=300.0)
    tau_t = st.number_input("Temperature coupling constant", value=0.1)
    pcoupl = st.selectbox("Barostat", ["Parrinello-Rahman", "Berendsen", "C-rescale"], index=0)
    ref_p = st.number_input("Pressure (bar)", value=1.0)
    tau_p = st.number_input("Pressure coupling constant", value=2.0)

    mdp_content = f"""integrator    = md
nsteps        = {nsteps}
dt            = {dt}
cutoff-scheme = Verlet
nstxout       = 1000
nstvout       = 1000
nstenergy     = 1000
nstlog        = 1000
continuation  = yes
tcoupl        = {tcoupl}
tc-grps       = System
tau_t         = {tau_t}
ref_t         = {ref_t}
pcoupl        = {pcoupl}
pcoupltype    = isotropic
tau_p         = {tau_p}
ref_p         = {ref_p}
compressibility = 4.5e-5
constraints   = h-bonds
pbc           = xyz
"""

    st.text_area("Preview MDP", mdp_content, height=280)
    st.download_button("📥 Download npt.mdp", data=mdp_content, file_name="npt.mdp", mime="text/plain")


# ───────────────────────────────────────────────
# PRODUCTION MD
# ───────────────────────────────────────────────

def generate_md_mdp():
    st.subheader("Production MD (.mdp)")

    nsteps = st.number_input("Number of Steps", value=2500000, step=100000)
    dt = st.number_input("Timestep (ps)", value=0.002)
    tcoupl = st.selectbox("Thermostat", ["V-rescale", "Berendsen", "Nose-Hoover"], index=0)
    ref_t = st.number_input("Temperature (K)", value=300.0)
    pcoupl = st.selectbox("Barostat", ["Parrinello-Rahman", "Berendsen", "C-rescale"], index=0)
    ref_p = st.number_input("Pressure (bar)", value=1.0)

    mdp_content = f"""integrator    = md
nsteps        = {nsteps}
dt            = {dt}
cutoff-scheme = Verlet
nstxout       = 1000
nstvout       = 1000
nstenergy     = 1000
nstlog        = 1000
continuation  = yes
tcoupl        = {tcoupl}
tc-grps       = System
tau_t         = 0.1
ref_t         = {ref_t}
pcoupl        = {pcoupl}
pcoupltype    = isotropic
tau_p         = 2.0
ref_p         = {ref_p}
compressibility = 4.5e-5
constraints   = h-bonds
pbc           = xyz
"""

    st.text_area("Preview MDP", mdp_content, height=280)
    st.download_button("📥 Download md.mdp", data=mdp_content, file_name="md.mdp", mime="text/plain")
