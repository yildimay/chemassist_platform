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

    st.text_area("Preview MDP", mdp_content, height=300)
    st.download_button("📥 Download nvt.mdp", mdp_content, file_name="nvt.mdp", mime="text/plain")


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

    st.text_area("Preview MDP", mdp_content, height=300)
    st.download_button("📥 Download npt.mdp", mdp_content, file_name="npt.mdp", mime="text/plain")


def generate_md_mdp():
    st.subheader("Production MD (.mdp)")

    nsteps = st.number_input("Number of Steps", value=2500000, step=100000)  # ~5 ns at 2 fs
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

    st.text_area("Preview MDP", mdp_content, height=300)
    st.download_button("📥 Download md.mdp", mdp_content, file_name="md.mdp", mime="text/plain")
