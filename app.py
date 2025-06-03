import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import streamlit as st
from gaussian.gaussian_input_ui import gaussian_input_ui
from gaussian.gaussian_fixer_ui import gaussian_fixer_ui
from smiles_to_mol import smiles_ui

# Set Streamlit page configuration
st.set_page_config(page_title="ChemAssist Platform", layout="wide")

# Sidebar module selection
st.sidebar.title("ChemAssist Modules")
selected_module = st.sidebar.selectbox("Choose a tool:", ["Gaussian", "Molecule Builder", "VASP (coming soon)", "GROMACS (coming soon)"])

if selected_module == "Gaussian":
    gaussian_tool = st.sidebar.radio("Gaussian Tool", ["Error Fixer", "Input Creator"])
    if gaussian_tool == "Error Fixer":
        gaussian_fixer_ui()
    elif gaussian_tool == "Input Creator":
        gaussian_input_ui()

elif selected_module == "Molecule Builder":
    smiles_ui()

elif selected_module == "VASP (coming soon)":
    st.markdown("### VASP module is under development.")

elif selected_module == "GROMACS (coming soon)":
    st.markdown("### GROMACS module is under development.")
    
import streamlit as st
from rdkit import Chem

st.sidebar.title("Modüller")
secim = st.sidebar.radio("Hangi modülü kullanmak istiyorsun?", 
                         ["Gaussian Fixer", "DFT Yardımcısı", "Çevresel Risk Değerlendirmesi", "..."])

if secim == "Çevresel Risk Değerlendirmesi":
    st.header("Molekülün Çevresel Tehlike Paneli")

    # 1. Molekül Girişi
    smiles_input = st.text_input("Molekül SMILES girin:", value="")
    file_input = st.file_uploader("Veya .mol/.sdf dosyası yükleyin:", type=["mol", "sdf"])
    
    # SMILES ve dosya kontrolü
    mol_obj = None
    if smiles_input:
        mol_obj = Chem.MolFromSmiles(smiles_input)
    elif file_input:
        mol_obj = Chem.MolFromMolBlock(file_input.read().decode("utf-8"))
    
    if mol_obj:
        st.subheader("Molekül Görünümü")
        from rdkit.Chem.Draw import rdMolDraw2D
        drawer = rdMolDraw2D.MolDraw2DSVG(300, 200)
        rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol_obj)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace('svg:','')
        st.write(f'<div>{svg}</div>', unsafe_allow_html=True)
        
        # 2. Fizikokimyasal Özellik Hesabı
        props = hesapla_fiziksel_kimyasal(smiles_input)
        st.subheader("Fizikokimyasal Özellikler")
        st.table(props)
        
        # 3. Ekotoksisite ve Biyobirikim Tahmini
        with st.spinner("Ekotoksisite parametreleri hesaplanıyor..."):
            # Örneğin, PubChem ECOSAR sorgusu veya T.E.S.T. çağrısı
            # Örnek placeholder:
            ecosar = get_ecosar_classification(cids=[...])  # önce CID almanız lazım
            st.subheader("ECOSAR Sınıflandırması")
            if ecosar:
                st.json(ecosar)
            else:
                st.write("ECOSAR verisi bulunamadı.")
            
            # EPI Suite tahminleri
            # episuite_results = run_episuitemodule(smiles_input, module_name="BCFBAF")
            # st.write("Biyobirikim Katsayısı (BCF):", episuite_results["BCF"])
        
        # 4. Risk Skoru & Renkli Uyarı
        st.subheader("Çevresel Risk Skoru")
        # Örnek: basitçe logP > 4 ve BCF > 5000 ise “Yüksek Risk”
        # Detaylı hesaplama buraya geliyor
        risk_skoru = 0  # placeholder
        uyarı_metni = ""
        if props["logP"] > 4:
            risk_skoru += 50
            uyarı_metni = "Yüksek biyobirikim potansiyeli"
        # ... diğer kurallar eklenebilir ...
        st.metric(label="Risk Skoru (0–100)", value=f"{risk_skoru}", delta=None)
        if risk_skoru > 75:
            st.error("Bu molekül çevre için yüksek risk taşıyor!")
        elif risk_skoru > 40:
            st.warning("Orta düzeyde çevresel risk.")
        else:
            st.success("Düşük çevresel risk.")
