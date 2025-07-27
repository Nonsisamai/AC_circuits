# AC/DC Educational Circuit Visualizer with Simulation, PDF Export, and Transient Analysis
# Streamlit app designed by an electrical engineer & educator

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import io
import graphviz
from math import cos, sin, radians, degrees, pi, sqrt, atan2, exp
from fpdf import FPDF

st.set_page_config(page_title="AC/DC Circuit Visualizer", layout="wide")

# Theme toggle
dark_mode = st.sidebar.toggle("🌙 Tmavý režim")
if dark_mode:
    plt.style.use('dark_background')
else:
    plt.style.use('default')

st.title("🔌 AC/DC Obvody – Vizualizácia, Schéma, Prechody a Výpočty")

st.markdown("""
Interaktívny elektroinžiniersky nástroj:
- ✅ AC/DC simulácie
- ✅ Reálne výpočty R, L, C
- ✅ Časové priebehy, výkon, fázory
- ✅ Interaktívna schéma
- ✅ Prechodové deje: **nabíjanie kondenzátora**, **prúd cievkou**
- ✅ Export výpočtov a schémy do PDF protokolu
""")

st.sidebar.header("🎛️ Parametre obvodu")

# Typ obvodu: AC, DC alebo prechodový
type_choice = st.sidebar.selectbox("Režim obvodu", ["AC", "DC", "DC - Prechodový dej (R-C / R-L)"])

# Zadanie RMS alebo peak
input_mode = st.sidebar.radio("Zadávaš hodnoty ako:", ["Efektívne (RMS)", "Maximálne (peak)"])
U_in = st.sidebar.number_input("Napätie [V]", value=230.0, step=0.1)
I_in = st.sidebar.number_input("Prúd [A]", value=5.0, step=0.1)

# Frekvencia a fázový posun
f = st.sidebar.number_input("Frekvencia [Hz]", value=50.0, step=1.0) if type_choice == "AC" else 0.0
phi_manual = st.sidebar.number_input("Fázový posun φ [°] (ak je známy)", value=0.0 if type_choice != "AC" else 30.0)
phi_manual_rad = radians(phi_manual)

# Súčiastky
st.sidebar.markdown("---")
st.sidebar.markdown("🧩 **Zadanie súčiastok**")
R = st.sidebar.number_input("Odpor R [Ω]", value=0.0, step=0.1)
L = st.sidebar.number_input("Indukčnosť L [H]", value=0.0, step=0.001)
C = st.sidebar.number_input("Kapacita C [F]", value=0.0, step=0.00001)

# Auto-zistenie typu záťaže
zataz_popis = []
if R > 0: zataz_popis.append("R")
if L > 0: zataz_popis.append("L")
if C > 0: zataz_popis.append("C")
zataz_type = "+".join(zataz_popis) if zataz_popis else "Žiadna záťaž"

omega = 2 * pi * f if f > 0 else 0
XL = omega * L if L > 0 else 0.0
XC = 1 / (omega * C) if (C > 0 and omega > 0) else 0.0
Z_complex = complex(R, XL - XC)
Z_abs = abs(Z_complex)
phi_calc_rad = atan2(Z_complex.imag, Z_complex.real) if Z_abs > 0 else 0.0
phi_calc_deg = degrees(phi_calc_rad)
cos_phi = cos(phi_calc_rad) if Z_abs > 0 else 1.0

if input_mode == "Efektívne (RMS)":
    Uef = U_in
    Ief = I_in
    Umax = Uef * sqrt(2)
    Imax = Ief * sqrt(2)
else:
    Umax = U_in
    Imax = I_in
    Uef = Umax / sqrt(2)
    Ief = Imax / sqrt(2)

S = Uef * Ief
P = S * cos_phi
Q = sqrt(abs(S ** 2 - P ** 2)) if type_choice == "AC" else 0.0

# Časová os
if type_choice == "AC":
    T = 1 / f if f > 0 else 1.0
    x = np.linspace(0, 2 * T, 1000)
    napatie = Umax * np.sin(omega * x)
    prud = Imax * np.sin(omega * x - phi_calc_rad)
elif type_choice == "DC":
    x = np.linspace(0, 0.1, 1000)
    napatie = np.full_like(x, Uef)
    prud = np.full_like(x, Ief)
else:
    # Prechodový dej (RC/RL nabíjanie/vybíjanie)
    x = np.linspace(0, 0.2, 1000)
    if C > 0 and R > 0:
        napatie = Uef * (1 - np.exp(-x / (R * C)))
        prud = (Uef / R) * np.exp(-x / (R * C))
    elif L > 0 and R > 0:
        napatie = np.full_like(x, Uef)
        prud = (Uef / R) * (1 - np.exp(-R * x / L))
    else:
        napatie = np.zeros_like(x)
        prud = np.zeros_like(x)

vykon = napatie * prud
vykon_avg = np.mean(vykon)

# Výsledky
st.subheader("📐 Výsledky")
st.markdown(f"""
- **Režim:** {type_choice}  
- **Záťaž:** {zataz_type}  
- **Z =** {Z_abs:.2f} Ω, **φ =** {phi_calc_deg:.2f}°, **cosφ =** {cos_phi:.3f}  
- **XL =** {XL:.2f} Ω, **XC =** {XC:.2f} Ω  
- **Uef =** {Uef:.2f} V, **Umax =** {Umax:.2f} V  
- **Ief =** {Ief:.2f} A, **Imax =** {Imax:.2f} A  
- **S =** {S:.2f} VA, **P =** {P:.2f} W, **Q =** {Q:.2f} VAR  
- **⟨P(t)⟩ =** {vykon_avg:.2f} W
""")

# Schéma
st.subheader("🔧 Schéma zapojenia")
g = graphviz.Digraph()
g.node("V", "Zdroj")
last = "V"
if R > 0:
    g.node("R", "R")
    g.edge(last, "R")
    last = "R"
if L > 0:
    g.node("L", "L")
    g.edge(last, "L")
    last = "L"
if C > 0:
    g.node("C", "C")
    g.edge(last, "C")
    last = "C"
g.edge(last, "Z")
g.node("Z", "Uzemnenie")
st.graphviz_chart(g)

# Grafy
st.subheader("📈 Časové priebehy")
fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
ax[0].plot(x * 1000, napatie, label='Napätie [V]', color='tab:blue')
ax[1].plot(x * 1000, prud, label='Prúd [A]', color='tab:orange')
ax[2].plot(x * 1000, vykon, label=f'Výkon [W] ⟨P⟩={vykon_avg:.2f}', color='tab:green')
for a in ax:
    a.legend()
    a.grid(True)
    a.set_ylabel("Hodnota")
ax[2].set_xlabel("Čas [ms]")
st.pyplot(fig)


# PDF export
def export_pdf():
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt="Protokol AC/DC obvodu", ln=True, align="C")
    pdf.ln(10)
    pdf.multi_cell(0, 10, txt=f"""
Režim: {type_choice}
Záťaž: {zataz_type}
Z = {Z_abs:.2f} Ω
XL = {XL:.2f} Ω
XC = {XC:.2f} Ω
φ = {phi_calc_deg:.2f}°
cosφ = {cos_phi:.3f}

Uef = {Uef:.2f} V
Ief = {Ief:.2f} A
Umax = {Umax:.2f} V
Imax = {Imax:.2f} A

Výkony:
S = {S:.2f} VA
P = {P:.2f} W
Q = {Q:.2f} VAR
⟨P⟩ = {vykon_avg:.2f} W

Výpočty podľa: Z = R + j(XL - XC), S = Uef·Ief, P = S·cosφ

Vytvorené pomocou elektroinžinierskeho nástroja v Streamlit
    """)
    return pdf.output(dest='S').encode('latin1')


st.subheader("📤 Export")
df = pd.DataFrame({
    "čas [s]": x,
    "napätie [V]": napatie,
    "prúd [A]": prud,
    "výkon [W]": vykon
})
csv = df.to_csv(index=False).encode('utf-8')
st.download_button("⬇️ CSV Export", csv, file_name="vysledky.csv", mime="text/csv")
pdf_bytes = export_pdf()
st.download_button("⬇️ PDF Protokol", data=pdf_bytes, file_name="protokol_obvodu.pdf", mime="application/pdf")
st.markdown("---")
st.markdown("👨Autor: Adrian Mahdon")