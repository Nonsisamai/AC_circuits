# AC/DC Educational Circuit Visualizer with Expert Logic
# Streamlit app designed by an electrical engineer & educator

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import io
from math import cos, sin, radians, pi, sqrt

# Page setup
st.set_page_config(page_title="AC/DC Circuit Visualizer", layout="wide")

# Theme toggle
dark_mode = st.sidebar.toggle("🌙 Tmavý režim")
if dark_mode:
    plt.style.use('dark_background')
else:
    plt.style.use('default')

st.title("🔌 AC/DC Obvody – Vizualizácia a Výpočty pre študentov aj odborníkov")

st.markdown("""
Tento nástroj slúži na pochopenie správania sa AC a DC obvodov s rôznymi záťažami (R, L, C, RLC):
- Zobrazí časové priebehy napätia, prúdu a výkonu
- Umožní výpočet výkonov a reaktancií
- Rozlišuje medzi AC a DC režimom
- Exportuje výsledky a vizualizácie
""")

st.sidebar.header("🎛️ Parametre obvodu")

# Typ obvodu: AC alebo DC
type_choice = st.sidebar.selectbox("Režim obvodu", ["AC", "DC"])

# Typ záťaže
if type_choice == "AC":
    zataz_type = st.sidebar.selectbox("Typ záťaže (AC)", ["R", "L", "C", "RL", "RC", "LC", "RLC"])
else:
    zataz_type = st.sidebar.selectbox("Typ záťaže (DC)", ["R", "RL (iba odpor)", "RC (iba odpor)", "LC (iba prechod)"])

# Zadanie: efektívne alebo maximálne hodnoty
input_mode = st.sidebar.radio("Zadávaš hodnoty ako:", ["Efektívne (RMS)", "Maximálne (peak)"])

# Zadávanie hodnôt
U_in = st.sidebar.number_input("Napätie [V]", value=230.0, step=0.1)
I_in = st.sidebar.number_input("Prúd [A]", value=5.0, step=0.1)
f = st.sidebar.number_input("Frekvencia [Hz]", value=50.0, step=1.0) if type_choice == "AC" else 0.0
phi_deg = st.sidebar.number_input("Fázový posun φ [°]", value=30.0 if type_choice == "AC" else 0.0)

# Hodnoty R, L, C
tab = st.sidebar.expander("⚙️ Parametre súčiastok")
with tab:
    R = st.number_input("Odpor R [Ω]", value=10.0, step=0.1)
    L = st.number_input("Indukčnosť L [H]", value=0.05, step=0.001)
    C = st.number_input("Kapacita C [F]", value=0.0001, step=0.00001)

# Výpočty základných veličín
phi_rad = radians(phi_deg)
omega = 2 * pi * f if type_choice == "AC" else 0.0

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

# Reaktancie
XL = omega * L if type_choice == "AC" else 0.0
XC = 1 / (omega * C) if (type_choice == "AC" and C > 0 and omega > 0) else 0.0
Z = R  # impedance
if type_choice == "AC":
    if zataz_type == "R":
        Z = R
    elif zataz_type == "L":
        Z = complex(R, XL)
    elif zataz_type == "C":
        Z = complex(R, -XC)
    elif zataz_type == "RL":
        Z = complex(R, XL)
    elif zataz_type == "RC":
        Z = complex(R, -XC)
    elif zataz_type == "LC":
        Z = complex(0, XL - XC)
    elif zataz_type == "RLC":
        Z = complex(R, XL - XC)

    absZ = abs(Z)
    phi_calc = np.angle(Z, deg=True)
    cos_phi = cos(radians(phi_calc))
else:
    absZ = R
    phi_calc = 0.0
    cos_phi = 1.0

# Výkony
S = Uef * Ief
P = S * cos_phi
Q = sqrt(abs(S**2 - P**2)) if type_choice == "AC" else 0.0

st.subheader("📐 Výsledky a vypočítané hodnoty")
st.markdown(f"""
- **Uef** = {Uef:.2f} V, **Umax** = {Umax:.2f} V  
- **Ief** = {Ief:.2f} A, **Imax** = {Imax:.2f} A  
- **Záťaž** = {zataz_type} → Z = {absZ:.2f} Ω  
- **Fázový posun φ** = {phi_calc:.2f}°, **cos φ** = {cos_phi:.3f}  

**Výkony:**  
- Zdanlivý výkon **S** = {S:.2f} VA  
- Činný výkon **P** = {P:.2f} W  
- Jalový výkon **Q** = {Q:.2f} VAR
""")

# Priebeh veličín v čase
T = 1 / f if f else 1.0
x = np.linspace(0, 2*T, 1000)
napatie = Umax * np.sin(omega * x)
prud = Imax * np.sin(omega * x - radians(phi_calc))
vykon = napatie * prud
vykon_avg = np.mean(vykon)

st.subheader("📈 Časové priebehy")
fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
ax[0].plot(x * 1000, napatie, color='tab:blue', label='Napätie')
ax[1].plot(x * 1000, prud, color='tab:orange', label='Prúd')
ax[2].plot(x * 1000, vykon, color='tab:green', label=f'Výkon (Pavg = {vykon_avg:.2f} W)')

for a in ax:
    a.legend()
    a.grid(True)
    a.set_ylabel("Hodnota")

ax[2].set_xlabel("Čas [ms]")
st.pyplot(fig)

# Export výsledkov
df = pd.DataFrame({
    "čas [s]": x,
    "napätie [V]": napatie,
    "prúd [A]": prud,
    "výkon [W]": vykon
})

st.subheader("📤 Export výpočtov")
col1, col2 = st.columns(2)

with col1:
    csv = df.to_csv(index=False).encode('utf-8')
    st.download_button("📥 CSV export", csv, file_name="vysledky.csv", mime="text/csv")

with col2:
    excel = io.BytesIO()
    with pd.ExcelWriter(excel, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False)
    st.download_button("📥 Excel export", excel.getvalue(), file_name="vysledky.xlsx")

# Fázorová animácia (len pre AC)
if type_choice == "AC":
    st.subheader("🧭 Fázorové zobrazenie")
    st.markdown("*Vizuálne znázornenie napätia a prúdu ako vektorov v čase.*")

    t_fazor = st.slider("Čas [ms]", 0.0, float(2*T*1000), step=1.0)
    t_rad = (t_fazor / 1000.0) * omega
    fig2, ax2 = plt.subplots(figsize=(4, 4))

    ax2.arrow(0, 0, np.cos(t_rad), np.sin(t_rad), head_width=0.05, color='tab:blue', label="U")
    ax2.arrow(0, 0, np.cos(t_rad - radians(phi_calc)), np.sin(t_rad - radians(phi_calc)), head_width=0.05, color='tab:orange', label="I")
    ax2.set_xlim(-1.2, 1.2)
    ax2.set_ylim(-1.2, 1.2)
    ax2.set_aspect('equal')
    ax2.grid(True)
    ax2.legend()
    st.pyplot(fig2)

st.markdown("---")
st.markdown("👨‍🏫 Tento nástroj navrhol elektroinžinier s cieľom poskytnúť nástroj vhodný pre výučbu, sebavzdelávanie aj odbornú analýzu AC a DC obvodov.")
st.markdown("---")
st.markdown("👨Autor: Adrian Mahdon")