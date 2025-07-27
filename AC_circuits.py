# AC/DC Educational Circuit Visualizer with Expert Logic
# Streamlit app designed by an electrical engineer & educator

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import io
from math import cos, sin, radians, degrees, pi, sqrt, atan2

# Page setup
st.set_page_config(page_title="AC/DC Circuit Visualizer", layout="wide")

# Theme toggle
dark_mode = st.sidebar.toggle("🌙 Tmavý režim")
if dark_mode:
    plt.style.use('dark_background')
else:
    plt.style.use('default')

st.title("🔌 AC/DC Obvody – Vizualizácia a Výpočty pre študentov")

st.markdown("""
Tento nástroj slúži na pochopenie správania sa AC a DC obvodov s rôznymi záťažami (R, L, C, RLC):
- Zobrazí časové priebehy napätia, prúdu a výkonu
- Umožní výpočet výkonov a reaktancií
- Rozlišuje medzi AC a DC režimom
- Zobrazí medzičísla, výpočtové vzorce a výsledky
- Exportuje výsledky a vizualizácie
""")

st.sidebar.header("🎛️ Parametre obvodu")

# Typ obvodu: AC alebo DC
type_choice = st.sidebar.selectbox("Režim obvodu", ["AC", "DC"])

# Zadávanie hodnôt
input_mode = st.sidebar.radio("Zadávaš hodnoty ako:", ["Efektívne (RMS)", "Maximálne (peak)"])

# Napätie a prúd
U_in = st.sidebar.number_input("Napätie [V]", value=230.0, step=0.1)
I_in = st.sidebar.number_input("Prúd [A]", value=5.0, step=0.1)

# Frekvencia a fázový posun (len pre AC)
f = st.sidebar.number_input("Frekvencia [Hz]", value=50.0, step=1.0) if type_choice == "AC" else 0.0
phi_deg = st.sidebar.number_input("Fázový posun φ [°]", value=30.0 if type_choice == "AC" else 0.0)
phi_rad = radians(phi_deg)
omega = 2 * pi * f if f > 0 else 0

# Hodnoty súčiastok
R = st.sidebar.number_input("Odpor R [Ω]", value=0.0, step=0.1)
L = st.sidebar.number_input("Indukčnosť L [H]", value=0.0, step=0.001)
C = st.sidebar.number_input("Kapacita C [F]", value=0.0, step=0.00001)

# Automatický výber typu záťaže
zataz_popis = []
if R > 0: zataz_popis.append("R")
if L > 0: zataz_popis.append("L")
if C > 0: zataz_popis.append("C")
zataz_type = "+".join(zataz_popis) if zataz_popis else "Žiadna záťaž"

# Výpočty základných hodnôt
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

# Reaktancie a impedancia
XL = omega * L if L > 0 else 0.0
XC = 1 / (omega * C) if (C > 0 and omega > 0) else 0.0
Z_complex = complex(R, XL - XC)
Z_abs = abs(Z_complex)
phi_calc_rad = atan2(Z_complex.imag, Z_complex.real)
phi_calc_deg = degrees(phi_calc_rad)
cos_phi = cos(phi_calc_rad) if Z_abs > 0 else 1.0

# Výpočty výkonov
S = Uef * Ief
P = S * cos_phi
Q = sqrt(abs(S**2 - P**2)) if type_choice == "AC" else 0.0

# Výpočet napätí, prúdov, výkonu v čase
t_duration = 2 / f if f > 0 else 0.1
x = np.linspace(0, t_duration, 1000)
napatie = Umax * np.sin(omega * x) if type_choice == "AC" else np.full_like(x, Uef)
prud = Imax * np.sin(omega * x - phi_calc_rad) if type_choice == "AC" else np.full_like(x, Ief)
vykon = napatie * prud
vykon_avg = np.mean(vykon)

# Zobrazenie výpočtov
st.subheader("📐 Vypočítané hodnoty")
st.markdown(f"""
**Zadané hodnoty:**  
- Režim: **{type_choice}**, Zadávanie: **{input_mode}**  
- Typ záťaže: **{zataz_type if zataz_type else 'žiadna'}**

**Základné výpočty:**  
- Uef = {Uef:.2f} V, Umax = {Umax:.2f} V  
- Ief = {Ief:.2f} A, Imax = {Imax:.2f} A  
- Z = {Z_abs:.2f} Ω, φ = {phi_calc_deg:.2f}°, cosφ = {cos_phi:.3f}  

**Výkony:**  
- S = {S:.2f} VA  
- P = {P:.2f} W  
- Q = {Q:.2f} VAR (len pre AC)  
- ⟨P(t)⟩ ≈ {vykon_avg:.2f} W  
""")

# Grafy
st.subheader("📈 Časové priebehy")
fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
ax[0].plot(x * 1000, napatie, label='Napätie [V]', color='tab:blue')
ax[1].plot(x * 1000, prud, label='Prúd [A]', color='tab:orange')
ax[2].plot(x * 1000, vykon, label=f'Výkon [W] (avg={vykon_avg:.2f})', color='tab:green')

for a in ax:
    a.legend()
    a.grid(True)
    a.set_ylabel("Hodnota")

ax[2].set_xlabel("Čas [ms]")
st.pyplot(fig)

# Fázorové zobrazenie (len AC)
if type_choice == "AC" and zataz_type != "":
    st.subheader("🧭 Fázorový diagram")
    t_fazor = st.slider("Fázorový čas [ms]", 0.0, float(t_duration * 1000), step=1.0)
    t_rad = (t_fazor / 1000.0) * omega
    fig2, ax2 = plt.subplots(figsize=(3.5, 3.5))
    ax2.arrow(0, 0, cos(t_rad), sin(t_rad), head_width=0.05, color='tab:blue', label="Napätie")
    ax2.arrow(0, 0, cos(t_rad - phi_calc_rad), sin(t_rad - phi_calc_rad), head_width=0.05, color='tab:orange', label="Prúd")
    ax2.set_xlim(-1.2, 1.2)
    ax2.set_ylim(-1.2, 1.2)
    ax2.set_aspect('equal')
    ax2.grid(True)
    ax2.legend()
    st.pyplot(fig2)

# Export
st.subheader("📤 Export")
df = pd.DataFrame({
    "čas [s]": x,
    "napätie [V]": napatie,
    "prúd [A]": prud,
    "výkon [W]": vykon
})
csv = df.to_csv(index=False).encode('utf-8')
st.download_button("📥 Export do CSV", csv, file_name="vypocet.csv", mime='text/csv')

st.markdown("---")
st.markdown("👨‍🏫 všetky výpočty vychádzajú zo základných fyzikálnych zákonov a reálnych elektrických modelov.")
st.markdown("---")
st.markdown("👨Autor: Adrian Mahdon")