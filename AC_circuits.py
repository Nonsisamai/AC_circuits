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

# Časové rozlíšenie pre prechodové javy
if type_choice == "DC - Prechodový dej (R-C / R-L)":
    st.sidebar.markdown("---")
    st.sidebar.markdown("⏱️ **Čas simulácie prechodu**")
    t_max = st.sidebar.number_input("Maximálny čas simulácie [s]", value=1.0, min_value=0.01, step=0.1)
    t_points = st.sidebar.number_input("Počet bodov", value=1000, step=100)
else:
    t_max = 0.1
    t_points = 1000

# Auto-zistenie typu záťaže
zataz_popis = []
if R > 0: zataz_popis.append("R")
if L > 0: zataz_popis.append("L")
if C > 0: zataz_popis.append("C")
zataz_type = "+".join(zataz_popis) if zataz_popis else "Žiadna záťaž"

st.subheader("📋 Prehľad zapojenia")
st.markdown(f"""
- **Zvolený režim:** {type_choice}  
- **Záťaž v obvode:** {zataz_type if zataz_type else "(žiadna)"}
""")

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
Q = sqrt(abs(S**2 - P**2)) if type_choice == "AC" else 0.0

# Časová os a simulácie
x = np.linspace(0, t_max, int(t_points))
annotation_time = None
if type_choice == "AC":
    T = 1 / f if f > 0 else 1.0
    x = np.linspace(0, 2*T, 1000)
    napatie = Umax * np.sin(omega * x)
    prud = Imax * np.sin(omega * x - phi_calc_rad)
elif type_choice == "DC":
    napatie = np.full_like(x, Uef)
    prud = np.full_like(x, Ief)
    tau = None
else:
    if C > 0 and R > 0:
        tau = R * C
        napatie = Uef * (1 - np.exp(-x / tau))
        prud = (Uef / R) * np.exp(-x / tau)
        annotation_time = 5 * tau
    elif L > 0 and R > 0:
        tau = L / R
        napatie = np.full_like(x, Uef)
        prud = (Uef / R) * (1 - np.exp(-x / tau))
        annotation_time = 5 * tau
    else:
        napatie = np.zeros_like(x)
        prud = np.zeros_like(x)
        tau = None

vykon = napatie * prud
vykon_avg = np.mean(vykon)

# Doplnková informácia o τ (časová konštanta)
if type_choice.startswith("DC") and tau is not None:
    st.markdown(f"**Časová konštanta τ =** {tau:.4f} s")

# Zobrazenie bodu, kedy sa kondenzátor nabije na 99 %
if annotation_time:
    st.markdown(f"⚡ **Prechod ustálený do:** {annotation_time:.3f} s (≈ 5τ)")
st.markdown("---")
st.markdown("👨Autor: Adrian Mahdon")