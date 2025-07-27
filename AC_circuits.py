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

# Interaktívna schéma
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

# 🔁 Interaktívny blok – voľba zobrazenia teórie a výpočtov
st.markdown("---")
teach_level = st.radio("Zobraziť výpočty na úrovni:", ["Stredoškolská úroveň", "Vysokoškolská úroveň"], horizontal=True)

if teach_level == "Stredoškolská úroveň":
    st.success("Zobrazujú sa zjednodušené vzorce a výsledky pre ľahšie pochopenie.")
else:
    st.info("Zobrazujú sa všetky kroky, kompletné vývody a technické detaily.")

if st.button("📐 Zobraziť všetky výpočtové kroky"):
    st.subheader("📘 Podrobné výpočty")

    if type_choice == "DC - Prechodový dej (R-C / R-L)":
        if C > 0:
            st.markdown("""
            ### Prechodový dej pre RC obvod
            - Rovnica pre napätie: \( u_C(t) = U \cdot (1 - e^{-t/RC}) \)
            - Rovnica pre prúd: \( i(t) = \frac{U}{R} \cdot e^{-t/RC} \)
            - Časová konštanta: \( \tau = RC \)
            - Po čase \( 5\tau \) sa kondenzátor považuje za plne nabitý (ustálený stav)

            **Správanie kondenzátora:**
            - Na začiatku: \( i \approx \frac{U}{R} \), kondenzátor sa správa ako skrat
            - Po čase: \( i \to 0 \), správa sa ako otvorený obvod
            """)
        elif L > 0:
            st.markdown("""
            ### Prechodový dej pre RL obvod
            - Rovnica pre prúd: \( i_L(t) = \frac{U}{R} \cdot (1 - e^{-Rt/L}) \)
            - Napätie na cievke: \( u_L(t) = U \cdot e^{-Rt/L} \)
            - Časová konštanta: \( \tau = \frac{L}{R} \)
            - Po čase \( 5\tau \) sa prúd ustáli a cievka sa správa ako vodič

            **Správanie cievky:**
            - Na začiatku: \( i = 0 \), správa sa ako otvorený obvod
            - Po čase: \( i = \frac{U}{R} \), správa sa ako skrat
            """)
    elif type_choice == "AC":
        st.markdown("""
        ### Ustálený stav v AC obvodoch
        - Impedancia: \( Z = R + j(X_L - X_C) \)
        - Reaktancie:
          - \( X_L = \omega L \)
          - \( X_C = \frac{1}{\omega C} \)
        - Fázový posun: \( \varphi = \arctan \left( \frac{X_L - X_C}{R} \right) \)
        - Výkony:
          - Zdanlivý: \( S = U_{ef} \cdot I_{ef} \)
          - Činný: \( P = S \cdot \cos\varphi \)
          - Jalový: \( Q = \sqrt{S^2 - P^2} \)
        """)
    else:
        st.markdown("""
        ### Jednoduchý výpočet DC obvodu
        - \( R = \frac{U}{I} \)
        - \( P = U \cdot I \)
        """)

    st.success("Výpočty boli úspešne zobrazené.")
st.markdown("---")
st.markdown("👨Autor: Adrian Mahdon")