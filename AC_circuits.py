# AC/DC Circuit Designer & Analyzer – Advanced Educational Simulator
# Autor: Adrian Mahdon / Elektroinzinier & Veduci Nastroj pre Vizualizaciu a Simulaciu

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import graphviz
from fpdf import FPDF
from streamlit_drawable_canvas import st_canvas
import io
from math import cos, sin, radians, degrees, pi, sqrt, atan2, exp

st.set_page_config(page_title="AC/DC Circuit Designer", layout="wide")

st.title("⚡ AC/DC Circuit Designer & Simulator")
st.markdown("""
Elektroinžinierstvo: **Komplexný nástroj pre návrh, simuláciu a analýzu obvodov**.
- Navrhni obvod interaktívne (Canvas)
- Automatické rozpoznanie topológie
- Reálne simulácie AC/DC prûdov, napätí, prechodov, fázorov
- Fyzikálne presné modely: R, L, C, τ, ϕ, Z, S, P, Q
""")

# ----------------------------- Sidebar vstupy ----------------------------------
st.sidebar.header("🎛️ Základné nastavenia")
type_choice = st.sidebar.selectbox("Režim obvodu", ["AC", "DC", "Prechodový (RC/RL)"])
input_mode = st.sidebar.radio("Zadávanie hodnôt", ["Efektívne (RMS)", "Maximálne (Peak)"])
U_in = st.sidebar.number_input("Napätie [V]", value=230.0, step=1.0)
I_in = st.sidebar.number_input("Prúd [A]", value=5.0, step=0.1)
f = st.sidebar.number_input("Frekvencia [Hz]", value=50.0 if type_choice=="AC" else 0.0)
phi_manual = radians(st.sidebar.number_input("Fázový posun [°] (ak známy)", value=0.0))

# ----------------------------- Interaktívny canvas ------------------------------
st.subheader("🖊️ Interaktívne kreslenie obvodu")
st.markdown("Nakresli obvod s komponentmi R, L, C a Zdrojom. Budú automaticky rozpoznané.")
canvas_result = st_canvas(
    fill_color="rgba(255,255,255,0.0)",
    stroke_width=3,
    background_color="#fff",
    drawing_mode="freedraw",
    height=300,
    width=800,
    key="circuit_canvas"
)

# Debug: vypísanie objektov z canvasu
if canvas_result.json_data is not None:
    objects = canvas_result.json_data.get("objects", [])
    st.write(f"🔍 Nájdené objekty: {len(objects)}")
    st.json(objects)
else:
    objects = []

# ----------------------------- Komponenty (alternatívne zadanie) ----------------
st.subheader("🧩 Konfigurovateľné komponenty")
components = []

col1, col2, col3 = st.columns(3)
with col1:
    R = st.number_input("Odpor R [Ω]", value=100.0)
with col2:
    L = st.number_input("Indukčnosť L [H]", value=0.1)
with col3:
    C = st.number_input("Kapacita C [F]", value=0.00001)

# Automatický výpočet impedance
omega = 2 * pi * f
Z_R = complex(R, 0)
Z_L = complex(0, omega*L)
Z_C = complex(0, -1/(omega*C)) if C > 0 and omega > 0 else complex(0, 0)
Z_total = Z_R + Z_L + Z_C
Z_abs = abs(Z_total)
phi_calc_rad = atan2(Z_total.imag, Z_total.real)
phi_calc_deg = degrees(phi_calc_rad)
cos_phi = cos(phi_calc_rad)

# Prepocet U, I ak su zadané ako RMS/peak
if input_mode == "Efektívne (RMS)":
    Uef, Ief = U_in, I_in
    Umax, Imax = Uef * sqrt(2), Ief * sqrt(2)
else:
    Umax, Imax = U_in, I_in
    Uef, Ief = Umax / sqrt(2), Imax / sqrt(2)

# ----------------------------- Simulácia ----------------------------------------
x = np.linspace(0, 0.1, 1000)
if type_choice == "AC":
    napatie = Umax * np.sin(omega * x)
    prud = Imax * np.sin(omega * x - phi_calc_rad)
elif type_choice == "DC":
    napatie = np.full_like(x, Uef)
    prud = np.full_like(x, Ief)
elif type_choice == "Prechodový (RC/RL)":
    if R > 0 and C > 0:
        tau = R*C
        napatie = Uef * (1 - np.exp(-x / tau))
        prud = (Uef/R) * np.exp(-x / tau)
    elif R > 0 and L > 0:
        tau = L/R
        napatie = np.full_like(x, Uef)
        prud = (Uef/R)*(1 - np.exp(-x / tau))
    else:
        napatie = np.zeros_like(x)
        prud = np.zeros_like(x)

vykon = napatie * prud
vykon_avg = np.mean(vykon)

# ----------------------------- Grafy ---------------------------------------------
st.subheader("📊 Priebeh veličí v čase")
fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
ax[0].plot(x, napatie, label="U(t) [V]", color='tab:blue')
ax[1].plot(x, prud, label="I(t) [A]", color='tab:orange')
ax[2].plot(x, vykon, label=f"P(t) [W], ⟨P⟩={vykon_avg:.2f}", color='tab:green')
for a in ax:
    a.legend()
    a.grid(True)
    a.set_ylabel("Hodnota")
ax[2].set_xlabel("Čas [s]")
st.pyplot(fig)

# ----------------------------- Fázorový diagram ----------------------------------
st.subheader("🧭 Fázorový diagram")
fig2, ax2 = plt.subplots()
ax2.quiver(0, 0, Umax*cos(0), Umax*sin(0), angles='xy', scale_units='xy', scale=1, color='b', label="U")
ax2.quiver(0, 0, Imax*cos(-phi_calc_rad), Imax*sin(-phi_calc_rad), angles='xy', scale_units='xy', scale=1, color='r', label="I")
ax2.set_xlim(-Umax, Umax)
ax2.set_ylim(-Umax, Umax)
ax2.set_aspect('equal')
ax2.grid(True)
ax2.legend()
st.pyplot(fig2)

# ----------------------------- Výsledky ------------------------------------------
S = Uef * Ief
P = S * cos_phi
Q = sqrt(abs(S**2 - P**2)) if type_choice == "AC" else 0.0

st.subheader("🧮 Výpočty")
st.markdown(f"""
- **Z =** {Z_total:.2f} Ω  
- **|Z| =** {Z_abs:.2f} Ω  
- **φ =** {phi_calc_deg:.2f}°  
- **cos(φ) =** {cos_phi:.3f}  
- **S =** {S:.2f} VA  
- **P =** {P:.2f} W  
- **Q =** {Q:.2f} VAR
""")

st.markdown("---")
st.markdown("👨‍🏫 *Autor: Adrian Mahdon – Interaktívny elektro-vzdelávací simulátor*")