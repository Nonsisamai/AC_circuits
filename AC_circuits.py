# Interactive AC circuit calculator with animated plots
# Streamlit-based educational tool

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import cos, sin, radians, pi

st.set_page_config(page_title="AC Circuit Interactive Tool", layout="wide")
st.title("🔌 Interaktívna vizualizácia AC obvodov")

st.markdown("""
Tento nástroj ti umožní:
- zadať efektívne alebo maximálne hodnoty napätia/prúdu
- pracovať s fázovým posunom, účinníkom
- zobraziť interaktívne grafy priebehov
- získať výpočty výkonov (činný, jalový, zdanlivý)
- pozorovať jednoduchú animáciu striedavých veličín
""")

st.sidebar.header("🎛️ Zadanie parametrov")

# Input mode
input_mode = st.sidebar.radio("Zadávaš efektívne alebo maximálne hodnoty?", ["Efektívne (RMS)", "Maximálne (peak)"])

# Universal float input (with step for decimals)
def float_input(label, default):
    return st.sidebar.number_input(label, value=float(default), step=0.01, format="%.3f")

# Values
U_input = float_input("Napätie [V]", 230)
I_input = float_input("Prúd [A]", 5)
cos_phi = float_input("Účinník (cos φ)", 0.9)
phi_deg = float_input("Fázový posun φ [stupne]", np.degrees(np.arccos(0.9)))
f = float_input("Frekvencia [Hz]", 50)

# Derived quantities
phi_rad = radians(phi_deg)
omega = 2 * pi * f

# Convert RMS to peak if needed
if input_mode == "Efektívne (RMS)":
    Umax = U_input * np.sqrt(2)
    Imax = I_input * np.sqrt(2)
else:
    Umax = U_input
    Imax = I_input

# Check for missing values or inconsistencies
if not (0 <= cos_phi <= 1):
    st.error("⚠️ Účinník musí byť medzi 0 a 1")
    st.stop()

# Power calculations
S = U_input * I_input                     # zdanlivý výkon [VA]
P = S * cos_phi                           # činný výkon [W]
Q = np.sqrt(S**2 - P**2)                 # jalový výkon [VAR]

st.subheader("📊 Výpočty výkonov")
st.markdown(f"""
- **Zdanlivý výkon S** = {S:.2f} VA  
- **Činný výkon P** = {P:.2f} W  
- **Jalový výkon Q** = {Q:.2f} VAR  
- **Fázový posun φ** = {phi_deg:.2f}°
""")

# Time vector
T = 1 / f
x = np.linspace(0, 2*T, 1000)

tension = Umax * np.sin(omega * x)
current = Imax * np.sin(omega * x - phi_rad)

# Plotting section
st.subheader("📈 Priebehy napätia a prúdu")

fig, ax = plt.subplots(figsize=(10, 4))
ax.plot(x * 1000, tension, label="Napätie [V]", color="tab:blue")
ax.plot(x * 1000, current, label="Prúd [A]", color="tab:orange")
ax.set_xlabel("Čas [ms]")
ax.set_ylabel("Hodnota")
ax.grid(True)
ax.legend()
st.pyplot(fig)

# Animation section
st.subheader("🎞️ Animácia fázového posunu")

fig2, ax2 = plt.subplots(figsize=(6, 6))

line1, = ax2.plot([], [], 'b-', lw=2, label="Napätie")
line2, = ax2.plot([], [], 'r-', lw=2, label="Prúd")
ax2.set_xlim(-1.2, 1.2)
ax2.set_ylim(-1.2, 1.2)
ax2.set_aspect('equal')
ax2.grid(True)
ax2.legend()

circle = plt.Circle((0, 0), 1, color='lightgray', fill=False)
ax2.add_artist(circle)

frames = 100
theta = np.linspace(0, 2 * pi, frames)

def animate(i):
    angle = theta[i % frames]
    line1.set_data([0, np.cos(angle)], [0, np.sin(angle)])
    line2.set_data([0, np.cos(angle - phi_rad)], [0, np.sin(angle - phi_rad)])
    return line1, line2

ani = animation.FuncAnimation(fig2, animate, frames=frames, interval=50, blit=True)
st.pyplot(fig2)

st.markdown("---")
st.markdown("👨‍🏫 Tento nástroj slúži na výuku a intuitívne pochopenie javov v striedavých obvodoch. Pre odborné výpočty použite verifikovaný software.Autor: Adrian Mahdon")
st.markdown("---")
st.markdown("👨Autor: Adrian Mahdon")