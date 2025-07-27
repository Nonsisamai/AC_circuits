import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
st.title("AC Obvod: Napätie, prúd, výkon")

# Autor: Adrian Mahdon
# Parametre od používateľa
Uef = st.slider("Napätie (Uef)", 100, 400, 230)
Ief = st.slider("Prúd (Ief)", 1, 20, 7)
phi_deg = st.slider("Fázový posun (°)", 0, 90, 30)
f = st.slider("Frekvencia (Hz)", 10, 60, 50)

# Výpočty
phi = np.deg2rad(phi_deg)
Um = Uef * np.sqrt(2)
Im = Ief * np.sqrt(2)
omega = 2 * np.pi * f
T = 1 / f
t = np.linspace(0, 2*T, 1000)
u_t = Um * np.sin(omega * t)
i_t = Im * np.sin(omega * t - phi)
p_t = u_t * i_t
P_avg = (Um * Im / 2) * np.cos(phi)

# Grafy
fig, axs = plt.subplots(3, 1, figsize=(10, 8))

axs[0].plot(t, u_t, label="Napätie u(t)")
axs[0].plot(t, i_t, label="Prúd i(t)")
axs[0].legend()
axs[0].set_title("Napätie a prúd")
axs[0].grid()

axs[1].plot(t, p_t, color="orange", label="Výkon p(t)")
axs[1].axhline(P_avg, color="red", linestyle="--", label="Pavg")
axs[1].legend()
axs[1].set_title("Okamžitý výkon")
axs[1].grid()

axs[2].plot(t, p_t - P_avg, color="green", label="p(t) - Pavg")
axs[2].legend()
axs[2].set_title("Oscilačná zložka výkonu")
axs[2].grid()

st.pyplot(fig)