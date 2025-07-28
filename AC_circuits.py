# AC Circuit Designer & Simulator – Profesionálny elektroinžiniersky nástroj



## 📦 AC_circuits.py – Hlavný skript aplikácie

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import networkx as nx
import math
from dataclasses import dataclass
# Vypocet prudu

# -- Komponenty
@dataclass
class Component:
    id: str
    type: str
    value: float
    unit: str
    node_from: str
    node_to: str

# -- Zoznam komponentov
components = []
nodes = set()

st.set_page_config(layout="wide")
st.title("🔌 AC Circuit Designer & Simulator")

# -- Panel na pridávanie komponentov
st.sidebar.header("🧱 Pridaj komponent")
comp_type = st.sidebar.selectbox("Typ komponentu", ["AC Zdroj", "Rezistor", "Kondenzátor", "Cievka", "LED"])
id_ = st.sidebar.text_input("ID (napr. R1)", value=f"{comp_type[0]}{len(components)+1}")
value = st.sidebar.number_input("Hodnota", min_value=0.0, value=100.0)
unit_map = {"AC Zdroj": "V", "Rezistor": "Ω", "Kondenzátor": "μF", "Cievka": "mH", "LED": "V"}
unit = unit_map.get(comp_type, "")
node_from = st.sidebar.text_input("Uzel z")
node_to = st.sidebar.text_input("Uzel do")

if st.sidebar.button("➕ Pridať komponent"):
    c = Component(id=id_, type=comp_type, value=value, unit=unit, node_from=node_from, node_to=node_to)
    components.append(c)
    nodes.add(node_from)
    nodes.add(node_to)

# -- Zobrazenie komponentov
st.subheader("🧮 Zoznam komponentov")
for c in components:
    st.write(f"**{c.id}** ({c.type}): {c.value} {c.unit} | {c.node_from} → {c.node_to}")

# -- Kreslenie obvodu
st.subheader("📐 Schéma obvodu")
G = nx.Graph()
for node in nodes:
    G.add_node(node)
for c in components:
    G.add_edge(c.node_from, c.node_to, label=f"{c.id}\n{c.value} {c.unit}")

pos = nx.spring_layout(G)
fig, ax = plt.subplots()
nx.draw(G, pos, with_labels=True, node_color='lightblue', node_size=2000, ax=ax)
edge_labels = nx.get_edge_attributes(G, 'label')
nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, ax=ax)
st.pyplot(fig)

# -- Simulácia
st.subheader("🔍 Výpočty a analýza")
f = st.number_input("Frekvencia f [Hz]", value=50.0)
omega = 2 * np.pi * f
# Vypocet prudu



Z_total = 0 + 0j
if Z_total == 0:
    st.error("Celková impedancia obvodu je nulová. Skontroluj zapojenie alebo chýbajúce komponenty.")
    st.stop()
for c in components:
    if c.type == "Rezistor":
        Z_total += c.value
    elif c.type == "Cievka":
        L = c.value / 1000  # mH to H
        Z_total += 1j * omega * L
    elif c.type == "Kondenzátor":
        C = c.value * 1e-6  # μF to F
        Z_total += 1 / (1j * omega * C)

U_m = next((c.value for c in components if c.type == "AC Zdroj"), 230)
I = U_m / Z_total
phi = np.angle(Z_total)
P = U_m * abs(I) * np.cos(phi)
Q = U_m * abs(I) * np.sin(phi)
S = U_m * abs(I)

st.write(f"**Celková impedancia** Z = {Z_total:.2f} Ω")
st.write(f"**Fázový posun** φ = {np.degrees(phi):.2f}°")
st.write(f"**Činný výkon** P = {P:.2f} W")
st.write(f"**Jalový výkon** Q = {Q:.2f} var")
st.write(f"**Zdanlivý výkon** S = {S:.2f} VA")

# -- Fázorový diagram
st.subheader("📊 Fázorový diagram")
fig2, ax2 = plt.subplots()
ax2.arrow(0, 0, U_m, 0, head_width=0.5, color='r', label='Napätie')
ax2.arrow(0, 0, abs(I)*np.cos(phi), abs(I)*np.sin(phi), head_width=0.5, color='b', label='Prúd')
ax2.set_xlim(-U_m, U_m)
ax2.set_ylim(-U_m, U_m)
ax2.grid(True)
ax2.set_aspect('equal')
ax2.legend()
st.pyplot(fig2)

# -- Časové priebehy
st.subheader("⏱️ Časový priebeh u(t) a i(t)")
t = np.linspace(0, 0.1, 1000)
u_t = U_m * np.sin(omega * t)
i_t = abs(I) * np.sin(omega * t + phi)
fig3, ax3 = plt.subplots()
ax3.plot(t, u_t, label="Napätie u(t)", color='r')
ax3.plot(t, i_t, label="Prúd i(t)", color='b')
ax3.set_xlabel("čas [s]")
ax3.set_ylabel("veľkosť")
ax3.grid(True)
ax3.legend()
st.pyplot(fig3)
