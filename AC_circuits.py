import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import graphviz as graphviz

st.set_page_config(layout="wide")
st.title("AC/DC Obvody s vizualizáciou prechodových javov")

# Jednotkové prefixy
unit_prefixes = {
    'p': 1e-12,
    'n': 1e-9,
    'μ': 1e-6,
    'm': 1e-3,
    '': 1,
    'k': 1e3,
    'M': 1e6
}

st.sidebar.markdown("### Nastavenie parametrov obvodu")

# Odpor R
unit_R = st.sidebar.selectbox("Predpona R (napr. k = kilo, M = mega)", list(unit_prefixes.keys()), index=4)
R = st.sidebar.number_input("Odpor R", value=100.0, help="Zadaj hodnotu bez prefixu, napr. 1 a vyber 'k' pre 1 kΩ") * unit_prefixes[unit_R]

# Indukčnosť L
unit_L = st.sidebar.selectbox("Predpona L (napr. m = mili, μ = mikro)", list(unit_prefixes.keys()), index=3)
L = st.sidebar.number_input("Indukčnosť L", value=1.0, help="Zadaj hodnotu bez prefixu, napr. 1 a vyber 'm' pre 1 mH") * unit_prefixes[unit_L]

# Kapacita C
unit_C = st.sidebar.selectbox("Predpona C (napr. μ = mikro, n = nano)", list(unit_prefixes.keys()), index=2)
C = st.sidebar.number_input("Kapacita C", value=1.0, help="Zadaj hodnotu bez prefixu, napr. 1 a vyber 'μ' pre 1 μF") * unit_prefixes[unit_C]

# Typ zdroja
source_type = st.sidebar.radio("Typ zdroja", ("AC", "DC", "DC – prechodový dej"))

# Napätie a frekvencia
U = st.sidebar.number_input("Amplitúda napätia [V]", value=10.0)
f = st.sidebar.number_input("Frekvencia [Hz] (iba pre AC)", value=50.0)

# Časová os
tau = max(L / R if R != 0 else 1, R * C if C != 0 else 1, np.sqrt(L * C) if L > 0 and C > 0 else 1)
t_max = 5 * tau
T_dynamic = st.sidebar.slider("Maximálny čas simulácie [s]", min_value=0.01, max_value=float(t_max), value=float(t_max))
t = np.linspace(0, T_dynamic, 1000)

# Výpočty napätia a prúdu
if source_type == "AC":
    omega = 2 * np.pi * f
    Z = complex(R, omega * L - 1 / (omega * C) if C != 0 else omega * L)
    I = U / abs(Z)
    phi = np.angle(Z)
    u_t = U * np.sin(omega * t)
    i_t = I * np.sin(omega * t - phi)

elif source_type == "DC":
    u_t = U * np.ones_like(t)
    if R > 0:
        i_t = U / R * np.ones_like(t)
    else:
        i_t = np.zeros_like(t)

elif source_type == "DC – prechodový dej":
    if L > 0 and C == 0:
        def model(i, t):
            return (U - R * i) / L
        i_t = odeint(model, 0, t).flatten()
        u_t = U * np.ones_like(t)

    elif C > 0 and L == 0:
        def model(u, t):
            return (U - u) / (R * C)
        u_t = odeint(model, 0, t).flatten()
        i_t = (U - u_t) / R

    elif L > 0 and C > 0:
        def model(x, t):
            q, i = x
            dqdt = i
            didt = (U - R * i - q / C) / L
            return [dqdt, didt]
        x = odeint(model, [0, 0], t)
        q_t, i_t = x[:, 0], x[:, 1]
        u_t = U * np.ones_like(t)

# Vizualizácia napätia a prúdu
fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
ax[0].plot(t, u_t, label="Napätie [V]", color="tab:red")
ax[0].legend()
ax[0].grid(True)
ax[1].plot(t, i_t, label="Prúd [A]", color="tab:blue")
ax[1].legend()
ax[1].grid(True)
ax[1].set_xlabel("Čas [s]")
st.pyplot(fig)

# Schéma obvodu
st.subheader("Schematické znázornenie obvodu")
g = graphviz.Digraph()
g.attr(rankdir='LR')
g.node("1", "Zdroj")
g.node("2", "R")
g.node("3", "L")
g.node("4", "C")
g.node("5", "Uzemnenie")
g.edges([("1", "2"), ("2", "3"), ("3", "4"), ("4", "5")])
st.graphviz_chart(g)

# Edukačné poznámky
st.markdown("""
### Poznámky:
- Použi predpony na rýchle nastavenie hodnôt (napr. `μ` = mikro = 1e-6).
- Prechodový dej sa prejaví len pri zmenách napätia v čase a v prítomnosti L alebo C.
- Pre AC sa počíta efektívny prúd aj fázový posun medzi U a I.
""")