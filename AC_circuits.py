import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import math

st.set_page_config(page_title="RLC Obvod Simulátor", layout="wide")

st.title("🔌 Interaktívna simulácia RLC obvodov")
st.markdown("""
Tento simulátor ti umožňuje analyzovať chovanie odporu (R), cievky (L) a kondenzátora (C) v rôznych typoch obvodov – **DC aj AC**, vrátane prechodových javov, výpočtov výkonov a fázorového správania.

---
""")

# Pomocná funkcia na aplikovanie SI prefixov
prefix_map = {
    "p" : 1e-12,
    "n" : 1e-9,
    "µ" : 1e-6,
    "m" : 1e-3,
    ""  : 1,
    "k" : 1e3,
    "M" : 1e6,
    "G" : 1e9
}

def apply_prefix(value, prefix):
    return value * prefix_map[prefix]

def format_value_si(value, unit):
    for prefix, factor in reversed(prefix_map.items()):
        if abs(value) >= factor:
            return f"{value / factor:.2f} {prefix}{unit}"
    return f"{value:.2e} {unit}"

# Výber typu zapojenia
type_of_circuit = st.selectbox("Vyber typ zapojenia:", ["DC - R", "DC - RL", "DC - RC", "DC - RLC (prechodový dej)", "AC - R", "AC - RL", "AC - RC", "AC - RLC"])

# Sekcia zadávania hodnôt
st.header("🧮 Zadaj hodnoty prvkov")
col1, col2, col3 = st.columns(3)

with col1:
    R_val = st.number_input("Odpor R", min_value=0.0, value=100.0)
    R_prefix = st.selectbox("Prefix R", list(prefix_map.keys()), index=4, key="Rprefix")
    R = apply_prefix(R_val, R_prefix)
    st.markdown(f"➡️ R = {format_value_si(R, 'Ω')}")

with col2:
    L_val = st.number_input("Indukčnosť L", min_value=0.0, value=10.0)
    L_prefix = st.selectbox("Prefix L", list(prefix_map.keys()), index=2, key="Lprefix")
    L = apply_prefix(L_val, L_prefix)
    st.markdown(f"➡️ L = {format_value_si(L, 'H')}")

with col3:
    C_val = st.number_input("Kapacita C", min_value=0.0, value=1.0)
    C_prefix = st.selectbox("Prefix C", list(prefix_map.keys()), index=2, key="Cprefix")
    C = apply_prefix(C_val, C_prefix)
    st.markdown(f"➡️ C = {format_value_si(C, 'F')}")

st.markdown("---")

st.header("⚙️ Nastavenie zdroja")
colu1, colu2 = st.columns(2)

with colu1:
    U_val = st.number_input("Napätie (U)", value=10.0)
    U_prefix = st.selectbox("Prefix napätia", list(prefix_map.keys()), index=4, key="Uprefix")
    U = apply_prefix(U_val, U_prefix)
    st.markdown(f"➡️ U = {format_value_si(U, 'V')}")

with colu2:
    if "AC" in type_of_circuit:
        f = st.number_input("Frekvencia [Hz]", value=50.0)
        omega = 2 * math.pi * f
        st.markdown(f"➡️ ω = {omega:.2f} rad/s")
    else:
        omega = 0
        f = 0

st.markdown("---")

# Výpočty impedancií
Z_R = complex(R, 0)
Z_L = complex(0, omega * L)
Z_C = complex(0, -1 / (omega * C)) if C > 0 and omega > 0 else complex(0, 0)

if "RLC" in type_of_circuit:
    Z_total = Z_R + Z_L + Z_C
elif "RL" in type_of_circuit:
    Z_total = Z_R + Z_L
elif "RC" in type_of_circuit:
    Z_total = Z_R + Z_C
else:
    Z_total = Z_R

# Výpočty
if Z_total != 0:
    I_complex = U / Z_total
else:
    I_complex = 0

I_mag = abs(I_complex)
phi = np.angle(I_complex, deg=True)

# Výkony
S = U * I_mag  # zdanlivý
P = U * I_mag * math.cos(np.radians(phi))  # činný
Q = U * I_mag * math.sin(np.radians(phi))  # jalový

# Zobrazenie výsledkov
st.header("📊 Výsledky")
colr1, colr2, colr3 = st.columns(3)
colr1.metric("Impedancia |Z|", f"{abs(Z_total):.2f} Ω")
colr2.metric("Prúd |I|", f"{I_mag:.3f} A")
colr3.metric("Fázový posun φ", f"{phi:.2f}°")

st.subheader("💡 Výkony")
st.markdown(f"- **Činný výkon (P):** {P:.2f} W")
st.markdown(f"- **Jalový výkon (Q):** {Q:.2f} VAR")
st.markdown(f"- **Zdanlivý výkon (S):** {S:.2f} VA")

# Fázorový diagram
if "AC" in type_of_circuit:
    st.subheader("📐 Fázorový diagram")
    fig, ax = plt.subplots()
    ax.quiver(0, 0, U, 0, angles='xy', scale_units='xy', scale=1, color='r', label='Napätie U')
    ax.quiver(0, 0, I_mag * math.cos(np.radians(phi)), I_mag * math.sin(np.radians(phi)),
              angles='xy', scale_units='xy', scale=1, color='b', label='Prúd I')
    ax.set_xlim(-U, U)
    ax.set_ylim(-U, U)
    ax.set_aspect('equal')
    ax.grid(True)
    ax.legend()
    st.pyplot(fig)

# Časový priebeh – len pre RL, RC, RLC
if "DC" in type_of_circuit and ("RL" in type_of_circuit or "RC" in type_of_circuit or "RLC" in type_of_circuit):
    st.subheader("⏱️ Časový priebeh napätia a prúdu (prechodový jav)")
    t = np.linspace(0, 5, 1000)
    if "RL" in type_of_circuit:
        tau = L / R if R > 0 else 0
        i_t = (U / R) * (1 - np.exp(-t / tau))
        u_L = U - R * i_t
        st.markdown(f"Časová konštanta τ = {tau:.3f} s")
    elif "RC" in type_of_circuit:
        tau = R * C if R > 0 else 0
        u_C = U * (1 - np.exp(-t / tau))
        i_t = (U / R) * np.exp(-t / tau)
        st.markdown(f"Časová konštanta τ = {tau:.3f} s")
    elif "RLC" in type_of_circuit:
        # Podobné ako tlmený oscilátor – zložitejší priebeh
        alpha = R / (2 * L)
        omega_0 = 1 / np.sqrt(L * C)
        if alpha > omega_0:
            # Preťažený (aperiodický)
            s1 = -alpha + np.sqrt(alpha**2 - omega_0**2)
            s2 = -alpha - np.sqrt(alpha**2 - omega_0**2)
            i_t = (U / L) * (np.exp(s1 * t) - np.exp(s2 * t))
        else:
            omega_d = np.sqrt(omega_0**2 - alpha**2)
            i_t = (U / (L * omega_d)) * np.exp(-alpha * t) * np.sin(omega_d * t)
    fig2, ax2 = plt.subplots()
    ax2.plot(t, i_t, label="Prúd i(t)", color="b")
    ax2.set_xlabel("Čas [s]")
    ax2.set_ylabel("Prúd [A]")
    ax2.grid(True)
    ax2.legend()
    st.pyplot(fig2)

# Diagnostika
st.markdown("---")
st.markdown("### ℹ️ Diagnostika a kontrola")
if R == 0 and L == 0 and C == 0:
    st.warning("Nie je zvolený žiadny prvok – prosím zadaj aspoň jeden.")
else:
    prvky = []
    if R > 0: prvky.append("R")
    if L > 0: prvky.append("L")
    if C > 0: prvky.append("C")
    st.success(f"Zvolené prvky: {', '.join(prvky)}")