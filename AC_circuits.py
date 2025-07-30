
import streamlit as st
import numpy as np
import matplotlib.pyplot as plt

# --- Constants ---
prefixes = {
    "p": 1e-12, "n": 1e-9, "µ": 1e-6, "m": 1e-3,
    "": 1,
    "k": 1e3, "M": 1e6
}

# --- Helper Functions ---
def apply_prefix(value, prefix):
    try:
        return float(value) * prefixes.get(prefix, 1)
    except:
        return 0

def calculate_impedance(R, L, C, omega):
    Z_R = R
    Z_L = 1j * omega * L if L > 0 else 0
    Z_C = -1j / (omega * C) if C > 0 else 0
    return Z_R + Z_L + Z_C

def calculate_current(U, Z):
    return U / Z if Z != 0 else 0

def calculate_powers(U, I, phi):
    S = U * I
    P = U * I * np.cos(phi)
    Q = U * I * np.sin(phi)
    return abs(S), abs(P), abs(Q)

def calculate_voltage_drop(I, Z):
    return I * Z

def temperature_adjust_resistance(R, temp, temp_coeff=0.004):  # Approx. for copper
    return R * (1 + temp_coeff * (temp - 20))

# --- Streamlit UI ---
st.title("💡 Interaktívna simulácia AC/DC obvodov VO VYVOJI viac verzii")
st.markdown("### Vyber typ simulácie")
sim_type = st.selectbox("Typ simulácie", ["AC", "DC", "DC - prechodový dej"])
st.markdown("### Zadaj základné parametre obvodu")

col1, col2, col3 = st.columns(3)
with col1:
    R_val = st.text_input("Odpor (R)", "100")
    R_pre = st.selectbox("Prefix R", list(prefixes.keys()), index=4)
with col2:
    L_val = st.text_input("Indukčnosť (L)", "10")
    L_pre = st.selectbox("Prefix L", list(prefixes.keys()), index=3)
with col3:
    C_val = st.text_input("Kapacita (C)", "100")
    C_pre = st.selectbox("Prefix C", list(prefixes.keys()), index=2)

col4, col5 = st.columns(2)
with col4:
    U_val = st.text_input("Napätie zdroja (U)", "10")
    U_pre = st.selectbox("Prefix U", list(prefixes.keys()), index=4)
with col5:
    f_val = st.text_input("Frekvencia (Hz)", "50")

st.markdown("### Rozšírené nastavenia")
col6, col7 = st.columns(2)
with col6:
    temp_env = st.slider("Teplota okolia (°C)", -20, 100, 20)
with col7:
    include_temp = st.checkbox("Zohľadniť vplyv teploty na odpor", value=True)

# --- Prepočty vstupov ---
R = apply_prefix(R_val, R_pre)
L = apply_prefix(L_val, L_pre)
C = apply_prefix(C_val, C_pre)
U = apply_prefix(U_val, U_pre)
f = float(f_val) if f_val else 0
omega = 2 * np.pi * f

if include_temp:
    R = temperature_adjust_resistance(R, temp_env)

# --- Výpočty ---
Z = calculate_impedance(R, L, C, omega if sim_type == "AC" else 0.001)
I = calculate_current(U, Z)
phi = np.angle(Z)
S, P, Q = calculate_powers(U, abs(I), phi)
UL = abs(calculate_voltage_drop(I, 1j * omega * L)) if L > 0 else 0
UC = abs(calculate_voltage_drop(I, -1j / (omega * C))) if C > 0 else 0
UR = abs(I * R)

# --- Výstup ---
st.markdown("## 💬 Výsledky simulácie")
st.write(f"**Celková impedancia Z:** {abs(Z):.2f} Ω")
st.write(f"**Fázový posun φ:** {np.degrees(phi):.2f}°")
st.write(f"**Prúd v obvode I:** {abs(I):.2f} A")
st.write(f"**Napätie na R:** {UR:.2f} V, na L: {UL:.2f} V, na C: {UC:.2f} V")
st.write(f"**Výkony:** Zdanlivý = {S:.2f} VA, Činný = {P:.2f} W, Jalový = {Q:.2f} var")

# --- Vzorce a vysvetlenia ---
st.markdown("### 📘 Edukačné vysvetlenie")
st.info("'Impedancia RLC obvodu Z = R + jωL - j/(ωC)'"        "Fázorový prúd: 'I = U / Z'"        "Výkony: 'S = UI', 'P = UI·cos(φ)', 'Q = UI·sin(φ)'"        "Zohľadnený vplyv teploty na odpor: 'R = R₀(1 + α·(T - 20))'")

# --- Grafy ---
t = np.linspace(0, 0.1, 1000)
u_t = U * np.sin(omega * t)
i_t = abs(I) * np.sin(omega * t - phi)

fig, ax = plt.subplots()
ax.plot(t * 1000, u_t, label="Napätie U(t)", color="blue")
ax.plot(t * 1000, i_t, label="Prúd I(t)", color="orange")
ax.set_xlabel("Čas [ms]")
ax.set_ylabel("Hodnota")
ax.set_title("🧭 Časové priebehy napätia a prúdu")
ax.grid(True)
ax.legend()
st.pyplot(fig)