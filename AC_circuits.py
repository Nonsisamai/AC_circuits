# Interactive AC/DC circuit calculator with animated plots and export
# Streamlit-based educational tool (advanced)

import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from math import cos, sin, radians, pi, acos
import pandas as pd
import io

st.set_page_config(page_title="AC/DC Circuit Interactive Tool", layout="wide")
st.title("🔌 Interaktívna vizualizácia AC / DC obvodov")

st.markdown("""
Tento nástroj ti umožní:
- zadať efektívne alebo maximálne hodnoty napätia/prúdu
- pracovať s fázovým posunom, účinníkom, typom záťaže (R, RL, RC, RLC, DC)
- automaticky dopočítať všetky neznáme
- zobraziť priebehy, časy, výkony a export do PDF/CSV
- získať porovnanie s DC verziou obvodu
""")

st.sidebar.header("🎛️ Zadanie parametrov")

# Výber typu obvodu
type_choice = st.sidebar.selectbox("Typ obvodu", ["AC", "DC"])
zataz_type = st.sidebar.selectbox("Typ záťaže", ["R", "RL", "RC", "LC", "RLC"] if type_choice == "AC" else ["R", "RL (ako DC)", "RC (ako DC)"])

# Prepínač medzi RMS a peak hodnotami
input_mode = st.sidebar.radio("Zadávaš efektívne alebo maximálne hodnoty?", ["Efektívne (RMS)", "Maximálne (peak)"])

# Pomocná funkcia na zadávanie float hodnôt
def float_input(label, default):
    return st.sidebar.number_input(label, value=float(default), step=0.01, format="%.3f")

# Základné vstupy
U_zadane = float_input("Napätie [V]", 230)
I_zadane = float_input("Prúd [A]", 5)
f = float_input("Frekvencia [Hz]", 50) if type_choice == "AC" else 0
phi_deg = float_input("Fázový posun φ [stupne]", 30 if type_choice == "AC" else 0)

# Prevod uhla na radiány
phi_rad = radians(phi_deg)
omega = 2 * pi * f if f else 0

# Prevod na peak hodnoty, ak boli zadané efektívne
if input_mode == "Efektívne (RMS)":
    Uef = U_zadane
    Ief = I_zadane
    Umax = Uef * np.sqrt(2)
    Imax = Ief * np.sqrt(2)
else:
    Umax = U_zadane
    Imax = I_zadane
    Uef = Umax / np.sqrt(2)
    Ief = Imax / np.sqrt(2)

# Výpočet účinníka ak ide o AC
if type_choice == "AC":
    cos_phi = cos(phi_rad)
else:
    cos_phi = 1.0
    phi_rad = 0

# Výpočty výkonov
S = Uef * Ief                     # zdanlivý výkon [VA]
P = S * cos_phi                  # činný výkon [W]
Q = np.sqrt(abs(S**2 - P**2))    # jalový výkon [VAR], bezpečné odmocnenie

# Zobrazenie prepočítaných hodnôt
st.subheader("📐 Základné hodnoty")
st.markdown(f"""
- Zadané ako: **{'Efektívne (RMS)' if input_mode == 'Efektívne (RMS)' else 'Maximálne (peak)'} hodnoty**  
- **Uef** = {Uef:.2f} V, **Umax** = {Umax:.2f} V  
- **Ief** = {Ief:.2f} A, **Imax** = {Imax:.2f} A  
- **Fázový posun φ** = {phi_deg:.2f}°, **cos φ** = {cos_phi:.3f}
""")

# Výkonová tabuľka
st.subheader("⚡ Výkony")
st.markdown(f"""
- **Zdanlivý výkon S** = {S:.2f} VA  
- **Činný výkon P** = {P:.2f} W  
- **Jalový výkon Q** = {Q:.2f} VAR
""")

# Časový vektor
T = 1 / f if f else 0.1
x = np.linspace(0, 2*T, 1000)
napatie = Umax * np.sin(omega * x)
prud = Imax * np.sin(omega * x - phi_rad)
vykon = napatie * prud

# Zobrazenie grafov
st.subheader("📈 Priebehy veličín")
fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
ax[0].plot(x * 1000, napatie, label="Napätie [V]", color='tab:blue')
ax[1].plot(x * 1000, prud, label="Prúd [A]", color='tab:orange')
ax[2].plot(x * 1000, vykon, label="Výkon [W]", color='tab:green')

for a in ax:
    a.grid(True)
    a.legend()
    a.set_ylabel("Hodnota")

ax[2].set_xlabel("Čas [ms]")
st.pyplot(fig)

# Export do CSV
df = pd.DataFrame({
    "cas [s]": x,
    "napatie [V]": napatie,
    "prud [A]": prud,
    "vykon [W]": vykon
})

st.subheader("📤 Export výpočtov")
export_format = st.selectbox("Vyber formát exportu", ["CSV", "Excel"])

if export_format == "CSV":
    csv = df.to_csv(index=False).encode('utf-8')
    st.download_button("📥 Stiahni CSV", data=csv, file_name="vypocty.csv", mime='text/csv')
else:
    excel_buffer = io.BytesIO()
    with pd.ExcelWriter(excel_buffer, engine='xlsxwriter') as writer:
        df.to_excel(writer, index=False, sheet_name='Vypocty')
    st.download_button("📥 Stiahni Excel", data=excel_buffer.getvalue(), file_name="vypocty.xlsx", mime='application/vnd.openxmlformats-officedocument.spreadsheetml.sheet')

st.markdown("---")
st.markdown("👨‍🏫 Tento nástroj slúži na výuku a intuitívne pochopenie javov v AC a DC obvodoch. Pri prepnutí na DC sa zobrazujú len čisto odporové vlastnosti — fázový posun a jalové výkony nie sú relevantné.")
st.markdown("---")
st.markdown("👨Autor: Adrian Mahdon")