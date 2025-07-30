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
from scipy.integrate import odeint

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
# === Súčiastky s jednotkovými násobkami a spätným zobrazením ===
st.sidebar.markdown("---")
st.sidebar.markdown("🧩 **Zadanie súčiastok**")

# Mapovanie predpôn na hodnoty
prefix_dict = {
    "p": ("piko", 1e-12),
    "n": ("nano", 1e-9),
    "µ": ("mikro", 1e-6),
    "m": ("mili", 1e-3),
    "": ("", 1),
    "k": ("kilo", 1e3),
    "M": ("mega", 1e6),
    "G": ("giga", 1e9),
    "T": ("tera", 1e12),
}

# Invertovaný slovník pre formátovanie
def find_prefix(val):
    if val == 0:
        return "0", ""
    abs_val = abs(val)
    for sym, (_, factor) in reversed(prefix_dict.items()):
        if abs_val >= factor:
            scaled = val / factor
            if scaled >= 0.1:
                return f"{scaled:.3g}", sym
    return f"{val:.3g}", ""

# Vstup s prefixom
def vstup_so_skalou(label, jednotka, default_val=0.0):
    col1, col2 = st.sidebar.columns([2, 1])
    hodnota = col1.number_input(f"{label} [{jednotka}]", value=default_val, step=0.1, key=label)
    prefix_list = [f"{sym} ({prefix_dict[sym][0]})" if sym else "(základná)" for sym in prefix_dict]
    symbol_map = {f"{sym} ({prefix_dict[sym][0]})" if sym else "(základná)": sym for sym in prefix_dict}
    default_prefix = "(základná)"
    col2_key = f"{label}_scale"
    prefix_label = col2.selectbox("×", prefix_list, index=prefix_list.index(default_prefix), key=col2_key)
    symbol = symbol_map[prefix_label]
    return hodnota * prefix_dict[symbol][1]

# Zadanie a výpočet R, L, C
R = vstup_so_skalou("Odpor R", "Ω")
L = vstup_so_skalou("Indukčnosť L", "H")
C = vstup_so_skalou("Kapacita C", "F")

# === Zobrazenie spätne formátovaných hodnôt ===
st.markdown("### 🔍 Zadali ste:")
r_str, r_prefix = find_prefix(R)
l_str, l_prefix = find_prefix(L)
c_str, c_prefix = find_prefix(C)
if R > 0 or L > 0 or C > 0:
    st.markdown(f"""
    - **Odpor R:** {r_str} {r_prefix}Ω  
    - **Indukčnosť L:** {l_str} {l_prefix}H  
    - **Kapacita C:** {c_str} {c_prefix}F  
    """)
else:
    st.info("Zatiaľ nebola zadaná žiadna súčiastka R, L ani C.")

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
    t_points = st.sidebar.number_input("Počet bodov", value=10000, step=10)
else:
    t_max = 0.001
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

if input_mode == "Efektívvne (RMS)":
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
#-----------------------------------------------------------------

#----------------------------------------
# Časová os a simulácie
x = np.linspace(0, t_max, int(t_points))
annotation_time = None
if type_choice == "AC":
    T = 1 / f if f > 0 else 1.0
    x = np.linspace(0, 2*T, 1000)
    u = Umax * np.sin(omega * x)
    i = Imax * np.sin(omega * x - phi_calc_rad)
elif type_choice == "DC":
    u = np.full_like(x, Uef)
    i = np.full_like(x, Ief)
    tau = None
elif type_choice == "DC - Prechodový dej (R-C / R-L)":
    # Vstupy
    st.sidebar.header("Parametre obvodu")
    #voltage_type = st.sidebar.selectbox("Typ napätia", ["AC", "DC", "DC - Prechodový dej"])

    U = st.sidebar.number_input("Napätie [V]", value=10.0)
    I = st.sidebar.number_input("Prúd [A] (len pre zobrazenie výkonu)", value=1.0)
    R = st.sidebar.number_input("Odpor R [Ω]", value=10.0)
    L = st.sidebar.number_input("Indukčnosť L [H]", value=0.0)
    C = st.sidebar.number_input("Kapacita C [F]", value=0.0)

    annotation_time = st.sidebar.number_input("Čas pre anotáciu výstupov [s]", value=1.0)


    # Časová os podľa tau:
    def calculate_dynamic_tmax(R, L, C):
        tau = 1.0
        if R > 0 and L > 0 and C == 0:
            tau = L / R
        elif R > 0 and C > 0 and L == 0:
            tau = R * C
        elif L > 0 and C > 0:
            tau = 2 * np.sqrt(L * C)
        return 5 * tau if tau > 0 else 5.0


    T_dynamic = calculate_dynamic_tmax(R, L, C)
    T_max = st.sidebar.slider("Čas simulácie [s]", min_value=0.1, max_value=float(T_dynamic * 2),
                              value=float(T_dynamic), step=0.1)
    t = np.linspace(0, T_max, 2000)

    # Výpočty a simulácia:
#if voltage_type == "DC - Prechodový dej":
    #   st.subheader("DC - Prechodový dej")

    if R > 0 and L > 0 and C == 0:
        tau = L / R


        def di_dt(i, t):
            return (U - R * i) / L


        i = odeint(lambda i, t: di_dt(i, t), 0, t).flatten()
        u_R = R * i
        u_L = U - u_R
        P_R = u_R * i
        P_L = u_L * i
        P_total = U * i

        st.write(f"\n**RL obvod:** τ = L/R = {tau:.4f} s")
        st.write(f"Prúd v čase {annotation_time}s: {np.interp(annotation_time, t, i):.4f} A")
        st.write(f"Napätie na R: {np.interp(annotation_time, t, u_R):.4f} V")
        st.write(f"Napätie na L: {np.interp(annotation_time, t, u_L):.4f} V")
        st.write(f"Výkon na R: {np.interp(annotation_time, t, P_R):.4f} W")
        st.write(f"Výkon na L: {np.interp(annotation_time, t, P_L):.4f} W")
        st.write(f"Celkový výkon: {np.interp(annotation_time, t, P_total):.4f} W")

        fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        ax[0].plot(t, i, label='Prúd [A]')
        ax[0].legend();
        ax[0].grid(True)
        ax[1].plot(t, u_R, label='Napätie na R [V]')
        ax[1].plot(t, u_L, label='Napätie na L [V]')
        ax[1].legend();
        ax[1].grid(True)
        st.pyplot(fig)

    elif R > 0 and C > 0 and L == 0:
        tau = R * C


        def uc_t(uC, t):
            return (U - uC) / (R * C)


        u_C = odeint(lambda uC, t: uc_t(uC, t), 0, t).flatten()
        i = (U - u_C) / R
        u_R = U - u_C
        P_R = u_R * i
        P_C = u_C * i
        P_total = U * i

        st.write(f"\n**RC obvod:** τ = RC = {tau:.4f} s")
        st.write(f"Prúd v čase {annotation_time}s: {np.interp(annotation_time, t, i):.4f} A")
        st.write(f"Napätie na R: {np.interp(annotation_time, t, u_R):.4f} V")
        st.write(f"Napätie na C: {np.interp(annotation_time, t, u_C):.4f} V")
        st.write(f"Výkon na R: {np.interp(annotation_time, t, P_R):.4f} W")
        st.write(f"Výkon na C: {np.interp(annotation_time, t, P_C):.4f} W")
        st.write(f"Celkový výkon: {np.interp(annotation_time, t, P_total):.4f} W")

        fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        ax[0].plot(t, i, label='Prúd [A]')
        ax[0].legend();
        ax[0].grid(True)
        ax[1].plot(t, u_R, label='Napätie na R [V]')
        ax[1].plot(t, u_C, label='Napätie na C [V]')
        ax[1].legend();
        ax[1].grid(True)
        st.pyplot(fig)

    elif L > 0 and C > 0:
        def rlc_series(y, t):
            uC, i = y
            duC_dt = i / C
            di_dt = (U - uC - R * i) / L
            return [duC_dt, di_dt]


        y0 = [0.0, 0.0]
        sol = odeint(rlc_series, y0, t)
        u_C = sol[:, 0]
        i = sol[:, 1]
        u_R = R * i
        u_L = U - u_C - u_R
        P_R = u_R * i
        P_C = u_C * i
        P_L = u_L * i
        P_total = U * i

        st.write("\n**RLC obvod - simulácia numerickým riešením 2. rádu**")
        st.write(f"Prúd v čase {annotation_time}s: {np.interp(annotation_time, t, i):.4f} A")
        st.write(f"Napätie na R: {np.interp(annotation_time, t, u_R):.4f} V")
        st.write(f"Napätie na L: {np.interp(annotation_time, t, u_L):.4f} V")
        st.write(f"Napätie na C: {np.interp(annotation_time, t, u_C):.4f} V")
        st.write(f"Celkový výkon: {np.interp(annotation_time, t, P_total):.4f} W")

        fig, ax = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        ax[0].plot(t, i, label='Prúd [A]')
        ax[0].legend();
        ax[0].grid(True)
        ax[1].plot(t, u_R, label='Napätie na R [V]')
        ax[1].plot(t, u_L, label='Napätie na L [V]')
        ax[1].plot(t, u_C, label='Napätie na C [V]')
        ax[1].legend();
        ax[1].grid(True)
        st.pyplot(fig)

    else:
        st.warning("Zadaj aspoň dve vhodné súčiastky (napr. R a L, R a C alebo L a C)")

else:
    st.info("Zvoľ režim \"DC - Prechodový dej\" pre simuláciu dynamiky zapínania obvodov.")

#vykon = u * i
vykon_avg = np.mean(P_total)

# Doplnková informácia o τ (časová konštanta)
if type_choice.startswith("DC") and tau is not None:
    st.markdown(f"**Časová konštanta τ =** {tau:.4f} s")

# Zobrazenie bodu, kedy sa kondenzátor nabije na 99 %
if annotation_time:
    st.markdown(f"⚡ **Prechod ustálený do:** {annotation_time:.3f} s (≈ 5τ)")

# Graf s anotáciou
fig, ax = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
ax[0].plot(x, u, label='Napätie [V]', color='tab:blue')
ax[1].plot(x, i, label='Prúd [A]', color='tab:orange')
ax[2].plot(x, vykon, label=f'Výkon [W] ⟨P⟩={vykon_avg:.2f}', color='tab:green')

# Pridanie anotácie pre čas 5τ
if annotation_time and annotation_time <= x[-1]:
    for a in ax:
        a.axvline(annotation_time, color='red', linestyle='--', alpha=0.5)
        a.text(annotation_time, a.get_ylim()[1]*0.8, '5τ', color='red')

for a in ax:
    a.legend()
    a.grid(True)
    a.set_ylabel("Hodnota")
ax[2].set_xlabel("Čas [s]")
st.subheader("📊 Priebeh veličín v čase")
st.pyplot(fig)

# Popis prechodového deja
if type_choice == "DC - Prechodový dej (R-C / R-L)":
    if C > 0:
        st.info("Kondenzátor sa nabíja exponenciálne podľa vzťahu: \n **U(t) = U(1 - e^(-t/RC))**. \n Prúd na začiatku prudko klesá, až dosiahne nulu v ustálenom stave.")
    elif L > 0:
        st.info("Cievka spôsobí oneskorený nábeh prúdu: \n **I(t) = (U/R)(1 - e^(-Rt/L))**. \n Prúd stúpa od nuly, až sa ustáli. Napätie na cievke počas prechodu klesá.")

# Výpočtové výsledky
st.subheader("🧮 Výpočty")
st.markdown(f"""
- **Zdanlivý výkon (S):** {S:.2f} VA  
- **Činný výkon (P):** {P:.2f} W  
- **Jalový výkon (Q):** {Q:.2f} VAR  
- **Fázový posun φ:** {phi_calc_deg:.2f}°  
- **Účinník (cosφ):** {cos_phi:.3f}  
- **Uef / Ief:** {Uef:.2f} V / {Ief:.2f} A  
- **Umax / Imax:** {Umax:.2f} V / {Imax:.2f} A
""")

st.markdown("---")
st.markdown("👨Autor: Adrian Mahdon")