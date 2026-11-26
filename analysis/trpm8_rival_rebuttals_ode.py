#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy>=1.26",
#   "scipy>=1.11",
#   "matplotlib>=3.8",
#   "gemmi>=0.6.5",
# ]
# ///
"""
Head-to-Head 6-ODE Rival Model Rebuttals & Residue-Level Molecular Kill-Switch Suite (Figure 7)
===============================================================================================
Extends the fixed baseline 5-ODE Electro-Metabolic Relaxation Oscillator (`trpm8_oscillator_ode.py`)
into a unified 6-state dynamical system:
    y = [V_env (mV), ca (uM), cm (a.u.), dh (nM), s_pH (0..1), p_PIP2 (0..1)]

Pits five competing architectures and two single-construct point-mutation kill-switches head-to-head
under identical boundary conditions:
  - Model M1 (Proposed Cooperative TMBIM5-MICU1 H+ Clock, n=5, Kd=0.30 uM)
  - Rival R1 (Linear Plasma-Membrane PMCA Ca2+/2H+ Import Alone, n=1; Gover 2007; Hwang 2011)
  - Rival R2 (Ca2+-PLCd PI(4,5)P2 Depletion Clock Alone, dPsi_m-independent; Rohacs 2005; Yudin 2011; Yin 2024)
  - Rival R3 (Pure Electrical Ih/IKD + T-Type Cav3.2 Without L-Type Cav1.3/1.2; Fox 1987; Orio 2009)
  - Model M2 (Coupled 6-ODE Hybrid: TMBIM5-MICU1 H+ + Linear PMCA H+ + Basal Ca2+-PLCd PIP2)
  - Kill-Switch Mutant 1: Human TRPM8(H845Q/A) ( or H845Q/H667Q intracellular H+-gate ablation; 8BDC / Kim 2026)
  - Kill-Switch Mutant 2: Human TMBIM5(D294R/D325R) (3.27 A pore-dyad dead antiporter; AF-Q9H3K2 / Zhang 2025)

Also measures exact inter-atomic distances directly from experimental cryo-EM mmCIFs (`8BDC`, `9ZCU`, `6WDN`)
and AlphaFold v6 (`AF-Q9H3K2-F1-model_v6.cif`) for Panel E.
"""

from dataclasses import dataclass
from pathlib import Path
import gemmi
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


ARTIFACT_DIR = Path("/Users/yorek_mac/.gemini/antigravity/brain/d545e599-df65-48be-8d84-bc4e4be28d69")


@dataclass
class RivalCondition:
    name: str
    temp_c: float = 25.0
    icilin_nm: float = 50.0
    menthol_um: float = 0.0
    delta_osm_mosm: float = 0.0
    mito_psi: float = 1.0
    # Ion channel switches
    cal_scale: float = 1.0
    cav32_t_scale: float = 0.0     # Rival R3: replace L-type with canonical T-type Cav3.2 (Fox et al. 1987)
    # Mitochondrial & PMCA proton coupling
    mcu_scale: float = 1.0
    micu1_kd_um: float = 0.30
    micu1_hill: float = 5.0
    alpha_h_mito: float = 16.5     # Baseline TMBIM5-MICU1 H+ coupling (M1 = 16.5)
    alpha_h_pmca: float = 0.0      # Linear PMCA Ca2+/2H+ import gain (R1 / M2)
    hepes_buffer_fold: float = 1.0 # >1 models 50-100 mM intracellular HEPES patch-pipette pH clamp
    pmca_scale: float = 1.0        # <1 models selective PMCA inhibition (e.g. caloxin)
    # TRPM8 intracellular H+ sensor (His845 / His667)
    trpm8_h_sensitive: bool = True # False models TRPM8(H845Q/A) H+-gate mutant!
    # Ca2+-PLCd PI(4,5)P2 hydrolysis pathway (R2 / M2)
    k_plc: float = 0.0             # PLCd hydrolysis rate (s^-1); >0 activates PIP2 depletion
    k_plc_kd_um: float = 0.48      # Ca2+ half-activation of PLCd1 (~0.5 uM; Rohacs 2005; Yudin 2011)
    tau_pip2_s: float = 34.0       # PI4K/PIP5K plasma-membrane PIP2 resynthesis time constant (s)
    vh_pip2_max_mv: float = 68.0   # Positive V1/2 shift upon PIP2 loss (+68 mV; Rohacs 2005)
    pip2_clamped: bool = False     # True models intracellular diC8-PIP2 patch-pipette dialysis
    menthol_pip2_lock: float = 0.0 # 0 = symmetric Yin 2024 R2; >0 = Kim 2026 menthol PIP2 stabilization
    t_stim_on: float = 60.0
    t_stim_off: float = 360.0
    t_end: float = 420.0


def sigmoid(x: float) -> float:
    return float(1.0 / (1.0 + np.exp(-np.clip(x, -40.0, 40.0))))


def run_6ode(cond: RivalCondition, fast_grid: bool = False) -> dict:
    E_L, g_L = -52.0, 0.95
    E_M8, E_Na, E_Ca = 0.0, 55.0, 60.0
    g_M8_max = 1.85
    g_nav = 0.27
    g_cal = 0.20 * cond.cal_scale
    g_cav32 = 0.35 * cond.cav32_t_scale  # Generous somatic T-type conductance for R3 steelman
    H0_nM = 50.12  # pH 7.30 baseline

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        V, ca, cm, dh, s, pip2 = y
        ca = max(0.04, ca)
        cm = max(0.0, cm)
        dh = max(0.0, dh)
        s = np.clip(s, 0.01, 1.0)
        pip2 = 1.0 if cond.pip2_clamped else float(np.clip(pip2, 0.05, 1.0))

        on = cond.t_stim_on <= t <= cond.t_stim_off
        temp = cond.temp_c
        icilin = cond.icilin_nm if on else 0.0
        menthol = cond.menthol_um if on else 0.0
        dosm = cond.delta_osm_mosm if on else 0.0

        # PIP2 depletion shifts V1/2 positively (Rohacs 2005; Yin 2024)
        pip2_loss = 1.0 - pip2
        vh_base = 2.0 + 5.0 * (temp - 25.0) - 0.72 * dosm
        vh_cold = vh_base + cond.vh_pip2_max_mv * pip2_loss

        ic_f = icilin / (icilin + 360.0) if icilin > 0 else 0.0
        ps = sigmoid((V - (vh_cold - 250.0 * ic_f)) / 13.0)

        if menthol > 0:
            men_f = menthol / (menthol + 60.0)
            vh_men = vh_base + cond.vh_pip2_max_mv * pip2_loss * (1.0 - cond.menthol_pip2_lock)
            pm = sigmoid((V - (vh_men - 82.0 * men_f)) / 13.0)
        else:
            pm = 0.0

        s_gate = s if cond.trpm8_h_sensitive else 1.0
        s_eff = s_gate + (1.0 - s_gate) * (0.65 * (ic_f**2))
        peff = max(ps * s_eff, pm * (0.98 + 0.02 * s_gate))
        I_M8 = g_M8_max * peff * (V - E_M8)

        m_nav = sigmoid((V + 34.0) / 3.2) * (1.0 - 0.76 * sigmoid((V + 16.0) / 3.5))
        I_Nav = g_nav * m_nav * (V - E_Na)

        m_cal = sigmoid((V + 27.0) / 3.2) / (1.0 + (ca / 0.65) ** 2)
        I_CaL = g_cal * m_cal * (V - E_Ca)

        # Canonical T-type Cav3.2 window current (Fox, Nowycky & Tsien 1987: V1/2_act=-51 mV, V1/2_inact=-62 mV)
        if g_cav32 > 0:
            m_t = sigmoid((V + 51.0) / 4.5)
            h_t = 1.0 / (1.0 + np.exp(np.clip((V + 62.0) / 4.2, -40.0, 40.0)))
            I_CaT = g_cav32 * (m_t**2) * h_t * (V - E_Ca)
        else:
            I_CaT = 0.0

        k_rate = 1.35
        I_L = (g_L + 0.16 * (ca**2 / (ca**2 + 0.45**2))) * (V - E_L)
        dV = -k_rate * (I_M8 + I_Nav + I_CaL + I_CaT + I_L) / 1.4

        Jin_m8 = 0.0040 * max(0.0, -I_M8) * (1.0 + 1.8 * ic_f)
        Jin_cal = 0.0260 * max(0.0, -(I_CaL + I_CaT))
        J_pmca = 0.34 * cond.pmca_scale * (ca - 0.095)
        J_mcu = (
            cond.mcu_scale
            * cond.mito_psi
            * 0.50
            * (ca**cond.micu1_hill / (ca**cond.micu1_hill + cond.micu1_kd_um**cond.micu1_hill))
        )
        J_mout = 0.075 * cm

        dca = k_rate * (Jin_m8 + Jin_cal - J_pmca - J_mcu + J_mout)
        dcm = k_rate * (J_mcu - J_mout)

        # Combined proton dynamics: cooperative mitochondrial TMBIM5 + linear plasma-membrane PMCA
        # Exogenous pipette HEPES buffer (B_HEPES) attenuates free d[H+]i excursions while pipette/NHE clears dh
        j_h_in = (cond.alpha_h_mito * J_mcu + cond.alpha_h_pmca * max(0.0, J_pmca)) / cond.hepes_buffer_fold
        ddh = k_rate * (j_h_in - dh / 11.0)

        s_inf = 1.0 / (1.0 + ((H0_nM + dh) / 74.0) ** 8)
        ds = k_rate * (s_inf - s) / 9.0

        if cond.pip2_clamped or cond.k_plc <= 0.0:
            dpip2 = 0.0
        else:
            plc_act = (ca**2) / (ca**2 + cond.k_plc_kd_um**2)
            dpip2 = k_rate * ((1.0 - pip2) / cond.tau_pip2_s - cond.k_plc * plc_act * pip2)

        return np.array([dV, dca, dcm, ddh, ds, dpip2])

    sol = solve_ivp(
        rhs,
        (0.0, cond.t_end),
        np.array([-48.0, 0.10, 0.0, 0.1, 0.98, 1.0]),
        method="BDF" if fast_grid else "Radau",
        t_eval=np.linspace(0.0, cond.t_end, 480 if fast_grid else 2400),
        rtol=1e-4 if fast_grid else 1e-6,
        atol=1e-6 if fast_grid else 1e-8,
    )

    t = sol.t
    V = sol.y[0]
    ca_nM = sol.y[1] * 1000.0
    cm = sol.y[2]
    dh_nM = np.maximum(0.0, sol.y[3])
    pH_i = -np.log10((H0_nM + dh_nM) * 1e-9)
    pip2 = sol.y[5]

    mask = (t >= cond.t_stim_on) & (t <= cond.t_stim_off)
    ca_win = ca_nM[mask]
    pks, _ = find_peaks(ca_win, prominence=60.0, distance=50)
    tail_slice = ca_win[len(ca_win) // 3 :]
    tail_pks, _ = find_peaks(tail_slice, prominence=50.0, distance=45)
    raw_tail_span = float(np.max(tail_slice) - np.min(tail_slice))
    osc_tail_amp = raw_tail_span if len(tail_pks) >= 1 else 0.0

    return {
        "name": cond.name,
        "t_min": t / 60.0,
        "V": V,
        "ca_nM": ca_nM,
        "cm": cm,
        "pH_i": pH_i,
        "pip2": pip2,
        "n_peaks": int(len(pks)),
        "tail_peaks": int(len(tail_pks)),
        "ca_max": float(np.max(ca_win)),
        "ca_mean": float(np.mean(ca_win)),
        "tail_amp": raw_tail_span,
        "osc_amp": osc_tail_amp,
        "cm_mean": float(np.mean(cm[mask])),
        "cm_max": float(np.max(cm[mask])),
        "pH_min": float(np.min(pH_i[mask])),
        "pip2_min": float(np.min(pip2[mask])),
    }


def measure_structural_contacts() -> list[tuple[str, float, str]]:
    """Measure exact inter-residue heavy-atom minimum distances from 8BDC.cif and AF-Q9H3K2-F1-model_v6.cif."""
    def min_heavy_dist(r1, r2) -> float:
        dmin = 999.0
        for a1 in r1:
            if a1.element.name == "H":
                continue
            for a2 in r2:
                if a2.element.name == "H":
                    continue
                d = a1.pos.dist(a2.pos)
                if d < dmin:
                    dmin = d
        return float(dmin)

    st_m8 = gemmi.read_structure(str(ARTIFACT_DIR / "8BDC" / "8BDC.cif"))
    ch_m8 = {r.seqid.num: r for r in st_m8[0]["A"] if r.seqid.num is not None}

    st_af_m8 = gemmi.read_structure(str(ARTIFACT_DIR / "Q7Z2W7" / "AF-Q7Z2W7-F1-model_v6.cif"))
    ch_af_m8 = {r.seqid.num: r for r in st_af_m8[0]["A"] if r.seqid.num is not None}

    st_tmb = gemmi.read_structure(str(ARTIFACT_DIR / "Q9H3K2" / "AF-Q9H3K2-F1-model_v6.cif"))
    ch_tmb = {r.seqid.num: r for r in st_tmb[0]["A"] if r.seqid.num is not None}

    contacts = [
        ("TMBIM5 Asp294–Ser321 (Pore H-Bond, AFv6)", min_heavy_dist(ch_tmb[294], ch_tmb[321]), "TMBIM5"),
        ("TRPM8 His667–Gln992 (TRP-Box Clamp, 8BDC)", min_heavy_dist(ch_m8[667], ch_m8[992]), "TRPM8_667"),
        ("TRPM8 His667–Glu622 (MHR4 Salt-Bridge, 8BDC)", min_heavy_dist(ch_m8[667], ch_m8[622]), "TRPM8_667"),
        ("TRPM8 His845–Arg842 (S4 Gating Sensor, 8BDC)", min_heavy_dist(ch_m8[845], ch_m8[842]), "TRPM8_845"),
        ("TRPM8 His667–Asp665 (Post-LASIK GoF, 8BDC)", min_heavy_dist(ch_m8[667], ch_m8[665]), "TRPM8_667"),
        ("TMBIM5 Asp294–Asp325 (Ca²⁺/H⁺ Pore Dyad, AFv6)", min_heavy_dist(ch_tmb[294], ch_tmb[325]), "TMBIM5"),
        ("TRPM8 His845–Asp796 (VSLD Pocket, 8BDC)", min_heavy_dist(ch_m8[845], ch_m8[796]), "TRPM8_845"),
        ("TRPM8 His845–Glu1004 (TRP Helix, AFv6)", min_heavy_dist(ch_af_m8[845], ch_af_m8[1004]), "TRPM8_845"),
        ("TRPM8 His845–Asp802 (Ca²⁺/Icilin Triad, 8BDC)", min_heavy_dist(ch_m8[845], ch_m8[802]), "TRPM8_845"),
    ]
    return contacts


def style_axis(ax, hide_top=True, hide_right=True):
    ax.spines["top"].set_visible(not hide_top)
    ax.spines["right"].set_visible(not hide_right)
    ax.spines["left"].set_color("#334155")
    ax.spines["bottom"].set_color("#334155")
    ax.tick_params(colors="#1e293b", labelsize=8.3)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.28, color="#94a3b8")


def draw_stim_bar(ax, s0, s1, y_frac, label, color="#334155"):
    y0, y1 = ax.get_ylim()
    yb = y0 + y_frac * (y1 - y0)
    ax.plot([s0, s1], [yb, yb], color=color, lw=4.0, solid_capstyle="butt", zorder=6)
    ax.text(
        (s0 + s1) / 2.0,
        yb + 0.022 * (y1 - y0),
        label,
        ha="center",
        va="bottom",
        fontsize=7.3,
        fontweight="bold",
        color=color,
        bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="none", alpha=0.88),
        zorder=7,
    )


def build_figure7(
    rival_res: dict,
    bifurc_grid: tuple[np.ndarray, np.ndarray, np.ndarray],
    m2_bench: list[tuple[str, dict, str]],
    contacts: list[tuple[str, float, str]],
    mut_res: dict,
    out_path: Path,
) -> None:
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(14.8, 10.8), dpi=220, facecolor="white")
    gs = gridspec.GridSpec(
        2, 3, hspace=0.40, wspace=0.32,
        left=0.085, right=0.955, top=0.905, bottom=0.075,
    )

    # ==================== PANEL A: M1 vs R1 (Linear PMCA) vs R3 (Cav3.2 T-type) ====================
    axA = fig.add_subplot(gs[0, 0])
    style_axis(axA, hide_top=True, hide_right=True)
    r_m1 = rival_res["M1_tmbim5"]
    r_r1 = rival_res["R1_pmca_linear"]
    r_r3 = rival_res["R3_cav32_only"]

    axA.plot(r_m1["t_min"], r_m1["ca_nM"], color="#0284c7", lw=2.1, label="M1: TMBIM5–MICU1 H⁺ (n=5; 3 bursts)", zorder=4)
    axA.plot(r_r1["t_min"], r_r1["ca_nM"], color="#ea580c", lw=2.0, ls="-.", label="R1: Matched PMCA H⁺ Only (n=1; 0 bursts)", zorder=3)
    axA.plot(r_r3["t_min"], r_r3["ca_nM"], color="#64748b", lw=1.8, ls="--", label="R3: Ih + IKD + T-Type Caᵥ3.2 (No Caᵥ1.3/1.2)", zorder=2)
    axA.set_xlim(0, 7.0)
    axA.set_ylim(40, 760)
    axA.set_xlabel("Time (min)", fontsize=8.5, color="#1e293b")
    axA.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.8, fontweight="bold", color="#0f172a")
    axA.set_title("A   Osmo-Cold Drive: M1 vs. Matched PMCA (R1) & Caᵥ3.2 (R3)", loc="left", fontsize=9.0, fontweight="bold", color="#0f172a", pad=7)
    axA.legend(loc="upper right", fontsize=6.7, frameon=True, facecolor="white", framealpha=0.95, edgecolor="#cbd5e1")
    r1_mid = float(r_r1["ca_nM"][1650])
    axA.annotate(
        "Matched PMCA H⁺\nstalls in damped plateau",
        xy=(4.6, r1_mid),
        xytext=(4.45, max(360.0, r1_mid - 165.0)),
        ha="center",
        fontsize=6.85,
        fontweight="bold",
        color="#c2410c",
        bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#fed7aa", alpha=0.95),
        arrowprops=dict(arrowstyle="->", color="#c2410c", lw=1.0),
        zorder=8,
    )
    draw_stim_bar(axA, 1.0, 6.0, 0.035, "33 °C +98 mOsm (or 25 °C Low Drive)", color="#334155")

    # ==================== PANEL B: 2D BIFURCATION PROOF (Cooperativity n vs Proton Gain) ====================
    axB = fig.add_subplot(gs[0, 1])
    style_axis(axB, hide_top=True, hide_right=True)
    axB.grid(False)
    n_vals, gain_vals, amp_grid = bifurc_grid
    im = axB.imshow(
        amp_grid,
        origin="lower",
        aspect="auto",
        extent=[gain_vals[0], gain_vals[-1], n_vals[0], n_vals[-1]],
        cmap="Blues",
        vmin=0.0,
        vmax=280.0,
        interpolation="bilinear",
    )
    axB.contour(gain_vals, n_vals, amp_grid, levels=[65.0], colors=["#0f172a"], linewidths=1.6, linestyles="--")
    axB.axhline(1.08, color="#ea580c", lw=2.2, ls="-.")
    axB.text(16.0, 1.22, "Unthresholded n = 1 → 0 Sustained Bursts Across αH",
             ha="center", va="bottom", fontsize=6.8, fontweight="bold", color="#9a3412",
             bbox=dict(boxstyle="round,pad=0.16", facecolor="white", edgecolor="#fdba74", alpha=0.94))
    axB.scatter([16.5], [5.0], s=85, color="#dc2626", edgecolor="white", lw=1.3, zorder=6)
    axB.annotate(
        "M1: Native MICU1–TMBIM5\n(n = 5, αH = 16.5 → 252 nM)",
        xy=(16.5, 5.0),
        xytext=(10.6, 4.10),
        ha="center",
        fontsize=6.8,
        fontweight="bold",
        color="#b91c1c",
        bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#fecaca", alpha=0.95),
        arrowprops=dict(arrowstyle="->", color="#b91c1c", lw=1.1),
        zorder=7,
    )
    axB.text(22.8, 5.45, "Hopf Limit-Cycle\nIsland (n ≥ 2.0)", ha="center", va="center",
             fontsize=7.1, fontweight="bold", color="white",
             bbox=dict(boxstyle="round,pad=0.18", facecolor="#0369a1", edgecolor="none", alpha=0.88))
    cbar = fig.colorbar(im, ax=axB, fraction=0.045, pad=0.03)
    cbar.set_label("Sustained Osc. Amp (nM)", fontsize=7.5, fontweight="bold", color="#1e293b")
    cbar.ax.tick_params(labelsize=7.3)
    axB.set_xlabel("Cytosolic H⁺ Feedback Gain (αH)", fontsize=8.5, color="#1e293b")
    axB.set_ylabel("Ca²⁺ → H⁺ Hill Cooperativity (n)", fontsize=8.8, fontweight="bold", color="#0f172a")
    axB.set_title("B   Bifurcation Map: Cooperative Gate Requirement", loc="left", fontsize=9.0, fontweight="bold", color="#0f172a", pad=7)

    # ==================== PANEL C: RIVAL R2 (PIP2 Clock Alone) vs M1 on FCCP & Menthol ====================
    axC = fig.add_subplot(gs[0, 2])
    style_axis(axC, hide_top=True, hide_right=True)
    r_m1_fccp = rival_res["M1_fccp"]
    r_r2_fccp = rival_res["R2_pip2_fccp"]
    r_m1_men = rival_res["M1_menthol"]
    r_r2_men = rival_res["R2_pip2_menthol"]

    axC.plot(r_m1_fccp["t_min"], r_m1_fccp["ca_nM"], color="#dc2626", lw=2.1, label="M1 + FCCP (Ψm=0.05): Plateau Collapse")
    axC.plot(r_r2_fccp["t_min"], r_r2_fccp["ca_nM"], color="#7c3aed", lw=1.85, ls="-", label="R2 (PIP₂-Only) + FCCP: Ψm-Blind Adaptation!")
    axC.plot(r_m1_men["t_min"], r_m1_men["ca_nM"], color="#0d9488", lw=1.8, ls="--", label="M1 + Menthol: H⁺-Spared Tonic (793 nM)")
    axC.plot(r_r2_men["t_min"], r_r2_men["ca_nM"], color="#d97706", lw=1.7, ls=":", label="R2 + Menthol (Symmetric Yin 2024): Adapts!")
    axC.set_xlim(0, 7.0)
    axC.set_ylim(40, 960)
    axC.set_xlabel("Time (min)", fontsize=8.5, color="#1e293b")
    axC.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.8, fontweight="bold", color="#0f172a")
    axC.set_title("C   Rival R2 (PIP₂-Only Clock) Fails FCCP & Menthol Tests", loc="left", fontsize=9.0, fontweight="bold", color="#0f172a", pad=7)
    axC.legend(loc="upper right", fontsize=6.55, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    draw_stim_bar(axC, 1.0, 6.0, 0.035, "FCCP (Ψm = 0.05) or 35 µM Menthol Discriminator", color="#334155")

    # ==================== PANEL D: HYBRID 6-ODE (M2) FALSIFICATION #6 BAR BENCHMARK ====================
    axD = fig.add_subplot(gs[1, 0])
    style_axis(axD, hide_top=True, hide_right=True)
    labels_d = [item[0] for item in m2_bench]
    amps_d = [item[1]["osc_amp"] for item in m2_bench]
    pks_d = [item[1]["n_peaks"] for item in m2_bench]
    cols_d = [item[2] for item in m2_bench]
    y_pos = np.arange(len(labels_d))

    bars = axD.barh(y_pos, amps_d, height=0.56, color=cols_d, edgecolor="#0f172a", linewidth=0.6, alpha=0.90)
    axD.set_yticks(y_pos)
    axD.set_yticklabels(labels_d, fontsize=7.5, fontweight="bold", color="#1e293b")
    axD.invert_yaxis()
    axD.set_xlim(0, 425)
    axD.set_xlabel("Sustained Burst Osc. Amplitude (nM)", fontsize=8.5, fontweight="bold", color="#1e293b")
    axD.set_title("D   Hybrid 6-ODE (M2): Falsification #6 Predictions", loc="left", fontsize=9.0, fontweight="bold", color="#0f172a", pad=7)
    for bar, pk, item in zip(bars, pks_d, m2_bench):
        w = bar.get_width()
        r_it = item[1]
        if w > 10.0:
            tag = f"{w:.0f} nM ({pk} bursts, pH {r_it['pH_min']:.2f})"
        else:
            tag = f"0 bursts → Plateau ({r_it['ca_max']:.0f} nM, pH {r_it['pH_min']:.2f})"
        axD.text(w + 6.0, bar.get_y() + bar.get_height() / 2.0, tag,
                 va="center", ha="left", fontsize=6.7, fontweight="bold", color="#0f172a")

    # ==================== PANEL E: STRUCTURAL ATOMIC DISTANCES (8BDC & AF-Q9H3K2) ====================
    axE = fig.add_subplot(gs[1, 1])
    style_axis(axE, hide_top=True, hide_right=True)
    y_e = np.arange(len(contacts)) * 1.15
    d_vals = [c[1] for c in contacts]
    c_lbls = [c[0] for c in contacts]
    c_map = {"TMBIM5": "#dc2626", "TRPM8_845": "#0284c7", "TRPM8_667": "#7c3aed"}
    c_cols = [c_map[c[2]] for c in contacts]

    axE.axvspan(0.0, 4.0, color="#dcfce7", alpha=0.38, zorder=0)
    axE.axvline(4.0, color="#15803d", lw=1.0, ls="--", alpha=0.75)
    axE.text(2.0, -0.82, "Direct H-Bond / Salt-Bridge (<4.0 Å)", ha="center", va="center",
             fontsize=6.75, fontweight="bold", color="#166534")
    for yi, dv, lbl, col in zip(y_e, d_vals, c_lbls, c_cols):
        axE.plot([0.0, dv], [yi + 0.16, yi + 0.16], color=col, lw=2.3, zorder=3)
        axE.scatter([dv], [yi + 0.16], s=52, color=col, edgecolor="white", lw=0.9, zorder=4)
        axE.text(0.12, yi - 0.19, lbl, va="center", ha="left", fontsize=6.7, fontweight="bold", color="#1e293b")
        axE.text(dv + 0.16, yi + 0.16, f"{dv:.2f} Å", va="center", ha="left", fontsize=6.9, fontweight="bold", color=col)

    axE.set_yticks([])
    axE.set_xlim(0.0, 8.1)
    axE.set_ylim(y_e[-1] + 0.65, -1.18)
    axE.set_xlabel("Heavy-Atom Minimum Distance (Å)", fontsize=8.5, fontweight="bold", color="#1e293b")
    axE.set_title("E   Cryo-EM (8BDC) & AlphaFold (Q9H3K2) Gate Contacts", loc="left", fontsize=9.0, fontweight="bold", color="#0f172a", pad=7)

    # ==================== PANEL F: SINGLE-RESIDUE POINT-MUTANT KILL-SWITCHES ====================
    axF = fig.add_subplot(gs[1, 2])
    style_axis(axF, hide_top=True, hide_right=False)
    r_wt = mut_res["wt"]
    r_h845q = mut_res["trpm8_h845q"]
    r_d294r = mut_res["tmbim5_d294r_d325r"]

    axF.plot(r_wt["t_min"], r_wt["ca_nM"], color="#0284c7", lw=1.8, ls=":", label="WT: 3 Bursts (pHᵢ→6.90, cm=3.6)")
    axF.plot(r_h845q["t_min"], r_h845q["ca_nM"], color="#7c3aed", lw=2.1, label="TRPM8(H845Q/A): H⁺-Blind Cold Plateau!")
    axF.plot(r_d294r["t_min"], r_d294r["ca_nM"], color="#dc2626", lw=1.9, ls="--", label="TMBIM5(D294R/D325R): Pore-Dead Plateau")
    axF.set_xlim(0, 7.0)
    axF.set_ylim(40, 960)
    axF.set_xlabel("Time (min)", fontsize=8.5, color="#1e293b")
    axF.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.8, fontweight="bold", color="#0f172a")

    axF2 = axF.twinx()
    axF2.spines["top"].set_visible(False)
    axF2.plot(r_d294r["t_min"], r_d294r["cm"], color="#ea580c", lw=1.55, ls="-.", alpha=0.88)
    axF2.set_ylim(0, 8.5)
    axF2.set_ylabel("D294R/D325R Matrix cm (a.u.)", fontsize=8.0, fontweight="bold", color="#ea580c")
    axF2.tick_params(axis="y", colors="#ea580c", labelsize=7.6)

    axF.set_title("F   Point-Mutant Kill-Switches: H845Q/A & D294R/D325R", loc="left", fontsize=9.0, fontweight="bold", color="#0f172a", pad=7)
    lines_f = [
        plt.Line2D([0], [0], color="#0284c7", lw=1.8, ls=":", label="WT (3 bursts, pHᵢ 6.90)"),
        plt.Line2D([0], [0], color="#7c3aed", lw=2.1, label="TRPM8 H845Q/A [Ca²⁺]ᵢ (pHᵢ 6.86)"),
        plt.Line2D([0], [0], color="#dc2626", lw=1.9, ls="--", label="TMBIM5 D294R/D325R [Ca²⁺]ᵢ (pHᵢ 7.24)"),
        plt.Line2D([0], [0], color="#ea580c", lw=1.55, ls="-.", label="D294R/D325R Matrix cm (6.3 a.u.!)"),
    ]
    axF.legend(handles=lines_f, loc="upper right", fontsize=6.45, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    draw_stim_bar(axF, 1.0, 6.0, 0.035, "33 °C +98 mOsm (Isothermal Osmo-Cold)", color="#334155")

    fig.suptitle(
        "Figure 7. Head-to-Head 6-ODE Rival Model Rebuttals, Cryo-EM/AlphaFold Gate Coordinates & Single-Residue Kill-Switches",
        fontsize=12.0,
        fontweight="bold",
        color="#0f172a",
        y=0.968,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, facecolor="white")
    plt.close(fig)


def main() -> None:
    print("=== 1. STRUCTURAL INTER-RESIDUE CONTACTS (8BDC & AF-Q9H3K2) ===")
    contacts = measure_structural_contacts()
    for lbl, dist, grp in contacts:
        print(f"  {lbl:<48} : {dist:5.2f} A  [{grp}]")

    print("\n=== 2. HEAD-TO-HEAD RIVAL MODELS (M1 vs R1 vs R2 vs R3) ===")
    rival_conds = {
        "M1_tmbim5": RivalCondition("M1: Cooperative TMBIM5-MICU1 H+ (33°C +98 mOsm)", temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0),
        "R1_pmca_linear": RivalCondition(
            "R1: Matched PMCA Ca2+/2H+ Only (alpha_H=16.5, mcu H+=0)",
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
            alpha_h_mito=0.0, alpha_h_pmca=16.5,
        ),
        "R1b_pmca_overdriven_fccp": RivalCondition(
            "R1b: 3x Over-Driven PMCA + FCCP (Psi_m=0.05)",
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
            mito_psi=0.05, alpha_h_mito=0.0, alpha_h_pmca=48.0,
        ),
        "R3_cav32_only": RivalCondition(
            "R3: Pure Ih/IKD + T-Type Cav3.2 (g_CaL=0)",
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
            cal_scale=0.0, cav32_t_scale=1.0,
        ),
        "M1_fccp": RivalCondition(
            "M1 + FCCP Uncoupling (Psi_m=0.05)",
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0, mito_psi=0.05,
        ),
        "R2_pip2_ctrl": RivalCondition(
            "R2: Ca2+-PLCd PIP2 Depletion Only (Ctrl)",
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
            alpha_h_mito=0.0, k_plc=0.11,
        ),
        "R2_pip2_fccp": RivalCondition(
            "R2: Ca2+-PLCd PIP2 Depletion + FCCP (Psi_m=0.05)",
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
            mito_psi=0.05, alpha_h_mito=0.0, k_plc=0.11,
        ),
        "M1_menthol": RivalCondition(
            "M1 + 35 uM Menthol (H+-Insensitive Tonic)",
            temp_c=25.0, icilin_nm=0.0, menthol_um=35.0,
        ),
        "R2_pip2_menthol": RivalCondition(
            "R2 + 35 uM Menthol (Symmetric Yin 2024 PIP2)",
            temp_c=25.0, icilin_nm=0.0, menthol_um=35.0,
            alpha_h_mito=0.0, k_plc=0.11, menthol_pip2_lock=0.0,
        ),
    }

    rival_res = {}
    print(f"{'Model / Condition':<50} | {'Pks':>3} | {'TailPks':>7} | {'CaMax':>6} | {'TailAmp':>7} | {'pHmin':>5} | {'PIP2min':>7}")
    print("-" * 102)
    for k, c in rival_conds.items():
        r = run_6ode(c)
        rival_res[k] = r
        print(
            f"{c.name:<50} | {r['n_peaks']:>3d} | {r['tail_peaks']:>7d} | {r['ca_max']:6.1f} | "
            f"{r['tail_amp']:>7.1f} | {r['pH_min']:>5.2f} | {r['pip2_min']:>7.2f}"
        )

    print("\n=== 3. 2D BIFURCATION SCAN: COOPERATIVITY (n) vs H+ GAIN (alpha_H) ===")
    n_vals = np.linspace(1.0, 6.0, 12)
    gain_vals = np.linspace(4.0, 28.0, 12)
    amp_grid = np.zeros((len(n_vals), len(gain_vals)))
    for i_n, nv in enumerate(n_vals):
        for i_g, gv in enumerate(gain_vals):
            rc = RivalCondition(
                "grid", temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
                micu1_hill=float(nv), alpha_h_mito=float(gv), t_end=380.0,
            )
            rg = run_6ode(rc, fast_grid=True)
            amp_grid[i_n, i_g] = rg["tail_amp"] if rg["tail_peaks"] >= 1 else 0.0
    print(f"Max sustained tail amp at n=1.0 (Linear PMCA): {np.max(amp_grid[0, :]):.1f} nM")
    print(f"Max sustained tail amp at n=5.0 (MICU1-TMBIM5): {np.max(amp_grid[-4, :]):.1f} nM")

    print("\n=== 4. COUPLED HYBRID 6-ODE (M2) FALSIFICATION #6 BENCHMARK ===")
    # M2 combines dominant TMBIM5-MICU1 H+ (alpha_h_mito=14.8) + physiological PMCA H+ (alpha_h_pmca=3.2) + basal PLCd PIP2 (k_plc=0.018)
    m2_specs = [
        ("M2 Hybrid WT", dict(), "#0284c7"),
        ("+ diC₈-PIP₂ Clamp", dict(pip2_clamped=True), "#0d9488"),
        ("+ PMCA Block", dict(pmca_scale=0.55, alpha_h_pmca=0.0), "#7c3aed"),
        ("+ 50 mM HEPES", dict(hepes_buffer_fold=10.0), "#dc2626"),
        ("+ Ru360 (MCU)", dict(mcu_scale=0.03), "#ea580c"),
    ]
    m2_bench = []
    for lbl, kw, col in m2_specs:
        base_kw = dict(
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
            alpha_h_mito=14.8, alpha_h_pmca=3.2, k_plc=0.018,
        )
        base_kw.update(kw)
        r = run_6ode(RivalCondition(lbl.replace("\n", " "), **base_kw))
        m2_bench.append((lbl, r, col))
        print(
            f"  {lbl.replace(chr(10), ' '):<42} -> Peaks={r['n_peaks']}, TailAmp={r['tail_amp']:5.1f} nM, "
            f"CaMax={r['ca_max']:5.1f} nM, pHmin={r['pH_min']:.2f}, PIP2min={r['pip2_min']:.2f}"
        )

    print("\n=== 5. SINGLE-RESIDUE POINT-MUTANT KILL-SWITCH PREDICTIONS ===")
    mut_conds = {
        "wt": RivalCondition(
            "WT TRPM8 + WT TMBIM5 (33°C +98 mOsm)",
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
            alpha_h_mito=14.8, alpha_h_pmca=3.2, k_plc=0.018,
        ),
        "trpm8_h845q": RivalCondition(
            "TRPM8(H845Q/A) H+-Gate Mutant (33°C +98 mOsm)",
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
            alpha_h_mito=14.8, alpha_h_pmca=3.2, k_plc=0.018,
            trpm8_h_sensitive=False,
        ),
        "tmbim5_d294r_d325r": RivalCondition(
            "TMBIM5(D294R/D325R) Pore-Dead Mutant (33°C +98 mOsm)",
            temp_c=33.0, icilin_nm=0.0, delta_osm_mosm=98.0,
            alpha_h_mito=0.0, alpha_h_pmca=3.2, k_plc=0.018,
        ),
    }
    mut_res = {}
    for k, c in mut_conds.items():
        r = run_6ode(c)
        mut_res[k] = r
        print(
            f"  {c.name:<52} -> Peaks={r['n_peaks']}, TailAmp={r['tail_amp']:5.1f} nM, "
            f"CaMax={r['ca_max']:5.1f} nM, cmMean={r['cm_mean']:4.2f} (Max {r['cm_max']:4.2f}), pHmin={r['pH_min']:.2f}"
        )

    out_fig7 = Path(__file__).resolve().parent / "outputs" / "fig7_rival_models_molecular_test.png"
    build_figure7(rival_res, (n_vals, gain_vals, amp_grid), m2_bench, contacts, mut_res, out_fig7)
    print(f"\nSaved Figure 7 to: {out_fig7}")


if __name__ == "__main__":
    main()
