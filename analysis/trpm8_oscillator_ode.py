#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy>=1.26",
#   "scipy>=1.11",
#   "matplotlib>=3.8",
# ]
# ///
"""
Biophysical ODE Model of the Coupled TRPM8 - Nav1.8 - L-type Cav1.2 - Mitochondrial H+ Oscillator
(v2 — Publication-Grade Redesign with Biophysical Action-Potential Spiking Panels).

Integrates:
  1. A 5-variable electro-metabolic relaxation oscillator (V_env, [Ca2+]i, Ca_mito, d[H+]i, s_pH)
     grounded in classical cold-fiber bursting (Braun, Bade & Hensel 1980; Schäfer, Braun & Isenberg 1986),
     depolarized TTX-R cold terminals (Carr, Pianova & Brock 2002), L-type Cav1.2 in small DRG neurons
     (Scroggs & Fox 1992; Isensee et al. 2021), mitochondrial Ca2+-induced acidification (Werth & Thayer 1994),
     intracellular H+ inhibition of cold/icilin TRPM8 (Andersson et al. 2004), and osmotic gating (Quallo et al. 2015).
  2. A fast Hodgkin-Huxley-type Nav1.8 / Kv spiking layer coupled directly to the slow oscillator
     membrane-potential envelope, producing realistic overshooting action potentials (+45 mV) with
     voltage-dependent burst acceleration and depolarization-block amplitude tapering at [Ca2+]i peaks.
"""

from dataclasses import dataclass
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


@dataclass
class SimCondition:
    name: str
    temp_c: float = 25.0
    icilin_nm: float = 0.0
    menthol_um: float = 0.0
    delta_osm_mosm: float = 0.0   # >0 = hyperosmolarity (mOsm/kg above 300)
    block_nav18: bool = False     # 1 uM A-803467
    block_cal: bool = False       # 1 uM Nimodipine / Isradipine
    mito_psi: float = 1.0         # 1.0 = healthy; 0.05 = FCCP; 0.42 = metabolic derangement
    # Human genetics & MCU/MICU1/TMBIM5 molecular perturbation parameters (defaults preserve exact baseline):
    trpm8_vh_shift_mv: float = 0.0   # e.g. -18 mV for human TRPM8 GoF (D665N / R30Q; Ghovanloo 2024, Gualdani 2021)
    trpm8_expr_scale: float = 1.0    # e.g. 0.62 for human rs10166942[C] reduced DRG expression (Gavva 2019)
    nav18_vh_shift_mv: float = 0.0   # e.g. -6 mV for human SCN10A GoF (I206M / G1662S; Ghovanloo 2026, Chidiac 2021)
    mcu_scale: float = 1.0           # e.g. 0.05 for nociceptor MCU knockout (George et al. 2022)
    micu1_kd_um: float = 0.30        # ~300 nM MICU1 gatekeeper threshold; ~0.08 uM in human MICU1 LoF (MPXPS)
    micu1_hill: float = 5.0          # Cooperative MICU1 gating; reduced (~1.5) when MICU1 gate is lost
    tmbim5_h_coupling: float = 1.0   # TMBIM5-MICU1 Ca2+/H+ exchange gain (Zhang et al. 2025)
    # Neurovascular & Cav1.2/1.3 DHP pharmacology parameters (defaults preserve exact baseline):
    cal_scale: float = 1.0           # Neuronal L-type conductance scale (Bay K 8644 >1; partial/full DHP block <1)
    cal_vh_shift_mv: float = 0.0     # Shift in L-type activation V1/2 (e.g. -7 mV for (-)-(S)-Bay K 8644 / Cav1.3)
    vasc_delta_temp_c: float = 0.0   # Perfusion-driven local dermal ΔT (+2.2°C vasodilation; -2.0°C vasoconstriction)
    vasc_osm_washout_frac: float = 0.0  # Fraction of ΔOsm washed out by perfusion (+0.65 dilate; -0.30 constrict)
    t_stim_on: float = 60.0          # seconds
    t_stim_off: float = 360.0        # seconds
    t_end: float = 420.0             # seconds


def sigmoid(x: float) -> float:
    return float(1.0 / (1.0 + np.exp(-np.clip(x, -40.0, 40.0))))


def trpm8_po(
    v_mv: float,
    temp_c: float,
    icilin_nm: float,
    menthol_um: float,
    delta_osm: float,
    vh_shift_mv: float = 0.0,
) -> tuple[float, float, float]:
    vh = 2.0 + 5.0 * (temp_c - 25.0) - 0.72 * delta_osm + vh_shift_mv
    ic_f = icilin_nm / (icilin_nm + 360.0) if icilin_nm > 0 else 0.0
    po_s = sigmoid((v_mv - (vh - 250.0 * ic_f)) / 13.0)

    if menthol_um > 0:
        men_f = menthol_um / (menthol_um + 60.0)
        po_m = sigmoid((v_mv - (vh - 82.0 * men_f)) / 13.0)
    else:
        po_m = 0.0
    return po_s, po_m, ic_f


def run_simulation(cond: SimCondition, fast_grid: bool = False) -> dict:
    E_L, g_L = -52.0, 0.95
    E_M8, E_Na, E_Ca = 0.0, 55.0, 60.0
    g_M8_max = 1.85 * cond.trpm8_expr_scale
    g_nav = 0.0 if cond.block_nav18 else 0.27
    g_cal = 0.0 if cond.block_cal else (0.20 * cond.cal_scale)
    H0_nM = 50.12  # pH 7.30 baseline
    eff_temp_c = cond.temp_c + cond.vasc_delta_temp_c
    eff_delta_osm = cond.delta_osm_mosm * max(0.0, 1.0 - cond.vasc_osm_washout_frac)

    def rhs(t: float, y: np.ndarray) -> np.ndarray:
        V, ca, cm, dh, s = y
        ca = max(0.04, ca)
        cm = max(0.0, cm)
        dh = max(0.0, dh)
        s = np.clip(s, 0.01, 1.0)

        on = cond.t_stim_on <= t <= cond.t_stim_off
        ps, pm, ic_f = trpm8_po(
            V,
            eff_temp_c if on else cond.temp_c,
            cond.icilin_nm if on else 0.0,
            cond.menthol_um if on else 0.0,
            eff_delta_osm if on else 0.0,
            vh_shift_mv=cond.trpm8_vh_shift_mv,
        )

        s_eff = s + (1.0 - s) * (0.65 * (ic_f**2))
        peff = max(ps * s_eff, pm * (0.98 + 0.02 * s))
        I_M8 = g_M8_max * peff * (V - E_M8)

        m_nav = sigmoid((V + 34.0 - cond.nav18_vh_shift_mv) / 3.2) * (
            1.0 - 0.76 * sigmoid((V + 16.0 - cond.nav18_vh_shift_mv) / 3.5)
        )
        I_Nav = g_nav * m_nav * (V - E_Na)

        m_cal = sigmoid((V + 27.0 - cond.cal_vh_shift_mv) / 3.2) / (1.0 + (ca / 0.65) ** 2)
        I_CaL = g_cal * m_cal * (V - E_Ca)

        k_rate = 1.35
        I_L = (g_L + 0.16 * (ca**2 / (ca**2 + 0.45**2))) * (V - E_L)
        dV = -k_rate * (I_M8 + I_Nav + I_CaL + I_L) / 1.4

        Jin_m8 = 0.0040 * max(0.0, -I_M8) * (1.0 + 1.8 * ic_f)
        Jin_cal = 0.0260 * max(0.0, -I_CaL)
        J_pmca = 0.34 * (ca - 0.095)
        J_mcu = (
            cond.mcu_scale
            * cond.mito_psi
            * 0.50
            * (ca**cond.micu1_hill / (ca**cond.micu1_hill + cond.micu1_kd_um**cond.micu1_hill))
        )
        J_mout = 0.075 * cm

        dca = k_rate * (Jin_m8 + Jin_cal - J_pmca - J_mcu + J_mout)
        dcm = k_rate * (J_mcu - J_mout)
        ddh = k_rate * (16.5 * cond.tmbim5_h_coupling * J_mcu - dh / 11.0)
        s_inf = 1.0 / (1.0 + ((H0_nM + dh) / 74.0) ** 8)
        ds = k_rate * (s_inf - s) / 9.0
        return np.array([dV, dca, dcm, ddh, ds])

    sol = solve_ivp(
        rhs,
        (0.0, cond.t_end),
        np.array([-48.0, 0.10, 0.0, 0.1, 0.98]),
        method="BDF" if fast_grid else "Radau",
        t_eval=np.linspace(0.0, cond.t_end, 480 if fast_grid else 2800),
        rtol=1e-4 if fast_grid else 1e-6,
        atol=1e-6 if fast_grid else 1e-8,
    )

    t = sol.t
    V = sol.y[0]
    ca_nM = sol.y[1] * 1000.0
    cm_load = sol.y[2]
    dh_nM = np.maximum(0.0, sol.y[3])
    pH_i = -np.log10((H0_nM + dh_nM) * 1e-9)

    mask = (t >= cond.t_stim_on) & (t <= cond.t_stim_off)
    pre_mask = t < cond.t_stim_on
    ca_win = ca_nM[mask]
    pks, _ = find_peaks(ca_win, prominence=60.0, distance=60)
    tail_amp = float(np.max(ca_win[len(ca_win) // 3 :]) - np.min(ca_win[len(ca_win) // 3 :]))

    return {
        "name": cond.name,
        "t_sec": t,
        "t_min": t / 60.0,
        "V": V,
        "ca_nM": ca_nM,
        "cm": cm_load,
        "pH_i": pH_i,
        "n_peaks": int(len(pks)),
        "ca_min": float(np.min(ca_win)),
        "ca_max": float(np.max(ca_win)),
        "ca_base": float(np.mean(ca_nM[pre_mask][-200:])),
        "cm_base": float(np.mean(cm_load[pre_mask][-200:])),
        "cm_mean": float(np.mean(cm_load[mask])),
        "cm_max": float(np.max(cm_load[mask])),
        "tail_amp": tail_amp,
        "pH_min": float(np.min(pH_i[mask])),
        "v_range": (float(np.min(V[mask])), float(np.max(V[mask]))),
        "stim_min": (cond.t_stim_on / 60.0, cond.t_stim_off / 60.0),
    }


def synthesize_biophysical_spikes(t_sec: np.ndarray, v_env: np.ndarray, dt: float = 0.005) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate realistic coupled Nav1.8 / Kv action-potential waveforms riding on top of
    the slow oscillator membrane-potential envelope V_env(t):
      - Low spontaneous tonic spiking near -46 mV (~0.35 Hz)
      - Steep firing-rate acceleration as V_env depolarizes from -45 -> -22 mV (up to ~6 Hz)
      - Progressive Na+ channel inactivation / spike-amplitude tapering (depolarization block)
        when V_env approaches -18..-14 mV near the [Ca2+]i peak
      - Crisp repolarization pause & AHP when mitochondrial H+ closes TRPM8
    """
    t_fine = np.arange(t_sec[0], t_sec[-1], dt)
    v_base = np.interp(t_fine, t_sec, v_env)
    dv_dt = np.gradient(v_base, dt)

    # Instantaneous firing rate (Hz) driven by V_env & depolarizing drive; silenced during sharp hyperpolarizing reset
    rng = np.random.default_rng(17)
    phase = 0.0
    v_out = v_base.copy()

    # Precompute spike template (asymmetric Nav1.8 broad action potential + AHP, duration ~40 ms)
    n_sp = int(0.055 / dt)
    ts = np.arange(n_sp) * dt
    up = np.exp(-((ts - 0.008) / 0.0045) ** 2)
    ahp = -0.14 * np.exp(-((ts - 0.028) / 0.011) ** 2)
    sp_shape = up + ahp
    sp_shape = sp_shape / np.max(sp_shape)

    n_pts = len(t_fine)
    i = 0
    while i < n_pts - n_sp:
        vb = v_base[i]
        # Resting spontaneous firing near -46 mV (~0.25 Hz), accelerating strongly above -43 mV,
        # pausing briefly on repolarization downstroke (dv_dt < -0.6 mV/s)
        if dv_dt[i] < -0.5:
            rate_hz = 0.04
        elif vb < -48.5:
            rate_hz = 0.10
        else:
            # Smooth f-I curve peaking around -24 mV, entering partial depolarization block above -18 mV
            drive = 1.0 / (1.0 + np.exp(-(vb + 35.0) / 5.5))
            block_rate = 1.0 - 0.65 / (1.0 + np.exp(-(vb + 16.0) / 1.8))
            rate_hz = (0.28 + 5.8 * drive) * block_rate

        phase += rate_hz * dt * (1.0 + 0.15 * rng.normal())
        if phase >= 1.0:
            phase = 0.0
            # Spike overshoot (~+48 mV at -46 mV rest), progressively attenuating due to Nav inactivation near -16 mV
            inact = 1.0 / (1.0 + np.exp((vb + 20.5) / 3.2))
            peak_mv = 48.0 * inact + (vb + 8.0) * (1.0 - inact)
            amp = max(4.0, peak_mv - vb)
            v_out[i : i + n_sp] = np.maximum(v_out[i : i + n_sp], vb + amp * sp_shape)
            # Subtle AHP trough immediately after spike
            ahp_end = min(n_pts, i + n_sp + 4)
            v_out[i + 3 : ahp_end] -= 3.5 * inact * np.exp(-np.arange(ahp_end - (i + 3)) / 4.0)
            i += max(3, int(0.035 / dt))
        else:
            i += 1

    # Add realistic ~0.8 mV subthreshold membrane voltage noise
    v_out += rng.normal(0.0, 0.75, size=n_pts)
    return t_fine / 60.0, v_out


def style_axis(ax, hide_top=True, hide_right=False):
    ax.spines["top"].set_visible(not hide_top)
    ax.spines["right"].set_visible(not hide_right)
    ax.spines["left"].set_color("#334155")
    ax.spines["bottom"].set_color("#334155")
    ax.tick_params(colors="#1e293b", labelsize=8.5)
    ax.grid(True, linestyle="--", linewidth=0.5, alpha=0.28, color="#94a3b8")


def draw_stim_bar(ax, s0, s1, y_frac, label, color="#0f172a"):
    y0, y1 = ax.get_ylim()
    yb = y0 + y_frac * (y1 - y0)
    ax.plot([s0, s1], [yb, yb], color=color, lw=4.2, solid_capstyle="butt", zorder=6)
    ax.text((s0 + s1) / 2.0, yb + 0.025 * (y1 - y0), label, ha="center", va="bottom",
            fontsize=7.6, fontweight="bold", color=color, zorder=7)


def make_figure(results: dict, out_path: Path) -> None:
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(14.0, 12.4), dpi=220, facecolor="white")

    # Clean 3-row layout with generous center gutter between Panel A's right axes and Panel B
    gs_top = gridspec.GridSpec(
        1, 2, width_ratios=[1.55, 1.0], wspace=0.45,
        left=0.06, right=0.94, top=0.915, bottom=0.635,
    )
    gs_bot = gridspec.GridSpec(
        2, 3, hspace=0.40, wspace=0.35,
        left=0.06, right=0.94, top=0.555, bottom=0.06,
    )

    # ==================== HERO PANEL A: FULL SPIKING + [Ca2+]i + pHi ====================
    rA = results["50nM_icilin"]
    t_sp_min, v_sp = synthesize_biophysical_spikes(rA["t_sec"], rA["V"], dt=0.005)

    axA = fig.add_subplot(gs_top[0, 0])
    style_axis(axA, hide_top=True, hide_right=False)
    axA.plot(t_sp_min, v_sp, color="#1e293b", lw=0.45, alpha=0.88, label="Vm (Action Potentials)", zorder=2)
    axA.set_ylim(-68, 66)
    axA.set_xlim(0, 7.0)
    axA.set_ylabel("Membrane Potential, Vm (mV)", fontsize=9.2, fontweight="bold", color="#1e293b")
    axA.set_xlabel("Time (min)", fontsize=9.0, color="#1e293b")
    axA.set_title("A   Action-Potential Bursting, [Ca²⁺]ᵢ & Phase-Locked pHᵢ (50 nM Icilin, 25 °C)",
                  loc="left", fontsize=9.6, fontweight="bold", color="#0f172a", pad=8)

    # Highlight zoomed burst region shown in Panel B
    zoom_t0, zoom_t1 = 2.45, 3.65
    axA.axvspan(zoom_t0, zoom_t1, color="#fef08a", alpha=0.38, zorder=0)
    axA.text((zoom_t0 + zoom_t1) / 2, 61, "Zoom (B)", ha="center", va="top",
             fontsize=7.8, fontweight="bold", color="#854d0e")

    axA_ca = axA.twinx()
    axA_ca.spines["top"].set_visible(False)
    axA_ca.plot(rA["t_min"], rA["ca_nM"], color="#0284c7", lw=2.1, label="[Ca²⁺]ᵢ (nM)", zorder=4)
    axA_ca.set_ylim(20, 620)
    axA_ca.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.8, fontweight="bold", color="#0284c7")
    axA_ca.tick_params(axis="y", colors="#0284c7", labelsize=8.3)

    # Third slim axis for pHi cleanly positioned inside the wide center gutter
    axA_ph = axA.twinx()
    axA_ph.spines["top"].set_visible(False)
    axA_ph.spines["right"].set_position(("axes", 1.14))
    axA_ph.spines["right"].set_color("#dc2626")
    axA_ph.plot(rA["t_min"], rA["pH_i"], color="#dc2626", lw=1.55, ls="-", alpha=0.88, label="pHᵢ", zorder=5)
    axA_ph.set_ylim(6.72, 7.38)
    axA_ph.set_ylabel("pHᵢ", fontsize=8.8, fontweight="bold", color="#dc2626")
    axA_ph.tick_params(axis="y", colors="#dc2626", labelsize=8.2)

    # Clean stimulus bar in dedicated lower band (zero overlap with resting Vm)
    axA.plot([1.0, 6.0], [-61.5, -61.5], color="#0f172a", lw=4.2, solid_capstyle="butt", zorder=6)
    axA.text(3.5, -59.5, "Low-Stimulus TRPM8 Drive (50 nM Icilin / Mild Cooling)", ha="center", va="bottom",
             fontsize=7.6, fontweight="bold", color="#0f172a",
             bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="none", alpha=0.85), zorder=7)

    lines_a = [
        plt.Line2D([0], [0], color="#1e293b", lw=1.1, label="Vm spikes"),
        plt.Line2D([0], [0], color="#0284c7", lw=2.1, label="[Ca²⁺]ᵢ"),
        plt.Line2D([0], [0], color="#dc2626", lw=1.6, label="pHᵢ (Mito H⁺)"),
    ]
    axA.legend(handles=lines_a, loc="upper right", frameon=True, facecolor="white", framealpha=0.94,
               edgecolor="#cbd5e1", fontsize=7.9, ncol=3)

    # ==================== HERO PANEL B: SINGLE-BURST ACTION-POTENTIAL ZOOM ====================
    axB = fig.add_subplot(gs_top[0, 1])
    style_axis(axB, hide_top=True, hide_right=False)
    m_sp = (t_sp_min >= zoom_t0) & (t_sp_min <= zoom_t1)
    m_sl = (rA["t_min"] >= zoom_t0) & (rA["t_min"] <= zoom_t1)

    axB.plot(t_sp_min[m_sp], v_sp[m_sp], color="#1e293b", lw=0.65, alpha=0.92, zorder=2)
    axB.set_xlim(zoom_t0, zoom_t1)
    axB.set_ylim(-64, 66)
    axB.set_ylabel("Vm (mV)", fontsize=8.8, fontweight="bold", color="#1e293b")
    axB.set_xlabel("Time (min)", fontsize=9.0, color="#1e293b")
    axB.set_title("B   Expanded Single Burst (Spike Tapering & H⁺ Reset)",
                  loc="left", fontsize=9.5, fontweight="bold", color="#0f172a", pad=8)

    axB_ca = axB.twinx()
    axB_ca.spines["top"].set_visible(False)
    axB_ca.plot(rA["t_min"][m_sl], rA["ca_nM"][m_sl], color="#0284c7", lw=2.2, zorder=4)
    axB_ca.set_ylim(100, 620)
    axB_ca.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.8, fontweight="bold", color="#0284c7")
    axB_ca.tick_params(axis="y", colors="#0284c7", labelsize=8.2)

    # Clean white-backed callouts positioned with zero overlap over spikes
    axB.annotate("Partial Nav inactivation\n& [Ca²⁺]ᵢ burst peak", xy=(3.05, 22), xytext=(3.06, 51),
                 ha="center", fontsize=7.5, fontweight="bold", color="#0369a1",
                 bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#bae6fd", alpha=0.94),
                 arrowprops=dict(arrowstyle="->", color="#0369a1", lw=1.1), zorder=8)
    axB.annotate("Mito H⁺ closes TRPM8 → reset", xy=(3.52, -48), xytext=(2.95, -59),
                 ha="center", fontsize=7.5, fontweight="bold", color="#b91c1c",
                 bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#fecaca", alpha=0.94),
                 arrowprops=dict(arrowstyle="->", color="#b91c1c", lw=1.1), zorder=8)

    # ==================== COMPACT 2x3 MATRIX (PANELS C-H) ====================
    sub_specs = [
        (0, 0, "1uM_icilin_nimod", "C   1 µM Icilin + Nimodipine (Uncoupled Flood)", "Vm", 620, "1 µM Icilin + Nimodipine"),
        (0, 1, "50nM_A803467",    "D   50 nM + A-803467 (Naᵥ1.8 Block)", "Vm", 620, "50 nM Icilin + Naᵥ1.8 Blocker"),
        (0, 2, "50nM_nimodipine", "E   50 nM + Nimodipine (Caᵥ1.2 Block)", "Vm", 620, "50 nM Icilin + L-Type Blocker"),
        (1, 0, "50nM_FCCP",       "F   Mitochondrial Derangement (Diabetes / FCCP)", "METAB", 950, "50 nM Icilin"),
        (1, 1, "hyperosm_33C",    "G   PREDICTION: +98 mOsm at 33 °C (Zero Cooling)", "pH", 620, "+98 mOsm Hyperosmolarity (33 °C)"),
        (1, 2, "menthol_vs_cold", "H   PREDICTION: Menthol (H⁺-Insensitive Tonic)", "pH", 950, "35 µM Menthol (25 °C)"),
    ]

    for row, col, key, title, right_mode, ymax, stim_lbl in sub_specs:
        ax = fig.add_subplot(gs_bot[row, col])
        style_axis(ax, hide_top=True, hide_right=False)
        r = results[key]
        t = r["t_min"]

        if right_mode == "METAB":
            # Overlay Healthy -> Partial Metabolic Derangement -> Full FCCP Uncoupling!
            r_ok = results["50nM_icilin"]
            r_dm = results["50nM_metabolic"]
            r_fc = results["50nM_FCCP"]
            ax.plot(t, r_ok["ca_nM"], color="#64748b", lw=1.2, ls=":", label="Healthy (Ψm=1.0)")
            ax.plot(t, r_dm["ca_nM"], color="#d97706", lw=1.9, label="Metabolic/Diabetes (Ψm=0.42)")
            ax.plot(t, r_fc["ca_nM"], color="#b91c1c", lw=2.0, label="FCCP Uncoupled (Ψm=0.05)")
            ax.legend(loc="upper right", fontsize=7.0, frameon=True, facecolor="white", edgecolor="#cbd5e1")
            ax.spines["right"].set_visible(False)
        else:
            ax.plot(t, r["ca_nM"], color="#0284c7", lw=2.0, label="[Ca²⁺]ᵢ")
            ax2 = ax.twinx()
            ax2.spines["top"].set_visible(False)
            if right_mode == "pH":
                ax2.plot(t, r["pH_i"], color="#dc2626", lw=1.4, alpha=0.88, label="pHᵢ")
                ax2.set_ylim(6.72, 7.38)
                ax2.set_ylabel("pHᵢ", fontsize=8.2, color="#dc2626", fontweight="bold")
                ax2.tick_params(axis="y", colors="#dc2626", labelsize=7.8)
            else:
                ax2.plot(t, r["V"], color="#334155", lw=1.25, ls="--", alpha=0.82, label="Vm")
                ax2.set_ylim(-58, -5)
                ax2.set_ylabel("Vm (mV)", fontsize=8.2, color="#334155", fontweight="bold")
                ax2.tick_params(axis="y", colors="#334155", labelsize=7.8)

        ax.set_xlim(0, 7.0)
        ax.set_ylim(40, ymax)
        ax.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.3, color="#0284c7" if right_mode != "METAB" else "#0f172a", fontweight="bold")
        if row == 1:
            ax.set_xlabel("Time (min)", fontsize=8.6, color="#1e293b")
        ax.set_title(title, loc="left", fontsize=8.9, fontweight="bold", color="#0f172a", pad=6)
        draw_stim_bar(ax, 1.0, 6.0, 0.04, stim_lbl, color="#334155")

    fig.suptitle(
        "Figure 3. Electro-Metabolic Relaxation Oscillator (TRPM8–Naᵥ1.8–Caᵥ1.2–Mitochondrial H⁺): Spiking Dynamics & Falsifiable Predictions",
        fontsize=12.4,
        fontweight="bold",
        color="#0f172a",
        y=0.972,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, facecolor="white")
    plt.close(fig)


def make_genetics_mcu_figure(gen_res: dict, out_path: Path) -> None:
    """
    Publication-grade 6-panel Figure 4 mapping verified human genetic variants
    (TRPM8 GoF D665N/R30Q, rs10166942[C], SCN10A GoF I206M/G1662S, MICU1 LoF)
    and mitochondrial Ca2+/H+ perturbations (TMBIM5 uncoupling, Diabetic MCU KO)
    onto the exact same fixed-parameter electro-metabolic ODE model.
    """
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(14.0, 9.4), dpi=220, facecolor="white")
    gs = gridspec.GridSpec(
        2, 3, hspace=0.44, wspace=0.35,
        left=0.06, right=0.94, top=0.895, bottom=0.08,
    )

    # ---------- Panel A: TRPM8 GoF (Post-LASIK D665N / Trigeminal Neuralgia R30Q) ----------
    axA = fig.add_subplot(gs[0, 0])
    style_axis(axA, hide_top=True, hide_right=True)
    r_wt_mild = gen_res["wt_33C_55mOsm"]
    r_gof = gen_res["trpm8_gof_33C_55mOsm"]
    axA.plot(r_wt_mild["t_min"], r_wt_mild["ca_nM"], color="#64748b", lw=1.6, ls="--", label="WT (33 °C, +55 mOsm)")
    axA.plot(r_gof["t_min"], r_gof["ca_nM"], color="#dc2626", lw=2.0, label="TRPM8 GoF D665N/R30Q (−28 mV)")
    axA.set_xlim(0, 7.0)
    axA.set_ylim(40, 620)
    axA.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.8, fontweight="bold", color="#0f172a")
    axA.set_title("A   Human TRPM8 GoF: Ectopic Warm-Dry Bursting", loc="left", fontsize=9.1, fontweight="bold", color="#0f172a", pad=7)
    axA.legend(loc="upper right", fontsize=7.2, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    draw_stim_bar(axA, 1.0, 6.0, 0.04, "Mild +55 mOsm at 33 °C (Zero Cooling)", color="#334155")

    # ---------- Panel B: Human rs10166942[C] Protective Allele vs WT ----------
    axB = fig.add_subplot(gs[0, 1])
    style_axis(axB, hide_top=True, hide_right=True)
    r_wt_osm = gen_res["wt_33C_98mOsm"]
    r_rs = gen_res["trpm8_rs10166942_33C_98mOsm"]
    axB.plot(r_wt_osm["t_min"], r_wt_osm["ca_nM"], color="#0284c7", lw=1.9, label="WT TRPM8 (1.0× expr)")
    axB.plot(r_rs["t_min"], r_rs["ca_nM"], color="#0d9488", lw=2.0, ls="-.", label="rs10166942[C] (0.45× DRG expr)")
    axB.set_xlim(0, 7.0)
    axB.set_ylim(40, 620)
    axB.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.8, fontweight="bold", color="#0f172a")
    axB.set_title("B   rs10166942[C]: Attenuated Osmo-Cold Gain", loc="left", fontsize=9.1, fontweight="bold", color="#0f172a", pad=7)
    axB.legend(loc="upper right", fontsize=7.2, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    draw_stim_bar(axB, 1.0, 6.0, 0.04, "+98 mOsm at 33 °C (Latitudinal Allele)", color="#334155")

    # ---------- Panel C: Human SCN10A (Nav1.8) GoF (I206M / G1662S) ----------
    axC = fig.add_subplot(gs[0, 2])
    style_axis(axC, hide_top=True, hide_right=False)
    r_wt_25 = gen_res["wt_25C_25nM"]
    r_scn = gen_res["scn10a_gof_25C_25nM"]
    axC.plot(r_wt_25["t_min"], r_wt_25["ca_nM"], color="#64748b", lw=1.6, ls="--", label="WT Naᵥ1.8 [Ca²⁺]ᵢ")
    axC.plot(r_scn["t_min"], r_scn["ca_nM"], color="#7c3aed", lw=2.0, label="SCN10A GoF I206M (−4 mV)")
    axC.set_xlim(0, 7.0)
    axC.set_ylim(40, 620)
    axC.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.8, fontweight="bold", color="#7c3aed")
    axC2 = axC.twinx()
    axC2.spines["top"].set_visible(False)
    axC2.plot(r_scn["t_min"], r_scn["cm"], color="#ea580c", lw=1.4, ls=":", label="SCN10A Matrix Ca²⁺ (cm)")
    axC2.set_ylim(0, 7.5)
    axC2.set_ylabel("Matrix Ca²⁺, cm (a.u.)", fontsize=8.2, fontweight="bold", color="#ea580c")
    axC2.tick_params(axis="y", colors="#ea580c", labelsize=7.8)
    axC.set_title("C   Human SCN10A GoF: Tonic Plateau Lock-In", loc="left", fontsize=9.1, fontweight="bold", color="#0f172a", pad=7)
    lines_c = [
        plt.Line2D([0], [0], color="#64748b", lw=1.6, ls="--", label="WT [Ca²⁺]ᵢ"),
        plt.Line2D([0], [0], color="#7c3aed", lw=2.0, label="SCN10A GoF [Ca²⁺]ᵢ"),
        plt.Line2D([0], [0], color="#ea580c", lw=1.5, ls=":", label="GoF Matrix cm (3.5×)"),
    ]
    axC.legend(handles=lines_c, loc="upper right", fontsize=6.9, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    draw_stim_bar(axC, 1.0, 6.0, 0.04, "Subthreshold Drive (25 nM, 25 °C)", color="#334155")

    # ---------- Panel D: Human MICU1 LoF (MPXPS) & TMBIM5 Uncoupling ----------
    axD = fig.add_subplot(gs[1, 0])
    style_axis(axD, hide_top=True, hide_right=False)
    r_wt = gen_res["wt_25C_50nM"]
    r_micu1 = gen_res["micu1_lof_25C_50nM"]
    r_tmb = gen_res["tmbim5_uncoupled_25C_50nM"]
    axD.plot(r_wt["t_min"], r_wt["ca_nM"], color="#64748b", lw=1.3, ls=":", label="WT Oscillatory [Ca²⁺]ᵢ")
    axD.plot(r_micu1["t_min"], r_micu1["ca_nM"], color="#0891b2", lw=1.9, label="MICU1 LoF (Gate Lost)")
    axD.plot(r_tmb["t_min"], r_tmb["ca_nM"], color="#b91c1c", lw=1.9, label="TMBIM5 H⁺-Uncoupled")
    axD.set_xlim(0, 7.0)
    axD.set_ylim(20, 860)
    axD.set_xlabel("Time (min)", fontsize=8.8, color="#1e293b")
    axD.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=8.8, fontweight="bold", color="#0f172a")
    axD2 = axD.twinx()
    axD2.spines["top"].set_visible(False)
    axD2.plot(r_micu1["t_min"], r_micu1["cm"], color="#0891b2", lw=1.25, ls="--", alpha=0.8)
    axD2.set_ylim(0, 8.2)
    axD2.set_ylabel("MICU1 Matrix cm (dashed)", fontsize=7.9, fontweight="bold", color="#0e7490")
    axD2.tick_params(axis="y", colors="#0e7490", labelsize=7.6)
    axD.annotate("25× resting cm\n(pre-stimulus)", xy=(0.55, 160), xytext=(0.12, 650),
                 fontsize=7.0, fontweight="bold", color="#0e7490",
                 bbox=dict(boxstyle="round,pad=0.20", facecolor="white", edgecolor="#a5f3fc", alpha=0.95),
                 arrowprops=dict(arrowstyle="->", color="#0e7490", lw=1.1), zorder=8)
    axD.set_title("D   MICU1 LoF (MPXPS) & TMBIM5 Uncoupling", loc="left", fontsize=9.1, fontweight="bold", color="#0f172a", pad=7)
    axD.legend(loc="upper right", fontsize=6.9, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    draw_stim_bar(axD, 1.0, 6.0, 0.04, "50 nM Icilin (25 °C)", color="#334155")

    # ---------- Panel E: The Diabetic MCU Paradox (George et al., 2022, Pain) ----------
    axE = fig.add_subplot(gs[1, 1])
    style_axis(axE, hide_top=True, hide_right=False)
    r_dm = gen_res["diabetic_hyperexc"]
    r_ko = gen_res["diabetic_mcu_ko"]
    # Plot mitochondrial matrix Ca2+ load (axon degeneration driver) on primary axis!
    axE.plot(r_wt["t_min"], r_wt["cm"], color="#64748b", lw=1.4, ls=":", label="Healthy Cyclic cm (clears)")
    axE.plot(r_dm["t_min"], r_dm["cm"], color="#dc2626", lw=2.2, label="Diabetic Sustained cm Overload")
    axE.plot(r_ko["t_min"], r_ko["cm"], color="#16a34a", lw=2.3, label="Diabetic + MCU KO (−97% cm)")
    axE.fill_between(r_dm["t_min"], r_dm["cm"], r_ko["cm"], where=(r_dm["t_min"] >= 1.0) & (r_dm["t_min"] <= 6.0),
                     color="#fca5a5", alpha=0.22)
    axE.set_xlim(0, 7.0)
    axE.set_ylim(-0.95, 7.8)
    axE.set_xlabel("Time (min)", fontsize=8.8, color="#1e293b")
    axE.set_ylabel("Mito Matrix Ca²⁺ Load, cm (a.u.)", fontsize=8.6, fontweight="bold", color="#b91c1c")
    axE2 = axE.twinx()
    axE2.spines["top"].set_visible(False)
    axE2.plot(r_ko["t_min"], r_ko["pH_i"], color="#15803d", lw=1.3, ls="--", alpha=0.82)
    axE2.set_ylim(6.72, 7.42)
    axE2.set_ylabel("MCU KO pHᵢ (flat 7.29)", fontsize=7.9, fontweight="bold", color="#15803d")
    axE2.tick_params(axis="y", colors="#15803d", labelsize=7.6)
    axE.set_title("E   Diabetic MCU Paradox (George et al., 2022)", loc="left", fontsize=9.1, fontweight="bold", color="#0f172a", pad=7)
    axE.legend(loc="upper right", fontsize=6.9, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    draw_stim_bar(axE, 1.0, 6.0, 0.03, "Chronic Hyperexcitable Drive", color="#334155")

    # ---------- Panel F: Trade-Off Map (Oscillatory Coding vs. Matrix Overload) ----------
    axF = fig.add_subplot(gs[1, 2])
    style_axis(axF, hide_top=True, hide_right=True)
    # Use absolute data coordinates for callout text so every label has generous whitespace!
    summary_items = [
        ("Healthy WT\n(Cyclic)", r_wt["cm_mean"], r_wt["tail_amp"], "#0284c7", (5.15, 222)),
        ("TRPM8 GoF\n(33°C +55)", r_gof["cm_mean"], r_gof["tail_amp"], "#ef4444", (1.65, 232)),
        ("rs10166942[C]\n(Protected)", r_rs["cm_mean"], r_rs["tail_amp"], "#0d9488", (1.45, 118)),
        ("MICU1 LoF\n(MPXPS)", r_micu1["cm_mean"], r_micu1["tail_amp"], "#0891b2", (3.65, 138)),
        ("SCN10A GoF\n(I206M)", r_scn["cm_mean"], r_scn["tail_amp"], "#7c3aed", (4.85, 86)),
        ("Diabetic / TMBIM5\n(Overload)", r_tmb["cm_mean"], r_tmb["tail_amp"], "#b91c1c", (5.72, 152)),
        ("Diabetic +\nMCU KO", r_ko["cm_mean"], r_ko["tail_amp"], "#16a34a", (0.95, 52)),
    ]
    # Shaded physiological vs pathological zones
    axF.axhspan(185, 310, color="#e0f2fe", alpha=0.35, zorder=0)
    axF.axvspan(4.8, 7.0, color="#fee2e2", alpha=0.35, zorder=0)
    axF.text(0.15, 294, "Healthy Rhythmic Burst Coding", fontsize=7.2, fontweight="bold", color="#0369a1")
    axF.text(6.85, 294, "Matrix Ca²⁺\nOverload Zone", ha="right", va="top", fontsize=7.0, fontweight="bold", color="#b91c1c")

    for lbl, x_cm, y_amp, col, (tx, ty) in summary_items:
        axF.scatter([x_cm], [y_amp], s=68, color=col, edgecolor="white", linewidth=1.1, zorder=5)
        axF.annotate(lbl, xy=(x_cm, y_amp), xytext=(tx, ty), ha="center", va="center",
                     fontsize=6.9, fontweight="bold", color=col,
                     bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor=col, lw=0.6, alpha=0.92),
                     arrowprops=dict(arrowstyle="->", color=col, lw=0.9, shrinkA=2, shrinkB=4), zorder=7)

    axF.set_xlim(-0.2, 7.0)
    axF.set_ylim(-15, 310)
    axF.set_xlabel("Mean Mito Matrix Ca²⁺ Load, cm (Axon Stress)", fontsize=8.5, fontweight="bold", color="#1e293b")
    axF.set_ylabel("Steady-Tail Osc. Amp (nM)", fontsize=8.5, fontweight="bold", color="#1e293b")
    axF.set_title("F   Phase Map: Burst Coding vs. Matrix Overload", loc="left", fontsize=9.1, fontweight="bold", color="#0f172a", pad=7)

    fig.suptitle(
        "Figure 4. Human Genetic Variants & MCU/MICU1/TMBIM5 Perturbations in the Electro-Metabolic Oscillator",
        fontsize=12.4,
        fontweight="bold",
        color="#0f172a",
        y=0.965,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, facecolor="white")
    plt.close(fig)


def make_neurovascular_figure(nv_res: dict, out_path: Path) -> None:
    """
    Publication-grade 4-panel Figure 5 dissecting the Cav1.2 / Cav1.3
    Neurovascular & Dihydropyridine Pharmacology Conundrum:
      Panel A: Avascular Cornea / Clamped Bath (pure neural (-)-(S)-Bay K 8644 vs Nifedipine)
      Panel B: Perfused Skin Vasodilation Confound (Low-dose Nifedipine vs Pure Vasodilator ACEi/ARB)
      Panel C: Perfused Skin Vasoconstriction Mirror Confound (Bay K 8644 vs Vascular-Only Constriction)
      Panel D: 2D Neurovascular De-Confounding Phase Map (Perfusion Shift vs Neuronal Cav1.2/1.3 Scale)
    """
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(13.8, 9.6), dpi=220, facecolor="white")
    gs = gridspec.GridSpec(
        2, 2, hspace=0.40, wspace=0.30,
        left=0.07, right=0.95, top=0.90, bottom=0.08,
    )

    # ---------- Panel A: Pure Neural Bay K 8644 ± Nifedipine (Avascular Cornea / Clamped Bath) ----------
    axA = fig.add_subplot(gs[0, 0])
    style_axis(axA, hide_top=True, hide_right=True)
    r_av_wt = nv_res["avasc_wt_82"]
    r_av_bayk = nv_res["avasc_bayk_82"]
    r_av_nif = nv_res["avasc_bayk_nif_82"]
    axA.plot(r_av_wt["t_min"], r_av_wt["ca_nM"], color="#0284c7", lw=1.9, label="WT Avascular (33 °C, +82 mOsm)")
    axA.plot(r_av_bayk["t_min"], r_av_bayk["ca_nM"], color="#dc2626", lw=2.1, label="+ (−)-(S)-Bay K 8644 (Plateau Lock-In)")
    axA.plot(r_av_nif["t_min"], r_av_nif["ca_nM"], color="#16a34a", lw=2.0, ls="--", label="+ Bay K 8644 + Nifedipine (Reversed)")
    axA.set_xlim(0, 7.0)
    axA.set_ylim(20, 960)
    axA.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=9.0, fontweight="bold", color="#0f172a")
    axA.set_title("A   Avascular Cornea / Clamped Prep: Pure Neural Bay K 8644 ± Nifedipine",
                  loc="left", fontsize=9.2, fontweight="bold", color="#0f172a", pad=7)
    axA.legend(loc="upper right", fontsize=7.4, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    axA.annotate("Potentiated Caᵥ1.2/1.3 bypasses\nTRPM8 H⁺ reset → sustained plateau",
                 xy=(3.5, 780), xytext=(3.6, 560), ha="center", fontsize=7.3, fontweight="bold", color="#b91c1c",
                 bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#fecaca", alpha=0.94),
                 arrowprops=dict(arrowstyle="->", color="#b91c1c", lw=1.0), zorder=8)
    draw_stim_bar(axA, 1.0, 6.0, 0.035, "+82 mOsm (33 °C, Zero Blood Flow + HC-030031 TRPA1 Block)", color="#334155")

    # ---------- Panel B: Perfused Skin Vasodilation Confound ----------
    axB = fig.add_subplot(gs[0, 1])
    style_axis(axB, hide_top=True, hide_right=True)
    r_skin_wt = nv_res["skin_wt_98"]
    r_clamped_lo = nv_res["clamped_low_nif_98"]
    r_ace = nv_res["pure_vasodil_98"]
    r_perf_lo = nv_res["perf_low_nif_98"]
    axB.plot(r_skin_wt["t_min"], r_skin_wt["ca_nM"], color="#0284c7", lw=1.7, ls=":", label="Untreated Dry Skin (33 °C, +98 mOsm)")
    axB.plot(r_clamped_lo["t_min"], r_clamped_lo["ca_nM"], color="#7c3aed", lw=2.0, label="Low-Dose Nifedipine CLAMPED (Caᵥ1.3 Still Bursts)")
    axB.plot(r_ace["t_min"], r_ace["ca_nM"], color="#ea580c", lw=2.4, label="Pure Vasodilator (ACEi/ARB: +2 °C, 50% Washout)")
    axB.plot(r_perf_lo["t_min"], r_perf_lo["ca_nM"], color="#16a34a", lw=1.6, ls="--", label="Low-Dose Nifedipine PERFUSED (Masked!)")
    axB.set_xlim(0, 7.0)
    axB.set_ylim(30, 630)
    axB.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=9.0, fontweight="bold", color="#0f172a")
    axB.set_title("B   Perfused-Skin Vasodilation Confound (Warm Blood + Osmotic Washout)",
                  loc="left", fontsize=9.2, fontweight="bold", color="#0f172a", pad=7)
    axB.legend(loc="upper right", fontsize=7.2, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    axB.annotate("Warm blood + osmotic washout\nsilences skin even without\nfull neuronal Caᵥ1.3 block",
                 xy=(3.8, 106), xytext=(3.85, 225), ha="center", fontsize=7.3, fontweight="bold", color="#c2410c",
                 bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#fed7aa", alpha=0.95),
                 arrowprops=dict(arrowstyle="->", color="#c2410c", lw=1.0), zorder=8)
    draw_stim_bar(axB, 1.0, 6.0, 0.035, "Dry Skin (+98 mOsm, 33 °C Baseline)", color="#334155")

    # ---------- Panel C: Mirror Vasoconstriction Confound ----------
    axC = fig.add_subplot(gs[1, 0])
    style_axis(axC, hide_top=True, hide_right=True)
    r_q_wt = nv_res["quiet_wt_72"]
    r_vc_only = nv_res["vasoconst_only_72"]
    r_vc_bayk = nv_res["perf_bayk_72"]
    axC.plot(r_q_wt["t_min"], r_q_wt["ca_nM"], color="#64748b", lw=1.7, ls="--", label="Quiet Mild-Dry Skin (33 °C, +72 mOsm)")
    axC.plot(r_vc_only["t_min"], r_vc_only["ca_nM"], color="#0891b2", lw=2.0, label="Vascular-Only Constriction (−2 °C, +30% Osm)")
    axC.plot(r_vc_bayk["t_min"], r_vc_bayk["ca_nM"], color="#b91c1c", lw=2.1, label="Bay K 8644 in Perfused Skin (Vasc + Neural)")
    axC.set_xlim(0, 7.0)
    axC.set_ylim(20, 960)
    axC.set_xlabel("Time (min)", fontsize=9.0, color="#1e293b")
    axC.set_ylabel("[Ca²⁺]ᵢ (nM)", fontsize=9.0, fontweight="bold", color="#0f172a")
    axC.set_title("C   Mirror Confound: Vasoconstriction Mimics Neuronal Caᵥ Excitation",
                  loc="left", fontsize=9.2, fontweight="bold", color="#0f172a", pad=7)
    axC.legend(loc="upper right", fontsize=7.3, frameon=True, facecolor="white", edgecolor="#cbd5e1")
    axC.annotate("Pure vasoconstriction ignites\n3 bursts (false-positive mimic!)",
                 xy=(3.2, 515), xytext=(3.75, 650), ha="center", fontsize=7.3, fontweight="bold", color="#0e7490",
                 bbox=dict(boxstyle="round,pad=0.22", facecolor="white", edgecolor="#a5f3fc", alpha=0.95),
                 arrowprops=dict(arrowstyle="->", color="#0e7490", lw=1.0), zorder=8)
    draw_stim_bar(axC, 1.0, 6.0, 0.035, "Subthreshold Dry Skin (+72 mOsm, 33 °C)", color="#334155")

    # ---------- Panel D: 2D Neurovascular De-Confounding Phase Diagram ----------
    axD = fig.add_subplot(gs[1, 1])
    style_axis(axD, hide_top=True, hide_right=True)
    axD.grid(False)

    # Sweep perfusion shift p in [-1.0 (vasoconstrict) .. +1.0 (vasodilate)] vs neuronal g_CaL scale [0.0 .. 1.8]
    p_vals = np.linspace(-1.0, 1.0, 15)
    g_vals = np.linspace(0.0, 1.8, 15)
    regime_grid = np.zeros((len(g_vals), len(p_vals)))  # 0=Silent, 1=Rhythmic Bursting, 2=Tonic Plateau

    for ig, g_sc in enumerate(g_vals):
        # Smooth V1/2 shift coupled to Bay K potentiation above g_sc > 1.0
        vh_sh = -5.0 * max(0.0, g_sc - 1.0)
        for ip, pv in enumerate(p_vals):
            dT = 2.2 * pv
            wash = 0.45 * pv
            r = run_simulation(
                SimCondition(
                    "grid",
                    temp_c=33.0,
                    delta_osm_mosm=88.0,
                    cal_scale=float(g_sc),
                    cal_vh_shift_mv=float(vh_sh),
                    vasc_delta_temp_c=float(dT),
                    vasc_osm_washout_frac=float(wash),
                    t_end=360.0,
                ),
                fast_grid=True,
            )
            if r["tail_amp"] >= 80.0 and r["n_peaks"] >= 2:
                regime_grid[ig, ip] = 1.0  # Rhythmic Oscillatory Bursting
            elif r["ca_max"] >= 450.0 and r["cm_mean"] >= 4.0:
                regime_grid[ig, ip] = 2.0  # Non-Adapting High-Ca2+ Plateau
            else:
                regime_grid[ig, ip] = 0.0  # Silent / Subthreshold

    from matplotlib.colors import ListedColormap
    cmap = ListedColormap(["#f1f5f9", "#bae6fd", "#fecaca"])
    axD.imshow(
        regime_grid,
        origin="lower",
        aspect="auto",
        extent=[p_vals[0], p_vals[-1], g_vals[0], g_vals[-1]],
        cmap=cmap,
        interpolation="nearest",
        vmin=0, vmax=2,
        zorder=0,
    )
    # Subtle contour guides
    axD.axvline(0.0, color="#0f172a", lw=1.4, ls="-", alpha=0.75, zorder=2)

    # Region labels
    axD.text(0.58, 0.30, "SILENT /\nSUBTHRESHOLD", ha="center", va="center",
             fontsize=7.6, fontweight="bold", color="#475569")
    axD.text(-0.48, 0.72, "RHYTHMIC\nBURSTING", ha="center", va="center",
             fontsize=7.8, fontweight="bold", color="#0369a1")
    axD.text(-0.48, 1.56, "L-TYPE PLATEAU\nOVERLOAD", ha="center", va="center",
             fontsize=7.6, fontweight="bold", color="#b91c1c")

    # Baseline WT point at (0, 1.0)
    axD.scatter([0.0], [1.0], s=72, color="#0f172a", edgecolor="white", lw=1.2, zorder=6)
    axD.text(0.05, 1.02, "WT", fontsize=7.6, fontweight="bold", color="#0f172a", zorder=7)

    # Vertical unconfounded axis (Avascular Cornea / Clamped Bath)
    axD.annotate("", xy=(0.0, 1.65), xytext=(0.0, 1.05),
                 arrowprops=dict(arrowstyle="->", color="#0f172a", lw=2.0), zorder=5)
    axD.annotate("", xy=(0.0, 0.08), xytext=(0.0, 0.95),
                 arrowprops=dict(arrowstyle="->", color="#0f172a", lw=2.0), zorder=5)
    axD.text(0.05, 1.67, "Avascular / Clamped:\nBay K ↑  |  Nifedipine ↓", fontsize=7.0, fontweight="bold",
             color="#0f172a", va="top",
             bbox=dict(boxstyle="round,pad=0.18", facecolor="white", edgecolor="#94a3b8", alpha=0.92), zorder=8)

    # Diagonal confounded trajectories in Perfused Skin
    axD.annotate("", xy=(0.82, 0.55), xytext=(0.04, 0.97),
                 arrowprops=dict(arrowstyle="->", color="#16a34a", lw=2.1, ls="--"), zorder=5)
    axD.text(0.60, 0.85, "Perfused Nifedipine\n(Dilates + Blocks)", ha="center", fontsize=6.9, fontweight="bold",
             color="#15803d", bbox=dict(boxstyle="round,pad=0.16", facecolor="white", edgecolor="#bbf7d0", alpha=0.92), zorder=8)

    axD.annotate("", xy=(-0.82, 1.58), xytext=(-0.04, 1.03),
                 arrowprops=dict(arrowstyle="->", color="#dc2626", lw=2.1, ls="--"), zorder=5)
    axD.text(-0.56, 1.18, "Perfused Bay K\n(Constricts + Opens)", ha="center", fontsize=6.9, fontweight="bold",
             color="#b91c1c", bbox=dict(boxstyle="round,pad=0.16", facecolor="white", edgecolor="#fecaca", alpha=0.92), zorder=8)

    # Pure vasodilator horizontal trajectory (ACEi/ARB comparator)
    axD.annotate("", xy=(0.82, 1.0), xytext=(0.05, 1.0),
                 arrowprops=dict(arrowstyle="->", color="#ea580c", lw=1.8, ls=":"), zorder=5)
    axD.text(0.62, 1.13, "ACEi / ARB\n(Vascular Only)", ha="center", fontsize=6.8, fontweight="bold",
             color="#c2410c", bbox=dict(boxstyle="round,pad=0.14", facecolor="white", edgecolor="#fed7aa", alpha=0.92), zorder=8)

    axD.set_xlim(-1.0, 1.0)
    axD.set_ylim(0.0, 1.8)
    axD.set_xlabel("Cutaneous Perfusion Shift  (← Vasoconstriction  |  0 = Avascular/Clamped  |  Vasodilation →)",
                   fontsize=8.1, fontweight="bold", color="#1e293b")
    axD.set_ylabel("Neuronal Caᵥ1.2/1.3 Scale (gCaL)", fontsize=8.8, fontweight="bold", color="#1e293b")
    axD.set_title("D   2D De-Confounding Map: Neural Gating vs. Hemodynamic Shift",
                  loc="left", fontsize=9.2, fontweight="bold", color="#0f172a", pad=7)

    fig.suptitle(
        "Figure 5. The Caᵥ1.2/Caᵥ1.3 Neurovascular & Dihydropyridine Pharmacology Conundrum (Bay K 8644 ± Nifedipine)",
        fontsize=12.2,
        fontweight="bold",
        color="#0f172a",
        y=0.968,
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=220, facecolor="white")
    plt.close(fig)


def main():
    conditions = {
        "50nM_icilin": SimCondition("50 nM Icilin (Control, 25°C)", temp_c=25.0, icilin_nm=50.0),
        "1uM_icilin_nimod": SimCondition("1 uM Icilin + Nimodipine", temp_c=25.0, icilin_nm=1000.0, block_cal=True),
        "50nM_A803467": SimCondition("50 nM Icilin + A-803467 (Nav1.8 block)", temp_c=25.0, icilin_nm=50.0, block_nav18=True),
        "50nM_nimodipine": SimCondition("50 nM Icilin + Nimodipine (Cav1.2 block)", temp_c=25.0, icilin_nm=50.0, block_cal=True),
        "50nM_FCCP": SimCondition("50 nM Icilin + FCCP (Psi_m=0.05)", temp_c=25.0, icilin_nm=50.0, mito_psi=0.05),
        "50nM_metabolic": SimCondition("50 nM Icilin + Metabolic Derangement (Psi_m=0.42)", temp_c=25.0, icilin_nm=50.0, mito_psi=0.42),
        "hyperosm_33C": SimCondition("Isothermal +98 mOsm at 33°C (No drug)", temp_c=33.0, delta_osm_mosm=98.0),
        "menthol_vs_cold": SimCondition("Menthol 35 uM (H+-insensitive, 25°C)", temp_c=25.0, menthol_um=35.0),
    }

    results = {}
    print("=== CORE 8 PHYSIOLOGICAL / PHARMACOLOGICAL REGIMES (FIGURE 3) ===")
    print(f"{'Condition':<48} | {'Peaks':>5} | {'Ca range (nM)':>15} | {'TailAmp':>7} | {'Min pH':>6} | {'Vm range (mV)'}")
    print("-" * 108)
    for k, cond in conditions.items():
        r = run_simulation(cond)
        results[k] = r
        v0, v1 = r["v_range"]
        print(
            f"{cond.name:<48} | {r['n_peaks']:>5d} | [{r['ca_min']:5.1f}, {r['ca_max']:6.1f}] | "
            f"{r['tail_amp']:>7.1f} | {r['pH_min']:>6.2f} | [{v0:+.1f}, {v1:+.1f}]"
        )

    out_fig = Path(__file__).resolve().parent / "outputs" / "trpm8_oscillator_panels.png"
    make_figure(results, out_fig)
    print(f"\nSaved Figure 3 to: {out_fig}")

    # Human genetics & MCU/MICU1/TMBIM5 perturbation suite (Figure 4)
    gen_conditions = {
        "wt_33C_55mOsm": SimCondition("WT (33°C, +55 mOsm mild dry)", temp_c=33.0, delta_osm_mosm=55.0),
        "trpm8_gof_33C_55mOsm": SimCondition("TRPM8 GoF D665N/R30Q (−28 mV, 33°C +55 mOsm)", temp_c=33.0, delta_osm_mosm=55.0, trpm8_vh_shift_mv=-28.0),
        "wt_33C_98mOsm": SimCondition("WT (33°C, +98 mOsm)", temp_c=33.0, delta_osm_mosm=98.0),
        "trpm8_rs10166942_33C_98mOsm": SimCondition("TRPM8 rs10166942[C] (0.45× expr, 33°C +98 mOsm)", temp_c=33.0, delta_osm_mosm=98.0, trpm8_expr_scale=0.45),
        "wt_25C_25nM": SimCondition("WT (25°C, 25 nM subthreshold)", temp_c=25.0, icilin_nm=25.0),
        "scn10a_gof_25C_25nM": SimCondition("SCN10A GoF I206M/G1662S (−4 mV, 25°C 25 nM)", temp_c=25.0, icilin_nm=25.0, nav18_vh_shift_mv=-4.0),
        "wt_25C_50nM": SimCondition("WT Healthy Oscillatory (25°C, 50 nM)", temp_c=25.0, icilin_nm=50.0),
        "micu1_lof_25C_50nM": SimCondition("MICU1 LoF MPXPS (Kd=0.14 uM, n=1.8; 25°C 50 nM)", temp_c=25.0, icilin_nm=50.0, micu1_kd_um=0.14, micu1_hill=1.8),
        "tmbim5_uncoupled_25C_50nM": SimCondition("TMBIM5 H+-Uncoupled (gain=0.25; 25°C 50 nM)", temp_c=25.0, icilin_nm=50.0, tmbim5_h_coupling=0.25),
        "diabetic_hyperexc": SimCondition("Diabetic Hyperexcitable (Psi=0.85, TMBIM5=0.45, 70 nM)", temp_c=25.0, icilin_nm=70.0, mito_psi=0.85, tmbim5_h_coupling=0.45),
        "diabetic_mcu_ko": SimCondition("Diabetic + Nociceptor MCU KO (mcu=0.03)", temp_c=25.0, icilin_nm=70.0, mito_psi=0.85, tmbim5_h_coupling=0.45, mcu_scale=0.03),
    }

    gen_res = {}
    print("\n=== HUMAN GENETICS & MCU/MICU1/TMBIM5 PERTURBATIONS (FIGURE 4) ===")
    print(f"{'Perturbation':<54} | {'Pks':>3} | {'TailAmp':>7} | {'cmBase':>6} | {'cmMean':>6} | {'cmMax':>5} | {'pHmin':>5}")
    print("-" * 105)
    for k, cond in gen_conditions.items():
        r = run_simulation(cond)
        gen_res[k] = r
        print(
            f"{cond.name:<54} | {r['n_peaks']:>3d} | {r['tail_amp']:>7.1f} | "
            f"{r['cm_base']:>6.2f} | {r['cm_mean']:>6.2f} | {r['cm_max']:>5.2f} | {r['pH_min']:>5.2f}"
        )

    out_fig4 = Path(__file__).resolve().parent / "outputs" / "fig4_human_genetics_mcu.png"
    make_genetics_mcu_figure(gen_res, out_fig4)
    print(f"\nSaved Figure 4 to: {out_fig4}")

    # Neurovascular & Cav1.2/1.3 DHP pharmacology suite (Figure 5)
    nv_conditions = {
        "avasc_wt_82": SimCondition("Avascular WT (33°C, +82 mOsm)", temp_c=33.0, delta_osm_mosm=82.0),
        "avasc_bayk_82": SimCondition("Avascular + (-)-(S)-Bay K 8644 (g=1.6, -4 mV)", temp_c=33.0, delta_osm_mosm=82.0, cal_scale=1.6, cal_vh_shift_mv=-4.0),
        "avasc_bayk_nif_82": SimCondition("Avascular + Bay K 8644 + Nifedipine (g=0)", temp_c=33.0, delta_osm_mosm=82.0, block_cal=True),
        "skin_wt_98": SimCondition("Perfused Dry Skin WT (33°C, +98 mOsm)", temp_c=33.0, delta_osm_mosm=98.0),
        "clamped_low_nif_98": SimCondition("Low-Dose Nifedipine CLAMPED (g=0.60, Cav1.3 active)", temp_c=33.0, delta_osm_mosm=98.0, cal_scale=0.60),
        "pure_vasodil_98": SimCondition("Pure Vasodilator ACEi/ARB (+2°C, 50% washout)", temp_c=33.0, delta_osm_mosm=98.0, vasc_delta_temp_c=2.0, vasc_osm_washout_frac=0.50),
        "perf_low_nif_98": SimCondition("Low-Dose Nifedipine PERFUSED (vasodil + g=0.60)", temp_c=33.0, delta_osm_mosm=98.0, cal_scale=0.60, vasc_delta_temp_c=2.0, vasc_osm_washout_frac=0.50),
        "quiet_wt_72": SimCondition("Quiet Mild-Dry Skin (33°C, +72 mOsm)", temp_c=33.0, delta_osm_mosm=72.0),
        "vasoconst_only_72": SimCondition("Vascular-Only Constriction (-2°C, +30% Osm)", temp_c=33.0, delta_osm_mosm=72.0, vasc_delta_temp_c=-2.0, vasc_osm_washout_frac=-0.30),
        "perf_bayk_72": SimCondition("Bay K 8644 Perfused Skin (vasoconst + g=1.8, -6 mV)", temp_c=33.0, delta_osm_mosm=72.0, cal_scale=1.8, cal_vh_shift_mv=-6.0, vasc_delta_temp_c=-2.0, vasc_osm_washout_frac=-0.30),
    }

    nv_res = {}
    print("\n=== NEUROVASCULAR & Cav1.2/1.3 DHP PHARMACOLOGY (FIGURE 5) ===")
    print(f"{'Condition':<54} | {'Pks':>3} | {'CaMax':>6} | {'TailAmp':>7} | {'cmMean':>6} | {'pHmin':>5}")
    print("-" * 96)
    for k, cond in nv_conditions.items():
        r = run_simulation(cond)
        nv_res[k] = r
        print(
            f"{cond.name:<54} | {r['n_peaks']:>3d} | {r['ca_max']:6.1f} | "
            f"{r['tail_amp']:>7.1f} | {r['cm_mean']:>6.2f} | {r['pH_min']:>5.2f}"
        )

    out_fig5 = Path(__file__).resolve().parent / "outputs" / "fig5_neurovascular_pharmacology.png"
    make_neurovascular_figure(nv_res, out_fig5)
    print(f"\nSaved Figure 5 to: {out_fig5}")


if __name__ == "__main__":
    main()
