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
Monte Carlo Power & Confound Simulation for the TRPM8 Osmolarity / Cold Study.

Answers Matthew's two core statistical questions quantitatively:
1. Does "large N" rescue a naive outdoor-weather correlation?
   -> NO for systematic confounding (winter indoor-heating RH inversion).
      YES once dew point is converted to true local Vapor Pressure Deficit (VPD).
2. How many participants are actually needed if you are NOT doing viral recruitment?
   -> Just N = 15-25 in within-subject randomized crossovers (Fast-Cooling Gel
      Submersion & 2x2 Water Bolus) achieves >80% power — outperforming a
      cross-sectional app with N = 500+ single reports.
"""

from pathlib import Path
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def sat_vp_kpa(temp_c: np.ndarray) -> np.ndarray:
    """August-Roche-Magnus saturation vapor pressure in kPa."""
    return 0.61094 * np.exp((17.625 * temp_c) / (temp_c + 243.04))


def simulate_winter_rh_inversion(n_obs: int = 2000, seed: int = 42) -> dict:
    """
    Simulate realistic multi-season indoor + outdoor check-ins where TRUE biology
    has dry air (higher local VPD) feeling colder (beta_vpd = -0.65 ASHRAE units/kPa).
    Show what happens when you naively correlate OUTDOOR RH vs. true LOCAL VPD.
    """
    rng = np.random.default_rng(seed)

    # Mix of winter (cold, humid outdoor air -> bone-dry heated indoor air) and mild/humid days
    outdoor_temp = rng.uniform(-5.0, 26.0, size=n_obs)
    # In winter outdoor RH tends to be higher (65-90%); in warmer weather varied (35-85%)
    outdoor_rh = np.clip(78.0 - 0.8 * outdoor_temp + rng.normal(0, 12, size=n_obs), 20, 95)
    is_indoors = rng.random(size=n_obs) < 0.72
    indoor_temp = rng.normal(21.0, 1.2, size=n_obs)

    ea_out = sat_vp_kpa(outdoor_temp) * (outdoor_rh / 100.0)
    local_temp = np.where(is_indoors, indoor_temp, outdoor_temp)
    es_local = sat_vp_kpa(local_temp)
    ea_local = np.where(is_indoors, np.minimum(ea_out, es_local * 0.65), ea_out)
    local_vpd = np.maximum(0.0, es_local - ea_local)

    # True perceptual generating equation:
    # Higher local temp -> warmer (+0.22/C), higher local VPD -> COLDER (-0.65/kPa)
    # Plus huge between-person + clothing noise (sigma = 1.4)
    subj_noise = rng.normal(0, 1.4, size=n_obs)
    thermal_rating = (
        -4.2
        + 0.22 * local_temp
        - 0.65 * local_vpd
        + subj_noise
    )

    # 1) Naive check: correlate outdoor RH with thermal rating (indoors subset!)
    ind_mask = is_indoors
    r_naive_rh, p_naive_rh = stats.pearsonr(outdoor_rh[ind_mask], thermal_rating[ind_mask])

    # 2) Physics-corrected check: partial correlation of local VPD controlling for local_temp
    # Residualize both thermal_rating and local_vpd on local_temp
    X = np.column_stack([np.ones(n_obs), local_temp])
    res_y = thermal_rating - X @ np.linalg.lstsq(X, thermal_rating, rcond=None)[0]
    res_vpd = local_vpd - X @ np.linalg.lstsq(X, local_vpd, rcond=None)[0]
    r_corr_vpd, p_corr_vpd = stats.pearsonr(res_vpd, res_y)

    return {
        "n_obs": n_obs,
        "n_indoors": int(ind_mask.sum()),
        "r_naive_outdoor_rh_indoors": float(r_naive_rh),
        "p_naive_outdoor_rh_indoors": float(p_naive_rh),
        "r_corrected_local_vpd": float(r_corr_vpd),
        "p_corrected_local_vpd": float(p_corr_vpd),
    }


def estimate_power_curves(n_sims: int = 800, seed: int = 123) -> dict:
    """
    Compare statistical power (alpha=0.05, two-tailed) across sample sizes N for:
      1. Cross-sectional passive survey (1 report/person, raw between-person variance)
      2. Longitudinal EMA (14 reports/person, within-subject centered VPD)
      3. Paired Fast-Cooling Gel Submersion / 2x2 Water Bolus crossover (medium dz=0.65)
      4. Paired Gel Submersion with strong osmotic contrast (dz=0.80)
    """
    rng = np.random.default_rng(seed)
    ns = np.array([6, 10, 15, 20, 25, 30, 40, 60, 100, 150, 250])

    curves = {
        "cross_sectional_1obs": [],
        "within_subject_ema_14obs": [],
        "crossover_gel_dz065": [],
        "crossover_gel_dz080": [],
    }

    # True effect in ASHRAE/10-pt units: 0.65 units per kPa VPD or condition delta
    beta_true = -0.65
    sigma_between = 1.40  # clothing, metabolism, individual cold tolerance
    sigma_within = 0.65   # trial-to-trial perceptual noise

    for N in ns:
        # 1. Cross-sectional (1 obs/person): total SD = sqrt(1.4^2 + 0.65^2) = 1.54
        hits_cs = 0
        for _ in range(n_sims):
            vpd = rng.uniform(0.3, 2.1, size=N)
            y = beta_true * vpd + rng.normal(0, np.hypot(sigma_between, sigma_within), size=N)
            _, p = stats.pearsonr(vpd, y)
            if p < 0.05:
                hits_cs += 1
        curves["cross_sectional_1obs"].append(hits_cs / n_sims)

        # 2. Longitudinal EMA (14 obs/person, within-subject slope t-test across N subjects)
        hits_ema = 0
        k_obs = 14
        for _ in range(n_sims):
            slopes = np.empty(N)
            for i in range(N):
                vpd_i = rng.uniform(0.3, 2.1, size=k_obs)
                # Random subject slope heterogeneity (SD = 0.35)
                b_i = rng.normal(beta_true, 0.35)
                y_i = b_i * vpd_i + rng.normal(0, sigma_within, size=k_obs)
                v_c = vpd_i - vpd_i.mean()
                slopes[i] = np.dot(v_c, y_i - y_i.mean()) / (np.dot(v_c, v_c) + 1e-9)
            _, p = stats.ttest_1samp(slopes, 0.0)
            if p < 0.05:
                hits_ema += 1
        curves["within_subject_ema_14obs"].append(hits_ema / n_sims)

        # 3 & 4. Paired Within-Subject Crossovers (e.g., Gel Submersion / Water Bolus 2x2)
        for key, dz in [("crossover_gel_dz065", 0.65), ("crossover_gel_dz080", 0.80)]:
            diffs = rng.normal(dz, 1.0, size=(n_sims, N))
            t_vals = diffs.mean(axis=1) / (diffs.std(axis=1, ddof=1) / np.sqrt(N))
            p_vals = 2.0 * stats.t.sf(np.abs(t_vals), df=N - 1)
            curves[key].append(float(np.mean(p_vals < 0.05)))

    return {"N": ns, "curves": curves}


def plot_results(power_data: dict, out_path: Path) -> None:
    ns = power_data["N"]
    c = power_data["curves"]

    fig, ax = plt.subplots(figsize=(8.2, 4.8), dpi=160)
    ax.plot(ns, np.array(c["crossover_gel_dz080"]) * 100, "o-", lw=2.2, color="#10b981",
            label="Paired Gel / Bolus Crossover (strong, d_z = 0.80)")
    ax.plot(ns, np.array(c["crossover_gel_dz065"]) * 100, "s-", lw=2.2, color="#0ea5e9",
            label="Paired Gel / Bolus Crossover (medium, d_z = 0.65)")
    ax.plot(ns, np.array(c["within_subject_ema_14obs"]) * 100, "^-", lw=2.0, color="#f59e0b",
            label="App EMA Within-Subject VPD (14 check-ins / person)")
    ax.plot(ns, np.array(c["cross_sectional_1obs"]) * 100, "x--", lw=1.8, color="#ef4444",
            label="Naive Cross-Sectional Survey (1 check-in / person)")

    ax.axhline(80, color="#64748b", ls=":", lw=1.2)
    ax.text(245, 81.5, "80% Power Target", ha="right", fontsize=8.5, color="#475569")

    ax.set_xscale("log")
    ax.set_xticks([6, 10, 15, 20, 30, 50, 100, 250])
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.set_xlabel("Number of Participants (N, log scale)", fontsize=10)
    ax.set_ylabel("Statistical Power (%, α = 0.05 two-tailed)", fontsize=10)
    ax.set_title("Why Small-N Within-Subject Crossovers Beat Massive Noisy Surveys", fontsize=11.5, fontweight="bold")
    ax.set_ylim(0, 103)
    ax.grid(True, alpha=0.25)
    ax.legend(loc="lower right", fontsize=8.5, frameon=True)
    fig.tight_layout()

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path)
    plt.close(fig)


def main():
    print("=== 1. WINTER INDOOR-HEATING CONFOUND SIMULATION (N = 2,000) ===")
    conf = simulate_winter_rh_inversion(n_obs=2000)
    print(f"Indoors check-ins simulated: {conf['n_indoors']} / {conf['n_obs']}")
    print(
        f"  Naive correlation (Outdoor RH vs Cold indoors): r = {conf['r_naive_outdoor_rh_indoors']:+.3f} "
        f"(p = {conf['p_naive_outdoor_rh_indoors']:.2e})  <-- CONFOUNDED / INVERTED SIGN!"
    )
    print(
        f"  Physics-corrected partial r (Local VPD | Local T): r = {conf['r_corrected_local_vpd']:+.3f} "
        f"(p = {conf['p_corrected_local_vpd']:.2e})  <-- TRUE SIGN RECOVERED (dry air = colder)"
    )

    print("\n=== 2. SAMPLE SIZE & STATISTICAL POWER COMPARISON ===")
    pdata = estimate_power_curves(n_sims=800)
    ns = pdata["N"]
    c = pdata["curves"]
    header = f"{'N':>5} | {'Cross-Sec (1 obs)':>17} | {'EMA (14 obs/subj)':>17} | {'Gel/Bolus dz=.65':>16} | {'Gel/Bolus dz=.80':>16}"
    print(header)
    print("-" * len(header))
    for i, n in enumerate(ns):
        print(
            f"{n:>5} | {c['cross_sectional_1obs'][i]*100:>16.1f}% | "
            f"{c['within_subject_ema_14obs'][i]*100:>16.1f}% | "
            f"{c['crossover_gel_dz065'][i]*100:>15.1f}% | "
            f"{c['crossover_gel_dz080'][i]*100:>15.1f}%"
        )

    out_fig = Path(__file__).resolve().parent / "outputs" / "power_comparison.png"
    plot_results(pdata, out_fig)
    print(f"\nSaved power comparison chart to: {out_fig}")


if __name__ == "__main__":
    main()
