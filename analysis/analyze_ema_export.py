#!/usr/bin/env python3
# /// script
# requires-python = ">=3.10"
# dependencies = [
#   "numpy>=1.26",
#   "pandas>=2.1",
#   "scipy>=1.11",
# ]
# ///
"""
Popperian Analysis Pipeline for TRPM8 Cold & Hydration Study Exports.

Ingests exported JSON (single or multiple participant files) or CSV from app/,
runs controlled statistical tests across H1–H5, and prints explicit
CORROBORATED / FALSIFIED / UNDERPOWERED verdicts for each mechanistic link.
"""

import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats


def generate_demo_dataset(seed: int = 7) -> dict:
    """Create a realistic synthetic multi-participant pilot dataset (N=20) to verify all tests."""
    rng = np.random.default_rng(seed)
    checkins = []
    trials = []

    for pid_idx in range(20):
        pid = f"P-DEMO{pid_idx:02d}"
        subj_offset = rng.normal(0, 0.9)

        # 12 EMA check-ins per participant
        for k in range(12):
            out_t = rng.uniform(0.0, 24.0)
            out_rh = np.clip(75 - 0.7 * out_t + rng.normal(0, 10), 25, 95)
            indoors = rng.random() < 0.65
            loc_t = rng.normal(21.0, 1.0) if indoors else out_t
            es_out = 0.61094 * np.exp(17.625 * out_t / (out_t + 243.04))
            ea_out = es_out * (out_rh / 100.0)
            es_loc = 0.61094 * np.exp(17.625 * loc_t / (loc_t + 243.04))
            ea_loc = min(ea_out, es_loc * 0.65) if indoors else ea_out
            loc_vpd = max(0.0, es_loc - ea_loc)

            water_recent = rng.random() < 0.30
            w_vol = rng.choice(["medium", "large"]) if water_recent else "none"
            w_rec = rng.choice(["lt10", "10_30"]) if water_recent else "none"
            w_typ = rng.choice(["plain", "electrolyte"]) if water_recent else "none"

            dry_skin = int(np.clip(round(1.5 + 3.2 * loc_vpd + rng.normal(0, 1.2)), 0, 10))
            dry_eyes = int(np.clip(round(1.0 + 3.0 * loc_vpd + rng.normal(0, 1.3)), 0, 10))

            # Synthetic true DGP: dry air feels colder (-0.65/kPa), plain recent water cools (-0.45)
            water_shift = -0.45 if (water_recent and w_typ == "plain") else 0.0
            therm = float(
                np.clip(
                    -4.0 + 0.21 * loc_t - 0.65 * loc_vpd + water_shift + subj_offset + rng.normal(0, 0.55),
                    -4.0,
                    4.0,
                )
            )
            checkins.append(
                {
                    "participantId": pid,
                    "locationMode": "indoors" if indoors else "outside",
                    "outdoorTempC": round(out_t, 2),
                    "outdoorRhPct": round(out_rh, 1),
                    "localTempC": round(loc_t, 2),
                    "localVpdKPa": round(loc_vpd, 3),
                    "thermalRating": round(therm * 2) / 2,
                    "drySkinRating": dry_skin,
                    "dryEyesRating": dry_eyes,
                    "wetnessRating": int(np.clip(round(1 + rng.normal(0, 1.2)), 0, 10)),
                    "waterVolume": w_vol,
                    "waterRecency": w_rec,
                    "waterType": w_typ,
                }
            )

        # Trial A: Fast-cooling gel submersion (1 per participant)
        ctrl_cold = rng.normal(5.0, 1.0)
        delta_gel = rng.normal(1.1, 1.1)  # true medium-strong non-evaporative osmotic effect
        test_cold = ctrl_cold + delta_gel
        fc = "test_colder" if delta_gel > 0.25 else ("ctrl_colder" if delta_gel < -0.25 else "same")
        trials.append(
            {
                "participantId": pid,
                "trialType": "gel_submersion",
                "mode": "precond_dry_vs_hydrated",
                "bathTempC": 22.0,
                "preSkinDeltaC": round(float(rng.normal(0, 0.08)), 2),
                "forcedChoice": fc,
                "testColdRating": round(float(np.clip(test_cold, 0, 10)), 1),
                "ctrlColdRating": round(float(np.clip(ctrl_cold, 0, 10)), 1),
            }
        )

        # Trial B: Water Bolus 2x2 (15 min post-bolus across 4 arms)
        for arm in ["warm_water", "warm_saline", "cool_saline", "cool_water"]:
            # Simulate TRPV4 vasoconstriction scenario: warm_water drops finger temp AND chills
            if arm == "warm_water":
                finger_dt = rng.normal(-1.8, 0.7)
                therm_15 = rng.normal(-0.8, 0.6)
                shiver = int(np.clip(round(rng.normal(3.2, 1.5)), 0, 10))
            elif arm == "warm_saline":
                finger_dt = rng.normal(-0.1, 0.5)
                therm_15 = rng.normal(0.0, 0.5)
                shiver = int(np.clip(round(rng.normal(0.7, 0.8)), 0, 10))
            elif arm == "cool_saline":
                finger_dt = rng.normal(-0.9, 0.6)
                therm_15 = rng.normal(-1.0, 0.6)
                shiver = int(np.clip(round(rng.normal(3.5, 1.4)), 0, 10))
            else:
                finger_dt = rng.normal(-2.3, 0.7)
                therm_15 = rng.normal(-1.6, 0.6)
                shiver = int(np.clip(round(rng.normal(5.0, 1.5)), 0, 10))

            trials.append(
                {
                    "participantId": pid,
                    "trialType": "water_bolus_2x2",
                    "arm": arm,
                    "timepointMin": 15,
                    "fingerIrTempC": round(float(33.0 + finger_dt), 2),
                    "thermalRating": round(float(np.clip(therm_15, -4, 4)), 1),
                    "shiverUrge": shiver,
                    "fingerColdness": int(np.clip(round(3 - finger_dt * 1.5), 0, 10)),
                }
            )

        # Trial D: Cornea blink check
        vpd_b = rng.uniform(0.4, 2.2)
        trials.append(
            {
                "participantId": pid,
                "trialType": "cornea_blink",
                "blinksPerMin": int(np.clip(round(10 + 7.5 * vpd_b + rng.normal(0, 3)), 4, 40)),
                "corneaDryRating": int(np.clip(round(1 + 3.4 * vpd_b + rng.normal(0, 1.2)), 0, 10)),
                "eyelidCoolRating": int(np.clip(round(3 + 0.4 * vpd_b + rng.normal(0, 1.5)), 0, 10)),
                "localVpdKPa": round(float(vpd_b), 2),
            }
        )

    return {"checkins": checkins, "trials": trials}


def load_inputs(paths: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    all_checkins = []
    all_trials = []
    for p in paths:
        fp = Path(p)
        if fp.suffix.lower() == ".json":
            data = json.loads(fp.read_text())
            all_checkins.extend(data.get("checkins", []))
            all_trials.extend(data.get("trials", []))
        elif fp.suffix.lower() == ".csv":
            df = pd.read_csv(fp)
            if "trialType" in df.columns:
                all_trials.extend(df.to_dict(orient="records"))
            else:
                all_checkins.extend(df.to_dict(orient="records"))
    return pd.DataFrame(all_checkins), pd.DataFrame(all_trials)


def ols_summary(y: np.ndarray, X: np.ndarray, names: list[str]) -> list[dict]:
    n, k = X.shape
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    s2 = float(np.sum(resid**2) / max(1, n - k))
    cov = s2 * np.linalg.pinv(X.T @ X)
    se = np.sqrt(np.maximum(np.diag(cov), 1e-12))
    t_vals = beta / se
    p_vals = 2.0 * stats.t.sf(np.abs(t_vals), df=max(1, n - k))
    return [
        {"name": names[i], "beta": float(beta[i]), "se": float(se[i]), "t": float(t_vals[i]), "p": float(p_vals[i])}
        for i in range(k)
    ]


def analyze_checkins(df: pd.DataFrame) -> None:
    print("\n================ H5: EMA CHECK-INS (VPD & WATER RECENCY) ================")
    if df.empty or len(df) < 8:
        print(f"Insufficient check-ins with weather (n={len(df)}). Collect ≥15 check-ins.")
        return

    valid = df.dropna(subset=["localTempC", "localVpdKPa", "thermalRating"]).copy()
    print(f"Valid check-ins: {len(valid)} across {valid['participantId'].nunique()} participant(s)")
    if len(valid) < 8:
        print("Need ≥8 weather-linked check-ins for regression.")
        return

    # Within-subject centering (eliminates between-person baseline variance!)
    for col in ["thermalRating", "localTempC", "localVpdKPa", "drySkinRating"]:
        if col in valid.columns:
            valid[f"{col}_ws"] = valid[col] - valid.groupby("participantId")[col].transform("mean")

    valid["recentPlainWater"] = (
        (valid.get("waterRecency", "").isin(["lt10", "10_30"]))
        & (valid.get("waterType", "") == "plain")
    ).astype(float)
    valid["recentPlainWater_ws"] = (
        valid["recentPlainWater"] - valid.groupby("participantId")["recentPlainWater"].transform("mean")
    )

    X = np.column_stack([
        np.ones(len(valid)),
        valid["localTempC_ws"].to_numpy(),
        valid["localVpdKPa_ws"].to_numpy(),
        valid["recentPlainWater_ws"].to_numpy(),
    ])
    res = ols_summary(
        valid["thermalRating_ws"].to_numpy(),
        X,
        ["Intercept", "LocalTemp (°C)", "LocalVPD (kPa)", "RecentPlainWater (<30m)"],
    )
    for r in res[1:]:
        print(f"  {r['name']:<26} beta = {r['beta']:+.3f} ± {r['se']:.3f}  (t = {r['t']:+.2f}, p = {r['p']:.4g})")

    vpd_row = res[2]
    if vpd_row["p"] < 0.05 and vpd_row["beta"] < 0:
        print("  -> VERDICT (H5-VPD): CORROBORATED — Higher drying stress (VPD) predicts colder perception at matched local T.")
    elif vpd_row["p"] < 0.05 and vpd_row["beta"] > 0:
        print("  -> VERDICT (H5-VPD): FALSIFIED (Opposite Sign) — Drier air felt warmer at matched local T.")
    else:
        print("  -> VERDICT (H5-VPD): NULL / UNDERPOWERED — No significant VPD effect yet.")


def analyze_trials(df: pd.DataFrame) -> None:
    if df.empty:
        print("\nNo N-of-1 trials recorded yet.")
        return

    # --- TRIAL A: FAST-COOLING GEL SUBMERSION (H2-Gel) ---
    gel = df[df["trialType"] == "gel_submersion"].copy()
    if not gel.empty:
        print(f"\n================ H2-GEL: FAST-COOLING GEL SUBMERSION (n={len(gel)}) ================")
        diffs = gel["testColdRating"].astype(float) - gel["ctrlColdRating"].astype(float)
        mean_d = float(diffs.mean())
        sd_d = float(diffs.std(ddof=1)) if len(diffs) > 1 else np.nan
        dz = mean_d / sd_d if (len(diffs) > 1 and sd_d > 0) else np.nan
        t_stat, p_val = stats.ttest_1samp(diffs, 0.0) if len(diffs) >= 3 else (np.nan, np.nan)

        n_test = int((gel["forcedChoice"] == "test_colder").sum())
        n_ctrl = int((gel["forcedChoice"] == "ctrl_colder").sum())
        n_decisive = n_test + n_ctrl
        p_binom = float(stats.binomtest(n_test, n_decisive, 0.5).pvalue) if n_decisive >= 3 else np.nan

        print(f"  Mean Coldness Diff (Dry/Hyper − Hydrated/Iso): {mean_d:+.2f} (d_z = {dz:.2f}, t = {t_stat:.2f}, p = {p_val:.4g})")
        print(f"  Forced-Choice Wins: Test Colder = {n_test}, Control Colder = {n_ctrl} (Binomial p = {p_binom:.4g})")
        if len(gel) < 10:
            print("  -> VERDICT (H2-Gel): PILOT STAGE (need n ≥ 15 counterbalanced trials for formal Popperian decision).")
        elif p_val < 0.05 and mean_d > 0:
            print("  -> VERDICT (H2-Gel): CORROBORATED — Non-evaporative (J_evap=0) osmotic/dry cold sensitization detected!")
        else:
            print("  -> VERDICT (H2-Gel): FALSIFIED — Eliminating evaporation removed the dry/hyperosmotic cold effect.")

    # --- TRIAL B: WATER BOLUS 2x2 (H4) ---
    bol = df[(df["trialType"] == "water_bolus_2x2") & (df["timepointMin"] == 15)].copy()
    if not bol.empty:
        print(f"\n================ H4: WATER BOLUS 2×2 AT t=15 MIN (n={len(bol)}) ================")
        summary = bol.groupby("arm")[["thermalRating", "shiverUrge", "fingerIrTempC"]].mean()
        print(summary.round(2).to_string())

        ww = bol[bol["arm"] == "warm_water"]
        ws = bol[bol["arm"] == "warm_saline"]
        if len(ww) >= 4 and len(ws) >= 4:
            t_ch, p_ch = stats.ttest_ind(ww["thermalRating"], ws["thermalRating"])
            t_ft, p_ft = stats.ttest_ind(ww["fingerIrTempC"].dropna(), ws["fingerIrTempC"].dropna())
            print(f"  37°C Water vs 37°C Saline Thermal diff: t = {t_ch:+.2f}, p = {p_ch:.4g}")
            print(f"  37°C Water vs 37°C Saline Finger IR Temp diff: t = {t_ft:+.2f}, p = {p_ft:.4g}")
            if p_ch < 0.05 and ww["thermalRating"].mean() < ws["thermalRating"].mean():
                if p_ft < 0.05 and ww["fingerIrTempC"].mean() < ws["fingerIrTempC"].mean():
                    print("  -> VERDICT (H4): HYPOTHESIS SPLIT — 37°C water chills, ACCOMPANIED by finger cooling (consistent with TRPV4 sympathetic vasoconstriction! Test skin-temp clamping).")
                else:
                    print("  -> VERDICT (H4): CORROBORATED — 37°C water chills WITHOUT fingertip vasoconstriction (supports direct osmotic effect).")
            else:
                print("  -> VERDICT (H4): FALSIFIED — Isothermal 37°C hypoosmotic water did not induce cold vs isotonic saline.")

    # --- TRIAL D: CORNEA BLINK vs EYELID (H1) ---
    blk = df[df["trialType"] == "cornea_blink"].copy()
    if len(blk) >= 5:
        print(f"\n================ H1: CORNEA BLINK vs EYELID QUALIA (n={len(blk)}) ================")
        r_blk, p_blk = stats.pearsonr(blk["localVpdKPa"], blk["blinksPerMin"])
        r_cd, p_cd = stats.pearsonr(blk["localVpdKPa"], blk["corneaDryRating"])
        r_lc, p_lc = stats.pearsonr(blk["localVpdKPa"], blk["eyelidCoolRating"])
        print(f"  VPD vs Blink Rate:       r = {r_blk:+.2f} (p = {p_blk:.4g})")
        print(f"  VPD vs Cornea Dryness:   r = {r_cd:+.2f} (p = {p_cd:.4g})")
        print(f"  VPD vs Eyelid Skin Cool: r = {r_lc:+.2f} (p = {p_lc:.4g})")


def main():
    ap = argparse.ArgumentParser(description="Analyze TRPM8 Cold & Hydration Study Exports")
    ap.add_argument("--input", nargs="*", help="One or more exported JSON or CSV files from app/")
    ap.add_argument("--demo", action="store_true", help="Run on synthetic N=20 pilot dataset to verify pipeline")
    args = ap.parse_args()

    if args.demo:
        print("Running on SYNTHETIC N=20 DEMO DATASET (for pipeline verification)...")
        data = generate_demo_dataset()
        df_c = pd.DataFrame(data["checkins"])
        df_t = pd.DataFrame(data["trials"])
    elif args.input:
        df_c, df_t = load_inputs(args.input)
    else:
        ap.error("Provide --input <file.json ...> or --demo")

    analyze_checkins(df_c)
    analyze_trials(df_t)


if __name__ == "__main__":
    main()
