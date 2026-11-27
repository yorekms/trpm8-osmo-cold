# Electro-Metabolic Oscillator Model of Cutaneous TRPM8 Osmosensation & Cold Coding (v2.4)

Open-source computational biophysics, OpenMM/PROPKA cryo-EM electrostatics, cryo-EM/AlphaFold molecular kill-switches, human-genetics modeling, and psychrometric mobile-application suite accompanying:

> **Yorek MS (2026).** *An Electro-Metabolic Relaxation Oscillator Model of Peripheral TRPM8 Osmosensation and Cold Coding: Human Genetic, Mitochondrial, and Neurovascular Predictions* (Open Chapter-Scale Preprint v2.4 • 8 Figures • 8 Tables • 193 Peer-Reviewed Citations).

- **Live Interactive Portal & 16-Preset ODE Cockpit:** https://yorekms.github.io/trpm8-osmo-cold/
- **Open Preprint PDF (v2.4):** https://yorekms.github.io/trpm8-osmo-cold/assets/TRPM8_Electro_Metabolic_Oscillator_Preprint_Yorek.pdf
- **Companion Dew-Point Local VPD & N-of-1 Psychophysics Web App:** https://yorekms.github.io/trpm8-osmo-cold/app/index.html

Built around a **100% independent, reproducible Human-Behavior-First → Molecular-Specialization** research program uniting:
- **Classical & corneal cold-fiber electrophysiology (1935–2019):** Endogenous $\text{Ca}^{2+}$-controlled oscillating inward currents, $\text{CO}_2$/acidosis inhibition, metabolic sensitivity, biphasic menthol paradox, and James A. Brock's guinea-pig and mouse corneal nerve-terminal impulse (NTI) corpus proving chronic terminal depolarization ($\ge -40\text{ mV}$) and mitochondria-packed TRPM8+ complex bulbar endings.
- **Cutaneous water-activity leverage & electro-metabolic dynamics (2002–2026):** Thermodynamic proof that 99.8% stratum corneum + SG2 tight-junction shielding still transmits $+56\text{ to }+89\text{ mOsm/kg}$ to intraepidermal TRPM8+ endings ($dc_{\text{osm}}/da_w \approx -55.8\text{ Osm/kg}$), coupled to $\text{Na}_\text{v}1.8$, low-voltage $\text{Ca}_\text{v}1.3/\text{Ca}_\text{v}1.2$, cooperative >300 nM `MICU1/MCU`, and physical `TMBIM5–MICU1` $\text{Ca}^{2+}/\text{H}^+$ exchange.
- **Full 6-paradigm dialectical rebuttals, 6-ODE rival benchmarks & empirical OpenMM/PROPKA $\text{p}K_a$ electrostatics (v2.4):** Complete dialectical steelmanning and quantitative rebuttals of the six rival paradigms (**Table 5**), head-to-head 6-ODE comparisons against linear PMCA $\text{H}^+$ import, $\text{Ca}^{2+}$-PLC$\delta$ $\text{PI(4,5)P}_2$ depletion, and pure $\text{Ca}_\text{v}3.2/I_h$ pacing (**Figure 7**), and empirical `PDBFixer` + **OpenMM 8.1** Amber14/OBC2 + **PROPKA 3.5** $\text{p}K_a$ & tetrameric Wyman proton-linkage analysis across **11 experimental cryo-EM TRPM8 conformations** (`8BDC`, `9P8Y`, `9ZCU`, `9ZF0`, `9PB5`, `9ZCV`, `9VJN`, `8E4L`, `9VJR`, `9VJQ`, `9VNA`) and **human TMBIM5** (`AF-Q9H3K2`) (**Figure 8**, **Tables 6–7**)—showing `His667` upshifts to $\text{p}K_a = 7.05$ in desensitized TRPM8 (matching $\text{pH}_{1/2} = 7.13$), menthol depresses `His845` $\text{p}K_a$ to `3.22–4.23`, and `TMBIM5 Asp294–Asp325` splits into a $\text{p}K_a = 4.84 / 9.60$ carboxyl proton cannon.

---

## Repository Map

```text
trpm8-osmo-cold/
├── index.html                                     # Interactive 16-Preset ODE Cockpit, SC Leverage Calculator, Figures 1–8 & Tables 6–8
├── assets/
│   └── TRPM8_Electro_Metabolic_Oscillator_Preprint_Yorek.pdf  # Complete chapter-scale preprint PDF (v2.4, 193 refs)
├── docs/
│   ├── manuscript_biorxiv_v2.md                   # Formal chapter-scale bioRxiv / journal manuscript (v2.4, 8 Figs, 8 Tables, 193 refs)
│   ├── 06_competing_alternatives_and_rebuttals.md # Full dialectical steelman, 6-ODE & OpenMM/PROPKA rebuttals to 6 rival paradigms
│   ├── 05_classical_oscillator_to_cmr1_synthesis.md # 90-year historical & corneal cold-fiber literature synthesis
│   ├── 04_manuscript_outline_and_data_map.md      # Manuscript architecture & reproducibility map
│   ├── 03_nsf_grant_strategy_and_aims.md          # Specific Aims & experimental planning notes
│   ├── 02_human_protocols_H1_to_H5.md             # Preregistered human falsification protocols (incl. Zero-Evaporation H2-Gel)
│   └── 01_theoretical_framework.md                # Formal biophysical & anti-labeled-line theory
├── analysis/
│   ├── trpm8_oscillator_ode.py                    # 5-ODE electro-metabolic oscillator + AP spiking + genetics + neurovascular module (Figs 3–5)
│   ├── trpm8_rival_rebuttals_ode.py               # 6-ODE rival model benchmarks + cryo-EM/AlphaFold point-mutant kill-switches (Fig 7)
│   ├── trpm8_propka_openmm_pka.py                 # PDBFixer + OpenMM Amber14/OBC2 + PROPKA 3.5 electrostatics across 11 cryo-EM states & TMBIM5 (Fig 8)
│   ├── build_drive_figures.py                     # 100% original copyright-free Figures 1 & 2 + STIX equation generator
│   ├── power_and_simulation.py                    # Monte Carlo psychrometric sign-inversion & power simulator (Fig 6)
│   ├── analyze_ema_export.py                      # Automated Popperian analysis & falsification pipeline for EMA exports
│   └── outputs/                                   # High-DPI publication figures (Figs 1–8 + CSV/JSON pKa tables + STIX equations)
├── app/                                           # Zero-dependency installable Mobile PWA (iOS & Android, offline retry sync queue)
│   ├── index.html                                 # 4-axis EMA Check-in, N-of-1 Micro-Trials, °F/°C toggle, Dashboard & Export
│   ├── styles.css
│   ├── psychrometrics.js                          # Outdoor dew-point -> Indoor RH & Local Vapor Pressure Deficit (VPD) physics engine
│   ├── app.js                                     # Open-Meteo geocoding/weather, timers, SVG plots, CSV/JSON export, cloud queue
│   ├── manifest.webmanifest
│   └── sw.js
└── cloud/
    ├── google_apps_script.gs                      # $0 Google Sheets collector webhook
    ├── supabase_schema.sql                        # Free-tier Supabase PostgreSQL schema (INSERT-only RLS)
    └── DEPLOY_GUIDE.md                            # Static hosting (GitHub Pages / Netlify Drop) + Sheets setup guide
```

---

## Reproducing All Manuscript Figures & Simulations

Each Python script embeds [PEP 723](https://peps.python.org/pep-0723/) dependency metadata (`numpy`, `scipy`, `matplotlib`, `gemmi`, `openmm`, `pdbfixer`, `propka`, `pillow`, `pandas`), so [`uv`](https://docs.astral.sh/uv/) runs them directly without creating a local `.venv`:

```bash
uv run analysis/build_drive_figures.py            # Figures 1, 2 & Equations 1–7
uv run analysis/trpm8_oscillator_ode.py           # Figures 3, 4 & 5 (29 physiological, genetic & neurovascular regimes)
uv run analysis/power_and_simulation.py           # Figure 6 (Monte Carlo psychrometrics & crossover power)
uv run analysis/trpm8_rival_rebuttals_ode.py      # Figure 7 (6-ODE rival models & TRPM8 H845Q / TMBIM5 D294R-D325R mutants)
uv run analysis/trpm8_propka_openmm_pka.py --plot-only  # Figure 8 from cached JSON (omit --plot-only for full OpenMM relaxation)
uv run analysis/analyze_ema_export.py --demo      # Popperian EMA analysis pipeline
```
