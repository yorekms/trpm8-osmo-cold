# Electro-Metabolic Oscillator Model of Cutaneous TRPM8 Osmosensation & Cold Coding (v2.3)

Open-source computational biophysics, cryo-EM/AlphaFold molecular kill-switches, human-genetics modeling, and psychrometric mobile-application suite accompanying:

> **Yorek MS (2026).** *An Electro-Metabolic Relaxation Oscillator Model of Peripheral TRPM8 Osmosensation and Cold Coding: Human Genetic, Mitochondrial, and Neurovascular Predictions* (Open Preprint v2.3 • 7 Figures • 6 Tables • 189 Peer-Reviewed Citations).

- **Live Interactive Portal & 16-Preset ODE Cockpit:** https://yorekms.github.io/trpm8-osmo-cold/
- **Open Preprint PDF (v2.3):** https://yorekms.github.io/trpm8-osmo-cold/assets/TRPM8_Electro_Metabolic_Oscillator_Preprint_Yorek.pdf
- **Companion Dew-Point Local VPD & N-of-1 Psychophysics Web App:** https://yorekms.github.io/trpm8-osmo-cold/app/index.html

Built around a **100% independent, reproducible Human-Behavior-First → Molecular-Specialization** research program uniting:
- **Classical & corneal cold-fiber electrophysiology (1935–2019):** Endogenous $\text{Ca}^{2+}$-controlled oscillating inward currents, $\text{CO}_2$/acidosis inhibition, metabolic sensitivity, biphasic menthol paradox, and James A. Brock's guinea-pig and mouse corneal nerve-terminal impulse (NTI) corpus proving chronic terminal depolarization ($\ge -40\text{ mV}$) and mitochondria-packed TRPM8+ complex bulbar endings.
- **Cutaneous water-activity leverage & electro-metabolic dynamics (2002–2026):** Thermodynamic proof that 99.8% stratum corneum + SG2 tight-junction shielding still transmits $+56\text{ to }+89\text{ mOsm/kg}$ to intraepidermal TRPM8+ endings ($dc_{\text{osm}}/da_w \approx -55.8\text{ Osm/kg}$), coupled to $\text{Na}_\text{v}1.8$, low-voltage $\text{Ca}_\text{v}1.3/\text{Ca}_\text{v}1.2$, cooperative >300 nM `MICU1/MCU`, and physical `TMBIM5–MICU1` $\text{Ca}^{2+}/\text{H}^+$ exchange.
- **Head-to-head 6-ODE rival benchmarks & residue-level kill-switches (v2.3):** Quantitative 6-ODE comparisons against linear PMCA $\text{H}^+$ import, $\text{Ca}^{2+}$-PLC$\delta$ $\text{PI(4,5)P}_2$ depletion, and pure $\text{Ca}_\text{v}3.2/I_h$ pacing (**Figure 7**), paired with cryo-EM (`8BDC`, `9ZCU`, `6WDN`) and AlphaFold (`AF-Q9H3K2`) point-mutation kill-switches (`TRPM8 His845Gln/Ala` & `His667Gln/Ala`; `TMBIM5 Asp294Arg/Asp325Arg`).

---

## Repository Map

```text
trpm8-osmo-cold/
├── index.html                                     # Interactive 16-Preset ODE Cockpit, SC Leverage Calculator, Figures 1–7 & Tables 5–6
├── assets/
│   └── TRPM8_Electro_Metabolic_Oscillator_Preprint_Yorek.pdf  # Complete typeset preprint PDF (v2.3, 189 refs)
├── docs/
│   ├── manuscript_biorxiv_v2.md                   # Formal IMRaD bioRxiv / journal submission manuscript (v2.3, 189 refs)
│   ├── 06_competing_alternatives_and_rebuttals.md # Dialectical steelman & quantitative/molecular rebuttals to 6 rival paradigms
│   ├── 05_classical_oscillator_to_cmr1_synthesis.md # 90-year historical & corneal cold-fiber literature synthesis
│   ├── 04_manuscript_outline_and_data_map.md      # Manuscript architecture & reproducibility map
│   ├── 03_nsf_grant_strategy_and_aims.md          # Specific Aims & experimental planning notes
│   ├── 02_human_protocols_H1_to_H5.md             # Preregistered human falsification protocols (incl. Zero-Evaporation H2-Gel)
│   └── 01_theoretical_framework.md                # Formal biophysical & anti-labeled-line theory
├── analysis/
│   ├── trpm8_oscillator_ode.py                    # 5-ODE electro-metabolic oscillator + AP spiking + genetics + neurovascular module (Figs 3–5)
│   ├── trpm8_rival_rebuttals_ode.py               # 6-ODE rival model benchmarks + cryo-EM/AlphaFold point-mutant kill-switches (Fig 7)
│   ├── build_drive_figures.py                     # 100% original copyright-free Figures 1 & 2 + STIX equation generator
│   ├── power_and_simulation.py                    # Monte Carlo psychrometric sign-inversion & power simulator (Fig 6)
│   ├── analyze_ema_export.py                      # Automated Popperian analysis & falsification pipeline for EMA exports
│   └── outputs/                                   # High-DPI publication figures (Figs 1–7 + STIX equation blocks)
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

Each Python script embeds [PEP 723](https://peps.python.org/pep-0723/) dependency metadata (`numpy`, `scipy`, `matplotlib`, `gemmi`, `pillow`, `pandas`), so [`uv`](https://docs.astral.sh/uv/) runs them directly without creating a local `.venv`:

```bash
uv run analysis/build_drive_figures.py        # Figures 1, 2 & Equations 1–7
uv run analysis/trpm8_oscillator_ode.py       # Figures 3, 4 & 5 (29 physiological, genetic & neurovascular regimes)
uv run analysis/power_and_simulation.py       # Figure 6 (Monte Carlo psychrometrics & crossover power)
uv run analysis/trpm8_rival_rebuttals_ode.py  # Figure 7 (6-ODE rival models & TRPM8 H845Q / TMBIM5 D294R-D325R mutants)
uv run analysis/analyze_ema_export.py --demo  # Popperian EMA analysis pipeline
```
