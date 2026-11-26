# Electro-Metabolic Oscillator Model of Cutaneous TRPM8 Osmosensation & Cold Coding

Open-source computational biophysics, human-genetics modeling, and psychrometric mobile-application suite accompanying:

> **Yorek MS.** *An Electro-Metabolic Relaxation Oscillator Model of Peripheral TRPM8 Osmosensation and Cold Coding: Human Genetic, Mitochondrial, and Neurovascular Predictions* (bioRxiv Submission Draft v2.0).

Built around a **100% independent, reproducible Human-Behavior-First → Molecular-Specialization** research program uniting:
- **Classical cold-fiber electrophysiology (1935–2002):** Endogenous $\text{Ca}^{2+}$-controlled oscillating inward currents, $\text{CO}_2$/acidosis inhibition, metabolic sensitivity, biphasic menthol paradox, and depolarized TTX-resistant cold terminals.
- **Molecular TRPM8 osmosensation & bioenergetics (2002–2026):** Corneal and somatic DRG hyperosmolar TRPM8 gating ([Quallo et al., 2015](https://europepmc.org/article/PMC/PMC4455064); [Yorek MS et al., *IOVS* 2016](https://europepmc.org/article/PMC/PMC5113984)), small-DRG L-type $\text{Ca}_\text{v}1.2/\text{Ca}_\text{v}1.3$ amplification, physical `TMBIM5–MICU1` mitochondrial $\text{Ca}^{2+}/\text{H}^+$ coupling, diabetic nociceptor `MCU` deletion protection, and cutaneous neurovascular dihydropyridine pharmacology.

---

## Repository Map

```text
trpm8-osmo-cold/
├── README.md
├── docs/
│   ├── manuscript_biorxiv_v2.md                   # Formal IMRaD bioRxiv / journal submission manuscript (v2.0)
│   ├── manuscript_draft.md                        # Working theoretical synthesis draft (v1.0)
│   ├── 01_theoretical_framework.md                # Formal biophysical & anti-labeled-line theory
│   ├── 02_human_protocols_H1_to_H5.md             # Preregistered human falsification protocols (incl. Fast-Cooling Gel Submersion)
│   ├── 03_nsf_grant_strategy_and_aims.md          # NSF PAPPG eligibility & Specific Aims planning notes
│   ├── 04_manuscript_outline_and_data_map.md      # Manuscript architecture & reproducibility map
│   └── 05_classical_oscillator_to_cmr1_synthesis.md # 90-year historical cold-fiber literature synthesis
├── analysis/
│   ├── trpm8_oscillator_ode.py                    # 5-ODE electro-metabolic oscillator + AP spiking + human genetics + neurovascular module (Figs 3–5)
│   ├── build_drive_figures.py                     # 100% original copyright-free Figures 1 & 2 generator
│   ├── power_and_simulation.py                    # Monte Carlo psychrometric sign-inversion & power simulator (Fig 6)
│   ├── analyze_ema_export.py                      # Automated Popperian analysis & falsification pipeline for EMA exports
│   └── outputs/                                   # High-DPI publication figures (Figs 1–6 + STIX equation blocks)
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

Requires Python 3.10+ with `numpy`, `scipy`, `matplotlib`, `pillow`, and `pandas` (see `requirements.txt`). Each script also embeds [PEP 723](https://peps.python.org/pep-0723/) dependency metadata, so [`uv`](https://docs.astral.sh/uv/) runs them directly from cache without creating a `.venv` inside Google Drive:

```bash
# Option A (Recommended with uv — zero local .venv needed):
uv run analysis/build_drive_figures.py
uv run analysis/trpm8_oscillator_ode.py
uv run analysis/power_and_simulation.py
uv run analysis/analyze_ema_export.py --demo

# Option B (Standard Python venv / pip):
python3 -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

### 0. Generate Figures 1, 2, and STIX Equation Blocks (Eqs. 1–7)
```bash
uv run analysis/build_drive_figures.py
```
Writes `analysis/outputs/fig1_corneal_nerve_interface.png`, `analysis/outputs/fig2_mechanistic_circuit.png`, and `analysis/outputs/eq_block{1,2,3}_*.png`.

### 1. Generate Figures 3, 4, and 5 (Core Oscillator + AP Spiking, Human Genetics/MCU, Neurovascular Pharmacology)
```bash
python3 analysis/trpm8_oscillator_ode.py
```
Outputs all quantitative tables to stdout and writes:
- `analysis/outputs/trpm8_oscillator_panels.png` (**Figure 3**)
- `analysis/outputs/fig4_human_genetics_mcu.png` (**Figure 4**)
- `analysis/outputs/fig5_neurovascular_pharmacology.png` (**Figure 5**)

### 2. Generate Figure 6 (Psychrometric Sign-Inversion & Statistical Power)
```bash
python3 analysis/power_and_simulation.py
```
Writes `analysis/outputs/power_comparison.png` (**Figure 6**).

### 3. Launch the Companion Psychrometric Mobile PWA Locally
```bash
python3 -m http.server 8765 --bind 127.0.0.1 -d app
```
Open `http://127.0.0.1:8765` in any browser, or host `app/` statically on GitHub Pages / Netlify Drop (HTTPS enables GPS + Home Screen installation; see `cloud/DEPLOY_GUIDE.md`).

### 4. Run the Automated Popperian Analysis Pipeline on Exported EMA Data
```bash
python3 analysis/analyze_ema_export.py --demo
python3 analysis/analyze_ema_export.py --input path/to/trpm8_export.json
```
