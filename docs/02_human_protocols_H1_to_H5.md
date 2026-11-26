# Phase 1 Preregistered Human Psychophysics Protocols (H1–H5) & Tier-1 Home Physiology Rig

**Strategic Order:** Establish reproducible, falsification-resistant **human behavioral phenotypes first** (low cost, zero wet lab required) before investing in Phase 2 molecular/channel-modification assays.

> **Safety & Ethics Note:** All Phase 1 protocols use non-invasive, everyday physiological stimuli (commercial preservative-free ophthalmic saline/artificial tears within safe osmolarity ranges, food-grade hydrogel/saline hand baths, drinking water/electrolytes, OTC petrolatum/menthol, and benign skin cooling ≥15 °C). Formal multi-participant publication requires IRB exemption/expedited approval; self-piloting ($N=1$) can begin immediately to calibrate effect sizes and timing.

---

## Protocol H2-Gel (Flagship Experiment): Fast-Cooling Gel Submersion — Eliminating Evaporative Cooling Completely ($J_{\text{evap}} = 0$)

### Why This Solves the Core Physics Confound
In air, dry skin or low humidity *always* accelerates evaporative heat loss ($\dot{q}_{\text{evap}} = L_v J_{\text{evap}}$). **Submerging the hand/forearm inside an aqueous or polymer hydrogel bath physically abolishes the air–skin vapor interface ($J_{\text{evap}} \equiv 0$).** Heat transfer is 100% conductive/convective, while **bath temperature ($T_{\text{bath}}$)** and **bath osmolarity ($\text{Osm}_{\text{bath}}$)** become completely independent, measurable dials!

### Materials (~$40 total)
- Two identical 2–4 L insulated containers (for paired simultaneous Left vs. Right hand/forearm immersion, eliminating memory drift!).
- Unscented carbomer / xanthan / aloe-free **hydrogel base** (or pure water/agarose 0.3% viscous gel — viscosity suppresses convective plumes and gives uniform rapid conductive contact!).
- **Osmolyte:** Food-grade NaCl (or non-ionic mannitol / glycerol control to rule out ionic strength / cutaneous electrical conductivity artifacts!! *Using mannitol vs. NaCl at matched osmolarity is an elite biophysical control*).
- Two calibrated digital thermocouple / probe thermometers (0.1 °C accuracy, verified side-by-side in the same water cup before every session).
- Consumer infrared (IR) forehead/surface thermometer (or phone thermal camera).

### 2×2 Within-Subject Crossover Design (Blinded to Osmolarity)
Have a helper (or coded bottles A/B) prepare **two thermally matched gel baths** within $\pm 0.1\,^\circ\text{C}$ of each other:
- **Condition 1 (Isothermal Room-Cool Fast Submersion, $22.0\,^\circ\text{C}$):**
  - Bath A: **Isoosmotic gel** (~300 mOsm/kg; ~0.9% NaCl or 5.5% w/v mannitol)
  - Bath B: **Hyperosmotic gel** (~900–1200 mOsm/kg; ~3.0–3.5% NaCl or ~18% mannitol; note stratum corneum diffuses slowly, so also test a **15-min pre-soak** or **lip/vermilion / mild tape-stripped volar forearm** patch!)
- **Condition 2 (Skin-Temperature Matched Control, $33.0\,^\circ\text{C}$):**
  - Identical osmotic contrast at thermally neutral skin temperature (33.0 °C) — isolates **pure non-thermal osmotic TRPM8 gating** with zero cooling AND zero evaporation!
- **Condition 3 (Pre-Hydrated vs. Dry-Conditioned Hand into Identical $22.0\,^\circ\text{C}$ Neutral Gel):**
  - Pre-condition Hand L (15 min sealed wet warm glove → hydrated stratum corneum) vs. Hand R (15 min warm dry fan → desiccated stratum corneum), equilibrate BOTH hands to identical skin surface temperature (verified by IR thermometer within 0.2 °C!), then **simultaneously plunge both hands into the SAME $22.0\,^\circ\text{C}$ fast-cooling gel bath**!
  - *Thermal conductivity control:* Hydrated stratum corneum actually has slightly **higher** thermal conductivity ($k$) than dry stratum corneum (meaning hydrated skin physically loses heat *faster* on contact!). Therefore, if the **dry-preconditioned hand still feels COLDER** upon simultaneous immersion into the identical gel bath despite lower conductive heat flux, **pure thermal physics is decisively falsified in favor of osmotic/channel sensitization!!**

### Measurements (at $t = 5\text{ s}, 30\text{ s}, 60\text{ s}, 180\text{ s}$)
- **Forced-choice comparison:** "Which hand feels colder?" (Left / Right) + magnitude difference (−5 to +5).
- **4-Axis Qualia (0–10 each hand):** `Coldness`, `Wetness/Clamminess`, `Dry/Tight/Prickly`, `Chill/Shiver urge`.
- Swap Left/Right assignments on Day 2 (counterbalanced!).

### Severe Popperian Falsification Criterion
- **Falsifies Osmotic Cold Theory:** If, once bath temperatures match ($\pm 0.1\,^\circ\text{C}$) and pre-immersion IR skin temperatures match ($\pm 0.2\,^\circ\text{C}$), **hyperosmotic gel** and **dry-preconditioned skin** feel *equal or warmer* than isoosmotic / hydrated skin across $N \ge 15$ counterbalanced trials (paired two-tailed $p > 0.05$ and Cohen's $d < 0.2$).
- **Corroborates Theory:** Dry-preconditioned or hyperosmotic-exposed skin feels significantly **colder** despite $J_{\text{evap}} = 0$ (and despite dry skin having lower thermal conductivity!).

---

## Protocol H1: Cornea vs. Palpebral Lid Margin vs. Eyelid Skin (Anti-Labeled-Line $V_1$ Mapping)

Extends your 2016 rat hyperosmolar blink assay ([Yorek MS et al., *IOVS* 2016](https://europepmc.org/article/PMC/PMC5113984)) to human qualia dissociation across adjacent trigeminal territories:
- **Stimuli (Strictly Isothermal, warmed in water bath to $34.0\,^\circ\text{C}$):**
  1. Sterile preservative-free isotonic saline (290–300 mOsm — standard Refresh/saline ampoule).
  2. Sterile preservative-free mildly hypertonic ophthalmic saline (450–500 mOsm, well within standard clinical corneal hyperosmolarity testing; *do NOT exceed safe ophthalmic range without IRB/clinical supervision—3% hypertonic ophthalmic drops such as OTC Muro 128 2%–5% are FDA-OTC approved specifically for human eyes, ~600–1500 mOsm, though 2% or diluted sterile 450 mOsm is gentler*).
  3. **Menthol vapor control** (non-contact cool vapor held 3 cm from closed eyelid vs. open eye).
- **Sites:** (A) Open cornea (drop/vapor), (B) Inner palpebral lid margin, (C) Closed upper eyelid skin.
- **Readouts:** 30-sec video-recorded **blink/squint rate** (built into the companion app's Trial D!) + 5-axis qualia (`Cool`, `Wet`, `Dry/Gritty`, `Fresh`, `Sting/Blink urge`).
- **Falsification Criterion:** If corneal, palpebral, and eyelid sites report an identical qualia profile (differing only in scalar intensity) to isothermal hyperosmolarity vs. menthol vapor, then peripheral tissue-specific coding across $V_1$ is falsified.

---

## Protocol H3: Breaking Scalar Labeled Lines via Temporal Pattern & $\text{H}^+$ Sensitivity (Menthol vs. Cold)

Directly translates your Usachev-lab **Slide 33–34** ($\text{H}^+$ feedback block) into human psychophysics:
- **Part A — Menthol ($\text{H}^+$-Insensitive) vs. Cold ($\text{H}^+$-Sensitive) Time-Course:**
  - Apply 1% topical menthol gel on Forearm A (occluded under Tegaderm to prevent evaporation!) vs. a clamped $24\,^\circ\text{C}$ dry Peltier/gel pack on Forearm B, matched for initial perceived cold intensity at $t = 30\text{ s}$.
  - Participants continuously track perceived **Cold**, **Wetness**, and **Waxing/Waning Pulse** every 20 s for 6 min.
  - *Optional metabolic $\text{pH}_i$ perturbation:* Repeat immediately after 90 s of rhythmic isometric forearm grip exercise (local tissue lactate/$\text{H}^+$ rise!).
  - **Prediction (Andersson 2004 + Yorek Slide 34):** Local acidosis suppresses **cold** intensity/oscillation while **sparing menthol**!
- **Part B — Oscillatory vs. Static Cooling Waveform (Tier-1 Peltier Rig):**
  - Deliver 4 minutes of **steady $24.0\,^\circ\text{C}$** vs. **slow sinusoidal $24.0 \pm 1.5\,^\circ\text{C}$ at $0.03\text{ Hz}$** (period ~33 s, matching Slide 17–18's burst-oscillation timescale; identical mean temperature!).
  - **Falsification Criterion:** If menthol and cold respond identically to local acidosis AND steady vs. $0.03\text{ Hz}$ oscillating cooling produce identical wetness/chill qualia, the $\text{H}^+$-clocked temporal-code model is falsified.

---

## Protocol H4: Fast Water Drinking 2×2 Crossover (Direct Osmolarity vs. TRPV4 Vasoconstriction vs. Thermal Sink)

Tests your water-drinking chill prediction cleanly while separating **all three competing physiological mechanisms**:

| Arm (500 mL consumed in <2 min, blinded opaque insulated cup + straw) | Hypoosmotic? | Cold Thermal Sink? | Prediction: Direct Osmotic TRPM8 | Prediction: Hepatic TRPV4 Vasoconstriction | Prediction: Gastric Thermal Sink |
| :--- | :---: | :---: | :--- | :--- | :--- |
| **1. 37 °C Plain Water** | **Yes** | No | Cold/shiver ↑ **without** (or before) skin cooling | Fingertip IR temp drops 2–4 °C **first**, then cold ↑ | No effect |
| **2. 37 °C Isotonic Oral Rehydration / 0.9% NaCl** | No | No | No effect (Baseline control) | No effect (Baseline control) | No effect |
| **3. 15 °C Isotonic Oral Rehydration / 0.9% NaCl** | No | **Yes** | Mild/no osmotic effect | Thermal only | Cold/shiver ↑ equally to Arm 4 |
| **4. 15 °C Plain Water** | **Yes** | **Yes** | Strongest cold/shiver (synergy) | Vasoconstriction + thermal | Equal to Arm 3 |

- **Timepoints ($t = -2, 0, 2, 5, 10, 15, 20, 30\text{ min}$):**
  1. **Fingertip + Volar Forearm IR Skin Temp** (critical discriminator for TRPV4 sympathetic vasoconstriction!),
  2. **Coldness (−4 to +4), Shivering/Goosebumps, Sudden Chill/Shudder**,
  3. **Standardized Cold Touch Test** (touch room-temperature steel block/spoon to volar wrist: rate 0–10 coldness).
- **Falsification Criterion:** If **Arm 1 (37 °C plain water)** causes zero change in cold sensitivity when fingertip IR temperature is controlled—or causes chills strictly downstream of a drop in fingertip IR skin temperature—direct peripheral hypoosmotic TRPM8 sensitization is **falsified** (and hepatic TRPV4 osmopressor vasoconstriction wins).

---

## Protocol H5: Distributed Smartphone EMA + Embedded N-of-1 Trials
Implemented directly in `app/`! Captures 4-axis qualia (`Cold`, `Wet/Clammy`, `Dry Skin`, `Dry Eyes`, `Shiver/Shudder`), **water volume + recency + temperature + tonicity**, and live **Open-Meteo weather → Dew-Point-Derived Indoor VPD**, PLUS built-in randomized N-of-1 timers for **H2-Gel**, **H4 Water Bolus**, **Forearm Balm Occlusion**, and **H1 Blink Count**.

---

## Tier-1 Home Human Physiology Rig (Total Cost: ~$150–$420)

| Component | Purpose | Approximate Cost |
| :--- | :--- | :---: |
| **Smartphone IR Thermal Camera** (e.g., Topdon TC001 / FLIR One / InfiRay P2 Pro) *or* 0.1 °C medical IR thermometer | Separates **TRPV4 cutaneous vasoconstriction** & **evaporative skin cooling** from direct osmotic/neural cold perception (critical for H2-Gel, H3, H4) | $35 (IR gun) – $250 (thermal camera) |
| **Dual Digital Thermocouple Meter** (e.g., K-type 2-channel, 0.1 °C resolution) | Verifies paired Gel Submersion baths match within $\pm 0.1\,^\circ\text{C}$ (H2-Gel) | $25–$45 |
| **ESP32 / Arduino + 40×40 mm TEC1-12706 Peltier + BTS7960 H-Bridge + DS18B20/NTC Thermistor + Aluminum Heatsink** | DIY programmable QST thermode delivering exact **steady** vs. **0.03 Hz sinusoidal oscillatory** cooling ramps (H3-Part B) | $40–$60 |
| **Precision Indoor Thermo-Hygrometer** (Sensirion SHT45 / Govee Bluetooth logger) | Ground-truth indoor VPD validation alongside the app's dew-point weather model (H5) | $15–$30 |
| **Hydrogel + NaCl + Food-Grade Mannitol + Tegaderm Film** | Non-evaporative gel submersion (H2-Gel) & occlusion splits (H2/H3) | $30–$40 |
