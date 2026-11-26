# Standalone Theoretical & Computational Paper Blueprint (100% Independent Ownership)

## 1. Strategic Advantage: 100% Independent, Clean Provenance

This manuscript is designed to stand **completely on its own**—requiring **zero** outreach to former labs and **zero** reliance on unpublished archival data. Every empirical premise comes from:
1. **90 years of published primary electrophysiology (1935–2025):** Hensel & Zotterman (1951a/b), Dodt & Zotterman (1952), Hensel (1953), Pierau, Torrey & Carpenter (1974), Braun, Bade & Hensel (1980), Schäfer, Braun & Isenberg (1986, `PMC2228855`), Swandulla et al. (1987), Werth & Thayer (1994), Huber & Braun (1998/2000), Reid & Flonta (2001), McKemy/Julius CMR1 (2002), Viana & Belmonte (2002), Carr, Pianova & Brock (2002, `PMC2229520`), Andersson, Chase & Bevan (2004), Zimmermann et al. (2007), Orio & Olivares (2009–2015, `PMC4591370`), Quallo et al. (2015, `PMC4455064`), Gkika et al. (2015, `PMC4284226`), Paricio-Montesinos et al. (2020, `PMC7272120`), and Selescu, Zimmermann & Babes (2025, `PMC12830821`).
2. **Your own published first-author work:** Hyperosmolar TRPM8 corneal reflex assay in diabetic peripheral neuropathy ([Yorek MS et al., *IOVS* 2016, PMC5113984](https://europepmc.org/article/PMC/PMC5113984)).
3. **Original open-source biophysical modeling (`analysis/trpm8_oscillator_ode.py`)** + **Statistical power analysis (`analysis/power_and_simulation.py`)**.
4. **A preregistered, falsifiable Phase 1 human behavioral program** (`app/` + `docs/02_human_protocols_H1_to_H5.md`, including Fast-Cooling Gel Submersion and Water Bolus 2×2) that can be published either as a **standalone Theoretical + Modeling Hypothesis Paper first** (bioRxiv / *Frontiers in Cellular Neuroscience* / *Journal of General Physiology* Perspectives / *eLife* Ideas & Speculation / *Trends in Neurosciences* / *Neuroscience & Biobehavioral Reviews*) followed immediately by **Paper 2 (Empirical Human Behavioral Falsifications)**.

---

## 2. Working Title & Target Venues

**Working Title:**
> **Beyond the Cold Labeled Line: Re-Unifying Classical Oscillating Cold Fibers, CMR1/TRPM8 Osmosensation, and Mitochondrial Proton Feedback into a Falsifiable Theory of Thermal and Wetness Perception**

**Best Venues for a Standalone Theoretical + Biophysical Modeling Paper (Zero Wet-Lab Gatekeeping):**
1. **`bioRxiv` (Neuroscience) Preprint (Immediate, Free, Citable DOI)** — Establishes priority immediately and gives you a clean, peer-ready link for NSF SBIR, micro-grants, or new academic collaborators.
2. **Peer-Reviewed Theory / Modeling Venues:**
   - *Journal of Theoretical Biology* or *PLOS Computational Biology* (computational/biophysical oscillator + falsifiable predictions).
   - *Frontiers in Cellular Neuroscience* / *Frontiers in Physiology* (Hypothesis & Theory article type — specifically built for mechanistic synthesis + ODE modeling + testable predictions!).
   - *eLife* (Reviewed Preprint) or *Journal of General Physiology* (Perspectives / Modeling).
   - *Medical Hypotheses* or *Neuroscience & Biobehavioral Reviews* (strong home for the metabolic/mitochondrial cold intolerance & diabetic neuropathy link).

---

## 3. Section-by-Section Paper Structure & Figures

### Section 1 — Abstract (~250 words)
- **The Problem:** Post-2002 molecular cloning of CMR1/TRPM8 fostered a static "cold labeled-line" model that cannot account for non-evaporative wetness reflexes ("wet dog shakes"), adjacent trigeminal qualia divergence (cornea vs. palpebra vs. skin), or clinical cold intolerance across metabolic/mitochondrial disorders.
- **The Synthesis:** Re-integrates pre-molecular oscillating inward-current cold-fiber physiology (Hensel 1951–1953; Pierau 1974; Braun & Schäfer 1980–1986; Huber–Braun 1998; Carr, Pianova & Brock 2002) with CMR1/TRPM8 osmotic gating (Quallo 2015; Yorek 2016) and intracellular $\text{H}^+$ inhibition (Andersson 2004).
- **The Model:** Formulates a 5-variable coupled electro-metabolic relaxation oscillator ($\text{TRPM8} - \text{Na}_\text{v}1.8 - \text{Ca}_\text{v}1.2 - \text{mitochondrial } \text{Ca}^{2+}/\text{H}^+$) demonstrating how tissue osmolarity, mitochondrial membrane potential ($\Delta\Psi_m$), and $\text{H}^+$-sensitive vs. $\text{H}^+$-insensitive agonists switch cold terminals between quiet tone, rhythmic bursting, and pathological $\text{Ca}^{2+}$ overload.
- **Popperian Program:** Proposes five decisive human behavioral falsifications—led by zero-evaporation **Fast-Cooling Gel Submersion ($J_{\text{evap}} = 0$)** and a **2×2 isothermal water-bolus crossover**—with open-source psychrometric PWA and power analysis.

### Section 2 — Introduction: Historical Amnesia After the Cloning of CMR1 (2002)
- Contrast post-2002 single-channel labeled-line dogma against the 90-year classical record (`docs/05_classical_oscillator_to_cmr1_synthesis.md`):
  - **Figure 1 (Historical & Conceptual Architecture):** Timeline & circuit map uniting Hensel (1951–1953), Pierau (1974), Braun & Schäfer (1980–1986), Carr & Brock (2002), McKemy/Peier CMR1/TRPM8 (2002), Andersson (2004), Quallo (2015), Yorek (2016), and Selescu/Zimmermann/Babes (2025).

### Section 3 — Biophysical Synthesis: Thermodynamic Co-Gating & The Electro-Metabolic Oscillator
- Formalize $\Delta G_{\text{total}}(V, T, \text{Osm}_{\text{ext}} - \text{Osm}_{\text{int}}, \gamma, \text{pH}_i)$ (`docs/01_theoretical_framework.md`).
- Explain why **Carr, Pianova & Brock (2002, `PMC2229520`)** found corneal cold nerve terminals chronically depolarized and packed with TTX-resistant $\text{Na}^+$ channels ($\text{Na}_\text{v}1.8$): basal TRPM8 inward current holds the terminal right on the regenerative $\text{Na}_\text{v}1.8$/L-type $\text{Ca}^{2+}$ activation knee.
- Show how **Werth & Thayer (1994)** mitochondrial $\text{Ca}^{2+}$ sequestration → cytosolic $\text{H}^+$ surge + **Andersson et al. (2004)** intracellular $\text{H}^+$ block of cold/icilin CMR1/TRPM8 supplies the missing **metabolic negative-feedback reset** to the classic Huber–Braun (1998) & Olivares (2015) cold-receptor oscillator—immediately explaining **Hensel's 1953 $\text{CO}_2$ suppression** AND **Schäfer, Braun & Isenberg's (1986, `PMC2228855`) biphasic menthol anomaly**!

### Section 4 — Computational Predictions (`analysis/trpm8_oscillator_ode.py`)
- **Figure 2 (`analysis/outputs/trpm8_oscillator_panels.png`):** 8-panel dynamical-systems proof from a single parameter set:
  - Low-stimulus oscillatory bursting vs. saturating uncoupled flood;
  - Requirement for $\text{Na}_\text{v}1.8$ and L-type $\text{Ca}_\text{v}1.2$ exclusively in the physiological oscillatory regime;
  - Collapse into sustained $\text{Ca}^{2+}$ plateau under mitochondrial uncoupling (FCCP) and partial metabolic derangement ($\Psi_m = 0.42$);
  - **Novel Prediction 1:** Isothermal hyperosmolarity (+98 mOsm at warm 33 °C, zero cooling) ignites rhythmic bursting;
  - **Novel Prediction 2:** Menthol ($\text{H}^+$-insensitive) diverges qualitatively from cold/osmolarity ($\text{H}^+$-sensitive) by bypassing the mitochondrial proton brake and locking into a tonic plateau.

### Section 5 — Clinical & Metabolic Implications: Why Bioenergetic Derangement Breaks Cold Perception
- Directly unify three clinical syndromes through Panel E–F of the ODE model:
  1. **Diabetic Small-Fiber Neuropathy:** Impaired terminal mitochondrial $\text{Ca}^{2+}$ buffering → loss of rhythmic $\text{H}^+$ reset → chronic $\text{Ca}^{2+}$ excitotoxicity → progressive corneal/epidermal C-fiber loss and blunted hyperosmolar TRPM8 squinting (Yorek MS et al., *IOVS* 2016, `PMC5113984`).
  2. **Oxaliplatin / Chemotherapy-Induced Cold Allodynia:** Acute DRG mitochondrial stress (Kono et al., 2012, `PMC3343300`) removes the $\text{H}^+$ brake, turning mild room cooling into non-adapting painful discharge.
  3. **Hypothyroidism & Systemic Thermoregulatory Dysfunction:** Reduced mitochondrial proton-motive reserve shifts peripheral TRPM8 burst gain and thermoregulatory shivering thresholds.

### Section 6 — Popperian Falsification Program (Human Behaviors First → Phase 2 Molecular Specialization)
- **Figure 3 (`analysis/outputs/power_comparison.png`):** Mathematical proof that naive outdoor-RH correlations invert indoors in winter, whereas small-$N$ within-subject crossovers ($N = 15\text{–}25$) achieve >80% power.
- Present the 5 falsifiable human protocols (`docs/02_human_protocols_H1_to_H5.md`): **H2-Gel (Fast-Cooling Gel Submersion, $J_{\text{evap}} = 0$ + mannitol control)**, **H1 (Cornea vs. Palpebra vs. Eyelid $V_1$ mapping)**, **H3 (Menthol vs. Cold acidosis & 0.03 Hz Peltier oscillation)**, **H4 (Water Bolus 2×2 separating direct osmolarity vs. TRPV4 vasoconstriction vs. gastric heat sink)**, and **H5 (Open-source Dew-Point VPD PWA)**.
- Close with **Phase 2 Molecular Candidates** explaining tissue-specific tuning across cornea, palpebra, and skin: **TCAF1/TCAF2** (`PMC4284226`, `PMC8387397`), **Asn934 N-glycosylation**, **4TM-TRPM8 splice isoforms** (`PMC5045378`), and mitochondrial $\text{Ca}^{2+}/\text{H}^+$ exchange (**TMBIM5/LETM1**).
