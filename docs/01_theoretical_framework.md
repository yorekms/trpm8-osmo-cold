# Theoretical Framework: Beyond the Cold Labeled Line — TRPM8 as an Osmotically Tuned, Metabolically Clocked Free-Energy Transducer

## 1. The Central Problem
Classical somatosensory doctrine treats cold perception as a two-part simplification:
1. **Static Labeled-Line Coding (Müller's Law):** TRPM8-expressing afferents act as passive thermometers whose firing rate maps monotonically onto the conscious sensation of "cold."
2. **Evaporative Wetness Reductionism:** Because mammals lack insect-like hygroreceptors, sensations of "wetness" or "dry-cold" are explained purely by evaporative heat loss (latent heat of vaporization cooling the skin) combined with tactile mechanoreception.

Four independent empirical observations contradict this simple picture:
1. **Non-evaporative wetness-expulsion motor programs ("Wet Dog Shakes"):** TRPM8 agonists such as icilin elicit full-blown wet-dog shaking across mammals and birds with **zero moisture and zero evaporative cooling** (Selescu, Zimmermann, Babes et al., 2025, *Commun Biol*, PMC12830821).
2. **Adjacent-Territory Qualia Divergence ($V_1$ Trigeminal):** Millimeters apart, corneal TRPM8 activation (via tear evaporation or isothermal hyperosmolar drops) yields *freshness → dryness → blink/squint → gritty sting* (Parra et al., 2010; Quallo et al., 2015, PMC4455064; Yorek MS et al., 2016, PMC5113984; Kovács et al., 2016, PMC4733818), whereas palpebral/cutaneous TRPM8 yields *coolness, damp/wet illusions, or chill/shivering*.
3. **TRPM8 Alone Cannot Oscillate (Cell-Intrinsic Partner Requirement):** Overexpressing TRPM8 in large DRG neurons fails to reconstitute the rhythmic electrical bursting and $[\text{Ca}^{2+}]_i$ oscillations seen in native cold-sensitive neurons at physiological low-stimulus drive (50 nM icilin, $EC_{50} = 360\text{ nM}$; Yorek MS, Usachev Lab 2009–2011).
4. **Cross-Fiber Thermal Coding:** Even human/rodent warmth perception is not a pure labeled line—it requires **TRPM8 silencing** alongside warm-fiber excitation (Paricio-Montesinos et al., 2020, *Neuron*, PMC7272120).

---

## 2. Three-Layer Mechanistic Model

### Layer 1: Thermodynamic Co-Gating of TRPM8 (Temperature × Osmolarity)
TRPM8 open probability $P_{\text{open}}$ is governed by a composite gating free energy integrating thermal enthalpy/entropy, voltage, transmembrane osmotic gradient, and membrane tension/lipid packing ($\text{PI(4,5)P}_2$ microdomains):

$$P_{\text{open}}(V, T, \text{Osm}, \text{pH}_i) = \left[ 1 + \exp\left( \frac{\Delta G_{\text{total}}}{RT} \right) \right]^{-1}$$

$$\Delta G_{\text{total}} = \Delta H(T) - T\Delta S(T) - zFV + \Delta G_{\text{osm}}(\text{Osm}_{\text{ext}} - \text{Osm}_{\text{int}}) + \Delta G_{\text{tens}}(\gamma) + \Delta G_{\text{H}^+}(\text{pH}_i)$$

- **External hyperosmolarity** ($\uparrow \text{Osm}_{\text{ext}}$, as in evaporating tears or dry/hypertonic skin) lowers $\Delta G_{\text{total}}$, shifting $V_{1/2}$ leftward so TRPM8 conducts tonically at physiological surface temperatures (30–35 °C).
- **Internal osmolarity** ($\text{Osm}_{\text{int}}$): Acute systemic water loading lowers $\text{Osm}_{\text{int}}$ transiently. Critically, **two competing physiological predictions** must be falsified head-to-head:
  - **H-Direct (Transmembrane Gradient / Rebound Model):** Lowering $\text{Osm}_{\text{int}}$ increases the relative outward osmotic gradient (or transiently silences basal TRPM8, clearing intracellular $\text{H}^+$ block and triggering sensitized rebound), increasing cutaneous cold/wet sensitivity directly.
  - **H-Rival (Hepatic TRPV4 Osmopressor Reflex):** Portal hypoosmolarity activates hepatic/spinal TRPV4 afferents (Lechner et al., 2011; May & Jordan, 2011, PMC3189227), causing sympathetic **cutaneous vasoconstriction** that physically cools peripheral skin into TRPM8's thermal band.

### Layer 2: The Native Peripheral Oscillator (Bench Data, Usachev Lab 2009–2011)
At physiological sub-saturating drive (mild cooling or 50 nM icilin—unlike saturating 1 µM icilin), native TRPM8+ sensory neurons operate as a **coupled electro-metabolic oscillator**:

1. **Trigger:** Basal/low-stimulus TRPM8 influx ($\text{Na}^+/\text{Ca}^{2+}$) holds resting $V_m$ near −44 mV (right at the subthreshold activation/window zone).
2. **Excitable Amplification:** Depolarization recruits **TTX-sensitive $\text{Na}_\text{v}$** and cold-resistant **$\text{Na}_\text{v}1.8$** (blocked by 1 µM A-803467, hyperpolarizing $V_m$ to −52 mV and abolishing oscillations), driving burst firing up toward depolarization block.
3. **L-Type $\text{Ca}^{2+}$ Amplification:** Spike bursts open **L-type $\text{Ca}^{2+}$ channels** (abolished by 1 µM nimodipine/isradipine; mimicked by Bay K 8644; **retained in $\text{Ca}_\text{v}1.3^{-/-}$ mice**, isolating **$\text{Ca}_\text{v}1.2$**; independent of N-type $\text{Ca}_\text{v}2.2$), surging $[\text{Ca}^{2+}]_i$ from ~100 nM to ~300–500 nM.
4. **Mitochondrial $\text{H}^+$ Brake:** Rapid mitochondrial $\text{Ca}^{2+}$ sequestration (Werth & Thayer, 1994; candidate $\text{Ca}^{2+}/\text{H}^+$ exchanger TMBIM5/LETM1) generates a **phase-locked cytosolic acidification** ($\downarrow\text{pH}_i$, measured directly via simultaneous Fura-2 + F480/440 pH imaging; abolished into a flat plateau by 1 µM FCCP).
5. **Channel Closure & Reset:** Intracellular $\text{H}^+$ directly inhibits **cold- and icilin-gated TRPM8** (Andersson et al., 2004, *J Neurosci*), repolarizing $V_m$ to −45 mV until $\text{Ca}^{2+}$ and $\text{H}^+$ clear and the next burst initiates.

> **Crucial Chemical Prediction (The Menthol vs. Cold/Osmolarity Wedge):** Because intracellular $\text{H}^+$ blocks **cold and icilin** activation of TRPM8 *but spares **menthol*** (Andersson et al., 2004), cold/osmotic stimuli engage the $\text{H}^+$-clocked oscillatory code, whereas menthol bypasses the proton brake and drives a tonic, unclocked percept.

### Layer 3: Peripheral Channel Specialization (Phase 2 Molecular Targets)
Why do cornea, palpebral conjunctiva/eyelid margin, lip mucosa, and hairy/glabrous skin decode TRPM8 differently? Beyond central projection anatomy, **four testable peripheral tuning mechanisms** differentiate TRPM8 behavior across tissues:
1. **Accessory Subunits (TCAF1 / TCAF2):** TCAF1 potentiates TRPM8 trafficking and open probability while TCAF2 inhibits it (Gkika et al., 2015, PMC4284226); humans carry uniquely duplicated *TCAF* loci under opposing selection (Hsieh, Eichler et al., 2021, PMC8387397).
2. **N-Glycosylation (Asn934) & Lipid Rafts:** Deglycosylation shifts TRPM8 $V_{1/2}$ by +60 mV and cold threshold by ~6 °C (Pertusa et al., 2012).
3. **Short Splice Isoforms (4TM-TRPM8 / sM8):** Modulate ER–mitochondria $\text{Ca}^{2+}$ coupling and full-length TRPM8 activity (Bidaux et al., 2016, PMC5045378).
4. **Partner Stoichiometry:** Local density of $\text{Na}_\text{v}1.8$, $\text{Ca}_\text{v}1.2$, $I_{\text{KD}}$ ($\text{K}_\text{v}1.1/1.2$), HCN1/2, and TMBIM5/LETM1 sets whether a terminal behaves as an unconscious tonic tear sensor, a bursting cold/wet detector, or a high-threshold nociceptor.

---

## 3. Summary of Core Popperian Claims Subject to Refutation
- **Claim 1 (Non-Evaporative Osmotic Cold):** Changing external osmolarity or skin hydration alters perceived coldness **even when evaporative mass transfer is strictly zero ($J_{\text{evap}} = 0$)** (tested via **Fast-Cooling Gel Submersion** and **Vapor-Occluded Forearm/Lip** protocols).
- **Claim 2 (Anti-Labeled-Line Regional & Temporal Coding):** Identical isothermal osmotic stimuli elicit qualitatively orthogonal percepts across adjacent $V_1$ territories (cornea vs. palpebra vs. eyelid skin), AND temporal pattern (slow oscillation vs. steady cooling; $\text{H}^+$-sensitive cold vs. $\text{H}^+$-insensitive menthol) alters perceived quality (*cool* vs. *wet* vs. *chill/shudder*) at matched mean thermal intensity.
- **Claim 3 (Systemic Hydration Sensitivity vs. Vasoconstriction):** Acute hypoosmotic water ingestion shifts human cold/shiver perception via a mechanism that can be cleanly decomposed into (a) visceral thermal sink, (b) TRPV4-mediated sympathetic skin cooling, and (c) residual direct osmotic sensitization.
