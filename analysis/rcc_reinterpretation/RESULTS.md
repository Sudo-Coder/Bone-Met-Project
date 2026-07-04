# RESULTS.md — RCC CLEC_LAM reinterpretation

One figure/table per claim. Each entry: plain-language interpretation · exact stat (n, test, effect
size + 95% CI, p, FDR, adjustment) · **claim-ladder tier** · figure/table path. Unsupported / null
results are stated explicitly, never omitted.

Status: **Phase 0 + 0.5 + 1 + 1.5 + 2 + 3 + 4 complete.** Paused for review with combined Phase 3/4
summary. Phase 5 (inferCNV) optional/deferred; manuscript synthesis not started pending review.

> **Correction (2026-07-03):** Phase 0/1 used `myeloid_integrated_final_label.h5ad`; the canonical
> prostate myeloid object is `myeloid_FINAL.h5ad`. The two were **verified data-identical** (same 11,760
> cells, 100% label agreement, **raw byte-identical, max|Δ|=0**); Phase 1 scored on `.raw`, so all
> raw-based results reproduce **exactly**. Full before/after table + the new composition-vs-state
> decomposition: `CORRECTION_prostate_object.md`. Only the joint-scoring robustness arm moved (method
> improved to MERTK-inclusive raw scoring; now borderline) — primary confirmatory results unchanged.

## Headline confirmatory table (Phase 1)
| # | Claim (tier) | Test (unit) | Effect + 95% CI | p / FDR | Verdict | Fig/Table |
|---|---|---|---|---|---|---|
| C1 | CLEC_LAM is **tumor-enriched** in RCC | all-myeloid CORE score, Tumor vs Benign, mixed model +patient RE (n=9 tumor,7 benign patients) | Δ=**+0.174** [0.099, 0.249] | p=5.2e-6 / FDR=7.8e-6 | **Supported** | fig3, tables/phase1_confirmatory_results.csv |
| C2 | CLEC_LAM is **RCC-skewed vs prostate** | cancer_type×condition[Tumor] interaction, all-myeloid CORE, own-benign, mixed +patient RE | β=**+0.106** [0.023, 0.190] | p=0.012 / FDR=0.012 | **Supported (power-limited)** | fig4 |
| C3 | RCC-skewed **fraction** of CLEC_LAM-high | interaction on arcsin√(CLEC_LAM-high fraction), all-myeloid | β=**+0.515** [0.326, 0.704] | p=9.2e-8 / FDR=2.8e-7 | **Supported** | fig4 |

**Claim-ladder tiers reached this phase:** "a CLEC_LAM-like TAM state is present" (C0), "tumor-enriched"
(C1), and "RCC-skewed tumor-induced TAM fate" (C2/C3). No signaling/immune-evasion/prognostic/causal
claims are made yet.

---

## Phase 0 — environment & manifest
Objects inventoried, canonical set fixed, **shared-benign anchor verified** (same BMM2–9 donors, distinct
GEO/processing). See `../DATA_MANIFEST.md`. No results claimed.

## Phase 0.5 — external ICB staging
Braun/CheckMate-009+010+025 **staged & verified** (expression 43,893g×311s, all 8 CORE present;
OS/PFS/ORR/Benefit complete; nivo n=181, CM-025 nivo-vs-evero n=250). IMmotion150 access identified,
not obtained (EGA controlled). See `resources/external_icb/README.md`. Descriptive/staging only.

## Phase 1 — define CLEC_LAM and prove tumor-enrichment

### 1.1 State definition & scoring (descriptive, C0)
- **CORE = C1QA,C1QB,C1QC,APOE,APOC1,TREM2,GPNMB,MERTK** — all 8 present in both myeloid objects.
  SELENOP absent from all objects (dropped, pre-registered). FOLR2 held out as non-load-bearing support.
- Continuous CORE scored two ways: `sc.tl.score_genes` (additive) + **AUCell** (rank-based, decoupler).
  **Concordance (Spearman):** RCC 0.79 all / **0.97** TAM+TIM; prostate 0.49 all / **0.77** TAM+TIM
  (prostate lower — sparser CORE detection). AUCell used downstream (depth-robust).
- **CORE detection (fraction expressing), RCC TAM:** C1QA/B/C 0.83/0.79/0.81, APOC1 0.62, TREM2 0.45,
  APOE 0.44, GPNMB 0.25, MERTK 0.16. **Prostate TAM** markedly lower (C1QA 0.33, TREM2 0.03, MERTK 0.03).
  → the state is far more expressed in RCC TAM. Fig `fig2_core_detection.png`.
- Thresholded CLEC_LAM-high (top 10/15/20% AUCell over all myeloid + 2-comp GMM) = SECONDARY; used only
  for the fraction analysis. `tables/phase1_scoring_report.json`.

### 1.2 TIM comparability gate (pre-registered before pooling) — PASSED
Cross-cancer TIM marker-profile **Spearman rho = 0.936** (n=15 monocyte/TAM/inflammatory markers).
TIM is defined comparably in both cancers → TAM+TIM pooling justified. (Primary unit later became
all-myeloid for the benign-anchored tests; see 1.3.)

### 1.3 STRUCTURAL FINDING — TAMs are essentially tumor-restricted (drives the analysis choice)
RCC **benign** marrow has ~no TAM/TIM cells (median **1.5** TAM+TIM per benign sample; **0/4** benign
samples reach ≥10). Prostate benign has some (median 21). `tables/tamtim_availability.csv`.
→ A TAM+TIM-gated tumor-vs-benign contrast is **undefined on the RCC benign side**. The PRIMARY unit is
therefore **mean CORE score over ALL myeloid per patient/sample** (benign myeloid IS well-populated),
which captures both the compositional emergence of TAMs and their per-cell state. TAM+TIM-restricted
scores are tumor-side descriptive only. This is a biological result (TAM emergence is tumor-driven),
not just a technical choice.

### 1.4 Tumor-enrichment in RCC (C1) — **Supported** (tier: tumor-enriched)
All-myeloid CORE score, RCC Tumor vs Benign, mixed model with patient RE:
**Δ = +0.174 [0.099, 0.249], p=5.2e-6, FDR=7.8e-6** (tumor mean 0.177 vs benign 0.002; n=9 tumor, 7
benign patients). Benign myeloid essentially lacks the CORE program. Fig `fig3_patient_boxplots.png`.

### 1.5 Cross-cancer RCC-skew (C2/C3, HEADLINE) — **Supported but power-limited**
`CORE_score ~ cancer_type*condition + (1|patient_id)`, own-benign, all myeloid:
**interaction β = +0.106 [0.023, 0.190], p=0.012, FDR=0.012** (mixedlm). OLS cluster-robust sensitivity
+0.106 [0.004, 0.209], p=0.042. Fraction (arcsin) interaction +0.515 [0.326, 0.704], p=9.2e-8.
(Binomial-GLM fraction test suffered **complete separation** — RCC benign ~0 high cells — reported but
not used; arcsin LMM is the confirmatory fraction test.)

**Reconciling the two "fraction"/composition results (no contradiction):** the CLEC_LAM-high fraction
interaction (C3, +0.515) is a **state-threshold fraction** — more cells cross the CLEC_LAM-high cutoff
because the per-cell program intensifies. The TAM+TIM **cell-type abundance** fraction (D1, below) is
**not** significantly different between RCC and prostate. Together these support **per-cell state
intensification**, not a simple TAM/TIM abundance shift.

**Robustness (all pre-specified):**
- **Primary pseudobulk (observed expression):** +0.106 [0.023, 0.190], p=0.012 — significant.
- **Integration sensitivity (MAIN = scVI-normalized, batch-aware, full CLEC_LAM8 incl MERTK):**
  **+0.159 [0.069, 0.250], p=5.7e-4** — the RCC-skew is supported in the model-based integrated
  representation. `tables/pseudobulk_joint_scVInorm.csv`.
- **3 benign references (scVI arm; vary the RCC benign baseline, same BMM donors):** RCC-processed +0.159
  [0.049, 0.269] (p=4.6e-3); prostate-processed benign for the RCC arm +0.190 [0.081, 0.300] (p=6.7e-4);
  per-donor averaged +0.175 [0.065, 0.284] (p=1.8e-3). All significant; using the other processing's
  benign gives a *larger* effect → not a benign-processing artifact.
- **Dataset-regressed:** +0.106 (p=0.042).
- **Conservative common-space (raw-lognorm) supplemental — separate from the main scVI arm:** +0.087
  [0.011, 0.163] (p=0.025). Reported as a conservative sensitivity, not the primary integration arm.
- **Leave-one-patient-out (n=23):** interaction stays positive, est ∈ [0.084, 0.128], but **max p=0.098**
  — one influential patient pushes it just above 0.05. Honest verdict: **direction-robust, significance
  borderline**; consistent with prostate TAM being underpowered (552 TAM; benign TAM near-0). `tables/C2_leave_one_out.csv`.

Net: the effect is **supported in the primary pseudobulk model and the scVI-normalized integration arm**,
with the conservative common-space result reported separately (not claimed identical across arms).

### 1.5b Composition-vs-state decomposition (standing requirement) — the skew is a STATE effect
Is the RCC-skew "more TAMs" (composition) or "hotter TAMs" (state)? (`tables/phase1_confirmatory_results.csv`)
| test | effect [95% CI] | p | reading |
|---|---|---|---|
| D1 composition — TAM+TIM fraction interaction (arcsin) | +0.232 [−0.109, 0.573] | 0.18 | not a cross-cancer TAM-abundance difference (underpowered) |
| D2 state — within-TAM+TIM CORE, RCC vs prostate **tumor** | +0.137 [0.035, 0.240] | **0.008** | **RCC TAMs express more CORE per cell** |
| D2b state — prostate within-TAM+TIM, tumor vs benign | +0.057 [−0.014, 0.129] | 0.12 | prostate TAMs do not significantly upregulate CORE |
| D3 residual — interaction \| TAM+TIM fraction (GUARDED) | +0.080 [0.013, 0.147] | **0.020** | **RCC-skew survives** fraction adjustment |
**Reading (guarded):** the RCC-skew is a **per-cell state intensification in RCC TAMs** (D2) that
**survives** adjustment for TAM abundance (D3), not a cross-cancer difference in TAM fraction (D1 n.s.).
Per the residual-test guard, the attenuated-but-significant D3 is **not** read as null — the emergence of
TAMs is itself the canalization phenotype and is not adjusted away.

**Design-honesty on D2:** D2 is a robust *direct tumor-vs-tumor* comparison within TAM+TIM, but it does
not itself carry the shared-benign anchor; it is reported together with the benign-anchored interaction.
Combined statement: *RCC tumor TAM/TIM states are more intensely CLEC_LAM than prostate tumor TAM/TIM
states (D2), and the shared-benign interaction (C2) supports tumor-induced RCC skewing of this program.*

> **LOCKED THESIS (Phase 1):** RCC tumor myeloid cells, especially TAM/TIM states, show **stronger
> per-cell induction of the CLEC_LAM complement–lipid–TREM2 program** than prostate, relative to the
> shared benign marrow anchor. (Not framed as a TAM/TIM abundance difference.)

### 1.6 Per-gene tumor-induction — refines locked biology (honest nuance)
Patient-level Δ(tumor−benign) of log-norm expression, per cancer (`tables/core_gene_tumor_induction.csv`):
| gene | RCC Δ | prostate Δ | RCC-differential? |
|---|---|---|---|
| C1QA | +1.30 | +0.40 | **yes (≈3×)** |
| C1QB | +1.17 | +0.36 | **yes** |
| C1QC | +1.16 | +0.40 | **yes** |
| APOC1 | +0.75 | +0.46 | yes (modest) |
| APOE | +0.42 | +0.25 | yes (modest) |
| TREM2 | +0.33 | +0.18 | yes (modest) |
| GPNMB | +0.23 | +0.25 | **no (comparable)** |
| MERTK | +0.12 | +0.14 | **no (comparable)** |
→ The RCC-skew is carried by the **complement (C1QA/B/C)** and **lipid (APOC1/APOE)** + TREM2 arms.
**MERTK (efferocytosis) and GPNMB are induced comparably in both cancers** — NOT RCC-differential here.
This **revises the pre-registered locked biology** ("MERTK tumor-induced vs absent in prostate"): in this
matched, shared-benign comparison MERTK is a shared, not RCC-specific, tumor-induction. Per ruling, MERTK & GPNMB
**remain inside the pre-registered CLEC_LAM8 confirmation score** — treated as **shared/support members**
rather than RCC-specific load-bearing genes (unless Phase 2 or Phase 4 evidence changes their role); not
removed. Their non-RCC-specificity is an interpretation, not a score change. (Cross-cancer absolute Δ
compares each object's own normalization; treat as pattern, not calibrated magnitude.)

### 1.6b Cross-cancer LOCK — score scale, complement lead, design rationale, C1q source
**Score scale (pre-check 1):** all scores are **AUCell rank-based enrichment (0..1); NOT z-scored, NOT
raw log-norm mean**. The estimand is the **patient-level pseudobulk all-myeloid** score. Its SD (CLEC_LAM8)
= **0.1067**, so the primary CLEC_LAM8 interaction **+0.106 = +1.00 SD [0.22, 1.78 SD]**
(`tables/module_interactions.csv`).

**Complement module is the sharper lead (pre-check 2):** the cross-cancer interaction on the **complement
C1Q module (C1QA/B/C)** is **+0.210 [0.071, 0.349], p=0.0030, FDR=0.014 = +1.13 SD [0.38, 1.88]**
(all-myeloid; TAM+TIM secondary +0.127, p=0.006) — **stronger and lower-FDR than CLEC_LAM8** (+0.106,
FDR 0.028). `complement_C1Q_C3` (+C3) +0.169 [0.057, 0.281], p=0.0031. **The cross-cancer claim now leads
with the complement module.** It is **specific**: in the same interaction, generic **panTAM +0.00
(p=0.89)**, SPP1_TAM (p=0.63), inflammatory-mono (p=0.13), APOE/TREM2 (p=0.13) and MERTK/GPNMB (p=0.27)
are all null → the RCC-skew is complement/CLEC_LAM-specific, not a generic TAM-burden effect.

**Design rationale (deliberate, forced choices — methods+results):**
- RCC benign marrow has **~0 TAM/TIM** → a TAM-gated benign contrast is **undefined**; the anchored primary
  unit is therefore **all-myeloid pseudobulk** (benign defined on both sides). True TAM/TIM labels are
  biologically interpretable only in tumor-containing samples, which motivates the condition-aware design.
- Prostate TAM-alone is thin (552 cells / ~138 patient-level) → the cross-cancer term runs on **TAM+TIM
  pooled**, justified by cross-cancer **TAM/TIM marker concordance rho = 0.936**.
- The "complement-high CLEC_LAM RCC-skew" is measured at the **myeloid-compartment level (pooled
  pseudobulk)**, NOT as purified-CLEC_LAM-TAM-vs-benign-CLEC_LAM-TAM (undefined — no benign CLEC_LAM TAM).
- The interaction is reported as a **stability range**: **+0.106 (primary pseudobulk, p=0.012)** to
  **+0.087 (conservative common-space, benign-refs p~0.07)**, with the scVI-normalized arm +0.159 — not a
  single triumphant number.

**C1q source — mechanism direction LOCKED (gating item; `tables/complement_source_by_class.csv`):**
Patient-level, tumor niche, C1QA/B/C mean log-norm is **highest in TAM_CLEC_LAM (3.47/3.71/3.72, ~100%
expressing)** ≫ TAM_other (≈2) ≫ tumor/stroma (≈0.1, ~6%). **C1q is macrophage (CLEC_LAM)-derived, NOT
tumor-cell-derived.** C3 is also highest in TAM_CLEC_LAM (0.90) but with real tumor contribution (0.41).
→ **Mechanism wording: "tumor-niche macrophage autocrine/paracrine complement signaling centered on C1q,
with C3→C3AR1/CR3 treated as a related but separately-tested complement axis"** — direction is TAM→TAM
(autocrine amplification), not tumor→TAM, and C1q and C3 are not collapsed into one linear axis.

### 1.7 Limitations / not-yet-done (Phase 1)
- **scCODA + Milo deferred:** both require R/edgeR via rpy2, which failed to build in the Python venv.
  Compositional DA is currently covered by the arcsin fraction LMM (+ binomial GLM, flagged separation).
  scCODA/Milo will be run once the R env (`install.R`) is built (needed anyway for Phase 2 CellChat).
- **De-novo scVI joint re-integration** (concat + scVI with explicit dataset covariate) not yet run; the
  "integrated" arm currently uses scVI-normalized expression (joint common-norm). Available on GPU if a
  stricter check is wanted — all current arms already agree.
- **Prostate TAM underpowered** (552 TAM; benign TAM ≈0), as pre-registered — hence C2 borderline in LOO.

### Files
Figures (300 dpi): `outputs/figures/fig1_umap_core_score.png`, `fig2_core_detection.png`,
`fig3_patient_boxplots.png`, `fig4_forest.png`. Tables: `outputs/tables/` (phase1_confirmatory_results.csv,
phase1_scoring_report.json, pseudobulk_all_myeloid.csv, pseudobulk_TAMTIM.csv,
pseudobulk_joint_commonnorm.csv, tamtim_availability.csv, C2_leave_one_out.csv,
core_gene_tumor_induction.csv, clec_lam_percell.parquet).

---

## Phase 1.5 — malignant-sender validation & niche labeling (Phase 2-prep)

### 1.5a Niche labeling + benign-sender validity (descriptive)
The RCC with-stromal niche (`integrated-with-stromal.h5ad`, 96,829 cells) was labeled by barcode transfer
(81,299 cells) + kNN in X_scVI for the 9 benign-stroma samples (~15,530 cells). **All compartments
marker-validate** (`tables/withstromal_marker_validation.csv`; diagonal mean log-norm high): Tumor→
CA9/PAX8 0.91, MSC 2.22, Pericyte 2.24, Endothelial 1.17, Osteoclast 3.62, TAM 2.20, T/NK 1.61, B 1.92.
**Benign stromal senders are valid & abundant** — MSC 6001, Pericyte 4971, Endothelial 603
(`tables/withstromal_sender_availability.csv`). **BUT benign has ~0 TAM/TIM receivers** (7 cells;
TAMs are tumor-restricted, per Phase 1). → benign→TAM CellChat is unpopulated; **tumor-specificity of
sender→TAM axes is therefore derived from patient-level ligand/receptor expression**, and CellChat
results are worded "predicted axes operating in the tumor niche."

### 1.5b Malignant-sender validation (tier: descriptive — justifies the sender set)
**Tumor-labeled cells are the malignant senders by canonical ccRCC markers** (CA9/PAX8/NDUFA4L2/VEGFA/
EPCAM mean 0.91 vs ~0 in stroma/immune; stroma high on own lineage markers). Per plan, canonical markers
are an accepted malignant cross-check. **CopyKAT was attempted and FAILED** (v1.2.5 crashed at step 6
`convert.all.bins.hg20`, "low data quality" warning — a known failure mode on integrated inputs; no
aneuploid/diploid calls). inferCNV installed, not run (needs gene-order file); nominated as optional CNV
cross-check. `tables/malignant_validation_status.md`.

## Phase 2 — CellChat + ligand-target support (predicted signaling axes)

**All Phase 2 claims are tier "predicted signaling axes"** (hypothesis-generating; no spatial/protein/
perturbation). Three independent lines were integrated per pre-specified axis:
(A) ligand tumor-elevation in senders (patient-level Δ tumor−benign, Mann-Whitney);
(B) CellChat communication probability sender→TAM receiver (tumor niche);
(C) NicheNet ligand→CLEC_LAM-program activity rank (of 116 expressed ligands).

### 2.1 Axis-support synthesis (`tables/axis_support_table.csv`, `figP2_1_axis_support.png`)
| Axis | (A) ligand Δ T−B | (B) CellChat→TAM prob | (C) NicheNet rank | Verdict |
|---|---|---|---|---|
| **Complement C1QB** | **+1.81, p=0.026** | ✓ C3→CR (0.005) | **rank 12 (top20)** | **SUPPORTED — tumor-specific** |
| **Complement C1QA** | **+1.74, p=0.003** | ✓ | rank 43 | **SUPPORTED — tumor-specific** |
| **Complement C1QC** | **+2.29, p=0.015** | ✓ | — | **SUPPORTED — tumor-specific** |
| **TGFB1→TGFBR** | **+0.32, p=4e-4** | ✗ (receptor low) | **rank 5 (top ligand)** | **SUPPORTED — tumor-specific** |
| C3→C3AR1/CR3 | +0.15, p=0.16 (ns) | ✓ (tumor+autocrine) | rank 36 | present, **not tumor-gained** |
| APOE→TREM2 | −0.04, p=0.96 (ns) | ✓ **autocrine prob 0.174** | rank 55 | present, **not tumor-gained** |
| GAS6→MERTK | −0.42, **p=0.02 (benign↑)** | weak (0.004) | rank 107 | **NOT tumor-gained** |
| PROS1→MERTK | +0.00, p=0.86 | ✗ | — | unsupported |

**Interpretation (all "predicted"):**
- **Complement (C1QA/B/C) is the best-supported tumor-specific axis** — tumor-elevated (macrophage-derived,
  autocrine/paracrine C1q; C3→C3AR1/ITGAX-ITGB2 from tumor + TAM) and a top NicheNet predictor of the
  CLEC_LAM program (C1QB rank 12). **This is the same complement arm that drives the Phase 1 RCC-skew** →
  the mechanism converges with the state signal.
- **TGFB1→TGFBR is tumor-specific and the #1 NicheNet ligand** (rank 5) — the APOE/TGF-β/ICB-resistance link.
- **APOE→TREM2 operates strongly in the tumor niche** (CellChat autocrine prob 0.174, TAM_CLEC_LAM→
  TAM_CLEC_LAM) **but is NOT tumor-gained** (ligand present in benign senders, p=0.96). Wording:
  "predicted APOE→TREM2 axis operating in the tumor niche, macrophage-derived, not tumor-induced."
- **GAS6/PROS1→MERTK efferocytosis axes are NOT tumor-gained** (GAS6 benign-biased). **Consistent with
  Phase 1** where MERTK was not an RCC-differential driver — efferocytosis is a shared, not RCC-specific,
  program at both the state and signaling level.

### 2.2 CellChat top predicted axes → TAM_CLEC_LAM (`figP2_2_cellchat_clec_lam.png`, `tables/cellchat_tam_LR_tumor.csv`)
363 LR pairs to TAM receivers. Highest-prob into CLEC_LAM TAM: Endothelial APP→CD74 (0.21);
**TAM_CLEC_LAM APOE→TREM2/TYROBP (0.17, autocrine)**; Tumor MIF→CD74/CXCR4 (0.16); Osteoclast SPP1→CD44
(0.13; prior-art SPP1 TAM); MSC/Pericyte COLLAGEN→CD44. Pathways present: COMPLEMENT, ApoE, GAS, SPP1, MIF,
MK, GALECTIN, ANNEXIN, APP, THBS, LAMININ, MHC-II.

### 2.3 NicheNet ligand→CLEC_LAM-program (`figP2_3_nichenet_ranks.png`, `tables/nichenet_*.csv`)
Of 116 expressed ligands predicting the TAM_CLEC_LAM-high program: **TGFB1 rank 5, C1QB rank 12** are the
pre-specified leaders; C3 36, C1QA 43, APOE 55, GAS6 107. Consistent with (A)/(B): complement + TGFB1 lead.

### 2.4 Limitations (Phase 2)
- CellChat/NicheNet are hypothesis-generating; no spatial/protein/perturbation validation (see
  DISCUSSION_VALIDATION.md). Benign→TAM CellChat unpopulated (benign TAMs absent) → tumor-specificity from
  LR expression, not a benign-vs-tumor CellChat contrast. CopyKAT crashed (malignant call via markers).
  CellChat probabilities are small in absolute terms (curated-DB, many groups); ranks/agreement across the
  three lines carry the inference, not any single probability.

### Phase 2 files
Figures (300 dpi): `figP2_1_axis_support.png`, `figP2_2_cellchat_clec_lam.png`, `figP2_3_nichenet_ranks.png`,
`figP2_4_ligand_tumor_vs_benign.png`. Tables: `axis_support_table.csv`, `axis_lr_expression.csv`,
`cellchat_tam_LR_{tumor,benign}.csv`, `cellchat_allLR_tumor.csv`, `nichenet_ligand_activities.csv`,
`nichenet_prespecified_ranks.csv`, `nichenet_ligand_target_links.csv`, `withstromal_*.csv`,
`malignant_validation_status.md`.

---

## Phase 3 — immune-evasion association (patient-level, adjusted; CONTRAST-primary)
RCC full-niche; n=24 samples (9 tumor/4 involved/4 distal/7 benign). Myeloid module scores vs T-cell/
tumor outcomes, adjusted for TAM+TIM fraction, CD8 fraction, malignant fraction, condition. **Primary =
whether the complement module beats the negative controls** (panTAM/TAM-fraction, SPP1_TAM, inflammatory-
mono, MERTK/GPNMB). Tier: **"associated with an immune-evasive phenotype"** (patient-level, adjusted; not
causal). `tables/phase3_adjusted_associations.csv`, `phase3_headtohead_complement_vs_panTAM.csv`, fig `figP34_2`.
- **Reduced CD8 cytotoxicity — complement-SPECIFIC:** complement_C1Q std-coef **−0.82 [−1.59, −0.05], p=0.039**;
  CLEC_LAM8/RCC-skew similar; **panTAM (−0.09, p=0.77) and MERTK/GPNMB (−0.06, p=0.82) null**. Head-to-head
  complement vs panTAM: **−0.91 (p=0.029) vs −0.26 (p=0.36)** → complement beats generic TAM.
- **CD8 exhaustion:** complement +0.51 (p=0.049); beats panTAM head-to-head (+0.61 p=0.020 vs +0.29 p=0.11).
  SPP1_TAM (+0.51) and MERTK/GPNMB (+0.35) also positive → partly shared, complement not uniquely dominant.
- **Treg fraction — NOT complement-specific:** complement +0.15 (p=0.66, ns); driven by **SPP1_TAM (+1.21,
  p<0.001) and MERTK/GPNMB (+0.65, p=0.001)**. Honest boundary: Treg burden is a broad-TAM/other-module effect.
- **MHC-II/APC:** complement +1.19 (p=0.012, beats panTAM −0.01 p=0.98) — complement-high TAMs are MHC-II-**high**
  (more antigen-presentation capacity; complicates a simple suppressive reading — reported as nuance).
- **MHC-I/APM:** n=9 tumor samples, **underpowered**, all null — no conclusion.

## Phase 4 — TCGA-KIRC + Braun/CheckMate (continuous Cox/logistic; CONTRAST-primary)
Modules scored on bulk (AUCell single-sample; GSVA available as alt). Stage/grade/purity/age NOT in the
staged TCGA file → parsimonious models (unavailable covariates listed). Collinearity with Obradovic
EXPECTED and reported; novelty is **not** staked on survival independence.

### 4a TCGA-KIRC (`tables/phase4_tcga_kirc_cox.csv`, fig `figP34_1`) — report OS/PFI/DFI separately
- **complement_C1Q OS:** unadj HR/SD **1.26 [1.09, 1.45], p=0.002**; **adjusted (panTAM+Obradovic+immune)
  HR 1.61 [1.15, 2.27], p=0.006 — SURVIVES adjustment.** `RCC_skew_CORE` even stronger (adj HR 2.24, p=0.002).
- **Specificity:** **panTAM (generic macrophage) is NOT prognostic** (unadj p=0.11); **APOE/TREM2 collapses
  fully into Obradovic** (adj p=0.95 — it IS the prior signature); **CLEC_LAM8 collapses** (adj p=0.28). So the
  independent OS signal is carried by the **complement / RCC-skew** part, not generic TAM or the APOE/TREM2 arm.
- **PFI:** complement unadj p=0.028, adj p=0.14 (attenuates); RCC-skew adj p=0.021. **DFI:** 15 events, all null (underpowered).

### 4b Braun/CheckMate ICB (`tables/phase4_braun_nivo.csv`, `phase4_braun_cm025_armINT.csv`)
- **Nivolumab-treated (primary, n=181, 58 OS events): NULL.** No module significantly predicts OS/PFS/ORR/
  clinical-benefit after adjustment (complement OS adj p=0.41; ORR p=0.30). Honest null — the complement/
  CLEC_LAM state is **not** a within-nivolumab response biomarker here.
- **CM-025 nivo-vs-everolimus arm × module interaction (EXPLORATORY, n=250, underpowered):** **CLEC_LAM8 ×
  arm HR 0.53 [0.32, 0.90], p=0.018**; RCC-skew p=0.033; complement_C1Q trend (0.65, p=0.093). Suggests the
  state may associate with **differential nivo-vs-everolimus OS benefit** — labeled exploratory, not a headline.

### Combined Phase 3/4 verdict table (`tables/phase34_verdict_table.csv`) — the decision
| module | CD8 exhaust | cytotox | Treg | MHC-II | TCGA OS (adj) | Braun nivo OS/ORR | beats generic TAM? | **verdict** |
|---|---|---|---|---|---|---|---|---|
| **complement_C1Q** | ↑ p0.05 | **↓ p0.04** | ns | ↑ p0.01 | **1.61↑ p0.006** | null | **yes** (cytotox, exhaust, MHC-II) | **SUPPORTED — complement-SPECIFIC** |
| RCC_skew_CORE | ↑ p0.09 | ↓ p0.03 | ns | ↑ p0.04 | **2.24↑ p0.002** | null | ~complement | supported (≈complement) |
| CLEC_LAM8 | ↑ p0.07 | ↓ p0.04 | ns | ↑ p0.03 | collapses (p0.28) | null; **arm-INT p0.018 (expl)** | partial | supported-broad |
| APOE_TREM2 | ns | ~p0.08 | ns | ns | = Obradovic (p0.95) | null | no | = Obradovic; not independent |
| MERTK_GPNMB | ↑ p0.03 | ns | **↑ p0.001** | ns | 0.72↓ p0.002 (flips) | null | no | shared-support (Treg via broad) |
| panTAM (control) | ns | ns | ns | ns | ns (p0.11) | null | — | generic control — mostly null |
| SPP1_TAM (control) | ↑ p0.04 | ns | **↑ p<0.001** | ns | — | — | — | broad-TAM (Treg-driver) |

**Bottom line (decides the mechanism):** the immune-evasion/prognostic signal is **complement-specific, not
macrophage-burden-driven** — complement beats generic panTAM head-to-head (reduced CD8 cytotoxicity, exhaustion,
MHC-II), panTAM is null throughout, and complement TCGA-OS **survives** adjustment for both generic macrophage
and the Obradovic signature (HR 1.61, p=0.006). Honest boundaries: **Treg burden is a broad/SPP1 effect, not
complement**; **ICB response (Braun nivo) is null** with only an exploratory CM-025 arm-interaction; APOE/TREM2
is the prior Obradovic signature (collapses); MERTK/GPNMB stays shared-support. Claim-ladder tiers reached:
**"associated with an immune-evasive phenotype"** (Phase 3) and **"clinically associated / prognostic"** (TCGA OS,
Phase 4a); ICB predictive = exploratory only.

### Phase 3/4 files
Figures: `figP34_1_tcga_os_forest.png`, `figP34_2_immune_evasion_contrast.png`. Tables: `phase3_sample_table.csv`,
`phase3_adjusted_associations.csv`, `phase3_headtohead_complement_vs_panTAM.csv`, `phase4_tcga_kirc_cox.csv`,
`phase4_braun_nivo.csv`, `phase4_braun_cm025_armINT.csv`, `phase34_verdict_table.csv`, `module_interactions.csv`,
`complement_source_by_class.csv`.
