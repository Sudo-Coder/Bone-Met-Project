# CLAIM_LEDGER.md — frozen, tier-checked claim set (RCC complement-high CLEC_LAM)

Populated **only** from committed results (RESULTS.md, `outputs/tables/*`, commit fb07e13). **No recompute.**
Stats quoted exactly as committed. Claim-ladder tiers per `CLAUDE.md`. Cross-ref: `phase34_verdict_table.csv`,
`CORRECTION_prostate_object.md`, `PRIOR_ART_AND_CLAIMS.md`.

Score scale for all module scores: **AUCell rank-based enrichment (0..1), not z-scored.** SD units are the
**patient-level pseudobulk all-myeloid score SD** (the estimand): CLEC_LAM8 SD=0.1067, complement_C1Q SD=0.1854.

Tiers: descriptive · tumor-enriched · RCC-skewed · predicted-axis · immune-evasion-assoc · prognostic · mechanism(causal).

---

## A. PRIMARY (lead) claims

| # | Claim (manuscript wording) | Exact stat | Model / unit | Tier | Role |
|---|---|---|---|---|---|
| A1 | A CLEC_LAM-like complement–lipid–TREM2 TAM state is present in RCC myeloid and is tumor-enriched. | Δ=+0.174 [0.099, 0.249], p=5.2e-6, FDR=7.8e-6 | all-myeloid pseudobulk CORE score, RCC Tumor vs Benign, mixed +patient RE | tumor-enriched | primary |
| A2 | RCC skews tumor myeloid toward this state more than prostate, relative to a shared benign marrow anchor. | +0.106 [0.023, 0.190], p=0.012, FDR=0.012 = **+1.00 SD [0.22, 1.78]** | cancer_type×condition interaction, all-myeloid pseudobulk, own-benign | RCC-skewed | primary |
| A3 | **The RCC skew is carried specifically by the complement (C1QA/B/C) module** — sharper than the composite and absent for generic TAM burden. | complement_C1Q **+0.210 [0.071, 0.349], p=0.003, FDR=0.014 = +1.13 SD [0.38, 1.88]**; vs CLEC_LAM8 +0.106 (=+1.00 SD); **panTAM control +0.00 [−0.029, 0.025], p=0.89** | same interaction, per module, all-myeloid pseudobulk | RCC-skewed | **primary (lead)** |
| A4 | The skew is a per-cell state intensification of RCC TAM/TIM, not merely more TAMs. | within-TAM+TIM RCC-vs-prostate tumor +0.137 [0.035, 0.240], **p=0.008**; residual after TAM+TIM-fraction adj +0.080 [0.013, 0.147], p=0.020; composition (TAM+TIM fraction) interaction n.s. (p=0.18) | D2/D3, tumor-side within-TAM+TIM & residual | RCC-skewed (state) | primary |
| A5 | C1q in the tumor niche is macrophage (CLEC_LAM)-derived, not tumor-cell-derived: "tumor-niche macrophage autocrine/paracrine complement centered on C1q, with C3→C3AR1/CR3 a related but separately tested axis." | C1QA TAM_CLEC_LAM mean 3.47 (99.8% expr) vs Tumor 0.12 (7.4%); C1QB 3.71 (99.4%) vs 0.15 (6.1%); C1QC 3.72 (100%) vs 0.07 (3.5%); C3 TAM_CLEC_LAM 0.90 (38.8%) vs Tumor 0.41 (31.6%) | patient-level mean log-norm by sender class, tumor niche | predicted-axis (source-resolved niche assoc; **not causal**) | primary |
| A6 | The complement-high state is **associated with an immune-evasive CD8 phenotype**, beating generic TAM burden. **Scope = CD8 cytotoxicity/exhaustion + MHC-II only; Tregs EXCLUDED (see C4).** | **Head-to-head vs panTAM (lead):** reduced CD8 cytotoxicity complement −0.909 (p=0.029) vs panTAM −0.264 (p=0.359); CD8 exhaustion +0.610 (p=0.020) vs +0.292 (p=0.114); MHC-II +1.189 (p=0.019) vs −0.007 (p=0.983). Marginal adjusted: cytotoxicity −0.817 [−1.586, −0.049] p=0.039 (FDR 0.087); exhaustion +0.509 [0.002, 1.016] p=0.049 (FDR 0.118); MHC-II +1.191 [0.299, 2.084] p=0.012 (FDR 0.041) | patient-level (n=24 RCC samples), adjusted for TAM/CD8/malignant fraction + condition | immune-evasion-assoc (**"associated," never "drives/suppresses"**) | primary |
| A7 | The complement/RCC-skew core is prognostic in TCGA-KIRC, independent of generic macrophage burden and the prior Obradovic signature. **Adjustment boundary: adjusted for macrophage burden + Obradovic + total-immune ONLY — NOT yet stage/grade/age/purity-adjusted** (those covariates absent from the staged TCGA file; see Gap D-purity). Manuscript must state this boundary, not imply full adjustment. Tumor purity is a known bulk-TCGA confounder (complement-macrophage signal could partly track immune/stromal content); purity/stage/grade-adjusted model is a one-line future addition. | complement_C1Q OS unadj HR/SD 1.257 [1.091, 1.448] p=0.002; **adjusted (panTAM+Obradovic+immune) HR 1.613 [1.146, 2.269] p=0.006** | continuous Cox per SD, TCGA-KIRC OS (n=532, 175 events) | prognostic (**partial adjustment — see boundary**) | primary |

## B. SUPPORTING claims

| # | Claim | Exact stat | Model / unit | Tier | Role |
|---|---|---|---|---|---|
| B1 | The RCC skew is robust across integration methods and benign references (stability range, not one number). | primary +0.106 (p=0.012); scVI-normalized arm +0.159 [0.069, 0.250] p=5.7e-4; conservative common-space +0.087 [0.011, 0.163] p=0.025; benign-refs p≈0.07 | interaction, 3 integration arms + 3 benign refs | RCC-skewed | supporting |
| B2 | The RCC-skewed fraction of CLEC_LAM-high cells is elevated (state-threshold fraction, not cell-type abundance). | +0.515 [0.326, 0.704], p=9.2e-8, FDR=2.8e-7 | arcsin√(CLEC_LAM-high fraction) interaction, all-myeloid | RCC-skewed | supporting |
| B3 | Complement is the best-supported tumor-specific predicted signaling axis into CLEC_LAM TAM. | ligand tumor-elevation C1QA +1.74 p=0.003 / C1QB +1.81 p=0.026 / C1QC +2.29 p=0.015; CellChat COMPLEMENT C3→C3AR1/CR3 to TAM; NicheNet C1QB rank 12/116 | LR-expr + CellChat + NicheNet integration (tumor niche) | predicted-axis | supporting |
| B4 | TGF-β is a supported tumor-specific co-axis (niche context, not headline). | ligand TGFB1 +0.32 p=4e-4; NicheNet rank 5/116; CellChat TGFb→TAM not significant (receptor low) | LR-expr + NicheNet | predicted-axis | supporting |
| B5 | The prognostic signal concentrates in the complement/RCC-skew core; RCC_skew_CORE is strongest. | RCC_skew_CORE OS adj HR 2.237 [1.358, 3.684] p=0.002 | continuous Cox per SD, TCGA-KIRC OS | prognostic | supporting |

## C. REQUIRED explanatory / boundary / null rows (credibility — do NOT omit)

| # | Claim / boundary | Exact stat | Why stated | Tier |
|---|---|---|---|---|
| C1 | The prognostic signal concentrates in the complement/RCC-skew core and is **diluted in the fuller CLEC_LAM8 composite** (dilution by shared-support genes MERTK/GPNMB/APOE-arm) — so "headline signature not prognostic" cannot be misread. | CLEC_LAM8 OS unadj HR 1.189 [1.034, 1.368] p=0.015 → **adjusted HR 1.543 [0.705, 3.376] p=0.278 (collapses)**; RCC_skew_CORE adjusted HR 2.237 p=0.002 (retains) | prevents mis-reading | prognostic (boundary) |
| C2 | APOE/TREM2 is **not independent** of the Obradovic signature — no independent APOE/TREM2 prognostic claim is made. | APOE_TREM2 OS unadj HR 1.165 p=0.033 → **adjusted HR 1.008 [0.782, 1.300] p=0.948** | collinearity honesty | prognostic (boundary) |
| C3 | The state is **prognostic, not predictive of ICB benefit** (hard ceiling). Braun nivolumab-treated is null; only an exploratory CM-025 arm-interaction. | Braun nivo NULL (complement OS adj p=0.41; ORR p=0.30). Exploratory CM-025 nivo-vs-evero arm×module: CLEC_LAM8 HR 0.533 [0.317, 0.898] p=0.018; RCC_skew 0.573 [0.343, 0.955] p=0.033; complement_C1Q 0.647 [0.389, 1.076] p=0.093 (n=250, underpowered) | **no predictive claim anywhere** | prognostic ceiling; ICB-predictive = exploratory only |
| C4 | Immune-evasion specificity is **CD8 cytotoxicity/exhaustion + MHC-II only — NOT Tregs**; Treg burden is a broad/SPP1 effect. **← this is the exclusion that scopes A6; the two must be stated together.** | complement Treg_fraction +0.154 p=0.665 (n.s.); SPP1_TAM Treg +1.211 p<0.001; MERTK_GPNMB Treg +0.654 p=0.001; head-to-head Treg complement +0.294 (p=0.408) < panTAM +0.404 (p=0.137) | scope the A6 specificity claim | immune-evasion-assoc (boundary) |
| C5 | MERTK/GPNMB is **shared-support across all phases**, not RCC-specific. | Phase-1 per-gene Δ RCC≈prostate (MERTK 0.12/0.14, GPNMB 0.23/0.25); interaction n.s. (p=0.27); Phase-2 GAS6/PROS1→MERTK not tumor-gained; Phase-4 TCGA OS flips protective when adjusted (HR 0.718 p=0.002) | keep in CLEC_LAM8 as shared member, not driver | descriptive (shared) |
| C6 | **No claim may state tumor-cell→TAM C1q signaling** — tumor/epithelial C1q is low (~6% expressing). | C1QA/B/C in Tumor: mean 0.07–0.15, frac 3.5–7.4% (vs TAM_CLEC_LAM ~100%) | direction of arrow is TAM→TAM | (constraint) |
| C7 | Cross-cancer term is power-limited; the **within-TAM state contrast (D2, p=0.008) is the firmer leg** than the benign-anchored interaction. | interaction LOO max p=0.098 (direction-stable); conservative common-space p=0.025; benign-refs p≈0.07; D2 p=0.008 | honest power framing | RCC-skewed (boundary) |
| C8 | We do **not** claim to discover C1Q/APOE/TREM2 TAMs in ccRCC (prior art: Obradovic 2021); novelty is the cross-cancer shared-benign architecture + complement-specific/source-resolved framing. | (qualitative — see `PRIOR_ART_AND_CLAIMS.md`) | novelty discipline | (constraint) |
| C9 | Malignant/epithelial senders are **CA9/PAX8 marker-validated**, NOT CNV-validated (CopyKAT crashed; inferCNV not run). | Tumor CA9/PAX8 panel mean 0.91 vs stroma/immune ~0 (`withstromal_marker_validation.csv`) | wording discipline | descriptive |

## D. Gaps (needed for a claim but NOT computed — list only, do not compute)

- **MHC-I/APM immune-evasion**: underpowered (n=9 tumor samples), all null (complement −1.175 p=0.223). No MHC-I claim; would need more tumor-containing samples.
- **TGF-β response module in Phase 3/4**: TGFB1 established as tumor-specific co-axis in Phase 2 (B4), but a TGF-β *response* score was not carried into the immune-evasion/survival models (ligand vs receptor vs response distinct). Gap — no Phase-3/4 TGF-β claim.
- **TCGA stage/grade/age/purity adjustment (has teeth — sits under headline claim A7)**: not in the staged TCGA
  file → the A7 Cox is adjusted for macrophage+Obradovic+immune ONLY. **Tumor purity is a known bulk-TCGA
  confounder**: a complement-high *macrophage* signal could partly track total immune/stromal content, so A7
  must be worded "adjusted for macrophage burden + Obradovic; not yet purity/stage/grade-adjusted." Any
  "independent of stage/grade/purity" claim is **unsupported today** — one-line future model if a reviewer asks.
- **scCODA / Milo compositional cross-check**: not installable (rpy2); compositional DA is the arcsin-fraction LMM only. No neighborhood-DA claim.
- **PFI/DFI prognosis**: PFI complement adj attenuates (p=0.14); DFI underpowered (15 events, null). Prognostic claim is **OS-specific**; PFI/DFI reported separately, not as support.
- **Spatial / protein / perturbation**: none — all Phase-2 axis claims capped at "predicted"; no causal "drives/contributes to" claim anywhere (see DISCUSSION_VALIDATION.md).

## E. Tier-reach audit (claims flagged if reaching above evidence)
- No claim reaches **mechanism (causal)** — highest tier used is prognostic (A7) and predicted-axis (A5/B3). ✔
- A5 (complement source) is capped at predicted-axis/niche-association, **not** mechanism — direction is descriptive co-localization of production, not a tested causal signal. ✔
- A6 worded "associated," never "drives/suppresses." ✔ Scope explicitly linked to C4 (Tregs excluded). ✔
- A7 tagged **partial adjustment** (macrophage+Obradovic+immune; not purity/stage/grade) — must not imply full adjustment. ✔
- C3 enforces "prognostic, not predictive" — no predictive claim exists in A/B. ✔
- Any stat not in a committed table above = **"unsupported — cut or soften"** (none currently in A/B; the D-gaps are pre-emptively excluded from claims).
