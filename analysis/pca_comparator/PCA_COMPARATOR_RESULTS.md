# PCA_COMPARATOR_RESULTS.md — prostate (PCa) mirror of the RCC complement-high CLEC_LAM analysis

Locked comparator to the RCC study (`../rcc_reinterpretation/`). Same patient-level pseudobulk
framework; scripts in `scripts/` (06/31/32); tables in `outputs/tables/`.

---

## PCa comparator (locked mirror analysis) — is the RCC complement program present / immune-associated / clinical in prostate?
Same patient-level pseudobulk framework applied to prostate myeloid (`myeloid_FINAL.h5ad`) + prostate
full-niche (`integrated_with_kfoury_labels.h5ad`, n=32 samples). Modules incl. ATF3/NF-κB.
Tables: `pca_vs_rcc_module_contrast.csv`, `pca_within_tumor_localization.csv`, `pca_phase3_adjusted_associations.csv`,
`pca_phase3_headtohead.csv`, `pca_vs_rcc_contrast_table.csv`.

### PCa-vs-RCC module contrast (the deliverable)
| module | RCC T-vs-B (SD) | PCa T-vs-B (SD) | interaction (RCC-skew) | PCa TAM/TIM localization | RCC cytotox assoc | PCa cytotox assoc | verdict |
|---|---|---|---|---|---|---|---|
| **complement_C1Q** | 1.41 | **1.00 (p=0.002)** | **+0.210 (p=0.003)** | **Δ+0.11, p=0.016 (localizes)** | **−0.82 (p=0.04)** | −0.02 (p=0.87) | **present+TAM-localized in PCa, RCC-SKEWED, RCC-functionally-specific** |
| CLEC_LAM8 | 1.37 | 1.03 (p=0.001) | +0.106 (p=0.012) | Δ+0.07, p=0.016 | −0.77 (p=0.04) | +0.01 (p=0.95) | same |
| RCC_skew_CORE | 1.37 | 1.03 (p<0.001) | +0.138 (p=0.010) | — | −0.78 (p=0.03) | +0.01 (p=0.96) | same |
| ATF3_NFkB | 1.54 | **1.52 (p<0.001)** | +0.079 (p=0.033, small) | Δ+0.10, p=0.016 | — | −0.31 (p=0.21) | **CONSERVED core (both strong)** |
| SPP1_TAM | 1.33 | 1.17 (p<0.001) | +0.020 (p=0.63, ns) | Δ+0.06, p=0.30 | +0.61 (p=0.12) | +0.04 (p=0.81) | shared (not RCC-skewed) |
| APOE_TREM2 | 1.27 | 0.92 (p<0.001) | +0.050 (p=0.13, ns) | Δ+0.04, p=0.016 | −0.55 (p=0.08) | +0.02 (p=0.87) | shared |
| MERTK_GPNMB | 1.21 | 0.96 (p<0.001) | +0.012 (p=0.27, ns) | — | −0.06 (p=0.82) | +0.03 (p=0.86) | shared |
| panTAM | 1.14 | 1.00 (p<0.001) | −0.002 (p=0.89, ns) | — | −0.09 (p=0.77) | −0.03 (p=0.87) | shared generic |

### Three-part answer to the locked comparator question
1. **Present?** YES. complement_C1Q / CLEC_LAM8 are **tumor-induced in PCa too** (≈+1.0 PCa-SD, p≤0.002), follow a
   benign→distal→involved→tumor **gradient** (complement slope p=0.033), and **localize to PCa TAM/TIM** vs monocytes
   (Wilcoxon p=0.016). The program is **not RCC-exclusive**.
2. **RCC-skewed?** YES (quantitative). The cross-cancer interaction is significant for complement (+0.210, p=0.003),
   CLEC_LAM8 (+0.106, p=0.012) and RCC-skew core (+0.138, p=0.010) — RCC induces the program **more strongly** —
   while shared programs (SPP1/APOE-TREM2/MERTK-GPNMB/panTAM) show no skew. The **ATF3/NF-κB core is comparably
   strong in both** (RCC 1.54 / PCa 1.52; only a small skew p=0.033) — confirming it as the **conserved** monocyte core.
3. **Immune-associated in PCa?** **NO — functionally decoupled.** In PCa (n=23–24 samples, comparable power to RCC's
   n=24), complement_C1Q is **inert across every immune outcome** (CD8 cytotoxicity coef −0.02 p=0.87; exhaustion
   p=0.83; Treg p=0.17; MHC-I p=0.68; MHC-II p=0.33; head-to-head vs panTAM null). The RCC cytotoxicity association
   (−0.82, p=0.04) is **absent** in PCa with a near-**zero point estimate** (not merely non-significant). → the same
   molecular program is coupled to CD8 suppression **only in the RCC niche**.

### Interpretation (sharpens, does not overturn, the thesis)
The complement-high CLEC_LAM program is **RCC-skewed and RCC-functionally-specific, not RCC-unique**: it exists and is
TAM-localized in prostate bone-met myeloid, but (i) RCC induces it more strongly and (ii) its coupling to an
immune-evasive CD8 phenotype is present only in RCC. The ATF3/NF-κB monocyte core is genuinely **conserved** across
both cancers. Manuscript wording: "RCC-skewed and functionally RCC-specific complement-high CLEC_LAM program on a
conserved ATF3/NF-κB monocyte core," never "RCC-unique / absent in prostate."

**Gap (clinical relevance in PCa):** not tested — no prostate bulk survival/ICB cohort is staged (the PCa data are
bone-met scRNA; TCGA-PRAD not staged). A PCa survival mirror of Phase 4 would need TCGA-PRAD; listed as a gap, not run.
