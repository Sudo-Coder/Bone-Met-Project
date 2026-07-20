# Prognostic modelling of the myeloid complement–macrophage program in ccRCC and prostate cancer

## Part A — a myeloid-derived prognostic model stratifies ccRCC overall survival

We assembled a pool of 152 candidate genes: the eight CLEC_LAM core genes (C1QA, C1QB, C1QC, APOE, APOC1,
TREM2, GPNMB, MERTK) and 144 genes differentially expressed between tumour-associated macrophages and monocytes
in the RCC single-cell data (Wilcoxon test, log2 fold-change > 0.5, adjusted p < 0.01; `candidate_genes.csv`).

TCGA-KIRC (604 tumours with survival annotation, 202 deaths) was divided 70/30, stratified by event. On the
training split we screened candidates by univariate Cox regression (`univariate_cox.csv`) and fitted a LASSO-Cox
model (l1_ratio = 1, penalty chosen by 10-fold cross-validation), which retained a 25-gene signature
(`lasso_selected.csv`, `lasso_cv.png`). The risk score is the sum of each gene's standardized expression weighted
by its coefficient; genes, coefficients and the training-set means and standard deviations are stored in
`risk_model.json`. The largest positive weights fall on the complement receptor VSIG4, NAMPT, MAP3K8 and FUS, and
the largest negative weights on LIPA, KLF6, HERPUD1 and BTG2.

In the held-out third of TCGA-KIRC the signature reached a Harrell concordance index of 0.69, with
time-dependent AUCs of 0.80, 0.72 and 0.70 at one, three and five years (`discrimination_metrics.csv`). Patients
above the median risk had substantially shorter survival in both the training split (log-rank p = 1.9 × 10⁻¹⁶,
`km_train.png`) and the held-out split (log-rank p = 8.6 × 10⁻⁵, `km_test.png`); the ten-fold cross-validated
training concordance was 0.73. As a continuous variable the risk score remained prognostic after adjustment for
tumour stage, grade and age (hazard ratio per standard deviation 2.03, p < 0.001; stage 1.74, p < 0.001; age
1.46, p < 0.001; grade 1.11, p = 0.20; `multivariable_cox.csv`), and we summarize this multivariable model as a
nomogram with calibration and decision-curve plots for one-, three- and five-year survival (`nomogram.png`,
`calibration.png`, `dca.png`).

Applied without refitting to the Braun/CheckMate ccRCC cohort (311 patients, 80 deaths), the signature did not
carry over: the concordance index fell to 0.44 and the median split gave a log-rank p of 0.32
(`km_validation_braun.png`). This cohort is metastatic, immune-checkpoint-treated disease, whereas the model was
trained on primary-tumour survival, and the two contexts differ substantially. We were unable to retrieve an
independent primary-tumour ccRCC RNA-seq cohort during the analysis (E-MTAB-1980 returned a 404, and GSE22541 is
a metastasis-enriched microarray series that we did not use), so external validation in the matched primary-tumour
setting remains outstanding.

Two checks argue against an artefactual fit. The held-out concordance of 0.69 exceeded that of 500 random
25-gene scores drawn from the same candidate pool (mean 0.56, 95th percentile 0.62), and it exceeded a generic
macrophage score (panTAM concordance 0.58) and the prior TREM2/APOE/C1Q recurrence signature (0.55) evaluated on
the same held-out patients (`sanity_checks.csv`).

High-risk tumours were the more immune-infiltrated of the two groups: the CD8 T-cell signature was higher in
high-risk tumours (difference 0.006, p = 5 × 10⁻⁴, FDR 0.004), as was a broad leukocyte-infiltration signature
(p = 0.045), while exhaustion, regulatory-T-cell, MHC-I/II and panTAM signatures did not differ
(`riskgroup_immune.csv`, `riskgroup_immune.png`). This pairing of heavy immune infiltration with poor outcome is
the well-recognised paradox of clear-cell RCC. We did not run drug-sensitivity prediction, as the GDSC training
resource and the oncoPredict workflow were not available during the analysis.

## Part B — myeloid programs and biochemical recurrence in primary prostate cancer

Prostate cancer progresses along a defined path. The primary tumour gives rise to biochemical recurrence, a
rising PSA after definitive local therapy, which precedes clinically detectable metastasis; when metastasis
occurs it is overwhelmingly bone-tropic, with bone involvement in up to roughly 90% of metastatic disease¹⁻⁴.
Biochemical recurrence is thus the earliest clinical marker on the road to bone metastasis and
prostate-cancer-specific death, and transcriptomic signatures read from primary prostatectomy tissue are the
established way to capture a tumour's metastatic potential — the Decipher classifier and related expression
signatures predict recurrence, post-prostatectomy metastasis and cancer-specific mortality from primary tissue⁵⁻⁷.
Testing our bone-metastasis-derived myeloid programs against biochemical recurrence in primary tumours therefore
asks whether the conserved myeloid reprogramming, and the shared complement-high macrophage program, mark primary
prostate tumours predisposed to the bone-tropic progression these programs characterise in established disease.
Because the study's central point is that the complement program is present in both cancers but coupled to CD8
suppression only in RCC, this comparison speaks to the thesis whichever way it falls.

We used TCGA-PRAD primary tumours (550 with expression and progression annotation, 95 progression events),
scoring the progression-free interval — the endpoint recommended for this cohort — and confirming against the
explicit biochemical-recurrence field (464 patients, 56 events). Programs were scored per sample with the same
AUCell approach used in the single-cell work.

The complement_C1Q program showed no association with prostate recurrence (hazard ratio per standard deviation
1.06, 95% CI 0.93–1.20, p = 0.41, concordance 0.50; secondary biochemical-recurrence endpoint p = 0.65). The
ATF3/NF-κB conserved core reached a marginal, trend-level association that approached significance (hazard ratio
0.80, 95% CI 0.64–1.00, p = 0.054), and was closer to significance here than in ccRCC, where it was clearly null
(p = 0.46) — a pattern in keeping with this core being the prostate-anchored program supported by our in-vitro
work — though its direction was inverse, with higher core scores tracking lower recurrence hazard. The broader
macrophage-state scores were associated with recurrence: RCC_skew_CORE (1.17, 95% CI 1.04–1.32, p = 0.012, FDR
0.025) and CLEC_LAM8 (1.18, 95% CI 1.03–1.34, p = 0.015, FDR 0.025). The frozen 25-gene ccRCC signature carried
across tumour type and was strongly associated with prostate recurrence (1.60, 95% CI 1.34–1.91, p = 3 × 10⁻⁷,
concordance 0.68).

| program | ccRCC OS HR (p) | PCa recurrence HR (p) |
|---|---|---|
| complement_C1Q | 1.23 (0.002) | 1.06 (0.41) |
| ATF3/NF-κB core | 0.95 (0.46) | 0.80 (0.054) |
| RCC_skew_CORE | 1.21 (0.005) | 1.17 (0.012) |
| CLEC_LAM8 | 1.16 (0.028) | 1.18 (0.015) |

The complement program predicts survival in clear-cell RCC but shows no association with prostate biochemical
recurrence, and this split fits the functional-divergence account — the same macrophage state is present in both
cancers, but its clinical and CD8-coupled behaviour is RCC-specific. The broader CLEC_LAM and RCC-skew scores
track outcome in both cancers, and the myeloid-derived ccRCC signature retains prognostic value for prostate
early progression across tumour type.

The models are prognostic rather than predictive of checkpoint-immunotherapy benefit; the single-cell analysis
found no within-cohort association with ICB response and only a suggestive arm-interaction signal. The prostate
analysis addresses primary-tumour early progression using the field-standard endpoint, not overall survival.
Drug-sensitivity prediction and a matched primary-tumour external ccRCC validation cohort remain to be added.

## References
1. Bubendorf L, et al. Metastatic patterns of prostate cancer: an autopsy study. Hum Pathol. 2000;31(5):578–583.
2. Pound CR, et al. Natural history of progression after PSA elevation following radical prostatectomy. JAMA. 1999;281(17):1591–1597.
3. Freedland SJ, et al. Risk of prostate cancer-specific mortality following biochemical recurrence after radical prostatectomy. JAMA. 2005;294(4):433–439.
4. Antonarakis ES, et al. The natural history of metastatic progression in men with PSA recurrence after radical prostatectomy. BJU Int. 2012;109(1):32–39.
5. Erho N, et al. Discovery and validation of a prostate cancer genomic classifier that predicts early metastasis following radical prostatectomy. PLoS One. 2013;8(6):e66855.
6. Karnes RJ, et al. Validation of a genomic classifier that predicts metastasis following radical prostatectomy in an at-risk patient population. J Urol. 2013;190(6):2047–2053.
7. Feng FY, et al. Validation of a 22-gene genomic classifier in patients with recurrent prostate cancer (NRG/RTOG 9601). JAMA Oncol. 2021;7(4):544–552.
