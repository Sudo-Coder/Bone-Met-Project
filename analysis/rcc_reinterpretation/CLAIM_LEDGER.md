# Summary of findings and their limits

Module scores are AUCell enrichment values (0–1). Effect sizes given in SD units use the standard deviation of
the patient-level pseudobulk all-myeloid score (CLEC_LAM8 SD 0.107, complement_C1Q SD 0.185).

## Main findings

A CLEC_LAM-like complement–lipid–TREM2 macrophage state is present in RCC myeloid cells and is enriched in
tumour relative to benign marrow (all-myeloid CORE score, tumour vs benign, Δ +0.174 [0.099, 0.249], p=5.2e-6,
FDR 7.8e-6). Against a shared benign anchor, RCC skews tumour myeloid toward this state more than prostate does
(cancer × condition interaction +0.106 [0.023, 0.190], p=0.012, FDR 0.012, equal to +1.00 SD [0.22, 1.78]).

The skew is carried specifically by the complement (C1QA/B/C) part of the program: the interaction is stronger
for the complement module (+0.210 [0.071, 0.349], p=0.003, FDR 0.014, +1.13 SD) than for the full CLEC_LAM8
score (+0.106), and a generic macrophage-burden score shows no skew at all (panTAM +0.00 [−0.029, 0.025],
p=0.89). It reflects a per-cell intensification of the RCC TAM/TIM state rather than a change in TAM abundance:
within TAM/TIM the RCC-versus-prostate tumour difference is +0.137 [0.035, 0.240] (p=0.008) and survives
adjustment for TAM/TIM fraction (+0.080 [0.013, 0.147], p=0.020), while the TAM/TIM cell-fraction interaction is
not significant (p=0.18).

C1q in the tumour niche comes from the macrophages themselves, not the tumour cells. C1QA/B/C are highest in
CLEC_LAM TAM (mean 3.47/3.71/3.72, expressed in ~100% of cells) and near-absent in tumour cells (0.07–0.15,
3.5–7.4%); C3 is highest in CLEC_LAM TAM (0.90) with a smaller tumour contribution (0.41). We describe this as
tumour-niche macrophage autocrine/paracrine complement centred on C1q, with C3→C3AR1/CR3 as a related but
separately tested axis, and we do not claim tumour-cell-to-TAM C1q signalling.

The complement-high state is associated with an immune-evasive CD8 phenotype at the patient level (24 RCC
samples, adjusted for TAM, CD8 and malignant fraction and condition). Compared head-to-head with generic
macrophage burden it tracks lower CD8 cytotoxicity (complement −0.909, p=0.029, versus panTAM −0.264, p=0.359),
higher exhaustion (+0.610, p=0.020 vs +0.292, p=0.114) and higher MHC-II (+1.189, p=0.019 vs −0.007, p=0.983);
the single-module adjusted effects are cytotoxicity −0.817 [−1.586, −0.049] (p=0.039), exhaustion +0.509
[0.002, 1.016] (p=0.049) and MHC-II +1.191 [0.299, 2.084] (p=0.012). This scope is CD8 cytotoxicity, exhaustion
and MHC-II; regulatory T cells are not part of it (see limits below).

In TCGA-KIRC the complement/RCC-skew core is prognostic for overall survival and remains so after adjusting for
generic macrophage burden and the prior Obradovic signature (complement_C1Q per-SD HR 1.26 [1.09, 1.45], p=0.002
unadjusted; 1.61 [1.15, 2.27], p=0.006 adjusted; 532 patients, 175 deaths). This adjustment covers macrophage
burden, the Obradovic signature and total immune content; it does not include stage, grade, age or tumour
purity, which were not in the clinical table used, so we do not claim independence from those factors.

## Supporting results

The RCC skew holds across scoring approaches and benign references, spanning +0.106 (primary pseudobulk,
p=0.012), +0.159 (scVI-normalized, p=5.7e-4) and +0.087 (raw common-space, p=0.025), with the benign-reference
variants around p=0.07. The share of CLEC_LAM-high cells is also skewed toward RCC (fraction interaction +0.515
[0.326, 0.704], p=9.2e-8). Complement is the best-supported tumour-specific signalling axis into CLEC_LAM TAM by
three independent readings — ligand elevation (C1QA +1.74 p=0.003, C1QB +1.81 p=0.026, C1QC +2.29 p=0.015),
CellChat complement signalling to TAM, and NicheNet (C1QB ranked 12 of 116). TGF-β is a supporting tumour-specific
co-axis (TGFB1 ligand +0.32 p=4e-4; NicheNet rank 5), not a headline. Among the modules, RCC_skew_CORE gives the
strongest adjusted survival signal (HR 2.24 [1.36, 3.68], p=0.002).

## Limits and boundaries

The survival signal sits in the complement/RCC-skew core and washes out in the fuller CLEC_LAM8 composite, which
also contains the shared-support genes: CLEC_LAM8 goes from HR 1.19 (p=0.015) unadjusted to 1.54 (p=0.28)
adjusted, while RCC_skew_CORE stays significant. APOE/TREM2 is not independent of the Obradovic signature (HR
1.17, p=0.033 unadjusted; 1.01, p=0.95 adjusted), so we make no separate APOE/TREM2 survival claim.

The state is prognostic, not a validated predictor of checkpoint-immunotherapy benefit. Within the Braun
nivolumab-treated patients there was no association (complement OS p=0.41, response p=0.30); the only signal was
an exploratory nivolumab-versus-everolimus arm interaction in CM-025 (CLEC_LAM8 HR 0.53 [0.32, 0.90], p=0.018;
RCC-skew 0.57, p=0.033; complement 0.65, p=0.093; 250 patients, underpowered).

The immune-evasion association is specific to CD8 cytotoxicity, exhaustion and MHC-II, not regulatory T cells:
complement shows no Treg association (+0.154, p=0.665), whereas SPP1 and MERTK/GPNMB do (SPP1 +1.211, p<0.001;
MERTK/GPNMB +0.654, p=0.001), so Treg burden is a broader effect. MERTK/GPNMB behaves as a shared-support
member across the study, not an RCC-specific driver (per-gene induction is similar in the two cancers, the
interaction is not significant, p=0.27, and its adjusted TCGA hazard turns protective, HR 0.72, p=0.002).

The cross-cancer interaction is power-limited: leave-one-patient-out keeps the direction but reaches p=0.098 at
worst, and the firmer leg is the within-TAM state contrast (p=0.008). We do not claim to discover
C1Q/APOE/TREM2 TAMs in ccRCC — that is prior art (Obradovic 2021); what is new is the cross-cancer shared-benign
design and the complement-specific, source-resolved framing. Tumour/epithelial cells are identified as senders
by CA9/PAX8 markers rather than by copy-number inference (CopyKAT did not complete and inferCNV was not run).

Not yet done: MHC-I antigen-presentation was underpowered (9 tumour samples, null); a TGF-β response score was
not carried into the survival or immune-evasion models; stage/grade/age/purity adjustment for the TCGA survival
model; a compositional neighbourhood analysis (scCODA/Milo); and any spatial, protein or perturbation validation
of the signalling axes, which stay descriptive.
