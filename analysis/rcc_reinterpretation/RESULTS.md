# Results

## A complement-high CLEC_LAM macrophage state is enriched in RCC bone-metastasis myeloid cells

We scored a CLEC_LAM core program (C1QA, C1QB, C1QC, APOE, APOC1, TREM2, GPNMB, MERTK) across myeloid cells from
RCC and prostate bone-metastasis marrow, using both an additive gene score and AUCell; the two agreed closely in
TAM/TIM (Spearman 0.97 in RCC, 0.77 in prostate). The core genes were far better expressed in RCC TAM than in
prostate TAM (for example C1QA in 83% versus 33% of cells, TREM2 in 45% versus 3%, MERTK in 16% versus 3%),
which pointed to a state that is much more developed in the RCC macrophage compartment.

A structural feature of the data shaped the analysis: RCC benign marrow contains almost no TAM or TIM cells
(median 1.5 per benign sample; none reached ten cells), because these tumour-associated states are essentially
tumour-restricted. A TAM-gated tumour-versus-benign contrast is therefore undefined on the benign side, so we
used the mean core score over all myeloid cells per sample as the anchored unit, which captures both the
emergence of TAMs and their per-cell state. On that unit the core score was strongly enriched in RCC tumour
relative to benign marrow (Δ +0.174 [0.099, 0.249], p=5.2e-6, FDR 7.8e-6; nine tumour and seven benign patients).

## RCC skews myeloid cells toward this state more than prostate

Against the shared benign donors that both cancers were profiled against, the state was skewed toward RCC: the
cancer-by-condition interaction on the all-myeloid core score was +0.106 [0.023, 0.190] (p=0.012, FDR 0.012),
which equals one standard deviation of the patient-level score. The share of CLEC_LAM-high cells was likewise
skewed (arcsin fraction interaction +0.515 [0.326, 0.704], p=9.2e-8). These are two views of the same shift:
more cells cross the high-state threshold because the per-cell program intensifies, while the TAM/TIM
cell-type fraction itself does not differ between the cancers (interaction p=0.18). The skew was a per-cell
state effect — within TAM/TIM the RCC-versus-prostate tumour difference was +0.137 [0.035, 0.240] (p=0.008) and
survived adjustment for TAM/TIM fraction (+0.080 [0.013, 0.147], p=0.020).

The skew was robust to how the cells were scored and to the choice of benign reference. It ranged from +0.106 in
the primary pseudobulk model (p=0.012) to +0.159 in the scVI-normalized representation (p=5.7e-4) and +0.087 in a
raw common-space scoring (p=0.025), with the benign-reference variants near p=0.07; leaving out one patient at a
time kept the direction but reached p=0.098 at worst, in keeping with the small prostate TAM compartment.

## Complement carries the skew, and the macrophages make it themselves

Scored as separate modules, the skew concentrated in the complement (C1QA/B/C) part of the program: its
interaction (+0.210 [0.071, 0.349], p=0.003, FDR 0.014, +1.13 SD) was stronger than the full CLEC_LAM8 score
(+0.106), while a generic macrophage-burden score showed no skew (panTAM +0.00, p=0.89), as did SPP1-TAM,
APOE/TREM2 and MERTK/GPNMB modules. Per gene, the complement genes were the most RCC-specific in their
tumour induction (C1QA Δ 1.30 in RCC versus 0.40 in prostate; C1QB 1.17 vs 0.36; C1QC 1.16 vs 0.40), followed by
the lipid and TREM2 genes, whereas MERTK and GPNMB were induced similarly in both cancers.

In the tumour niche, C1q was produced by the macrophages, not the tumour cells. C1QA/B/C were highest in
CLEC_LAM TAM (mean 3.47/3.71/3.72, expressed in essentially all such cells) and near-absent in tumour cells
(0.07–0.15, 3.5–7.4%); C3 was highest in CLEC_LAM TAM (0.90) with a smaller tumour contribution (0.41). We
therefore treat this as tumour-niche macrophage autocrine and paracrine complement centred on C1q, with
C3→C3AR1/CR3 as a related but separately tested axis, and make no claim of tumour-cell-to-macrophage C1q
signalling.

## Predicted signalling into the CLEC_LAM state

Because the benign niche has essentially no TAM to receive signals, we established the tumour-specificity of
candidate axes from patient-level ligand and receptor expression rather than from a benign-versus-tumour
communication contrast, and treat the CellChat and NicheNet results as predicted rather than causal. Three
independent readings agreed that complement is the best-supported tumour-specific axis into CLEC_LAM TAM: the
ligands were elevated in tumour (C1QA +1.74 p=0.003, C1QB +1.81 p=0.026, C1QC +2.29 p=0.015), CellChat placed
complement signalling (C3 to C3AR1 and the CR3 integrin) into the TAM, and NicheNet ranked C1QB 12th of 116
ligands for the CLEC_LAM program. TGF-β was a supporting tumour-specific co-axis (TGFB1 ligand +0.32, p=4e-4;
NicheNet rank 5). The APOE→TREM2 and GAS6/PROS1→MERTK efferocytosis axes were active in the niche but not
tumour-gained — their ligands were not elevated in tumour (APOE p=0.96) or were higher in benign (GAS6, p=0.02) —
consistent with MERTK and efferocytosis being shared rather than RCC-specific.

The tumour cells acting as senders were identified by the clear-cell RCC epithelial markers CA9 and PAX8 (mean
0.91 versus about zero in stroma and immune cells); a copy-number confirmation was attempted with CopyKAT but did
not complete, and inferCNV was not run.

## The state is associated with an immune-evasive CD8 phenotype

At the patient level (24 RCC samples, adjusted for TAM, CD8 and malignant fraction and condition), the
complement-high state tracked an immune-evasive CD8 phenotype, and did so more specifically than generic
macrophage burden. Head to head against panTAM it was associated with lower CD8 cytotoxicity (complement −0.909,
p=0.029, versus panTAM −0.264, p=0.359), higher exhaustion (+0.610, p=0.020 vs +0.292, p=0.114) and higher
MHC-II (+1.189, p=0.019 vs −0.007, p=0.983). The single-module adjusted associations were cytotoxicity −0.817
[−1.586, −0.049] (p=0.039), exhaustion +0.509 [0.002, 1.016] (p=0.049) and MHC-II +1.191 [0.299, 2.084]
(p=0.012). This specificity applies to CD8 cytotoxicity, exhaustion and MHC-II; regulatory-T-cell burden was not
associated with complement (+0.154, p=0.665) but was with SPP1-TAM (+1.211, p<0.001) and MERTK/GPNMB (+0.654,
p=0.001), and so reflects a broader macrophage effect. MHC-I antigen presentation could not be assessed with
only nine tumour samples (null).

## The complement/RCC-skew core is prognostic in ccRCC

In TCGA-KIRC (532 patients, 175 deaths) the complement module predicted overall survival and remained
prognostic after adjusting for generic macrophage burden, the prior Obradovic signature and total immune content
(per-SD HR 1.26 [1.09, 1.45], p=0.002 unadjusted; 1.61 [1.15, 2.27], p=0.006 adjusted); the RCC-skew core gave
the strongest adjusted signal (HR 2.24 [1.36, 3.68], p=0.002). Generic macrophage burden was not prognostic
(panTAM p=0.11), and the APOE/TREM2 module was not independent of the Obradovic signature (HR 1.17 p=0.033
unadjusted, 1.01 p=0.95 adjusted), so the survival signal is carried by the complement/RCC-skew part rather than
by macrophage abundance or the prior signature. This adjustment did not include stage, grade, age or tumour
purity, which were absent from the clinical table used, so we do not claim independence from those factors. The
signal was specific to overall survival; the fuller CLEC_LAM8 composite, which also carries the shared-support
genes, lost significance on adjustment (HR 1.19 p=0.015 unadjusted, 1.54 p=0.28 adjusted).

The state was prognostic but not a validated predictor of checkpoint-immunotherapy benefit. In the Braun
nivolumab-treated patients there was no association with survival or response (complement OS p=0.41, response
p=0.30); the only signal was an exploratory nivolumab-versus-everolimus arm interaction in CM-025 (CLEC_LAM8 HR
0.53 [0.32, 0.90], p=0.018, underpowered at 250 patients).

## Relation to prior work

We do not claim to discover C1Q/APOE/TREM2 macrophages in ccRCC — those were described by Obradovic and
colleagues in 2021. What is new here is the cross-cancer, shared-benign design and the finding that the skew is
carried specifically by a macrophage-autocrine complement program, that this program is coupled to an
immune-evasive CD8 phenotype, and that it is prognostic for overall survival beyond generic macrophage burden
and the prior signature.
