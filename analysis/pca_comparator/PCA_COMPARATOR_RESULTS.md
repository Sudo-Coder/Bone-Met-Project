# Prostate comparator: is the RCC complement-high CLEC_LAM program present, immune-associated, and clinically relevant in prostate cancer?

We applied the same patient-level pseudobulk framework used in the RCC study (`../rcc_reinterpretation/`) to
prostate myeloid cells (`myeloid_FINAL.h5ad`) and the prostate full niche (`integrated.h5ad`, 32 samples,
75,089 cells, with cell labels joined from `Data/cell-annotations.csv`), adding an ATF3/NF-κB module. Scripts
are in `scripts/` (06/31/32) and result tables in `outputs/tables/`.

## Module contrast between the two cancers

| module | RCC tumour-vs-benign (SD) | PCa tumour-vs-benign (SD) | cross-cancer interaction | PCa TAM/TIM localization | RCC cytotoxicity | PCa cytotoxicity |
|---|---|---|---|---|---|---|
| complement_C1Q | 1.41 | 1.00 (p=0.002) | +0.210 (p=0.003) | Δ+0.11, p=0.016 | −0.82 (p=0.04) | +0.00 (p=0.99) |
| CLEC_LAM8 | 1.37 | 1.03 (p=0.001) | +0.106 (p=0.012) | Δ+0.07, p=0.016 | −0.77 (p=0.04) | +0.04 (p=0.81) |
| RCC_skew_CORE | 1.37 | 1.03 (p<0.001) | +0.138 (p=0.010) | — | −0.78 (p=0.03) | +0.04 (p=0.82) |
| ATF3_NFkB | 1.54 | 1.52 (p<0.001) | +0.079 (p=0.033) | Δ+0.10, p=0.016 | — | −0.30 (p=0.25) |
| SPP1_TAM | 1.33 | 1.17 (p<0.001) | +0.020 (p=0.63) | Δ+0.06, p=0.30 | +0.61 (p=0.12) | +0.11 (p=0.53) |
| APOE_TREM2 | 1.27 | 0.92 (p<0.001) | +0.050 (p=0.13) | Δ+0.04, p=0.016 | −0.55 (p=0.08) | +0.03 (p=0.82) |
| MERTK_GPNMB | 1.21 | 0.96 (p<0.001) | +0.012 (p=0.27) | — | −0.06 (p=0.82) | +0.06 (p=0.72) |
| panTAM | 1.14 | 1.00 (p<0.001) | −0.002 (p=0.89) | — | −0.09 (p=0.77) | −0.01 (p=0.97) |

The complement_C1Q and CLEC_LAM8 programs are induced in prostate tumours as well (about +1.0 prostate-SD,
p ≤ 0.002), rise along the benign-to-tumour gradient (complement slope p=0.033), and are enriched in prostate
TAM/TIM over monocytes (Wilcoxon p=0.016). The program is therefore not confined to RCC. It is, however, induced
more strongly in RCC: the cross-cancer interaction is significant for complement (+0.210, p=0.003), CLEC_LAM8
(+0.106, p=0.012) and the RCC-skew core (+0.138, p=0.010), whereas SPP1, APOE/TREM2, MERTK/GPNMB and panTAM show
no skew. The ATF3/NF-κB core is close to equal in the two cancers (RCC 1.54, PCa 1.52, small skew p=0.033),
consistent with it being the conserved part of the monocyte response.

Where the two cancers part ways is function. In prostate (31–32 samples across all four conditions, against 24
in RCC) complement_C1Q showed no relationship with the CD8 readouts: cytotoxicity coefficient +0.003 (p=0.99),
exhaustion +0.038 (p=0.85), regulatory T cells p=0.245, MHC-I p=0.669 (n=14). The head-to-head against panTAM
was likewise null for CD8 (exhaustion: complement +0.121 p=0.693 versus panTAM −0.121 p=0.714; cytotoxicity
+0.014 p=0.951). The reduced-cytotoxicity association seen in RCC (−0.82, p=0.04) is essentially zero in
prostate, so the same molecular program is coupled to CD8 suppression only in the RCC setting.

The one prostate readout that is not flat is MHC-II. Complement_C1Q alone was not associated (+0.198, p=0.290),
but it reached nominal significance head-to-head against panTAM (+0.540, p=0.049), as did APOE_TREM2 (+0.333,
p=0.038) and MERTK_GPNMB (+0.339, p=0.044) in single-module models. These are uncorrected and marginal, and we
do not treat them as an immune association in prostate; they are noted so the "no immune relationship" wording
is confined to the CD8 and Treg readouts.

Taken together, the complement-high CLEC_LAM program is skewed toward RCC and functionally specific to RCC, but
not unique to it: it is present and TAM-localized in prostate bone-metastasis myeloid cells, yet RCC drives it
harder and only in RCC is it tied to an immune-evasive CD8 phenotype. The ATF3/NF-κB monocyte core is genuinely
shared between the two cancers.

We later tested clinical relevance in prostate directly, using TCGA-PRAD biochemical recurrence (see
`../rcc_reinterpretation/outputs/model/MODEL_RESULTS.md`): the complement program was not associated with
recurrence, mirroring the functional decoupling seen here at the single-cell level.
