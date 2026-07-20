# Prostate comparator: is the RCC complement-high CLEC_LAM program present, immune-associated, and clinically relevant in prostate cancer?

We applied the same patient-level pseudobulk framework used in the RCC study (`../rcc_reinterpretation/`) to
prostate myeloid cells (`myeloid_FINAL.h5ad`) and the prostate full niche
(`integrated_with_kfoury_labels.h5ad`, 32 samples), adding an ATF3/NF-κB module. Scripts are in `scripts/`
(06/31/32) and result tables in `outputs/tables/`.

## Module contrast between the two cancers

| module | RCC tumour-vs-benign (SD) | PCa tumour-vs-benign (SD) | cross-cancer interaction | PCa TAM/TIM localization | RCC cytotoxicity | PCa cytotoxicity |
|---|---|---|---|---|---|---|
| complement_C1Q | 1.41 | 1.00 (p=0.002) | +0.210 (p=0.003) | Δ+0.11, p=0.016 | −0.82 (p=0.04) | −0.02 (p=0.87) |
| CLEC_LAM8 | 1.37 | 1.03 (p=0.001) | +0.106 (p=0.012) | Δ+0.07, p=0.016 | −0.77 (p=0.04) | +0.01 (p=0.95) |
| RCC_skew_CORE | 1.37 | 1.03 (p<0.001) | +0.138 (p=0.010) | — | −0.78 (p=0.03) | +0.01 (p=0.96) |
| ATF3_NFkB | 1.54 | 1.52 (p<0.001) | +0.079 (p=0.033) | Δ+0.10, p=0.016 | — | −0.31 (p=0.21) |
| SPP1_TAM | 1.33 | 1.17 (p<0.001) | +0.020 (p=0.63) | Δ+0.06, p=0.30 | +0.61 (p=0.12) | +0.04 (p=0.81) |
| APOE_TREM2 | 1.27 | 0.92 (p<0.001) | +0.050 (p=0.13) | Δ+0.04, p=0.016 | −0.55 (p=0.08) | +0.02 (p=0.87) |
| MERTK_GPNMB | 1.21 | 0.96 (p<0.001) | +0.012 (p=0.27) | — | −0.06 (p=0.82) | +0.03 (p=0.86) |
| panTAM | 1.14 | 1.00 (p<0.001) | −0.002 (p=0.89) | — | −0.09 (p=0.77) | −0.03 (p=0.87) |

The complement_C1Q and CLEC_LAM8 programs are induced in prostate tumours as well (about +1.0 prostate-SD,
p ≤ 0.002), rise along the benign-to-tumour gradient (complement slope p=0.033), and are enriched in prostate
TAM/TIM over monocytes (Wilcoxon p=0.016). The program is therefore not confined to RCC. It is, however, induced
more strongly in RCC: the cross-cancer interaction is significant for complement (+0.210, p=0.003), CLEC_LAM8
(+0.106, p=0.012) and the RCC-skew core (+0.138, p=0.010), whereas SPP1, APOE/TREM2, MERTK/GPNMB and panTAM show
no skew. The ATF3/NF-κB core is close to equal in the two cancers (RCC 1.54, PCa 1.52, small skew p=0.033),
consistent with it being the conserved part of the monocyte response.

Where the two cancers part ways is function. In prostate (23–24 samples, comparable in size to the RCC set of
24) complement_C1Q showed no relationship with any immune readout: CD8 cytotoxicity coefficient −0.02 (p=0.87),
exhaustion p=0.83, regulatory T cells p=0.17, MHC-I p=0.68, MHC-II p=0.33, and nothing in the head-to-head
against panTAM. The reduced-cytotoxicity association seen in RCC (−0.82, p=0.04) is essentially zero in prostate,
so the same molecular program is coupled to CD8 suppression only in the RCC setting.

Taken together, the complement-high CLEC_LAM program is skewed toward RCC and functionally specific to RCC, but
not unique to it: it is present and TAM-localized in prostate bone-metastasis myeloid cells, yet RCC drives it
harder and only in RCC is it tied to an immune-evasive CD8 phenotype. The ATF3/NF-κB monocyte core is genuinely
shared between the two cancers.

We later tested clinical relevance in prostate directly, using TCGA-PRAD biochemical recurrence (see
`../rcc_reinterpretation/outputs/model/MODEL_RESULTS.md`): the complement program was not associated with
recurrence, mirroring the functional decoupling seen here at the single-cell level.
