# Publication figures (Figure 4 onward)

All figures are generated from committed data by `build_figures.py` on top of the shared style
module `style.py`. Nothing is screenshotted or re-photographed; every panel is regenerated as vector
art (PDF + SVG) and the multi-panel figures are assembled at exact journal width. Seeds are fixed at 0.

```
envs/rcc_reinterp_venv/bin/python figures/build_figures.py            # all figures
envs/rcc_reinterp_venv/bin/python figures/build_figures.py 4 5        # selected
```

Deliverables live in `figures/` (assembled `FigureN.pdf` + 600-dpi `FigureN.png` preview) and
`figures/panels/` (individual `*.pdf` + `*.svg` for downstream refinement in Illustrator/Inkscape).
`figures/cache/` holds regenerated per-cell/per-patient tables so rebuilds do not reload large objects.

## Style
`style.py` — Okabe-Ito / ColorBrewer colorblind-safe palette used identically everywhere: PCa `#E69F00`,
ccRCC `#0072B2`, ATF3/NF-κB `#7B3294`, complement_C1Q `#1B9E77`, panTAM/control `#999999`; a single
benign→tumor sequential gradient (`#FDE0C8 → #7F2704`) for all continuous scores. Fonts embedded as
editable vector text (`pdf.fonttype = 42`, `svg.fonttype = 'none'`); panel letters 15 pt bold, axis
titles 8 pt, ticks/legends 7 pt, nothing below 6 pt. Top/right spines removed, ticks outward.

Font note: Arial/Helvetica/Liberation Sans are not installed on this host; the specified
Helvetica-metric fallback (Nimbus Sans) is embedded instead. Swap `font.family` in `style.py` if Arial
is available at typesetting time.

## Figures

| Figure | Width | Panels | Source data |
|---|---|---|---|
| Figure 4 — cross-cancer module comparison | 183 mm | A tumor induction PCa vs ccRCC; B cancer-type×condition interaction (RCC skew); C PCa within-tumor TAM localization | `pca_comparator/.../pca_vs_rcc_contrast_table.csv`, `pca_within_tumor_localization.csv` |
| Figure 5 — C1q source and predicted complement signaling | 183 mm | A C1QA/B/C source by cell class; B CellChat ApoE/complement/GAS axes to CLEC_LAM; C NicheNet ligand ranking; D ccRCC myeloid complement_C1Q UMAP | `complement_source_by_class.csv`, `cellchat_allLR_tumor.csv`, `nichenet_prespecified_ranks.csv`, `module_scores_percell.parquet` + `myeloid_FINAL_labels.h5ad` |
| Figure 6 — CD8 dysfunction and survival | 183 mm | A complement_C1Q vs CD8 cytotoxicity (RCC vs PCa); B ccRCC immune context (exhaustion, MHC-II); C TCGA-KIRC OS forest (adjusted, per SD); D TCGA-KIRC KM; E CPTAC external KM | `pca_vs_rcc_contrast_table.csv`, `phase3_adjusted_associations.csv`, `phase4_tcga_kirc_cox.csv`, frozen `risk_model.json` on TCGA-KIRC + CPTAC |
| Figure 7 — myeloid-derived 25-gene signature | 183 mm | A LASSO-Cox coefficients; B discrimination (C-index + time-AUC, TCGA + CPTAC); C benchmark vs random/panTAM/Obradovic; D cross-tumor prostate BCR | `lasso_selected.csv`, `discrimination_metrics.csv`, `sanity_checks.csv`, `pca_univariate_bcr.csv` |
| Figure S1 — model construction and calibration | 183 mm | A LASSO-Cox regularization path (10-fold CV, seed 0); B calibration (3-yr); C decision-curve analysis; D nomogram contributions | frozen `risk_model.json` + `multivariable_cox.csv` on TCGA-KIRC |
| Figure S2 — risk-group landscape | 183 mm | A predicted drug sensitivity (oncoPredict); B TIDE and components; C functional enrichment of TAM-vs-monocyte DEGs; D risk-group immune infiltration | `riskgroup_drug.csv`, `riskgroup_tide.csv`, `enrichment_degs.csv`, `riskgroup_immune.csv` |

Figures 1–3 are not produced here. Significance stars: `*` p<0.05, `**` p<0.01, `***` p<0.001,
`****` p<0.0001; FDR-adjusted where more than one test. Non-significant comparisons are left
unmarked. Manuscript-ready legends for every panel are in `FIGURE_CAPTIONS.md`.
