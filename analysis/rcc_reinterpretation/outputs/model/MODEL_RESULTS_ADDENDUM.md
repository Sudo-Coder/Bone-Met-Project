# Prognostic-model addendum: therapeutic landscape, external validation and signature biology

These analyses reuse the frozen 25-gene ccRCC risk model without refitting. The risk score was recomputed
deterministically from `risk_model.json` on the TCGA-KIRC expression matrix and split at the previous median,
recovering the identical high- and low-risk groups (302 patients each).

## External validation in an independent primary-tumour cohort

Applied without refitting to CPTAC clear-cell RCC (110 primary tumours with survival, 21 deaths), the signature
reached a Harrell concordance of 0.69 and time-dependent AUCs of 0.73 and 0.68 at one and three years, matching
its performance in the held-out TCGA-KIRC set. The median split separated the survival curves in the expected
direction but did not reach significance at this event count (log-rank p = 0.096). Because CPTAC is RNA-seq and
the model was trained on the TCGA microarray-scale matrix, the signature genes were standardized within the CPTAC
cohort before applying the frozen coefficients. This independent primary-tumour validation removes the
single-cohort limitation of the original discrimination results; the metastatic, checkpoint-treated Braun cohort
remained the one setting where the model did not transfer (concordance 0.44). Results are added to
`discrimination_metrics.csv`, `km_validation_external.png` and `roc_validation_external.png`.

## Predicted drug sensitivity of the risk groups

We predicted response to ccRCC agents for the TCGA-KIRC tumours with oncoPredict, training on the GDSC1 cell-line
expression and IC50 matrices (GDSC2, the default oncoPredict training set, contains only axitinib among these
agents, so GDSC1 was used; everolimus is absent from GDSC1 and its mTOR-inhibitor class is represented by
temsirolimus). High-risk tumours were predicted to be more sensitive than low-risk tumours to axitinib (Δ log-IC50
−0.12, p = 7e-4, FDR = 1e-3), cabozantinib (−0.15, p = 3e-5, FDR = 4e-5) and temsirolimus (−0.15, p = 4e-5, FDR =
1e-4), with no difference for sunitinib (p = 0.99) or pazopanib (p = 0.16). The high-risk, complement-macrophage
group is thus predicted to retain, and in places exceed, sensitivity to the VEGFR/MET tyrosine-kinase inhibitors
and mTOR inhibition used in ccRCC. Values are in `riskgroup_drug.csv` and `riskgroup_drug.png`.

## Biology of the signature and its source DEGs

Over-representation analysis (Enrichr, BH-adjusted) of the 144 tumour-associated-macrophage-versus-monocyte DEGs
that seeded the model was dominated by antigen processing and presentation (adjusted p = 2e-25), the phagosome
(2e-16) and related immune and macrophage pathways, placing the candidate pool squarely in antigen-presenting
macrophage biology. The compact 25-gene signature itself did not reach adjusted significance for any single term;
its leading themes were cytokine-mediated signalling and the unfolded-protein/stress response, consistent with a
small predictor drawn from across the macrophage program rather than one pathway. Tables are `enrichment_signature.csv`
and `enrichment_degs.csv`, with `enrichment.png`.

## Immune-response prediction reinforces the prognostic-not-predictive boundary

TIDE scores computed for the TCGA-KIRC tumours were higher in the high-risk group (0.30 versus 0.06, p < 1e-4),
driven by greater T-cell dysfunction (0.35 higher, p < 1e-4) rather than T-cell exclusion (p = 0.18); predicted
MDSC infiltration was also higher (p < 1e-4). High-risk tumours therefore carry a more immune-dysfunctional,
checkpoint-resistant phenotype, which fits the prognostic-not-ICB-predictive reading of the state: the high-risk
group is immune-infiltrated yet dysfunctional. Results are in `riskgroup_tide.csv` and `riskgroup_tide.png`.
