#!/usr/bin/env python
"""Assemble additional_analysis.ipynb from the computed outputs (valid nbformat v4 JSON)."""
import json, os
ROOT="/autofs/projects-t3/hussain/scProj"
NB=f"{ROOT}/kidney-cancer/additional_analysis/additional_analysis.ipynb"

def md(*lines): return {"cell_type":"markdown","metadata":{},"source":[l+"\n" for l in lines]}
def code(src): return {"cell_type":"code","metadata":{},"execution_count":None,"outputs":[],
                       "source":[l+"\n" for l in src.strip("\n").split("\n")]}

cells=[]
cells.append(md(
 "# Additional analysis — LA-TAM / NR4A mechanism (RCC myeloid niche)",
 "",
 "Follow-up on the CellRank/celloracle myeloid trajectory (Mono1 → terminal TAM/mDC/osteoclast).",
 "Goal: test whether **NR4A2/NR4A3** regulate the LA-TAM program, what **upstream niche signals**",
 "activate them, and which **pathways/TFs** are active along the TAM trajectory. All results are",
 "hypothesis-generating. Full narrative + claim mapping: `outputs/SUMMARY.md`.",
 "",
 "**Object:** `Cleaned_Data/myeloid_FINAL_labels.h5ad` (16,881 myeloid cells; RCC bone-met + benign BM).",
 "TAM is tumor-restricted (2534/2548 TAM cells in the Tumor faction).",
 "",
 "> Tasks were executed by the scripts in `analysis/` on SLURM; this notebook just loads the results",
 "> from `outputs/`. It does **not** modify any pre-existing file."))

cells.append(code(
 "import pandas as pd, os\n"
 "from IPython.display import Image, display, Markdown\n"
 "pd.set_option('display.width',180); pd.set_option('display.max_columns',30)\n"
 "ROOT='/autofs/projects-t3/hussain/scProj/kidney-cancer/additional_analysis'; T=f'{ROOT}/outputs/tables'; F=f'{ROOT}/outputs/figures'\n"
 "def show(p, **kw):\n"
 "    return pd.read_csv(p, **kw)\n"
 "def fig(name):\n"
 "    display(Image(filename=f'{F}/{name}'))"))

# Task 1
cells.append(md("## Task 1 — pySCENIC regulons + AUCell",
 "262 regulons. **No NR4A2/NR4A3 regulon survived cisTarget motif-pruning.** NFKB1/HIF1A/NFE2L2",
 "regulon AUCell is elevated in TAM-Tumor; NR4A→LA-TAM link unsupported even at the GRNBoost",
 "co-expression level."))
cells.append(code(
 "print('Regulon target summary (focus TFs):')\n"
 "display(show(f'{T}/scenic_regulon_targets.csv'))\n"
 "print('\\nNR4A/NFKB regulon vs LA-TAM hypergeometric (NR4A regulons are EMPTY -> undefined):')\n"
 "display(show(f'{T}/scenic_nr4a_latam_hypergeom.csv'))\n"
 "print('\\nGRNBoost2 co-expression (pre-pruning) NR4A/NFKB/HIF1A vs LA-TAM:')\n"
 "display(show(f'{T}/scenic_grnboost_nr4a_latam_coexpr.csv'))\n"
 "print('\\nAUCell TAM-Tumor specificity (Mann-Whitney, TAM-Tumor vs rest):')\n"
 "display(show(f'{T}/scenic_tam_tumor_specificity.csv'))"))
cells.append(code("fig('scenic_aucell_heatmap.png')"))

# Task 2
cells.append(md("## Task 2 — NicheNet (niche → TAM)",
 "Top ligand driving the LA-TAM program = **TGFB1**; NR4A2's top predicted upstream ligand = **VEGFA**."))
cells.append(code(
 "print('Top 15 ligands by activity (Pearson) for the LA-TAM program:')\n"
 "display(show(f'{T}/nichenet_ligand_activity.csv').head(15))\n"
 "print('\\nTop ligands predicted to induce NR4A2:')\n"
 "display(show(f'{T}/nichenet_top_ligands_for_NR4A2.csv'))"))
cells.append(code("fig('nichenet_ligand_activity.png'); fig('nichenet_ligand_target_heatmap.png')"))

# Task 3
cells.append(md("## Task 3 — decoupleR PROGENy (pathways) + CollecTRI (TF activity)",
 "NF-κB & hypoxia rise strongly into the Tumor faction; NR4A2 activity elevated in tumor TAM/TIM."))
cells.append(code(
 "print('CollecTRI TF activity by faction:')\n"
 "display(show(f'{T}/decoupler_collectri_by_condition.csv', index_col=0)[['NR4A2','NR4A3','NFKB1','RELA','HIF1A']].round(3))\n"
 "print('PROGENy pathway activity by faction:')\n"
 "display(show(f'{T}/decoupler_progeny_by_condition.csv', index_col=0)[['NFkB','Hypoxia','JAK-STAT','TNFa']].round(3))"))
cells.append(code("fig('decoupler_collectri_heatmap.png'); fig('decoupler_progeny_heatmap.png')\n"
 "fig('decoupler_collectri_along_pseudotime.png'); fig('decoupler_progeny_along_pseudotime.png')"))

# Task 4
cells.append(md("## Task 4 — pseudotime gene-trend ordering (Mono1 → TAM)",
 "NR4A2 onset **precedes** the LA-TAM lipid genes (C1Q/APOE/CD9/GPNMB/LIPA), though peaks are",
 "near-simultaneous — a modest, onset-level lead."))
cells.append(code(
 "display(show(f'{T}/genetrend_ordering_summary.csv'))\n"
 "display(show(f'{T}/genetrend_peaktimes.csv'))"))
cells.append(code("fig('genetrends_nr4a_vs_latam.png')"))

# Task 6
cells.append(md("## Task 6 (optional) — TCGA-KIRC external validation (n=531)",
 "Bulk RNA-seq is essentially **silent** on the myeloid-intrinsic NR4A↔LA-TAM link."))
cells.append(code(
 "display(show(f'{T}/tcga_kirc_survival.csv'))\n"
 "display(show(f'{T}/tcga_kirc_latam_nr4a2_corr.csv'))"))
cells.append(code("fig('tcga_kirc_km_latam.png')"))

# Summary
cells.append(md("## Summary → claims",
 "See `outputs/SUMMARY.md` for the full mapping. In brief:",
 "- **(a)** LA-TAM identity — **supported**.",
 "- **(b)** NF-κB conserved regulator — **supported** (NFKB1 robust across SCENIC + decoupleR).",
 "- **(c)** NR4A2/NR4A3 focal tumor-restricted; NFE2L2/HIF1A — **mixed**: NR4A2 is niche-induced",
 "  (VEGFA/hypoxia), early, tumor-associated, but **not** a direct de-novo driver of the LA-TAM",
 "  lipid module; NR4A3 largely silent; NFE2L2/HIF1A regulons elevated in TAM.",
 "- **(d)** ATF3 (PCa) — **untestable here** (absent from this RCC object)."))

nb={"cells":cells,"metadata":{"kernelspec":{"display_name":"Python 3","language":"python","name":"python3"},
    "language_info":{"name":"python","version":"3.10"}},"nbformat":4,"nbformat_minor":5}
with open(NB,"w") as f: json.dump(nb,f,indent=1)
print("wrote", NB, "cells:", len(cells))
