#!/usr/bin/env python

import anndata as ad
import sys, os

paths = [

    "kidney-cancer/Cleaned_Data/combined.h5ad",
    "kidney-cancer/Cleaned_Data/combined-with-stromal.h5ad",
    "kidney-cancer/Cleaned_Data/integrated.h5ad",
    "kidney-cancer/Cleaned_Data/integrated-with-stromal.h5ad",
    "kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad",
    "kidney-cancer/Cleaned_Data/myeloid_FINAL_labels_old.h5ad",

    "prostate-cancer/Cleaned_Data/benign_combined.h5ad",
    "prostate-cancer/Cleaned_Data/combined.h5ad",
    "prostate-cancer/Cleaned_Data/integrated.h5ad",
    "prostate-cancer/Cleaned_Data/integrated_with_kfoury_labels.h5ad",
    "prostate-cancer/Cleaned_Data/kfoury_myeloid.h5ad",
    "prostate-cancer/Cleaned_Data/myeloid_FINAL.h5ad",
    "prostate-cancer/Cleaned_Data/myeloid_FINAL_kfoury_annotated.h5ad",
    "prostate-cancer/Cleaned_Data/myeloid_integrated_final_label.h5ad",
    "prostate-cancer/integrated.h5ad",
]

CHECK_GENES = ["ATF3","C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK","FOLR2","SELENOP","ABCA1","CH25H","C3","SERPING1","CA9","PAX8","CD68","CD163"]

for p in paths:
    print("="*90)
    print("FILE:", p)
    if not os.path.exists(p):
        print("  MISSING")
        continue
    try:
        a = ad.read_h5ad(p, backed="r")
    except Exception as e:
        print("  ERROR reading:", e)
        continue
    print(f"  shape: n_obs={a.n_obs}  n_vars={a.n_vars}")
    print("  obs columns:", list(a.obs.columns))

    for col in a.obs.columns:
        try:
            vc = a.obs[col]
            nun = vc.nunique(dropna=True)
            if nun <= 40:
                counts = vc.value_counts(dropna=False)
                print(f"    [{col}] ({nun} uniq): " + ", ".join(f"{k}={v}" for k,v in counts.items()))
        except Exception as e:
            pass
    print("  obsm:", list(a.obsm.keys()))
    print("  layers:", list(a.layers.keys()))
    print("  raw:", None if a.raw is None else f"raw n_vars={a.raw.n_vars}")
    vn = set(a.var_names)
    present = [g for g in CHECK_GENES if g in vn]
    absent = [g for g in CHECK_GENES if g not in vn]
    print("  genes present:", present)
    print("  genes ABSENT:", absent)
    a.file.close()
