#!/usr/bin/env python
"""00c_stage_braun_icb.py — Phase 0.5: stage the Braun 2020 (CheckMate-009/010/025) ccRCC ICB cohort.

Source (downloaded on IGS login node, internet OK):
  Braun DA et al., Nat Med 2020; 26:909-918. DOI 10.1038/s41591-020-0839-y
  Supplementary Table file MOESM2:
  https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-020-0839-y/MediaObjects/41591_2020_839_MOESM2_ESM.xlsx
  raw/ md5: aea91c06c171f090100eb3c0141e3428 (see raw/CHECKSUMS.md5). RAW LEFT UNTOUCHED.

Parses MOESM2 sheets into processed/ TSVs. Verifies all 8 CLEC_LAM CORE genes are present.
If the raw XLSX is absent, prints the required drop path and exits 2 (does NOT fabricate).

Run: envs/celloracle_env/bin/python (has openpyxl); no internet needed once raw/ is staged.
"""
import os, sys, hashlib
import pandas as pd

BASE = os.path.join(os.path.dirname(__file__), "..", "..", "..",
                    "resources", "external_icb", "braun_checkmate")
BASE = os.path.abspath(BASE)
RAW  = os.path.join(BASE, "raw", "41591_2020_839_MOESM2_ESM.xlsx")
OUT  = os.path.join(BASE, "processed")
os.makedirs(OUT, exist_ok=True)

CORE = ["C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK"]

if not os.path.exists(RAW):
    print(f"[BLOCKED] Braun XLSX not found. Drop it here:\n  {RAW}\n"
          f"Download: https://static-content.springer.com/esm/art%3A10.1038%2Fs41591-020-0839-y/"
          f"MediaObjects/41591_2020_839_MOESM2_ESM.xlsx", file=sys.stderr)
    sys.exit(2)

md5 = hashlib.md5(open(RAW,"rb").read()).hexdigest()
print(f"raw md5: {md5}")

xl = pd.ExcelFile(RAW, engine="openpyxl")
print("sheets:", xl.sheet_names)

# --- S1: clinical + immune (row 1 is a title; real header is row 2 -> skiprows=1) ---
s1 = xl.parse("S1_Clinical_and_Immune_Data", skiprows=1)
s1 = s1.rename(columns={s1.columns[0]: "SUBJID"})
s1.to_csv(os.path.join(OUT, "clinical.tsv"), sep="\t", index=False)

treatment = s1[["SUBJID","Cohort","Arm"]].copy()
treatment.to_csv(os.path.join(OUT, "treatment.tsv"), sep="\t", index=False)

resp_cols = ["SUBJID","Tumor_Shrinkage","ORR","Benefit","ExtremeResponder","irORR"]
s1[[c for c in resp_cols if c in s1.columns]].to_csv(
    os.path.join(OUT, "response.tsv"), sep="\t", index=False)

surv_cols = ["SUBJID","OS","OS_CNSR","PFS","PFS_CNSR","irPFS","irPFS_CNSR"]
s1[[c for c in surv_cols if c in s1.columns]].to_csv(
    os.path.join(OUT, "survival.tsv"), sep="\t", index=False)

# adjustment covariates that Phase 4 will use (kept as an explicit convenience table)
covar_cols = ["SUBJID","Cohort","Arm","Age","Sex","MSKCC","IMDC","Sarc_or_Rhab",
              "Received_Prior_Therapy","Number_of_Prior_Therapies","SampleType",
              "Tumor_Sample_Primary_or_Metastasis","Site_of_Metastasis",
              "Purity","Ploidy","TMB_Counts","WGII","ImmunoPhenotype","TM_CD8_Density",
              "PBRM1","Deletion_9p21.3","Angio","Teff","Myeloid","Javelin","Merck18"]
s1[[c for c in covar_cols if c in s1.columns]].to_csv(
    os.path.join(OUT, "covariates.tsv"), sep="\t", index=False)

# --- S4A: normalized RNA expression (gene x sample). Row1 title -> skiprows=1 ---
expr = xl.parse("S4A_RNA_Expression", skiprows=1)
expr = expr.rename(columns={expr.columns[0]: "gene_name"}).set_index("gene_name")
expr.to_csv(os.path.join(OUT, "expression_normalized.tsv"), sep="\t")

# --- S4C: CIBERSORTx immune deconvolution ---
deconv = xl.parse("S4C_CIBERSORTx", skiprows=1)
deconv = deconv.rename(columns={deconv.columns[0]: "cell_type"}).set_index("cell_type")
deconv.to_csv(os.path.join(OUT, "immune_deconv.tsv"), sep="\t")

# --- sample_manifest: link expression columns (RNA_ID) <-> SUBJID/clinical ---
expr_samples = list(expr.columns)
rna_map = s1.set_index("RNA_ID")[["SUBJID","Cohort","Arm"]] if "RNA_ID" in s1.columns else None
rows = []
for e in expr_samples:
    if rna_map is not None and e in rna_map.index:
        r = rna_map.loc[e]
        r = r.iloc[0] if hasattr(r, "iloc") and getattr(r, "ndim", 1) > 1 else r
        rows.append({"RNA_ID": e, "SUBJID": r["SUBJID"], "Cohort": r["Cohort"], "Arm": r["Arm"], "has_rna": True})
    else:
        rows.append({"RNA_ID": e, "SUBJID": None, "Cohort": None, "Arm": None, "has_rna": True})
man = pd.DataFrame(rows)
man.to_csv(os.path.join(OUT, "sample_manifest.tsv"), sep="\t", index=False)

# --- verification report ---
present = [g for g in CORE if g in expr.index]
absent  = [g for g in CORE if g not in expr.index]
n_link = int(man["SUBJID"].notna().sum())
print("\n=== STAGING REPORT ===")
print(f"clinical:   {s1.shape[0]} subjects x {s1.shape[1]} cols")
print(f"expression: {expr.shape[0]} genes x {expr.shape[1]} RNA samples")
print(f"deconv:     {deconv.shape[0]} cell types x {deconv.shape[1]} samples")
print(f"manifest:   {man.shape[0]} RNA samples, {n_link} linked to a SUBJID")
print(f"cohorts:    {s1['Cohort'].value_counts().to_dict()}")
print(f"arms:       {s1['Arm'].value_counts().to_dict()}")
print(f"CORE present ({len(present)}/8): {present}")
print(f"CORE ABSENT: {absent}")
for ep in ["OS","OS_CNSR","PFS","PFS_CNSR","ORR","Benefit"]:
    if ep in s1.columns:
        nn = int(s1[ep].notna().sum())
        print(f"  endpoint {ep}: {nn} non-null")
print("wrote:", sorted(os.listdir(OUT)))
