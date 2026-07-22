#!/usr/bin/env python
# Does the survival signal survive adjustment for the standard clinical prognosticators?
# Scoring is identical to 40_tcga_kirc.py; the only addition is stage / grade / age.
#   stage: GDC clinical (ajcc_pathologic_stage), collapsed to I/II vs III/IV
#   grade: cBioPortal kirc_tcga GRADE (Fuhrman) -- GDC exports it as "Not Reported"
#   age:   GDC age_at_index
import os, json, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd, anndata as ad, decoupler as dc
from lifelines import CoxPHFitter
np.random.seed(0)
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
TAB = os.path.join(ROOT, "analysis", "rcc_reinterpretation", "outputs", "tables")
MOD = os.path.join(ROOT, "analysis", "rcc_reinterpretation", "outputs", "model")
RES = os.path.join(ROOT, "resources", "tcga")

MODULES = {"complement_C1Q": ["C1QA", "C1QB", "C1QC"],
           "CLEC_LAM8": ["C1QA", "C1QB", "C1QC", "APOE", "APOC1", "TREM2", "GPNMB", "MERTK"],
           "RCC_skew_CORE": ["C1QA", "C1QB", "C1QC", "APOE", "APOC1", "TREM2"],
           "APOE_TREM2": ["APOE", "APOC1", "TREM2", "TYROBP"],
           "MERTK_GPNMB": ["MERTK", "GPNMB"],
           "panTAM": ["CD68", "CD163", "MRC1", "CSF1R", "LYZ", "AIF1", "FCGR3A"],
           "Obradovic_TREM2": ["TREM2", "APOE", "APOC1", "C1QA", "C1QB", "C1QC", "GPNMB", "FOLR2", "SPP1", "CTSD", "CD68"],
           "total_immune": ["PTPRC", "CD3D", "CD8A", "CD4", "CD19", "NKG7", "CD68"]}

expr = pd.read_csv(os.path.join(RES, "KIRC_HiSeqV2.gz"), sep = "\t", index_col = 0)
expr.index = [str(g).upper() for g in expr.index]
samples = [c for c in expr.columns if c[-2:] in ("01", "05")]
X = expr[samples].T
a = ad.AnnData(X.values.astype(float), obs = pd.DataFrame(index = X.index), var = pd.DataFrame(index = X.columns))
net = pd.concat([pd.DataFrame({"source": k, "target": [g for g in v if g in a.var_names], "weight": 1.0})
                 for k, v in MODULES.items()])
dc.mt.aucell(a, net, tmin = 2, verbose = False)
S = a.obsm["score_aucell"]
sc = pd.DataFrame({k: S[k].values for k in MODULES if k in S.columns}, index = X.index)

# frozen 25-gene risk score, same standardisation as the packaged model
rm = json.load(open(os.path.join(MOD, "risk_model.json")))
g, mu, sd, co = np.array(rm["genes"]), np.array(rm["mean"]), np.array(rm["sd"]), np.array(rm["coef"])
have = [i for i, q in enumerate(g) if q in X.columns]
sc["risk25"] = ((X[list(g[have])] - mu[have]) / sd[have]).values @ co[have]

sv = pd.read_csv(os.path.join(RES, "KIRC_survival.txt"), sep = "\t")
sv["bcr"] = sv["sample"].str[:15]
sc["bcr"] = [s[:15] for s in sc.index]
sc["patient"] = [s[:12] for s in sc.index]
D = sc.merge(sv, on = "bcr", how = "inner")

cl = pd.read_csv(os.path.join(RES, "TCGA-KIRC.clinical.tsv"), sep = "\t", low_memory = False)
cl = cl[["submitter_id", "ajcc_pathologic_stage.diagnoses", "age_at_index.demographic"]].dropna(subset = ["submitter_id"])
cl.columns = ["patient", "stage_raw", "age"]
cl = cl.drop_duplicates("patient")
gr = pd.DataFrame(json.load(open(os.path.join(RES, "grade.json"))))[["patientId", "value"]]
gr.columns = ["patient", "grade_raw"]
D = D.merge(cl, on = "patient", how = "left").merge(gr, on = "patient", how = "left")

D["stage_hi"] = D["stage_raw"].map({"Stage I": 0, "Stage II": 0, "Stage III": 1, "Stage IV": 1})
D["grade_hi"] = D["grade_raw"].map({"G1": 0, "G2": 0, "G3": 1, "G4": 1})
D["age"] = pd.to_numeric(D["age"], errors = "coerce")
print("merged n=%d | stage %d | grade %d | age %d"
      % (len(D), D.stage_hi.notna().sum(), D.grade_hi.notna().sum(), D.age.notna().sum()))

def zc(x):
    x = pd.to_numeric(x, errors = "coerce")
    return (x - x.mean()) / x.std()

TEST = ["complement_C1Q", "RCC_skew_CORE", "CLEC_LAM8", "APOE_TREM2", "MERTK_GPNMB",
        "panTAM", "Obradovic_TREM2", "risk25"]
MODELS = {
    "unadjusted":                       [],
    "molecular (panTAM+Obradovic+imm)":  ["z_pan", "z_obr", "z_imm"],
    "clinical (stage+grade+age)":        ["stage_hi", "grade_hi", "age_z"],
    "full (molecular+clinical)":         ["z_pan", "z_obr", "z_imm", "stage_hi", "grade_hi", "age_z"],
}
rows = []
for m in TEST:
    for name, cov in MODELS.items():
        d = D.dropna(subset = ["OS", "OS.time", m]).copy()
        d = d[d["OS.time"] > 0]
        d["z"] = zc(d[m])
        d["z_pan"] = zc(d["panTAM"])
        d["z_obr"] = zc(d["Obradovic_TREM2"])
        d["z_imm"] = zc(d["total_immune"])
        d["age_z"] = zc(d["age"])
        use = [c for c in cov if not (m in ("panTAM", "Obradovic_TREM2")
                                      and c in ("z_pan", "z_obr"))]
        d = d.dropna(subset = use) if use else d
        try:
            fit = CoxPHFitter().fit(d[["z"] + use + ["OS.time", "OS"]]
                                    .rename(columns = {"OS.time": "T", "OS": "E"}), "T", "E")
            r = fit.summary.loc["z"]
            rows.append(dict(module = m, model = name, HR = np.exp(r["coef"]),
                             ci_low = np.exp(r["coef lower 95%"]), ci_high = np.exp(r["coef upper 95%"]),
                             p = r["p"], n = len(d), events = int(d["OS"].sum())))
        except Exception as e:
            rows.append(dict(module = m, model = name, HR = np.nan, ci_low = np.nan,
                             ci_high = np.nan, p = np.nan, n = len(d), events = 0, note = str(e)[:40]))
R = pd.DataFrame(rows)
R["verdict"] = np.where(R.p.isna(), "n/a", np.where(R.p < 0.05, "survives", "collapses"))
R.to_csv(os.path.join(TAB, "phase4_tcga_clinical_adjusted.csv"), index = False)
pd.set_option("display.width", 220)
print("\n=== TCGA-KIRC overall survival, HR per SD ===")
print(R[["module", "model", "HR", "ci_low", "ci_high", "p", "n", "events", "verdict"]].round(3).to_string(index = False))

# what the clinical variables themselves do, for reference
d = D.dropna(subset = ["OS", "OS.time", "stage_hi", "grade_hi", "age"]).copy()
d = d[d["OS.time"] > 0]
d["age_z"] = zc(d["age"])
fit = CoxPHFitter().fit(d[["stage_hi", "grade_hi", "age_z", "OS.time", "OS"]]
                        .rename(columns = {"OS.time": "T", "OS": "E"}), "T", "E")
print("\n=== clinical variables alone (n=%d, %d events) ===" % (len(d), int(d.OS.sum())))
print(fit.summary[["exp(coef)", "exp(coef) lower 95%", "exp(coef) upper 95%", "p"]].round(3).to_string())
