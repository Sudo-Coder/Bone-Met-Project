#!/usr/bin/env python
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
TAB = os.path.join(ROOT, "analysis", "rcc_reinterpretation", "outputs", "tables")
p3 = pd.read_csv(os.path.join(TAB, "phase3_adjusted_associations.csv"))
h2h = pd.read_csv(os.path.join(TAB, "phase3_headtohead_complement_vs_panTAM.csv"))
tcga = pd.read_csv(os.path.join(TAB, "phase4_tcga_kirc_cox.csv"))
braun = pd.read_csv(os.path.join(TAB, "phase4_braun_nivo.csv"))
MODS = ["complement_C1Q", "CLEC_LAM8", "RCC_skew_CORE", "APOE_TREM2", "MERTK_GPNMB", "panTAM", "SPP1_TAM", "inflammatory_mono"]
def sig(df, **k):
    d = df.copy()
    for c, v in k.items(): d = d[d[c]==v]
    if len(d)==0 or pd.isna(d.iloc[0].get("p")): return "—"
    r = d.iloc[0]
    est = r.get("coef", r.get("HR_per_SD", np.nan))
    pv = r["p"]
    arrow = "↑" if est>(0 if "coef" in d.columns else 1) else "↓"
    star = "*" if pv<0.05 else ("~" if pv<0.10 else "")
    val = f"{est:.2f}{arrow}{star}(p{pv:.2f})"
    return val
rows = []
for m in MODS:
    r = {"module": m}
    r["CD8_exhaustion"] = sig(p3, outcome = "CD8_exhaustion", module = m)
    r["cytotoxicity"] = sig(p3, outcome = "cytotoxicity", module = m)
    r["Treg_fraction"] = sig(p3, outcome = "Treg_fraction", module = m)
    r["MHC_I_APM(n9)"] = sig(p3, outcome = "MHC_I_APM", module = m)
    r["MHC_II_APC"] = sig(p3, outcome = "MHC_II_APC", module = m)
    r["TCGA_OS_adj"] = sig(tcga, endpoint = "OS", module = m, model = "adj(panTAM+Obradovic+immune)") if m not in("panTAM", "SPP1_TAM", "inflammatory_mono") else sig(tcga, endpoint = "OS", module = m, model = "unadj")
    r["TCGA_PFI_adj"] = sig(tcga, endpoint = "PFI", module = m, model = "adj(panTAM+Obradovic+immune)") if m not in("panTAM", "SPP1_TAM", "inflammatory_mono") else sig(tcga, endpoint = "PFI", module = m, model = "unadj")
    r["Braun_nivo_OS_adj"] = sig(braun, endpoint = "OS", module = m, model = "adj(Myeloid+Purity)")
    r["Braun_nivo_ORR"] = sig(braun, endpoint = "responder", module = m, model = "logit(OR/SD)")
    rows.append(r)
V = pd.DataFrame(rows)

def beats(m):

    if m=="complement_C1Q":
        hh = h2h.set_index("outcome")
        wins = []
        for oc in ["CD8_exhaustion", "cytotoxicity", "MHC_II_APC"]:
            if oc in hh.index and hh.loc[oc, "complement_p"]<0.05 and hh.loc[oc, "panTAM_p"]>=0.05: wins.append(oc)
        return "yes ("+",".join(wins)+")" if wins else "partial"
    return ""
def verdict(m):
    if m in ("panTAM", "inflammatory_mono"): return "generic control (mostly null)"
    if m=="SPP1_TAM": return "broad-TAM (Treg/exhaustion, not specific)"
    if m=="complement_C1Q": return "SUPPORTED-specific (cytotox/exhaustion beat panTAM; TCGA-OS survives Obradovic)"
    if m=="CLEC_LAM8": return "supported-broad (TCGA-OS collapses vs Obradovic; ICB arm-INT exploratory)"
    if m=="RCC_skew_CORE": return "supported (TCGA-OS survives Obradovic; ~complement)"
    if m=="APOE_TREM2": return "= Obradovic (collapses); not independent"
    if m=="MERTK_GPNMB": return "shared-support (Treg via broad; OS flips protective when adj)"
    return "ambiguous"
V["beats_generic_controls"] = [beats(m) for m in V.module]
V["verdict"] = [verdict(m) for m in V.module]
V.to_csv(os.path.join(TAB, "phase34_verdict_table.csv"), index = False)
pd.set_option("display.width", 280, "display.max_colwidth", 40)
print(V.to_string(index = False))
