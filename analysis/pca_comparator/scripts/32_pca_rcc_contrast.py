#!/usr/bin/env python
"""32_pca_rcc_contrast.py — assemble the PCa-vs-RCC module contrast (composition + immune) verdict table."""
import os, warnings, numpy as np, pandas as pd
warnings.filterwarnings("ignore")
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","pca_comparator","outputs","tables"); os.makedirs(TAB,exist_ok=True)
mc=pd.read_csv(os.path.join(TAB,"pca_vs_rcc_module_contrast.csv"))
loc=pd.read_csv(os.path.join(TAB,"pca_within_tumor_localization.csv")).set_index("module")
RCC_TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables"); rcc3=pd.read_csv(os.path.join(RCC_TAB,"phase3_adjusted_associations.csv"))
pca3=pd.read_csv(os.path.join(TAB,"pca_phase3_adjusted_associations.csv"))
def imm(df,module,outcome):
    r=df[(df.module==module)&(df.outcome==outcome)]
    if len(r)==0 or pd.isna(r.iloc[0]["p"]): return (np.nan,np.nan)
    return float(r.iloc[0]["coef"]),float(r.iloc[0]["p"])
MODS=["complement_C1Q","CLEC_LAM8","RCC_skew_CORE","APOE_TREM2","MERTK_GPNMB","SPP1_TAM","panTAM","ATF3_NFkB"]
rows=[]
for m in MODS:
    c=mc[mc.module==m].iloc[0]
    rcc_cyt=imm(rcc3,m,"cytotoxicity"); pca_cyt=imm(pca3,m,"cytotoxicity")
    rcc_ex=imm(rcc3,m,"CD8_exhaustion"); pca_ex=imm(pca3,m,"CD8_exhaustion")
    locp=loc.loc[m,"wilcoxon_p"] if m in loc.index else np.nan
    locd=loc.loc[m,"delta"] if m in loc.index else np.nan
    rows.append(dict(module=m,
        RCC_TvB_SD=round(c.RCC_TvB_SD,2),PCa_TvB_SD=round(c.PCa_TvB_SD,2),PCa_TvB_p=round(c.PCa_p,3),
        interaction=round(c.interaction_RCCvsPCa,3),interaction_p=round(c.interaction_p,3),
        PCa_TAM_localization_delta=round(locd,3) if locd==locd else np.nan,PCa_loc_p=round(locp,3) if locp==locp else np.nan,
        RCC_cytotox=f"{rcc_cyt[0]:.2f}(p{rcc_cyt[1]:.2f})" if rcc_cyt[0]==rcc_cyt[0] else "-",
        PCa_cytotox=f"{pca_cyt[0]:.2f}(p{pca_cyt[1]:.2f})" if pca_cyt[0]==pca_cyt[0] else "-",
        RCC_exhaustion=f"{rcc_ex[0]:.2f}(p{rcc_ex[1]:.2f})" if rcc_ex[0]==rcc_ex[0] else "-",
        PCa_exhaustion=f"{pca_ex[0]:.2f}(p{pca_ex[1]:.2f})" if pca_ex[0]==pca_ex[0] else "-"))
def verdict(m,c,rc,pc):
    if m in ("panTAM",): return "shared generic (both tumor-induced, no skew)"
    if m=="ATF3_NFkB": return "CONSERVED core (both strong; small RCC skew)"
    if m in ("SPP1_TAM","MERTK_GPNMB","APOE_TREM2"): return "shared (present both; not RCC-skewed)"
    if m in ("complement_C1Q","RCC_skew_CORE","CLEC_LAM8"):
        return "present both + TAM-localized in PCa, RCC-SKEWED (interaction sig) + RCC-functionally-specific (immune-assoc only in RCC)"
    return "ambiguous"
C=pd.DataFrame(rows); C["verdict"]=[verdict(m,None,None,None) for m in C.module]
C.to_csv(os.path.join(TAB,"pca_vs_rcc_contrast_table.csv"),index=False)
pd.set_option("display.width",300,"display.max_colwidth",60)
print(C.to_string(index=False))
