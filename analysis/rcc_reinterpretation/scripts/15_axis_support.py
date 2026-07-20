#!/usr/bin/env python
import os, ast, warnings
warnings.filterwarnings("ignore")
import numpy as np, pandas as pd
ROOT=os.path.abspath(os.path.join(os.path.dirname(__file__),"..","..",".."))
TAB=os.path.join(ROOT,"analysis","rcc_reinterpretation","outputs","tables")
def rd(f):
    p=os.path.join(TAB,f); return pd.read_csv(p) if os.path.exists(p) else None
lr=rd("axis_lr_expression.csv"); cc=rd("cellchat_tam_LR_tumor.csv"); nn=rd("nichenet_prespecified_ranks.csv")

AX={"complement_C1QA":("C1QA","COMPLEMENT"),"complement_C1QB":("C1QB","COMPLEMENT"),
    "complement_C1QC":("C1QC","COMPLEMENT"),"complement_C3":("C3","COMPLEMENT"),
    "APOE_TREM2":("APOE","ApoE"),"GAS6_MERTK":("GAS6","GAS"),"PROS1_MERTK":("PROS1","PROS"),
    "TGFB1_TGFBR":("TGFB1","TGFb")}
rows=[]
for axis,(lig,pw) in AX.items():
    r={"axis":axis,"ligand":lig}

    if lr is not None and axis in set(lr.axis):
        a=lr[lr.axis==axis].iloc[0]
        r["lig_delta_TvsB"]=a["ligand_delta_TvsB"]; r["lig_p"]=a["ligand_p"]
        r["A_tumor_elevated"]=bool(a["ligand_delta_TvsB"]>0 and a["ligand_p"]<0.05)
        try: rec=ast.literal_eval(a["receptor_detect_in_tumorTAM"]); r["receptor_max_detect"]=max(rec.values())
        except Exception: r["receptor_max_detect"]=np.nan

    if cc is not None and pw is not None and "pathway_name" in cc.columns:
        sub=cc[(cc.pathway_name==pw)]
        r["B_cellchat_toTAM_pairs"]=int(len(sub)); r["B_cellchat_maxprob"]=float(sub["prob"].max()) if len(sub) else 0.0
        r["B_cellchat"]=bool(len(sub)>0)
    else:
        r["B_cellchat"]=None if pw is None else (False if cc is not None else None)

    if nn is not None and lig in set(nn.test_ligand):
        nr=nn[nn.test_ligand==lig].iloc[0]; r["C_nichenet_rank"]=int(nr["rank"]); r["C_nichenet_aupr"]=float(nr["aupr_corrected"])
        r["C_top20"]=bool(nr["rank"]<=20)
    rows.append(r)
df=pd.DataFrame(rows)
def verdict(x):

    A=x.get("A_tumor_elevated"); B=x.get("B_cellchat"); C=x.get("C_top20")
    delta=x.get("lig_delta_TvsB",0) or 0; lp=x.get("lig_p",1)
    present = (B is True) or (C is True)
    benign_biased = (A is False) and (delta < 0) and (lp is not None and lp < 0.05)
    if A is True and present:  return "SUPPORTED_TUMOR_SPECIFIC"
    if A is True:              return "AMBIGUOUS"
    if benign_biased:          return "NOT_TUMOR_GAINED"
    if present:                return "PRESENT_NOT_TUMOR_GAINED"
    return "UNSUPPORTED"
df["verdict"]=df.apply(verdict,axis=1)
df["lines_available"]=f"A={lr is not None},B={cc is not None},C={nn is not None}"
df.to_csv(os.path.join(TAB,"axis_support_table.csv"),index=False)
pd.set_option("display.width",240,"display.max_colwidth",22)
print(df.to_string(index=False))
print("\nlines: LR-expr(A)=%s CellChat(B)=%s NicheNet(C)=%s"%(lr is not None,cc is not None,nn is not None))
print("wrote axis_support_table.csv")
