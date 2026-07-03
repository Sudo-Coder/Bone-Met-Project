#!/usr/bin/env python
"""Task 3/4 bridge — overlay decoupleR pathway/TF activity along the Mono1->TAM pseudotime.
Reuses saved decoupleR per-cell matrices + Task-4 pseudotime (no recompute).
Run (after task 4): envs/mechanism_env/bin/python analysis/09_pseudotime_overlay.py
"""
import os, numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
ROOT="/autofs/projects-t3/hussain/scProj"
PROJ = ROOT + "/kidney-cancer"   # project data/outputs (relocated 2026-07-02)
DEC=f"{PROJ}/outputs/decoupler"; CR=f"{PROJ}/outputs/cellrank"; FIG=f"{PROJ}/outputs/figures"
pt=pd.read_csv(f"{CR}/pseudotime.csv", index_col=0)
ps=pt["dpt_pseudotime"]
for tag,cols in [("progeny",["NFkB","Hypoxia","JAK-STAT","TNFa"]),
                 ("collectri",["NR4A2","NR4A3","NFKB1","HIF1A"])]:
    act=pd.read_parquet(f"{DEC}/{tag}_percell.parquet")
    cols=[c for c in cols if c in act.columns]
    j=act[cols].join(ps,how="inner").dropna(subset=["dpt_pseudotime"]).sort_values("dpt_pseudotime")
    plt.figure(figsize=(7,4))
    for c in cols:
        y=j[c].rolling(300,min_periods=30,center=True).mean()
        plt.plot(j["dpt_pseudotime"],y,label=c,lw=2)
    plt.legend(); plt.xlabel("DPT pseudotime (Mono1->TAM)"); plt.ylabel("activity (rolling mean)")
    plt.title(f"{tag} activity along pseudotime")
    plt.tight_layout(); plt.savefig(f"{FIG}/decoupler_{tag}_along_pseudotime.png",dpi=150); plt.close()
    print("wrote overlay", tag, flush=True)
print("OVERLAY DONE", flush=True)
