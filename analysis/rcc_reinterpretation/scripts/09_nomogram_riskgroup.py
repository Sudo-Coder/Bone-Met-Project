import os
import json
import numpy as np
import pandas as pd

root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
out = os.path.join(root, 'analysis/rcc_reinterpretation/outputs/model')
np.random.seed(0)

model = json.load(open(os.path.join(out, 'risk_model.json')))
genes = np.array(model['genes']); mm = np.array(model['mean']); ss = np.array(model['sd']); cc = np.array(model['coef'])

kirc = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_HiSeqV2.gz'), sep = '\t', index_col = 0)
kirc.index = [g.upper() for g in kirc.index]
surv = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_survival.txt'), sep = '\t')
surv['bcr'] = surv['sample'].str[:15]
cl = pd.read_csv(os.path.join(root, 'resources/model_data/kirc_clinical/KIRC_clinicalMatrix'), sep = '\t', index_col = 0)

x = kirc.loc[[g for g in genes if g in kirc.index]].T
z = (x[list(genes)] - mm) / ss
risk = pd.Series(z.values @ cc, index = x.index)

df = pd.DataFrame({'risk': risk})
df['bcr'] = [s[:15] for s in df.index]
df = df.merge(surv[['bcr','OS','OS.time']], on = 'bcr')
df = df[df['OS.time'] > 0].dropna(subset = ['OS','OS.time'])

cc2 = cl[['pathologic_stage','neoplasm_histologic_grade','age_at_initial_pathologic_diagnosis']].copy()
cc2['bcr'] = cc2.index

def stage_num(x):
    m = {'Stage I': 1, 'Stage II': 2, 'Stage III': 3, 'Stage IV': 4}
    return m.get(x, np.nan)

def grade_num(x):
    if isinstance(x, str) and x.startswith('G') and x[1:].isdigit():
        return int(x[1:])
    return np.nan

cc2['stage'] = cc2['pathologic_stage'].map(stage_num)
cc2['grade'] = cc2['neoplasm_histologic_grade'].map(grade_num)
cc2['age'] = pd.to_numeric(cc2['age_at_initial_pathologic_diagnosis'], errors = 'coerce')
df = df.merge(cc2[['bcr','stage','grade','age']], on = 'bcr', how = 'left')

from lifelines import CoxPHFitter
mv = df[['risk','stage','grade','age','OS.time','OS']].dropna().copy()
for c in ['risk','stage','grade','age']:
    mv[c + '_z'] = (mv[c] - mv[c].mean()) / mv[c].std()
cox = CoxPHFitter().fit(mv[['risk_z','stage_z','grade_z','age_z','OS.time','OS']], 'OS.time', 'OS')
cox.summary[['coef','exp(coef)','p']].to_csv(os.path.join(out, 'multivariable_cox.csv'))
print(cox.summary[['exp(coef)','coef lower 95%','coef upper 95%','p']].round(4))

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

terms = ['risk','stage','grade','age']
beta = {t: cox.params_[t + '_z'] for t in terms}
ranges = {}
for t in terms:
    v = mv[t]
    ranges[t] = (v.min(), v.max(), (v - mv[t].mean()) / mv[t].std())
pts = {}
maxeff = max(abs(beta[t]) * ((ranges[t][1] - mv[t].mean()) / mv[t].std() - (ranges[t][0] - mv[t].mean()) / mv[t].std()) for t in terms)
plt.figure(figsize = (9,5))
yy = np.arange(len(terms))[::-1]
for i, t in enumerate(terms):
    lo_z = (ranges[t][0] - mv[t].mean()) / mv[t].std()
    hi_z = (ranges[t][1] - mv[t].mean()) / mv[t].std()
    p_lo = beta[t] * lo_z / maxeff * 100
    p_hi = beta[t] * hi_z / maxeff * 100
    plt.plot([min(p_lo, p_hi), max(p_lo, p_hi)], [yy[i], yy[i]], marker = '|', ms = 15, lw = 2)
    plt.text(-5, yy[i], t, ha = 'right', va = 'center')
plt.yticks([])
plt.xlabel('points')
plt.title('ccRCC OS nomogram (multivariable Cox: risk + stage + grade + age)')
plt.tight_layout(); plt.savefig(os.path.join(out, 'nomogram.png'), dpi = 200); plt.close()

from lifelines import KaplanMeierFitter
horizons = [365, 1095, 1825]
plt.figure(figsize = (6,6))
for t in horizons:
    sf = cox.predict_survival_function(mv, times = [t]).T[t]
    q = pd.qcut(sf, 4, labels = False, duplicates = 'drop')
    obs, pred = [], []
    kmf = KaplanMeierFitter()
    for g in sorted(pd.unique(q)):
        m = q == g
        kmf.fit(mv['OS.time'][m], mv['OS'][m])
        obs.append(kmf.predict(t))
        pred.append(sf[m].mean())
    plt.plot(pred, obs, marker = 'o', label = str(round(t / 365)) + '-yr')
plt.plot([0,1],[0,1],'k--', lw = 1)
plt.xlabel('predicted survival'); plt.ylabel('observed (KM)'); plt.legend(); plt.title('Calibration')
plt.tight_layout(); plt.savefig(os.path.join(out, 'calibration.png'), dpi = 200); plt.close()

t = 1825
p_death = 1 - cox.predict_survival_function(mv, times = [t]).T[t].values
time = mv['OS.time'].values; ev = mv['OS'].values
n = len(mv)
thr = np.linspace(0.01, 0.6, 60)
nb_model, nb_all = [], []
for pt in thr:
    pos = p_death >= pt
    tp = np.sum(pos & (ev == 1) & (time <= t)) / n
    fp = np.sum(pos & ~((ev == 1) & (time <= t))) / n
    nb_model.append(tp - fp * (pt / (1 - pt)))
    ev_rate = np.mean((ev == 1) & (time <= t))
    nb_all.append(ev_rate - (1 - ev_rate) * (pt / (1 - pt)))
plt.figure(figsize = (7,5))
plt.plot(thr, nb_model, label = 'nomogram')
plt.plot(thr, nb_all, label = 'treat all', ls = '--')
plt.axhline(0, color = 'k', lw = 0.8, label = 'treat none')
plt.ylim(-0.05, max(nb_model) + 0.02); plt.xlabel('threshold probability'); plt.ylabel('net benefit')
plt.legend(); plt.title('Decision curve analysis (5-yr OS)')
plt.tight_layout(); plt.savefig(os.path.join(out, 'dca.png'), dpi = 200); plt.close()

import anndata as ad
import decoupler as dc
imm = {'CD8_Tcell': ['CD8A','CD8B','GZMK','CD3D'], 'cytotoxicity': ['GZMB','PRF1','NKG7','GNLY','IFNG'],
       'exhaustion': ['TOX','PDCD1','HAVCR2','LAG3','TIGIT','CTLA4','ENTPD1'], 'Treg': ['FOXP3','IL2RA','IKZF2'],
       'MHC_I': ['HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2'], 'MHC_II': ['HLA-DRA','HLA-DRB1','CD74','CIITA'],
       'panTAM': ['CD68','CD163','MRC1','CSF1R','LYZ'], 'immune_infiltration': ['PTPRC','CD3D','CD2','CD53','LCK']}
kx = kirc.T
a = ad.AnnData(kx.values.astype(float), obs = pd.DataFrame(index = kx.index), var = pd.DataFrame(index = kx.columns))
net = pd.concat([pd.DataFrame({'source': k, 'target': [g for g in v if g in a.var_names], 'weight': 1.0}) for k, v in imm.items()])
dc.mt.aucell(a, net, tmin = 2, verbose = False)
isc = a.obsm['score_aucell']; isc['bcr'] = [s[:15] for s in kx.index]
rg = df[['bcr','risk']].copy()
rg['grp'] = np.where(rg['risk'] >= rg['risk'].median(), 'high', 'low')
isc = isc.merge(rg[['bcr','grp']], on = 'bcr')

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
rows = []
for k in imm:
    hi = isc[isc.grp == 'high'][k]; lo = isc[isc.grp == 'low'][k]
    u, p = mannwhitneyu(hi, lo)
    rows.append([k, hi.mean(), lo.mean(), hi.mean() - lo.mean(), p])
rgi = pd.DataFrame(rows, columns = ['signature','high_mean','low_mean','delta','p'])
rgi['fdr'] = multipletests(rgi['p'], method = 'fdr_bh')[1]
rgi.to_csv(os.path.join(out, 'riskgroup_immune.csv'), index = False)
print(rgi.round(4).to_string(index = False))

melt = isc.melt(id_vars = 'grp', value_vars = list(imm), var_name = 'signature', value_name = 'score')
plt.figure(figsize = (10,4))
import seaborn as sns
sns.boxplot(data = melt, x = 'signature', y = 'score', hue = 'grp', showfliers = False)
plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')
plt.title('Immune infiltration by ccRCC risk group')
plt.tight_layout(); plt.savefig(os.path.join(out, 'riskgroup_immune.png'), dpi = 200); plt.close()
