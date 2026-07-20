import os
import json
import numpy as np
import pandas as pd

root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
out = os.path.join(root, 'analysis/rcc_reinterpretation/outputs/model')
np.random.seed(0)

model = json.load(open(os.path.join(out, 'risk_model.json')))
genes = np.array(model['genes']); mu = np.array(model['mean']); sd = np.array(model['sd']); coef = np.array(model['coef'])

kirc = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_HiSeqV2.gz'), sep = '\t', index_col = 0)
kirc.index = [g.upper() for g in kirc.index]
surv = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_survival.txt'), sep = '\t')
surv['bcr'] = surv['sample'].str[:15]

x = kirc.loc[[g for g in genes if g in kirc.index]].T
z = (x[list(genes)] - mu) / sd
risk = pd.Series(z.values @ coef, index = x.index)
df = pd.DataFrame({'risk': risk})
df['bcr'] = [s[:15] for s in df.index]
df = df.merge(surv[['bcr','OS','OS.time']], on = 'bcr')
df = df[df['OS.time'] > 0].dropna(subset = ['OS','OS.time'])
df['grp'] = np.where(df['risk'] >= df['risk'].median(), 'high', 'low')
df['grp'].value_counts()

import gseapy as gp

sig_genes = list(genes)
cand = pd.read_csv(os.path.join(out, 'candidate_genes.csv'))
deg_genes = cand[cand['source'] == 'TAM_DEG']['gene'].tolist()

def enrich(gene_list, tag):
    res = gp.enrichr(gene_list = gene_list, gene_sets = ['GO_Biological_Process_2021','KEGG_2021_Human'], outdir = None)
    r = res.results.copy()
    r = r[['Gene_set','Term','Overlap','P-value','Adjusted P-value','Genes']]
    r = r.sort_values('Adjusted P-value')
    r.to_csv(os.path.join(out, 'enrichment_' + tag + '.csv'), index = False)
    return r

en_sig = enrich(sig_genes, 'signature')
en_deg = enrich(deg_genes, 'degs')
en_sig.head(10)

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

fig, ax = plt.subplots(1, 2, figsize = (13,5))
for i, (r, title) in enumerate([(en_sig, '25-gene signature'), (en_deg, 'TAM-vs-monocyte DEGs')]):
    top = r.head(10).copy()
    top['nlp'] = -np.log10(top['Adjusted P-value'])
    top['Term'] = top['Term'].str.replace(r' \(GO:\d+\)', '', regex = True).str[:45]
    sns.barplot(data = top, y = 'Term', x = 'nlp', ax = ax[i], color = '#4C72B0')
    ax[i].set_xlabel('-log10 adjusted p'); ax[i].set_ylabel(''); ax[i].set_title(title)
plt.tight_layout()
plt.savefig(os.path.join(out, 'enrichment.png'), dpi = 200)
plt.close()

print('risk groups: high', (df.grp=='high').sum(), 'low', (df.grp=='low').sum())
print('top signature term:', en_sig.iloc[0]['Term'], en_sig.iloc[0]['Adjusted P-value'])
print('top DEG term:', en_deg.iloc[0]['Term'], en_deg.iloc[0]['Adjusted P-value'])

import cptac
cc = cptac.Ccrcc()
tx = cc.get_transcriptomics('broad')
tx.columns = tx.columns.get_level_values(0)
tx = tx.groupby(level = 0, axis = 1).mean()
cl = cc.get_clinical('mssm')

have = [g for g in genes if g in tx.columns]
e = tx.loc[[i for i in tx.index if i in cl.index], have]
ze = (e - e.mean()) / e.std()
w = np.array([coef[list(genes).index(g)] for g in have])
ext = pd.DataFrame({'risk': ze.values @ w}, index = e.index)
ext['event'] = cl.loc[ext.index, 'vital_status_at_date_of_last_contact'].eq('Deceased').astype(int)
d2d = pd.to_numeric(cl.loc[ext.index, 'number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_death'], errors = 'coerce')
fu = cl.loc[ext.index, 'follow_up_period'].astype(str).str.extract(r'(\d+)')[0].astype(float) * 30.44
ext['time'] = np.where(ext['event'] == 1, d2d, fu)
ext = ext.dropna(subset = ['time'])
ext = ext[ext['time'] > 0]

from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored, cumulative_dynamic_auc
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

ci_ext = concordance_index_censored(ext['event'].astype(bool).values, ext['time'].values, ext['risk'].values)[0]
ext['g2'] = np.where(ext['risk'] >= ext['risk'].median(), 'high', 'low')
lr_ext = logrank_test(ext['time'][ext.g2=='high'], ext['time'][ext.g2=='low'], ext['event'][ext.g2=='high'], ext['event'][ext.g2=='low'])
kmf = KaplanMeierFitter()
plt.figure(figsize = (6,5))
for lev in ['low','high']:
    m = ext['g2'] == lev
    kmf.fit(ext['time'][m].values, ext['event'][m].values, label = lev + ' (n=' + str(int(m.sum())) + ')')
    kmf.plot_survival_function()
plt.title('CPTAC ccRCC OS by risk (logrank p=' + format(lr_ext.p_value, '.3f') + ')'); plt.xlabel('days')
plt.tight_layout(); plt.savefig(os.path.join(out, 'km_validation_external.png'), dpi = 200); plt.close()

yext = Surv.from_arrays(ext['event'].astype(bool), ext['time'])
etimes = [t for t in [365,1095,1825] if t < ext['time'][ext['event']==1].max()]
eauc, _ = cumulative_dynamic_auc(yext, yext, ext['risk'].values, etimes)
plt.figure(figsize = (6,5))
plt.plot([365,1095,1825][:len(eauc)], eauc, marker = 'o')
plt.ylim(0.4, 0.9); plt.xlabel('days'); plt.ylabel('time-dependent AUC'); plt.title('CPTAC ccRCC time-dependent AUC')
plt.tight_layout(); plt.savefig(os.path.join(out, 'roc_validation_external.png'), dpi = 200); plt.close()

auc_map = {365: np.nan, 1095: np.nan, 1825: np.nan}
for t, a in zip(etimes, eauc):
    auc_map[t] = a
disc = pd.read_csv(os.path.join(out, 'discrimination_metrics.csv'))
row = pd.DataFrame([{'cohort': 'CPTAC_ccRCC_primary', 'n': len(ext), 'events': int(ext['event'].sum()),
                     'cindex': ci_ext, 'auc_1yr': auc_map[365], 'auc_3yr': auc_map[1095], 'auc_5yr': auc_map[1825]}])
if 'CPTAC_ccRCC_primary' not in set(disc['cohort']):
    disc = pd.concat([disc, row], ignore_index = True)
disc.to_csv(os.path.join(out, 'discrimination_metrics.csv'), index = False)

from tidepy.pred import TIDE
tum = [c for c in kirc.columns if c[-2:] in ('01','05')]
tide = TIDE(kirc[tum], cancer = 'Other', force_normalize = True)
tide.index = [s[:15] for s in tide.index]
g = df.set_index('bcr')['grp']
tide['grp'] = g.reindex(tide.index)
tide = tide.dropna(subset = ['grp'])

from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
metrics = ['TIDE','Dysfunction','Exclusion','MDSC','CAF','TAM M2','CTL']
rows = []
for mt in metrics:
    hi = tide[tide.grp=='high'][mt]; lo = tide[tide.grp=='low'][mt]
    u, p = mannwhitneyu(hi, lo)
    rows.append([mt, hi.mean(), lo.mean(), hi.mean()-lo.mean(), p])
tide_res = pd.DataFrame(rows, columns = ['metric','high_mean','low_mean','delta','p'])
tide_res['fdr'] = multipletests(tide_res['p'], method = 'fdr_bh')[1]
tide_res.to_csv(os.path.join(out, 'riskgroup_tide.csv'), index = False)

melt = tide.melt(id_vars = 'grp', value_vars = metrics, var_name = 'metric', value_name = 'score')
plt.figure(figsize = (10,4))
sns.boxplot(data = melt, x = 'metric', y = 'score', hue = 'grp', showfliers = False)
plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')
plt.title('TIDE metrics by ccRCC risk group')
plt.tight_layout(); plt.savefig(os.path.join(out, 'riskgroup_tide.png'), dpi = 200); plt.close()

print('CPTAC external: n', len(ext), 'events', int(ext['event'].sum()), 'C-index', round(ci_ext,3), 'KM p', format(lr_ext.p_value,'.3f'))
print(tide_res.round(4).to_string(index = False))

pred_path = os.path.join(out, 'drug_out', 'DrugPredictions.csv')
if os.path.exists(pred_path):
    pred = pd.read_csv(pred_path, index_col = 0)
    pred.index = [s[:15] for s in pred.index]
    agents = ['sunitinib','pazopanib','axitinib','cabozantinib','everolimus','temsirolimus']
    rg = df.set_index('bcr')['grp']
    rows = []
    for a in agents:
        cols = [c for c in pred.columns if c.lower().split('_')[0] == a]
        if not cols:
            continue
        v = np.log(pred[cols[0]].clip(lower = 1e-6))
        v = v.groupby(level = 0).mean()
        gg = rg.reindex(v.index).dropna()
        hi = v.loc[gg[gg=='high'].index]; lo = v.loc[gg[gg=='low'].index]
        u, p = mannwhitneyu(hi, lo)
        rows.append([a, hi.mean(), lo.mean(), hi.mean()-lo.mean(), p])
    drug = pd.DataFrame(rows, columns = ['agent','high_mean_logIC50','low_mean_logIC50','delta','p'])
    drug['fdr'] = multipletests(drug['p'], method = 'fdr_bh')[1]
    drug.to_csv(os.path.join(out, 'riskgroup_drug.csv'), index = False)

    long = []
    for a in agents:
        cols = [c for c in pred.columns if c.lower().split('_')[0] == a]
        if not cols:
            continue
        v = np.log(pred[cols[0]].clip(lower = 1e-6)).groupby(level = 0).mean()
        gg = rg.reindex(v.index).dropna()
        for s in gg.index:
            long.append([a, gg[s], v[s]])
    ld = pd.DataFrame(long, columns = ['agent','grp','logIC50'])
    plt.figure(figsize = (10,4))
    sns.boxplot(data = ld, x = 'agent', y = 'logIC50', hue = 'grp', showfliers = False)
    plt.xticks(rotation = 35, rotation_mode = 'anchor', ha = 'right')
    plt.title('Predicted drug sensitivity by ccRCC risk group')
    plt.tight_layout(); plt.savefig(os.path.join(out, 'riskgroup_drug.png'), dpi = 200); plt.close()
    print(drug.round(4).to_string(index = False))
