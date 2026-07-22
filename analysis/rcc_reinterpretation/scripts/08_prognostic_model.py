import os
import json
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
out = os.path.join(root, 'analysis/rcc_reinterpretation/outputs/model')
os.makedirs(out, exist_ok = True)
np.random.seed(0)

core = ['C1QA', 'C1QB', 'C1QC', 'APOE', 'APOC1', 'TREM2', 'GPNMB', 'MERTK']

adata = ad.read_h5ad(os.path.join(root, 'kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad'))
mye = adata.raw.to_adata()
mye.obs = adata.obs.copy()

def map_group(x):
    if x == 'TAM':
        return 'TAM'
    elif x in ['Mono1', 'Mono2', 'Mono3']:
        return 'Mono'
    else:
        return 'other'

mye.obs['grp'] = mye.obs['final_label'].map(map_group)
sub = mye[mye.obs['grp'].isin(['TAM', 'Mono'])].copy()
sc.tl.rank_genes_groups(sub, 'grp', groups = ['TAM'], reference = 'Mono', method = 'wilcoxon')
deg = sc.get.rank_genes_groups_df(sub, group = 'TAM')
deg = deg[(deg['logfoldchanges'] > 0.5) & (deg['pvals_adj'] < 0.01)].sort_values('scores', ascending = False).head(150)

cand, seen = [], set()
for g in core:
    cand.append((g, 'core_module'))
    seen.add(g)
for g in deg['names']:
    if g not in seen:
        cand.append((g, 'TAM_DEG'))
        seen.add(g)
cand_df = pd.DataFrame(cand, columns = ['gene', 'source'])
cand_df.to_csv(os.path.join(out, 'candidate_genes.csv'), index = False)

kirc = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_HiSeqV2.gz'), sep = '\t', index_col = 0)
kirc.index = [g.upper() for g in kirc.index]
surv = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_survival.txt'), sep = '\t')
surv['bcr'] = surv['sample'].str[:15]

pool = [g for g in cand_df['gene'] if g in kirc.index]
expr = kirc.loc[pool].T
expr['bcr'] = [s[:15] for s in expr.index]
clin = surv[['bcr', 'OS', 'OS.time']].dropna()
clin = clin[clin['OS.time'] > 0]
df = expr.merge(clin, on = 'bcr')

from sklearn.model_selection import train_test_split
tr_idx, te_idx = train_test_split(np.arange(len(df)), test_size = 0.3, random_state = 0, stratify = df['OS'])
train = df.iloc[tr_idx].reset_index(drop = True)
test = df.iloc[te_idx].reset_index(drop = True)

from lifelines import CoxPHFitter
from statsmodels.stats.multitest import multipletests

uni = []
for g in pool:
    d = train[[g, 'OS', 'OS.time']].copy()
    d['z'] = (d[g] - d[g].mean()) / d[g].std()
    cox = CoxPHFitter().fit(d[['z', 'OS.time', 'OS']], 'OS.time', 'OS')
    r = cox.summary.loc['z']
    uni.append([g, np.exp(r['coef']), np.exp(r['coef lower 95%']), np.exp(r['coef upper 95%']), r['p']])
uni = pd.DataFrame(uni, columns = ['gene', 'HR', 'ci_low', 'ci_high', 'p'])
uni['fdr'] = multipletests(uni['p'], method = 'fdr_bh')[1]
uni = uni.sort_values('p')
uni.to_csv(os.path.join(out, 'univariate_cox.csv'), index = False)

sig = uni[uni['p'] < 0.05]['gene'].tolist()

from sksurv.linear_model import CoxnetSurvivalAnalysis
from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored, cumulative_dynamic_auc
from sklearn.model_selection import KFold
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

mu = train[sig].mean()
sd = train[sig].std()
xz = ((train[sig] - mu) / sd).values
y_tr = Surv.from_arrays(train['OS'].astype(bool), train['OS.time'])

path = CoxnetSurvivalAnalysis(l1_ratio = 1.0, alpha_min_ratio = 0.01, max_iter = 100000).fit(xz, y_tr)
alphas = path.alphas_

kf = KFold(n_splits = 10, shuffle = True, random_state = 0)
cv_mean = []
for a in alphas:
    fold = []
    for i_tr, i_te in kf.split(xz):
        try:
            m = CoxnetSurvivalAnalysis(l1_ratio = 1.0, alphas = [a], max_iter = 100000).fit(xz[i_tr], y_tr[i_tr])
            risk = m.predict(xz[i_te])
            fold.append(concordance_index_censored(train['OS'].astype(bool).values[i_te], train['OS.time'].values[i_te], risk)[0])
        except Exception:
            fold.append(np.nan)
    cv_mean.append(np.nanmean(fold))
cv_mean = np.array(cv_mean)
best = alphas[np.nanargmax(cv_mean)]
cv_cindex = float(np.nanmax(cv_mean))

plt.figure(figsize = (7, 4))
plt.plot(np.log10(alphas), cv_mean, marker = 'o', ms = 3)
plt.axvline(np.log10(best), color = 'red', ls = '--')
plt.xlabel('log10(alpha)')
plt.ylabel('cv concordance')
plt.title('LASSO-Cox 10-fold CV (train)')
plt.tight_layout()
plt.savefig(os.path.join(out, 'lasso_cv.png'), dpi = 200)
plt.close()

fin = CoxnetSurvivalAnalysis(l1_ratio = 1.0, alphas = [best], fit_baseline_model = True, max_iter = 100000).fit(xz, y_tr)
coef = pd.Series(fin.coef_.ravel(), index = sig)
coef = coef[coef != 0]
coef.to_frame('coef').to_csv(os.path.join(out, 'lasso_selected.csv'))

model = {'genes': coef.index.tolist(), 'coef': coef.values.tolist(),
         'mean': mu[coef.index].values.tolist(), 'sd': sd[coef.index].values.tolist(),
         'alpha': float(best), 'seed': 0, 'l1_ratio': 1.0,
         'training_cohort': 'TCGA-KIRC (70% split)', 'endpoint': 'OS',
         'cv_cindex_train': cv_cindex,
         'formula': 'risk = sum(coef_i * (expr_i - mean_i)/sd_i)'}
with open(os.path.join(out, 'risk_model.json'), 'w') as f:
    json.dump(model, f, indent = 2)

genes = np.array(model['genes'])
mm = np.array(model['mean'])
ss = np.array(model['sd'])
cc = np.array(model['coef'])

def risk_of(mat):
    have = [g for g in genes if g in mat.columns]
    idx = np.isin(genes, have)
    z = (mat[list(genes[idx])].values - mm[idx]) / ss[idx]
    return z @ cc[idx]

for d in (train, test):
    d['risk'] = risk_of(d)

from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
kmf = KaplanMeierFitter()

def km_plot(d, title, fname):
    d = d.copy()
    d['OS.time'] = pd.to_numeric(d['OS.time'])
    d['OS'] = pd.to_numeric(d['OS']).astype(int)
    cut = d['risk'].median()
    d['g2'] = np.where(d['risk'] >= cut, 'high', 'low')
    plt.figure(figsize = (6, 5))
    for lev in ['low', 'high']:
        m = d['g2'] == lev
        if m.sum() == 0:
            continue
        kmf.fit(d['OS.time'][m].values, d['OS'][m].values, label = lev + ' (n=' + str(int(m.sum())) + ')')
        kmf.plot_survival_function()
    hi, lo = d['g2'] == 'high', d['g2'] == 'low'
    lr = logrank_test(d['OS.time'][hi].values, d['OS.time'][lo].values, d['OS'][hi].values, d['OS'][lo].values)
    plt.title(title + ' (logrank p=' + format(lr.p_value, '.2e') + ')')
    plt.xlabel('days')
    plt.tight_layout()
    plt.savefig(os.path.join(out, fname), dpi = 200)
    plt.close()
    return lr.p_value

p_train = km_plot(train, 'TCGA-KIRC train OS by risk', 'km_train.png')
p_test = km_plot(test, 'TCGA-KIRC held-out OS by risk', 'km_test.png')

def disc(d, tag):
    ev = d['OS'].astype(bool).values
    c = concordance_index_censored(ev, d['OS.time'].values, d['risk'].values)[0]
    yrs = [365, 1095, 1825]
    tmax = d['OS.time'][d['OS'] == 1].max()
    times = [t for t in yrs if t < tmax]
    auc = {365: np.nan, 1095: np.nan, 1825: np.nan}
    try:
        a, _ = cumulative_dynamic_auc(y_tr, Surv.from_arrays(ev, d['OS.time']), d['risk'].values, times)
        for i, t in enumerate(times):
            auc[t] = a[i]
    except Exception:
        pass
    return {'cohort': tag, 'n': len(d), 'events': int(d['OS'].sum()), 'cindex': c,
            'auc_1yr': auc[365], 'auc_3yr': auc[1095], 'auc_5yr': auc[1825]}

metrics = [disc(train, 'TCGA-KIRC_train'), disc(test, 'TCGA-KIRC_heldout')]

braun_e = pd.read_csv(os.path.join(root, 'resources/external_icb/braun_checkmate/processed/expression_normalized.tsv'), sep = '\t', index_col = 0)
braun_e.index = [g.upper() for g in braun_e.index]
braun_s = pd.read_csv(os.path.join(root, 'resources/external_icb/braun_checkmate/processed/survival.tsv'), sep = '\t')
braun_m = pd.read_csv(os.path.join(root, 'resources/external_icb/braun_checkmate/processed/sample_manifest.tsv'), sep = '\t')
have = [g for g in genes if g in braun_e.index]
bt = braun_e.loc[have].T
bt['RNA_ID'] = bt.index
bt = bt.merge(braun_m[['RNA_ID', 'SUBJID']], on = 'RNA_ID').merge(braun_s[['SUBJID', 'OS', 'OS_CNSR']], on = 'SUBJID')
bt = bt.dropna(subset = ['OS', 'OS_CNSR'])
bt = bt[bt['OS'] > 0]
bt['ev'] = (bt['OS_CNSR'] == 0).astype(int)
bt = bt.rename(columns = {'OS': 'OS.time', 'ev': 'OS'})
bt['OS.time'] = pd.to_numeric(bt['OS.time'])
bt['OS'] = bt['OS'].astype(int)
bt['risk'] = risk_of(bt)
p_braun = km_plot(bt, 'Braun/CheckMate (metastatic, ICB) OS by risk', 'km_validation_braun.png')
metrics.append(disc(bt, 'Braun_CheckMate_metastatic_ICB'))

metrics = pd.DataFrame(metrics)
metrics.to_csv(os.path.join(out, 'discrimination_metrics.csv'), index = False)
print(metrics.to_string(index = False))

rng = np.random.default_rng(0)
ev_te = test['OS'].astype(bool).values
t_te = test['OS.time'].values
rand_c = []
for _ in range(500):
    gg = list(rng.choice(pool, size = len(model['genes']), replace = False))
    z = (test[gg] - train[gg].mean()) / train[gg].std()
    rc = z.values @ rng.normal(size = len(gg))
    rand_c.append(concordance_index_censored(ev_te, t_te, rc)[0])
rand_c = np.abs(np.array(rand_c) - 0.5) + 0.5

kz = kirc.copy()
kz.columns = [s[:15] for s in kz.columns]
ctrl = {'panTAM': ['CD68', 'CD163', 'MRC1', 'CSF1R', 'LYZ', 'AIF1', 'FCGR3A'],
        'Obradovic': ['TREM2', 'APOE', 'APOC1', 'C1QA', 'C1QB', 'C1QC', 'GPNMB', 'FOLR2', 'SPP1', 'CTSD', 'CD68']}
ctrl_c = {}
for k, gs in ctrl.items():
    gg = [g for g in gs if g in kz.index]
    s = ((kz.loc[gg].T - kz.loc[gg].T.mean()) / kz.loc[gg].T.std()).mean(axis = 1)
    v = test['bcr'].map(s)
    ctrl_c[k] = concordance_index_censored(ev_te, t_te, v.values)[0]

sig_c = float(metrics.loc[metrics.cohort == 'TCGA-KIRC_heldout', 'cindex'].iloc[0])
sanity = pd.DataFrame([{'metric': 'signature_heldout_cindex', 'value': sig_c},
                       {'metric': 'random_gene_cindex_mean', 'value': float(np.mean(rand_c))},
                       {'metric': 'random_gene_cindex_p95', 'value': float(np.percentile(rand_c, 95))},
                       {'metric': 'panTAM_cindex', 'value': ctrl_c['panTAM']},
                       {'metric': 'Obradovic_cindex', 'value': ctrl_c['Obradovic']}])
sanity.to_csv(os.path.join(out, 'sanity_checks.csv'), index = False)
print(sanity.to_string(index = False))

led = [['ccRCC held-out KM (test)', 'logrank_p', np.nan, np.nan, p_test, len(test), 'TCGA-KIRC_heldout', '08_prognostic_model.py', 'km_test'],
       ['ccRCC train KM', 'logrank_p', np.nan, np.nan, p_train, len(train), 'TCGA-KIRC_train', '08_prognostic_model.py', 'km_train'],
       ['ccRCC Braun external KM', 'logrank_p', np.nan, np.nan, p_braun, len(bt), 'Braun_metastatic_ICB', '08_prognostic_model.py', 'km_braun'],
       ['ccRCC 10-fold CV C-index (train)', cv_cindex, np.nan, np.nan, np.nan, len(train), 'TCGA-KIRC_train', '08_prognostic_model.py', 'cv']]
for _, r in metrics.iterrows():
    led.append(['C-index ' + r['cohort'], r['cindex'], np.nan, np.nan, np.nan, r['n'], r['cohort'], '08_prognostic_model.py', 'disc'])
    for yy in ['auc_1yr', 'auc_3yr', 'auc_5yr']:
        led.append([yy + ' ' + r['cohort'], r[yy], np.nan, np.nan, np.nan, r['n'], r['cohort'], '08_prognostic_model.py', 'disc'])
pd.DataFrame(led, columns = ['claim', 'value', 'ci_low', 'ci_high', 'p', 'n', 'cohort', 'script', 'line']).to_csv(os.path.join(out, 'model_claim_ledger_partA.csv'), index = False)

print('genes:', model['genes'])
print('cv C-index (train):', round(cv_cindex, 3))
print('held-out logrank p', format(p_test, '.2e'), '| braun logrank p', format(p_braun, '.2e'))
print('random-gene test C-index mean', round(float(np.mean(rand_c)), 3), '95pct', round(float(np.percentile(rand_c, 95)), 3))
