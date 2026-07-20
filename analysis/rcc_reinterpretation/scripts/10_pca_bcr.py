import os
import json
import numpy as np
import pandas as pd
import anndata as ad
import decoupler as dc

root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', '..'))
out = os.path.join(root, 'analysis/rcc_reinterpretation/outputs/model')
np.random.seed(0)

mods = {'ATF3_NFkB': ['ATF3','NFKB1','RELA','JUN','JUNB','DUSP1','TNFAIP3','NFKBIA'],
        'complement_C1Q': ['C1QA','C1QB','C1QC'],
        'RCC_skew_CORE': ['C1QA','C1QB','C1QC','APOE','APOC1','TREM2'],
        'CLEC_LAM8': ['C1QA','C1QB','C1QC','APOE','APOC1','TREM2','GPNMB','MERTK']}

def aucell_bulk(expr_ts):
    a = ad.AnnData(expr_ts.values.astype(float), obs = pd.DataFrame(index = expr_ts.index), var = pd.DataFrame(index = expr_ts.columns))
    net = pd.concat([pd.DataFrame({'source': k, 'target': [g for g in v if g in a.var_names], 'weight': 1.0}) for k, v in mods.items()])
    dc.mt.aucell(a, net, tmin = 2, verbose = False)
    return a.obsm['score_aucell']

prad = pd.read_csv(os.path.join(root, 'resources/model_data/prad_tcga/HiSeqV2.gz'), sep = '\t', index_col = 0)
prad.index = [g.upper() for g in prad.index]
sv = pd.read_csv(os.path.join(root, 'resources/model_data/prad_tcga/PRAD_survival.txt'), sep = '\t')
cl = pd.read_csv(os.path.join(root, 'resources/model_data/prad_tcga/PRAD_clinicalMatrix', ), sep = '\t', index_col = 0)

x = prad.T
sc = aucell_bulk(x)
sc.index = x.index

sv['id'] = sv['sample']
d = sc.copy()
d['id'] = d.index
d = d.merge(sv[['id','PFI','PFI.time']], on = 'id')
d = d[d['PFI.time'] > 0].dropna(subset = ['PFI','PFI.time'])

bcr = cl[['biochemical_recurrence','days_to_first_biochemical_recurrence','days_to_last_followup']].copy()
bcr['id'] = bcr.index
bcr['bcr_event'] = bcr['biochemical_recurrence'].map({'YES': 1, 'NO': 0})
bcr['bcr_time'] = np.where(bcr['bcr_event'] == 1, bcr['days_to_first_biochemical_recurrence'], bcr['days_to_last_followup'])
d = d.merge(bcr[['id','bcr_event','bcr_time']], on = 'id', how = 'left')
for c in ['PFI','PFI.time','bcr_event','bcr_time']:
    d[c] = pd.to_numeric(d[c], errors = 'coerce')
d.shape

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from statsmodels.stats.multitest import multipletests
from sksurv.util import Surv
from sksurv.metrics import concordance_index_censored, cumulative_dynamic_auc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

model = json.load(open(os.path.join(out, 'risk_model.json')))
genes = np.array(model['genes']); mm = np.array(model['mean']); ss = np.array(model['sd']); cc = np.array(model['coef'])
have = [g for g in genes if g in prad.index]
idx = np.isin(genes, have)
z = (x[list(genes[idx])].values - mm[idx]) / ss[idx]
sig_risk = pd.Series(z @ cc[idx], index = x.index)
d['partA_signature'] = d['id'].map(sig_risk)

def cox_endpoint(dat, cols, ev, tm):
    rows = []
    for c in cols:
        s = dat[[c, ev, tm]].dropna().copy()
        s = s[s[tm] > 0]
        s['zz'] = (s[c] - s[c].mean()) / s[c].std()
        fit = CoxPHFitter().fit(s[['zz', tm, ev]], tm, ev)
        r = fit.summary.loc['zz']
        ci = concordance_index_censored(s[ev].astype(bool).values, s[tm].values, s['zz'].values)[0]
        rows.append([c, np.exp(r['coef']), np.exp(r['coef lower 95%']), np.exp(r['coef upper 95%']), r['p'], ci, len(s), int(s[ev].sum())])
    out_df = pd.DataFrame(rows, columns = ['program','HR','ci_low','ci_high','p','cindex','n','events'])
    out_df['fdr'] = multipletests(out_df['p'], method = 'fdr_bh')[1]
    return out_df

cols = list(mods) + ['partA_signature']
pfi = cox_endpoint(d, cols, 'PFI', 'PFI.time')
pfi['endpoint'] = 'PFI'
bcr_res = cox_endpoint(d.rename(columns = {'bcr_event': 'BCR', 'bcr_time': 'BCR.time'}), cols, 'BCR', 'BCR.time')
bcr_res['endpoint'] = 'biochemical_recurrence'
pca_uni = pd.concat([pfi, bcr_res], ignore_index = True)
pca_uni.to_csv(os.path.join(out, 'pca_univariate_bcr.csv'), index = False)
print(pca_uni[['program','endpoint','HR','ci_low','ci_high','p','fdr','cindex','n','events']].round(3).to_string(index = False))

kmf = KaplanMeierFitter()
for c in ['ATF3_NFkB','complement_C1Q']:
    s = d[[c,'PFI','PFI.time']].dropna().copy()
    s['PFI.time'] = pd.to_numeric(s['PFI.time']); s['PFI'] = pd.to_numeric(s['PFI']).astype(int)
    s['g2'] = np.where(s[c] > s[c].median(), 'high', 'low')
    plt.figure(figsize = (6,5))
    for lev in ['low','high']:
        m = (s['g2'] == lev).values
        if m.sum() == 0:
            continue
        kmf.fit(s['PFI.time'].values[m], s['PFI'].values[m], label = lev + ' (n=' + str(int(m.sum())) + ')')
        kmf.plot_survival_function()
    lr = logrank_test(s['PFI.time'][s.g2=='high'].values, s['PFI.time'][s.g2=='low'].values, s['PFI'][s.g2=='high'].values, s['PFI'][s.g2=='low'].values)
    plt.title('TCGA-PRAD PFI by ' + c + ' (logrank p=' + format(lr.p_value, '.2e') + ')'); plt.xlabel('days')
    plt.tight_layout(); plt.savefig(os.path.join(out, 'pca_km_prad_' + c + '.png'), dpi = 200); plt.close()

disc_rows = []
yy = Surv.from_arrays(d['PFI'].astype(bool), d['PFI.time'])
for c in cols:
    s = d[[c,'PFI','PFI.time']].dropna().copy()
    s['zz'] = (s[c] - s[c].mean()) / s[c].std()
    ci = concordance_index_censored(s['PFI'].astype(bool).values, s['PFI.time'].values, s['zz'].values)[0]
    times = [t for t in [365,1095,1825] if t < s['PFI.time'][s['PFI']==1].max()]
    auc = {365: np.nan, 1095: np.nan, 1825: np.nan}
    try:
        a, _ = cumulative_dynamic_auc(Surv.from_arrays(s['PFI'].astype(bool), s['PFI.time']), Surv.from_arrays(s['PFI'].astype(bool), s['PFI.time']), s['zz'].values, times)
        for i, t in enumerate(times):
            auc[t] = a[i]
    except Exception:
        pass
    disc_rows.append({'program': c, 'cindex': ci, 'auc_1yr': auc[365], 'auc_3yr': auc[1095], 'auc_5yr': auc[1825], 'n': len(s), 'events': int(s['PFI'].sum())})
pca_disc = pd.DataFrame(disc_rows)
pca_disc.to_csv(os.path.join(out, 'pca_discrimination.csv'), index = False)

kirc = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_HiSeqV2.gz'), sep = '\t', index_col = 0)
kirc.index = [g.upper() for g in kirc.index]
ksv = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_survival.txt'), sep = '\t')
ksv['id'] = ksv['sample']
kx = kirc.T
ksc = aucell_bulk(kx); ksc.index = kx.index; ksc['id'] = ksc.index
kd = ksc.merge(ksv[['id','OS','OS.time']], on = 'id')
kd = kd[kd['OS.time'] > 0].dropna(subset = ['OS','OS.time'])
rcc_os = cox_endpoint(kd, list(mods), 'OS', 'OS.time')

cont = []
for m in mods:
    rr = rcc_os[rcc_os.program == m].iloc[0]
    pp = pfi[pfi.program == m].iloc[0]
    cont.append({'program': m,
                 'ccRCC_OS_HR': round(rr['HR'],3), 'ccRCC_OS_p': round(rr['p'],4), 'ccRCC_OS_cindex': round(rr['cindex'],3),
                 'PCa_BCR_HR': round(pp['HR'],3), 'PCa_BCR_p': round(pp['p'],4), 'PCa_BCR_cindex': round(pp['cindex'],3),
                 'ccRCC_sig': 'yes' if rr['p'] < 0.05 else 'no', 'PCa_sig': 'yes' if pp['p'] < 0.05 else 'no'})
cross = pd.DataFrame(cont)
cross.to_csv(os.path.join(out, 'cross_tumor_prognostic_table.csv'), index = False)
print()
print(cross.to_string(index = False))
