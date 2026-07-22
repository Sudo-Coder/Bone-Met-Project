import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

sys.path.insert(0, os.path.dirname(__file__))
import style
from style import PCA, RCC, ATF3, C1Q, CLEC, GRAY, GRAD_CMAP, mm, W1, W15, W2
from style import panel_label, style_ax, save_panel, stars, mlabel, MODULE_HIDE

np.random.seed(0)
root = '/autofs/projects-t3/hussain/scProj'
out = os.path.join(root, 'analysis/rcc_reinterpretation/outputs')
mod = os.path.join(out, 'model')
tab = os.path.join(out, 'tables')
pca = os.path.join(root, 'analysis/pca_comparator/outputs/tables')
fig_dir = os.path.dirname(__file__)
cache = os.path.join(fig_dir, 'cache')
os.makedirs(cache, exist_ok = True)

def savefig(fig, name, w):
    fig.savefig(os.path.join(fig_dir, name + '.pdf'))
    fig.savefig(os.path.join(fig_dir, name + '.png'), dpi = 600)
    plt.close(fig)

def load(p):
    return pd.read_csv(p)

def prep_km():
    tcga = os.path.join(cache, 'km_tcga.csv')
    if not os.path.exists(tcga):
        m = json.load(open(os.path.join(mod, 'risk_model.json')))
        g = np.array(m['genes']); mu = np.array(m['mean']); sd = np.array(m['sd']); co = np.array(m['coef'])
        k = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_HiSeqV2.gz'), sep = '\t', index_col = 0)
        k.index = [x.upper() for x in k.index]
        s = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_survival.txt'), sep = '\t')
        s['bcr'] = s['sample'].str[:15]
        x = k.loc[[q for q in g if q in k.index]].T
        z = (x[list(g)] - mu) / sd
        d = pd.DataFrame({'risk': z.values @ co}, index = x.index)
        d['bcr'] = [i[:15] for i in d.index]
        d = d.merge(s[['bcr','OS','OS.time']], on = 'bcr')
        d = d[d['OS.time'] > 0].dropna(subset = ['OS','OS.time'])
        d['grp'] = np.where(d['risk'] >= d['risk'].median(), 'high', 'low')
        d[['risk','OS','OS.time','grp']].to_csv(tcga, index = False)
    cptac_path = os.path.join(cache, 'km_cptac.csv')
    if not os.path.exists(cptac_path):
        try:
            import warnings; warnings.filterwarnings('ignore')
            import cptac
            m = json.load(open(os.path.join(mod, 'risk_model.json')))
            g = m['genes']; co = np.array(m['coef'])
            cc = cptac.Ccrcc()
            tx = cc.get_transcriptomics('broad'); tx.columns = tx.columns.get_level_values(0)
            tx = tx.groupby(level = 0, axis = 1).mean()
            cl = cc.get_clinical('mssm')
            have = [q for q in g if q in tx.columns]
            e = tx.loc[[i for i in tx.index if i in cl.index], have]
            ze = (e - e.mean()) / e.std()
            w = np.array([co[g.index(q)] for q in have])
            d = pd.DataFrame({'risk': ze.values @ w}, index = e.index)
            d['event'] = cl.loc[d.index, 'vital_status_at_date_of_last_contact'].eq('Deceased').astype(int)
            d2 = pd.to_numeric(cl.loc[d.index, 'number_of_days_from_date_of_initial_pathologic_diagnosis_to_date_of_death'], errors = 'coerce')
            fu = cl.loc[d.index, 'follow_up_period'].astype(str).str.extract(r'(\d+)')[0].astype(float) * 30.44
            d['time'] = np.where(d['event'] == 1, d2, fu)
            d = d.dropna(subset = ['time']); d = d[d['time'] > 0]
            d['grp'] = np.where(d['risk'] >= d['risk'].median(), 'high', 'low')
            d[['risk','event','time','grp']].to_csv(cptac_path, index = False)
        except Exception as e:
            print('CPTAC KM cache failed:', e)

def km_panel(ax, df, tcol, ecol, title, p, risk_times):
    # Reporting convention for survival curves (JAMA/NEJM): log-rank p, 95% CI bands, and a
    # number-at-risk table under the x-axis. Hazard ratios belong in the text, not the plot.
    from lifelines import KaplanMeierFitter
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    col = {'low': RCC, 'high': PCA}
    yrs = df[tcol].values / 365.25
    for lev in ['low','high']:
        mk = (df['grp'] == lev).values
        kmf = KaplanMeierFitter().fit(yrs[mk], df[ecol].values[mk], label = lev + ' risk')
        kmf.plot_survival_function(ax = ax, ci_show = True, ci_alpha = 0.12,
                                   color = col[lev], lw = 1.2)
    style_ax(ax)
    ax.set_ylabel('overall survival'); ax.set_ylim(0, 1.02)
    ax.set_xlim(0, max(risk_times)); ax.set_xticks(risk_times)
    ax.set_xlabel('')
    ax.tick_params(labelbottom = False)          # the at-risk table below carries the x axis
    ax.set_title('%s\nlog-rank p = %s' % (title, '%.3f' % p if p >= 1e-3 else '%.1e' % p), fontsize = 8)
    ax.legend(loc = 'lower left', fontsize = 6.5, handlelength = 1.2)

    # number at risk = still under observation at that time point. y=1 is the top row, so
    # 'high' sits on top; labels, colours and counts are all driven off the same mapping.
    tab = make_axes_locatable(ax).append_axes('bottom', size = '20%', pad = 0.42, sharex = ax)
    for sp in tab.spines.values():
        sp.set_visible(False)
    from matplotlib.transforms import blended_transform_factory
    ypos = {'high': 1, 'low': 0}
    tab.set_yticks([])                                 # row labels drawn as text instead, so
    tab.tick_params(axis = 'x', length = 2.5, width = 0.75)   # they share va with the counts
    tab.set_ylim(-0.7, 1.7)
    lab_tr = blended_transform_factory(tab.transAxes, tab.transData)
    for lev, yi in ypos.items():
        tab.text(-0.03, yi, lev + ' risk', transform = lab_tr, ha = 'right', va = 'center',
                 fontsize = 6, fontweight = 'bold', color = col[lev])
        t_lev = yrs[(df['grp'] == lev).values]
        for t in risk_times:
            tab.text(t, yi, str(int((t_lev >= t).sum())),
                     ha = 'left' if t == risk_times[0] else 'center', va = 'center',
                     fontsize = 6, color = col[lev])
    tab.set_xlabel('years')
    tab.text(0.0, 1.10, 'Number at risk', transform = tab.transAxes, ha = 'left', va = 'bottom',
             fontsize = 6, fontweight = 'bold', color = '#555555')


def draw4a(ax):
    d = load(os.path.join(pca, 'pca_vs_rcc_contrast_table.csv'))
    order = ['ATF3_NFkB','complement_C1Q','CLEC_LAM8','APOE_TREM2','MERTK_GPNMB','SPP1_TAM','panTAM']
    d = d.set_index('module').loc[order].reset_index()
    y = np.arange(len(d))[::-1]
    h = 0.38
    ax.barh(y + h/2, d['RCC_TvB_SD'], height = h, color = RCC, label = 'ccRCC')
    ax.barh(y - h/2, d['PCa_TvB_SD'], height = h, color = PCA, label = 'PCa')
    ax.set_yticks(y); ax.set_yticklabels([mlabel(m) for m in d['module']])
    ax.set_xlabel('tumor induction (SD, tumor vs shared benign)')
    ax.legend(loc = 'lower right', fontsize = 7)
    ax.set_xlim(0, max(d['RCC_TvB_SD'].max(), d['PCa_TvB_SD'].max()) * 1.12)
    style_ax(ax)

def draw4b(ax):
    d = load(os.path.join(pca, 'pca_vs_rcc_contrast_table.csv'))
    d = d[d['interaction'].notna() & ~d['module'].isin(MODULE_HIDE)].copy()
    order = d.sort_values('interaction')
    y = np.arange(len(order))
    cols = [RCC if p < 0.05 else GRAY for p in order['interaction_p']]
    ax.barh(y, order['interaction'], color = cols, height = 0.62)
    ax.set_yticks(y); ax.set_yticklabels([mlabel(m) for m in order['module']])
    ax.axvline(0, color = 'k', lw = 0.6)
    ax.set_xlabel('cancer-type × condition interaction (RCC skew)')
    for yi, (v, p) in enumerate(zip(order['interaction'], order['interaction_p'])):
        sig_label(ax, v + 0.004 if v >= 0 else v - 0.004, yi, stars(p), ha = 'left' if v >= 0 else 'right')
    ax.set_xlim(min(order['interaction'].min(), 0) - 0.03, order['interaction'].max() * 1.25)
    h = [plt.Rectangle((0,0),1,1,color=RCC), plt.Rectangle((0,0),1,1,color=GRAY)]
    ax.legend(h, ['RCC-skewed (p<0.05)','not skewed'], loc = 'lower right', fontsize = 6.5)
    style_ax(ax)

def draw4c(ax):
    d = load(os.path.join(pca, 'pca_vs_rcc_contrast_table.csv')).set_index('module')
    r = d.loc['ATF3_NFkB']
    ax.bar([0, 1], [r['RCC_TvB_SD'], r['PCa_TvB_SD']], color = [RCC, PCA], width = 0.6)
    ax.set_xticks([0, 1]); ax.set_xticklabels(['ccRCC', 'PCa'])
    ax.set_ylabel('ATF3/NF-κB induction (SD)')
    ax.set_title('conserved core', fontsize = 8, color = ATF3)
    for xi, v in enumerate([r['RCC_TvB_SD'], r['PCa_TvB_SD']]):
        ax.text(xi, v + 0.03, '%.2f' % v, ha = 'center', fontsize = 7)
    ax.set_ylim(0, 1.75)
    style_ax(ax)

def draw4d(ax):
    d = load(os.path.join(pca, 'pca_within_tumor_localization.csv'))
    order = ['ATF3_NFkB','complement_C1Q','CLEC_LAM8','APOE_TREM2','SPP1_TAM']
    d = d.set_index('module').loc[order].reset_index()
    y = np.arange(len(d))[::-1]
    h = 0.38
    ax.barh(y + h/2, d['TAM_TIM_mean'], height = h, color = C1Q, label = 'TAM/TIM')
    ax.barh(y - h/2, d['Mono_mean'], height = h, color = GRAY, label = 'monocyte')
    ax.set_yticks(y); ax.set_yticklabels([mlabel(m) for m in d['module']])
    ax.set_xlabel('module score (PCa within-tumor)')
    xmax = d['TAM_TIM_mean'].max() * 1.2
    for i in range(len(d)):
        sig_label(ax, d['TAM_TIM_mean'].iloc[i] + 0.006, y[i] + h / 2, stars(d['wilcoxon_p'].iloc[i]), ha = 'left')
    ax.set_xlim(0, xmax)
    ax.legend(loc = 'lower right', bbox_to_anchor = (1.04, -0.04), fontsize = 7)
    ax.set_title('PCa TAM localization', fontsize = 8)
    style_ax(ax)

def fig4_star_levels():
    # grades drawn in panels B (interaction) and C (PCa TAM localization)
    ps = list(load(os.path.join(pca, 'pca_vs_rcc_contrast_table.csv'))['interaction_p'].dropna())
    ps += list(load(os.path.join(pca, 'pca_within_tumor_localization.csv'))['wilcoxon_p'].dropna())
    used = {stars(p) for p in ps} - {'ns'}
    return [k for k in ['*','**','***','****'] if k in used]


def figure4():
    fig = plt.figure(figsize = (W2, mm(128)))
    gs = fig.add_gridspec(2, 2, hspace = 0.55, wspace = 0.42,
                          left = 0.12, right = 0.97, top = 0.9, bottom = 0.17)
    a = fig.add_subplot(gs[:, 0]); draw4a(a); panel_label(a, 'A', x = -0.22, y = 1.0)
    b = fig.add_subplot(gs[0, 1]); draw4b(b); panel_label(b, 'B')
    c = fig.add_subplot(gs[1, 1]); draw4d(c); panel_label(c, 'C')
    star_key(fig, levels = fig4_star_levels())
    savefig(fig, 'Figure4', W2)
    for nm, fn, w, h in [('F4A_induction', draw4a, W1, mm(55)), ('F4B_interaction', draw4b, W1, mm(55)),
                         ('F4C_localization', draw4d, W1, mm(52))]:
        f = plt.figure(figsize = (w, h)); ax = f.add_subplot(111); fn(ax); save_panel(f, nm)


def draw5a(ax):
    d = load(os.path.join(tab, 'complement_source_by_class.csv'))
    genes = ['C1QA','C1QB','C1QC']
    classes = ['TAM_CLEC_LAM','TAM_other','TIM','Tumor','Endothelial','MSC','Pericyte','Osteoclast']
    x = np.arange(len(classes))
    w = 0.26
    cols = [style.GRAD[4], style.GRAD[2], style.GRAD[1]]
    for i, g in enumerate(genes):
        sub = d[d['gene'] == g].set_index('compartment').reindex(classes)
        ax.bar(x + (i - 1) * w, sub['mean_lognorm'].values, width = w, color = cols[i], label = g)
    ax.set_xticks(x); ax.set_xticklabels([c.replace('_',' ') for c in classes], rotation = 35,
                                          rotation_mode = 'anchor', ha = 'right', fontsize = 6.5)
    ax.set_ylabel('mean log-norm expression')
    ax.set_title('C1q source by cell class', fontsize = 8)
    ax.legend(fontsize = 7, title = None)
    style_ax(ax)

def draw5b(ax):
    d = load(os.path.join(tab, 'cellchat_allLR_tumor.csv'))
    keep = d[(d['target'] == 'TAM_CLEC_LAM') & (d['pathway_name'].isin(['ApoE','COMPLEMENT','GAS']))]
    pairs = ['APOE - TREM2','C3 - (ITGAX+ITGB2)','C3 - C3AR1','GAS6 - AXL','GAS6 - MERTK']
    keep = keep[keep['interaction_name_2'].isin(pairs)]
    srcs = ['Tumor','TAM_CLEC_LAM','Endothelial','Osteoclast']
    yi = {p: i for i, p in enumerate(pairs[::-1])}
    xi = {s: i for i, s in enumerate(srcs)}
    pm = keep['prob'].max()
    def ssize(p): return 20 + 240 * p / pm
    for _, r in keep.iterrows():
        if r['source'] not in xi: continue
        ax.scatter(xi[r['source']], yi[r['interaction_name_2']], s = ssize(r['prob']),
                   c = [r['prob']], cmap = GRAD_CMAP, vmin = 0, vmax = pm, edgecolors = 'k', linewidths = 0.5)
    ax.set_xticks(range(len(srcs))); ax.set_xticklabels(srcs, rotation = 30, ha = 'right', fontsize = 6.5)
    ax.set_yticks(range(len(pairs))); ax.set_yticklabels(pairs[::-1], fontsize = 6.5)
    ax.set_xlim(-0.7, len(srcs) - 0.3); ax.set_ylim(-0.9, len(pairs) - 0.1)
    ax.set_title('predicted signaling to CLEC_LAM', fontsize = 8)
    ax.set_xlabel('sender', labelpad = 0)
    sm = plt.cm.ScalarMappable(cmap = GRAD_CMAP, norm = plt.Normalize(0, pm))
    cb = ax.figure.colorbar(sm, ax = ax, fraction = 0.045, pad = 0.03)
    cb.set_label('communication probability', fontsize = 6); cb.ax.tick_params(labelsize = 6)
    for pv in [0.03, 0.09, 0.15]:
        ax.scatter([], [], s = ssize(pv), facecolor = '#d9d9d9', edgecolors = 'k',
                   linewidths = 0.5, label = '%.2f' % pv)
    leg = ax.legend(title = 'communication probability\n', loc = 'upper center', bbox_to_anchor = (0.5, -0.27),
                    ncol = 3, fontsize = 7, title_fontsize = 7, columnspacing = 1.6,
                    handletextpad = 0.6, borderpad = 0.7, framealpha = 1.0, edgecolor = '#777777')
    leg.get_frame().set_linewidth(0.8)
    style_ax(ax)

def draw5c(ax):
    d = load(os.path.join(tab, 'nichenet_prespecified_ranks.csv')).sort_values('aupr_corrected')
    y = np.arange(len(d))
    comp = {'C1QA','C1QB','C1QC','C3'}
    cols = [C1Q if g in comp else GRAY for g in d['test_ligand']]
    ax.barh(y, d['aupr_corrected'], color = cols, height = 0.62)
    ax.set_yticks(y); ax.set_yticklabels(d['test_ligand'])
    for yi, (v, rk) in enumerate(zip(d['aupr_corrected'], d['rank'])):
        ax.text(v + 0.0004, yi, 'rank %d/%d' % (rk, d['n_ligands'].iloc[0]), va = 'center', fontsize = 6.5)
    ax.set_xlabel('NicheNet ligand activity (AUPR corrected)')
    ax.set_title('ligand prioritization', fontsize = 8)
    ax.set_xlim(0, d['aupr_corrected'].max() * 1.35)
    h = [plt.Rectangle((0, 0), 1, 1, color = C1Q), plt.Rectangle((0, 0), 1, 1, color = GRAY)]
    ax.legend(h, ['complement family', 'other axes'], loc = 'lower right', fontsize = 6.5)
    style_ax(ax)

def draw5d(ax):
    cf = os.path.join(cache, 'umap_c1q.parquet')
    if os.path.exists(cf):
        dfu = pd.read_parquet(cf)
    else:
        import anndata as ad
        a = ad.read_h5ad(os.path.join(root, 'kidney-cancer/Cleaned_Data/myeloid_FINAL_labels.h5ad'), backed = 'r')
        um = a.obsm['X_umap']
        sc = pd.read_parquet(os.path.join(tab, 'module_scores_percell.parquet'))
        sc = sc[sc['cancer_type'] == 'RCC']
        idx = a.obs_names
        common = idx.intersection(sc.index)
        score = pd.Series(np.nan, index = idx)
        score.loc[common] = sc.loc[common, 'complement_C1Q'].values
        dfu = pd.DataFrame({'x': um[:, 0], 'y': um[:, 1], 's': score.values}).dropna()
        dfu.to_parquet(cf)
    o = dfu['s'].argsort().values
    sca = ax.scatter(dfu['x'].values[o], dfu['y'].values[o], c = dfu['s'].values[o], s = 2,
                     cmap = GRAD_CMAP, vmin = 0, vmax = np.percentile(dfu['s'], 99), linewidths = 0, rasterized = True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel('UMAP-1'); ax.set_ylabel('UMAP-2')
    ax.set_title('ccRCC myeloid: Complement C1q score', fontsize = 8)
    for sp in ax.spines.values(): sp.set_visible(True); sp.set_linewidth(0.75)
    cb = ax.figure.colorbar(sca, ax = ax, fraction = 0.045, pad = 0.03)
    cb.set_label('score', fontsize = 6.5); cb.ax.tick_params(labelsize = 6)

def figure5():
    fig = plt.figure(figsize = (W2, mm(150)))
    gs = fig.add_gridspec(2, 2, hspace = 0.8, wspace = 0.5,
                          left = 0.1, right = 0.93, top = 0.94, bottom = 0.12)
    a = fig.add_subplot(gs[0, 0]); draw5a(a); panel_label(a, 'A')
    b = fig.add_subplot(gs[0, 1]); draw5b(b); panel_label(b, 'B')
    c = fig.add_subplot(gs[1, 0]); draw5c(c); panel_label(c, 'C')
    d = fig.add_subplot(gs[1, 1]); draw5d(d); panel_label(d, 'D')
    savefig(fig, 'Figure5', W2)
    for nm, fn, w, h in [('F5A_c1q_source', draw5a, W1, mm(52)), ('F5B_cellchat', draw5b, W1, mm(52)),
                         ('F5C_nichenet', draw5c, W1, mm(50)), ('F5D_umap', draw5d, W1, mm(52))]:
        f = plt.figure(figsize = (w, h)); ax = f.add_subplot(111); fn(ax); save_panel(f, nm)


def sig_label(ax, x, y, s, ha = 'left', fontsize = 9):
    if not s or s == 'ns':
        return
    dy = -0.32 * fontsize if set(s) <= {'*'} else 0
    ax.annotate(s, (x, y), textcoords = 'offset points', xytext = (0, dy), ha = ha, va = 'center',
                fontsize = fontsize, fontweight = 'bold')


def star_key(fig, y = 0.02, x = 0.5, ha = 'center', levels = None):
    # `levels` limits the key to the star grades actually drawn in that figure
    txt = {'*': '*  p<0.05', '**': '**  p<0.01', '***': '***  p<0.001', '****': '****  p<0.0001'}
    keys = levels if levels else list(txt)
    fig.text(x, y, '        '.join(txt[k] for k in keys),
             ha = ha, va = 'bottom', fontsize = 9, fontweight = 'bold',
             bbox = dict(boxstyle = 'round,pad=0.5', facecolor = '#f2f2f2', edgecolor = '#777777', linewidth = 0.8))


def hstars(ax, y, vals, ps, pad = 0.035, margin = 0.22):
    vals = np.asarray(vals, float)
    lo = min(0, vals.min()); hi = max(0, vals.max()); span = hi - lo
    for yi, v, p in zip(y, vals, ps):
        if v >= 0:
            sig_label(ax, v + span * pad, yi, stars(p), ha = 'left')
        else:
            sig_label(ax, v - span * pad, yi, stars(p), ha = 'right')
    ax.set_xlim(lo - span * margin, hi + span * margin)


def forest(ax, labels, hr, lo, hi, title, xlabel = 'HR per SD (95% CI)', ref = 1.0, logx = True, colors = None):
    y = np.arange(len(labels))[::-1]
    if colors is None:
        colors = [RCC if (l > ref and h > ref) or (l < ref and h < ref) else GRAY for l, h in zip(lo, hi)]
    for yi, l, h, c in zip(y, lo, hi, colors):
        ax.plot([l, h], [yi, yi], color = c, lw = 1.1)
    ax.scatter(hr, y, s = 20, c = colors, zorder = 3, edgecolors = 'k', linewidths = 0.4)
    ax.axvline(ref, color = 'k', lw = 0.6, ls = '--')
    ax.set_yticks(y); ax.set_yticklabels(labels)
    if logx:
        ax.set_xscale('log')
        ax.set_xticks([0.5, 1, 2, 4]); ax.set_xticklabels(['0.5','1','2','4'])
    ax.set_xlabel(xlabel); ax.set_title(title, fontsize = 8)
    style_ax(ax)


def draw6a(ax):
    # Estimate + 95% CI rather than bare bars: the PCa coefficient is +0.003, which has no
    # drawable bar height and would read as missing data instead of a null. Values come from
    # the association tables, not the 2-dp display strings in pca_vs_rcc_contrast_table.csv.
    def coef(path, module = 'complement_C1Q', outcome = 'cytotoxicity'):
        t = load(path)
        return t[(t['module'] == module) & (t['outcome'] == outcome)].iloc[0]
    r = coef(os.path.join(tab, 'phase3_adjusted_associations.csv'))
    p_ = coef(os.path.join(pca, 'pca_phase3_adjusted_associations.csv'))
    rows = [('ccRCC', r, RCC), ('PCa', p_, PCA)]
    y = np.arange(len(rows))[::-1]
    lo = min(x[1]['ci_low'] for x in rows); hi = max(x[1]['ci_high'] for x in rows)
    span = hi - lo
    ax.axvline(0, color = 'k', lw = 0.6, ls = '--', zorder = 0)
    for yi, (nm, v, c) in zip(y, rows):
        sig = v['p'] < 0.05
        # plain CI line + marker, matching panels B/C (no end caps)
        ax.plot([v['ci_low'], v['ci_high']], [yi, yi], color = c, lw = 1.1, zorder = 2)
        ax.scatter(v['coef'], yi, s = 26, c = c, zorder = 3, edgecolors = 'k',
                   linewidths = 0.4, alpha = 1.0 if sig else 0.55)
        fmt = '%+.3f' if abs(v['coef']) < 0.01 else '%+.2f'
        # centred on the estimate; white halo so it stays legible over the zero line.
        # U+2212 to match the axis tick labels
        ax.text(v['coef'], yi + 0.15, r'$\beta$ = ' + (fmt % v['coef']).replace('-', '−'),
                va = 'bottom', ha = 'center', fontsize = 6.5, zorder = 4,
                bbox = dict(boxstyle = 'square,pad=0.15', facecolor = 'white', edgecolor = 'none'))
        # stars in a fixed right-hand column, as in panels B and C; sig_label applies the
        # baseline shift that vertically centres asterisk glyphs
        if sig:
            sig_label(ax, hi + span * 0.08, yi, stars(v['p']), ha = 'left')
    ax.set_yticks(y)
    ax.set_yticklabels(['%s\n(n=%d)' % (nm, int(v['n'])) for nm, v, _ in rows])
    ax.set_ylim(-0.6, len(rows) - 0.25)
    ax.set_xlim(lo - span * 0.08, hi + span * 0.22)
    ax.set_xlabel(r'adjusted $\beta$ (95% CI)')
    ax.set_title('Complement C1q vs CD8 cytotoxicity', fontsize = 8)
    style_ax(ax)

def draw6b(ax):
    d = load(os.path.join(tab, 'phase3_adjusted_associations.csv'))
    d = d[d['module'] == 'complement_C1Q'].set_index('outcome')
    rows = ['cytotoxicity','CD8_exhaustion','MHC_II_APC']
    lab = ['CD8 cytotoxicity','CD8 exhaustion','MHC-II APC']
    y = np.arange(len(rows))[::-1]
    hi = d.loc[rows, 'ci_high'].max(); lo = d.loc[rows, 'ci_low'].min(); span = hi - lo
    sx = hi + span * 0.08
    for yi, o in zip(y, rows):
        r = d.loc[o]
        c = RCC if (r['ci_low'] > 0 and r['ci_high'] > 0) or (r['ci_low'] < 0 and r['ci_high'] < 0) else GRAY
        ax.plot([r['ci_low'], r['ci_high']], [yi, yi], color = c, lw = 1.1)
        ax.scatter(r['coef'], yi, s = 20, c = c, zorder = 3, edgecolors = 'k', linewidths = 0.4)
        sig_label(ax, sx, yi, stars(r['p']), ha = 'left')
        ax.text(r['coef'], yi + 0.15, r'$\beta$ = ' + ('%+.2f' % r['coef']).replace('-', '\u2212'),
                va = 'bottom', ha = 'center', fontsize = 6.5, zorder = 4,
                bbox = dict(boxstyle = 'square,pad=0.15', facecolor = 'white', edgecolor = 'none'))
    ax.axvline(0, color = 'k', lw = 0.6, ls = '--')
    ax.set_xlim(lo - span * 0.1, hi + span * 0.22)
    ax.set_ylim(-0.55, len(rows) - 0.35)
    ax.set_yticks(y); ax.set_yticklabels(lab)
    ax.set_xlabel(r'adjusted $\beta$ (95% CI)')
    ax.set_title('ccRCC immune context', fontsize = 8)
    style_ax(ax)

def draw6c(ax):
    d = load(os.path.join(tab, 'phase4_tcga_kirc_cox.csv'))
    d = d[(d['endpoint'] == 'OS') & (d['model'].str.startswith('adj'))]
    order = ['complement_C1Q','CLEC_LAM8','APOE_TREM2','MERTK_GPNMB','panTAM','Obradovic_TREM2']
    d = d.set_index('module').reindex(order).dropna(subset = ['HR_per_SD'])
    forest(ax, [mlabel(m) for m in d.index], d['HR_per_SD'].values, d['ci_low'].values, d['ci_high'].values,
           'TCGA-KIRC OS (adjusted)', xlabel = 'HR per SD (95% CI)')
    for yi, (hr, h, p) in zip(np.arange(len(d))[::-1], zip(d['HR_per_SD'], d['ci_high'], d['p'])):
        sig_label(ax, 4.2, yi, stars(p), ha = 'left')
        ax.text(hr, yi + 0.16, 'HR = %.2f' % hr, va = 'bottom', ha = 'center', fontsize = 6.5,
                zorder = 4, bbox = dict(boxstyle = 'square,pad=0.15', facecolor = 'white',
                                        edgecolor = 'none'))
    ax.set_xlim(0.45, 5); ax.set_ylim(-0.55, len(d) - 0.30)

def draw6d(ax):
    from lifelines.statistics import logrank_test
    d = pd.read_csv(os.path.join(cache, 'km_tcga.csv'))
    lr = logrank_test(d['OS.time'][d.grp=='high'], d['OS.time'][d.grp=='low'], d['OS'][d.grp=='high'], d['OS'][d.grp=='low'])
    km_panel(ax, d, 'OS.time', 'OS', 'TCGA-KIRC (n=%d)' % len(d), lr.p_value,
             risk_times = [0, 3, 6, 9, 12])

def draw6e(ax):
    cf = os.path.join(cache, 'km_cptac.csv')
    if not os.path.exists(cf):
        ax.text(0.5, 0.5, 'CPTAC data unavailable', ha = 'center', va = 'center', transform = ax.transAxes)
        style_ax(ax); return
    from lifelines.statistics import logrank_test
    d = pd.read_csv(cf)
    lr = logrank_test(d['time'][d.grp=='high'], d['time'][d.grp=='low'], d['event'][d.grp=='high'], d['event'][d.grp=='low'])
    km_panel(ax, d, 'time', 'event', 'CPTAC ccRCC (n=%d)' % len(d), lr.p_value,
             risk_times = [0, 1, 2, 3, 4, 5])

def fig6_star_levels():
    # the star grades actually drawn in panels A-C, so the key lists only those
    ps = []
    rcc = load(os.path.join(tab, 'phase3_adjusted_associations.csv'))
    rcc = rcc[rcc['module'] == 'complement_C1Q']
    ps += list(rcc[rcc['outcome'].isin(['cytotoxicity','CD8_exhaustion','MHC_II_APC'])]['p'])
    pc = load(os.path.join(pca, 'pca_phase3_adjusted_associations.csv'))
    ps += list(pc[(pc['module'] == 'complement_C1Q') & (pc['outcome'] == 'cytotoxicity')]['p'])
    cox = load(os.path.join(tab, 'phase4_tcga_kirc_cox.csv'))
    cox = cox[(cox['endpoint'] == 'OS') & (cox['model'].str.startswith('adj'))]
    ps += list(cox[cox['module'].isin(['complement_C1Q','RCC_skew_CORE','CLEC_LAM8',
                                       'MERTK_GPNMB','panTAM','Obradovic_TREM2'])]['p'])
    used = {stars(p) for p in ps if pd.notna(p)} - {'ns'}
    return [k for k in ['*','**','***','****'] if k in used]


def figure6():
    fig = plt.figure(figsize = (W2, mm(152)))
    # Two gridspecs rather than one 2x3: the top row keeps the original geometry and the
    # two-panel bottom row is centred under it. Spanning columns of a single grid would
    # make each panel swallow an internal gap and shift the panel labels off the canvas.
    # Row height h solves 2h + 0.62h = (0.93 - 0.16); column width w solves 3w + 1.24w = 0.89.
    gs_top = fig.add_gridspec(1, 3, wspace = 0.62, left = 0.08, right = 0.97,
                              top = 0.93, bottom = 0.631)
    gs_bot = fig.add_gridspec(1, 2, wspace = 0.62, left = 0.250, right = 0.800,
                              top = 0.459, bottom = 0.16)
    a = fig.add_subplot(gs_top[0, 0]); draw6a(a); panel_label(a, 'A', x = -0.28)
    b = fig.add_subplot(gs_top[0, 1]); draw6b(b); panel_label(b, 'B')
    c = fig.add_subplot(gs_top[0, 2]); draw6c(c); panel_label(c, 'C')
    d = fig.add_subplot(gs_bot[0, 0]); draw6d(d); panel_label(d, 'D')
    e = fig.add_subplot(gs_bot[0, 1]); draw6e(e); panel_label(e, 'E')
    star_key(fig, y = 0.05, levels = fig6_star_levels())
    savefig(fig, 'Figure6', W2)
    for nm, fn, w, h in [('F6A_cytotox', draw6a, mm(45), mm(52)), ('F6B_context', draw6b, W1, mm(50)),
                         ('F6C_os_forest', draw6c, W1, mm(50)), ('F6D_km_tcga', draw6d, W1, mm(52)),
                         ('F6E_km_cptac', draw6e, W1, mm(52))]:
        f = plt.figure(figsize = (w, h)); ax = f.add_subplot(111); fn(ax); save_panel(f, nm)


def draw7a(ax):
    d = load(os.path.join(mod, 'lasso_selected.csv')).rename(columns = {d0: 'gene' for d0 in ['Unnamed: 0']})
    d = d.sort_values('coef')
    y = np.arange(len(d))
    cols = [RCC if c > 0 else PCA for c in d['coef']]
    ax.barh(y, d['coef'], color = cols, height = 0.72)
    ax.set_yticks(y); ax.set_yticklabels(d['gene'], fontsize = 6)
    ax.axvline(0, color = 'k', lw = 0.6)
    ax.set_xlabel('LASSO-Cox coefficient')
    ax.set_title('25-gene signature', fontsize = 8)
    h = [plt.Rectangle((0,0),1,1,color=RCC), plt.Rectangle((0,0),1,1,color=PCA)]
    ax.legend(h, ['risk ↑','risk ↓'], loc = 'lower right', fontsize = 6.5)
    style_ax(ax)

def draw7b(ax):
    d = load(os.path.join(mod, 'discrimination_metrics.csv')).set_index('cohort')
    nm = load(os.path.join(mod, 'nomogram_discrimination.csv')).set_index('cohort')
    coh = ['TCGA-KIRC_train','TCGA-KIRC_heldout','CPTAC_ccRCC_primary','NOMO']
    lab = ['TCGA\ntraining','TCGA\nheld-out','CPTAC\nexternal','TCGA held-out\n+ clinical']
    vals = list(d.loc[coh[:3], 'cindex'].values) + [nm.loc['TCGA-KIRC_heldout', 'cindex']]
    cols = [RCC, RCC, RCC, ATF3]        # purple marks the combined model, not the signature alone
    x = np.arange(len(coh))
    dark = style.GRAD[3]
    ax.bar(x, vals, color = cols, width = 0.5)
    for xi, v in zip(x, vals):
        ax.text(xi, v - 0.065, '%.2f' % v, ha = 'center', va = 'center', color = 'white', fontsize = 6.5)
    auc_src = [d.loc[c, ['auc_1yr','auc_3yr','auc_5yr']].values.astype(float) for c in coh[:3]]
    auc_src.append(nm.loc['TCGA-KIRC_heldout', ['auc_1yr','auc_3yr','auc_5yr']].values.astype(float))
    for xi, aucs in zip(x, auc_src):
        tt = np.array([1, 3, 5]); m = ~np.isnan(aucs)
        ax.plot(xi + (tt[m] - 3) * 0.13, aucs[m], marker = 'o', ms = 3.5, color = dark, lw = 0.9, zorder = 3)
    ax.axhline(0.5, color = 'k', lw = 0.6, ls = '--')
    ax.set_xticks(x); ax.set_xticklabels(lab, fontsize = 6)
    ax.set_ylabel('C-index  /  time-AUC'); ax.set_ylim(0.4, 1.08)
    ax.set_title('discrimination', fontsize = 8)
    h = [plt.Rectangle((0, 0), 1, 1, color = RCC), plt.Rectangle((0, 0), 1, 1, color = ATF3),
         Line2D([0], [0], marker = 'o', color = dark, lw = 0.9, ms = 3.5)]
    ax.legend(h, ['signature alone', '+ stage/grade/age', 'time-AUC 1/3/5 yr'],
              loc = 'upper right', fontsize = 5.5, handlelength = 1.1,
              borderpad = 0.3, labelspacing = 0.3)
    style_ax(ax)

def draw7c(ax):
    b = load(os.path.join(mod, 'fair_benchmark.csv')).set_index('comparator')
    rnd = load(os.path.join(mod, 'fair_benchmark_random_draws.csv'))['cindex'].values
    names = ['25-gene signature','Obradovic','CLEC_LAM','panTAM']
    vals = [b.loc['25-gene signature (frozen)', 'heldout_cindex'],
            b.loc['Obradovic TREM2 genes (refit, same pipeline)', 'heldout_cindex'],
            b.loc['CLEC_LAM8 core (refit, same pipeline)', 'heldout_cindex'],
            b.loc['panTAM genes (refit, same pipeline)', 'heldout_cindex']]
    cols = [RCC, GRAY, GRAY, GRAY]
    x = np.arange(len(names))
    ax.bar(x, vals, color = cols, width = 0.6)
    sc = load(os.path.join(mod, 'sanity_checks.csv')).set_index('metric')['value']
    p95 = float(sc['random_gene_cindex_p95'])
    ax.axhspan(0.5, p95, color = GRAY, alpha = 0.2)
    ax.axhline(p95, color = GRAY, lw = 0.8, ls = ':')
    ax.text(len(names) - 0.4, p95 + 0.005, '95th percentile, random gene sets', ha = 'right',
            fontsize = 5.5, color = '#555')
    for xi, v in zip(x, vals):
        ax.text(xi, v + 0.006, '%.2f' % v, ha = 'center', fontsize = 6.5)
    ax.set_xticks(x); ax.set_xticklabels(names, rotation = 25, ha = 'right', fontsize = 6.5)
    ax.set_ylabel('held-out C-index'); ax.set_ylim(0.45, 0.75)
    ax.set_title('benchmark vs prior signatures', fontsize = 8)
    style_ax(ax)

def draw7d(ax):
    d = load(os.path.join(mod, 'pca_univariate_bcr.csv'))
    d = d[d['program'] == 'partA_signature']
    lab = {'PFI': 'progression-free\ninterval', 'biochemical_recurrence': 'biochemical\nrecurrence'}
    d = d.set_index('endpoint')
    order = ['PFI','biochemical_recurrence']
    y = np.arange(len(order))[::-1]
    for yi, o in zip(y, order):
        r = d.loc[o]
        ax.plot([r['ci_low'], r['ci_high']], [yi, yi], color = RCC, lw = 1.1)
        ax.scatter(r['HR'], yi, s = 22, c = RCC, zorder = 3, edgecolors = 'k', linewidths = 0.4)
        ax.text(2.02, yi, 'HR %.2f' % r['HR'], va = 'center', ha = 'left', fontsize = 6.5)
        sig_label(ax, 2.42, yi, stars(r['p']), ha = 'left')
    ax.axvline(1, color = 'k', lw = 0.6, ls = '--')
    ax.set_yticks(y); ax.set_yticklabels([lab[o] for o in order])
    ax.set_xticks([1.0, 1.5, 2.0])
    ax.set_xlabel('HR per SD (95% CI)')
    ax.set_xlim(0.9, 2.6); ax.set_ylim(-0.6, 1.6)
    ax.set_title('cross-tumor validation: TCGA-PRAD', fontsize = 8)
    style_ax(ax)

def fig7_star_levels():
    d = load(os.path.join(mod, 'pca_univariate_bcr.csv'))
    ps = list(d[d['program'] == 'partA_signature']['p'].dropna())
    used = {stars(p) for p in ps} - {'ns'}
    return [k for k in ['*','**','***','****'] if k in used]


def figure7():
    fig = plt.figure(figsize = (W2, mm(156)))
    gs = fig.add_gridspec(3, 2, hspace = 0.7, wspace = 0.42,
                          left = 0.12, right = 0.95, top = 0.94, bottom = 0.15,
                          width_ratios = [1, 1.2], height_ratios = [1, 1, 0.62])
    a = fig.add_subplot(gs[:, 0]); draw7a(a); panel_label(a, 'A', x = -0.22, y = 1.02)
    b = fig.add_subplot(gs[0, 1]); draw7b(b); panel_label(b, 'B')
    c = fig.add_subplot(gs[1, 1]); draw7c(c); panel_label(c, 'C')
    d = fig.add_subplot(gs[2, 1]); draw7d(d); panel_label(d, 'D')
    star_key(fig, levels = fig7_star_levels())
    savefig(fig, 'Figure7', W2)
    for nm, fn, w, h in [('F7A_lasso', draw7a, mm(60), mm(95)), ('F7B_discrimination', draw7b, W1, mm(52)),
                         ('F7C_benchmark', draw7c, W1, mm(52)), ('F7D_prostate', draw7d, W1, mm(40))]:
        f = plt.figure(figsize = (w, h)); ax = f.add_subplot(111); fn(ax); save_panel(f, nm)


def prep_model_diag():
    cf = os.path.join(cache, 'model_diag.npz')
    if os.path.exists(cf):
        return np.load(cf, allow_pickle = True)
    from sksurv.util import Surv
    from sksurv.linear_model import CoxnetSurvivalAnalysis
    from sklearn.model_selection import KFold
    m = json.load(open(os.path.join(mod, 'risk_model.json')))
    g = np.array(m['genes']); mu = np.array(m['mean']); sd = np.array(m['sd']); co = np.array(m['coef'])
    k = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_HiSeqV2.gz'), sep = '\t', index_col = 0)
    k.index = [x.upper() for x in k.index]
    s = pd.read_csv(os.path.join(root, 'resources/tcga/KIRC_survival.txt'), sep = '\t')
    s['bcr'] = s['sample'].str[:15]
    x = k.loc[[q for q in g if q in k.index]].T
    z = (x[list(g)] - mu) / sd
    d = pd.DataFrame({'risk': z.values @ co, 'sid': x.index}, index = x.index)
    d['bcr'] = [i[:15] for i in d.index]
    d = d.merge(s[['bcr','OS','OS.time']], on = 'bcr')
    d = d[d['OS.time'] > 0].dropna(subset = ['OS','OS.time'])
    Zd = z.loc[d['sid']].values
    y = Surv.from_arrays(d['OS'].astype(bool), d['OS.time'])
    cv = CoxnetSurvivalAnalysis(l1_ratio = 1.0, alpha_min_ratio = 0.01, n_alphas = 60, max_iter = 100000)
    cv.fit(Zd, y)
    alphas = cv.alphas_
    coefs = cv.coef_
    nnz = (np.abs(coefs) > 0).sum(axis = 0)
    from sksurv.metrics import concordance_index_censored
    kf = KFold(5, shuffle = True, random_state = 0)
    cidx = np.zeros((len(alphas),))
    counts = np.zeros((len(alphas),))
    Z = Zd; ev = d['OS'].astype(bool).values; ti = d['OS.time'].values
    for tr, te in kf.split(Z):
        mdl = CoxnetSurvivalAnalysis(l1_ratio = 1.0, alphas = alphas, max_iter = 100000)
        try:
            mdl.fit(Z[tr], Surv.from_arrays(ev[tr], ti[tr]))
        except Exception:
            continue
        for j, al in enumerate(alphas):
            try:
                pr = mdl.predict(Z[te], alpha = al)
                cidx[j] += concordance_index_censored(ev[te], ti[te], pr)[0]; counts[j] += 1
            except Exception:
                pass
    cidx = np.divide(cidx, counts, out = np.full_like(cidx, np.nan), where = counts > 0)
    np.savez(cf, alphas = alphas, nnz = nnz, cidx = cidx,
             risk = d['risk'].values, OS = d['OS'].values, time = d['OS.time'].values)
    return np.load(cf, allow_pickle = True)

def draw_s1_lasso(ax):
    dd = prep_model_diag()
    al = dd['alphas']; cidx = dd['cidx']
    ax.plot(np.log10(al), cidx, marker = 'o', ms = 2.5, color = RCC, lw = 0.9)
    j = int(np.nanargmax(cidx))
    ax.axvline(np.log10(al[j]), color = 'k', ls = '--', lw = 0.7)
    ax.set_xlabel('log10 λ'); ax.set_ylabel('10-fold CV C-index')
    ax.set_title('LASSO-Cox regularization path', fontsize = 8)
    ax2 = ax.twinx()
    ax2.plot(np.log10(al), dd['nnz'], color = GRAY, lw = 0.8)
    ax2.set_ylabel('non-zero genes', color = GRAY, fontsize = 7)
    ax2.tick_params(axis = 'y', labelcolor = GRAY, labelsize = 6.5)
    ax2.spines['top'].set_visible(False)
    style_ax(ax)

def draw_s1_calib(ax):
    from lifelines import KaplanMeierFitter
    dd = prep_model_diag()
    d = pd.DataFrame({'risk': dd['risk'], 'e': dd['OS'], 't': dd['time']})
    d['q'] = pd.qcut(d['risk'], 4, labels = False)
    yrs = 3 * 365.25
    obs = []; pos = []
    for qi in range(4):
        sub = d[d['q'] == qi]
        kmf = KaplanMeierFitter().fit(sub['t'], sub['e'])
        obs.append(1 - kmf.predict(yrs))
        pos.append(qi)
    pred = [np.mean(d[d['q'] == qi]['risk']) for qi in range(4)]
    pr = (np.array(pred) - min(pred)) / (max(pred) - min(pred))
    ax.plot(pr, obs, marker = 'o', ms = 4, color = RCC, lw = 1.0)
    ax.plot([0, 1], [min(obs), max(obs)], color = GRAY, ls = '--', lw = 0.7)
    ax.set_xlabel('relative predicted risk (quartile)'); ax.set_ylabel('observed 3-yr mortality')
    ax.set_title('calibration', fontsize = 8)
    style_ax(ax)

def draw_s1_dca(ax):
    from lifelines import CoxPHFitter
    dd = prep_model_diag()
    d = pd.DataFrame({'risk': dd['risk'], 'e': dd['OS'].astype(int), 't': dd['time']})
    tt = 3 * 365.25
    cph = CoxPHFitter().fit(d, 't', 'e')
    surv = cph.predict_survival_function(d, times = [tt]).T.values.ravel()
    risk = 1 - surv
    thr = np.linspace(0.01, 0.6, 60)
    ev = ((d['t'] <= tt) & (d['e'] == 1)).values
    n = len(d)
    nb_model = []; nb_all = []
    prev = ev.mean()
    for pt in thr:
        tp = ((risk >= pt) & ev).sum() / n
        fp = ((risk >= pt) & ~ev).sum() / n
        nb_model.append(tp - fp * (pt / (1 - pt)))
        nb_all.append(prev - (1 - prev) * (pt / (1 - pt)))
    ax.plot(thr, nb_model, color = RCC, lw = 1.1, label = 'risk model')
    ax.plot(thr, nb_all, color = GRAY, lw = 0.8, label = 'treat all')
    ax.axhline(0, color = 'k', lw = 0.7, label = 'treat none')
    ax.set_ylim(-0.05, max(nb_model) * 1.3 + 0.01)
    ax.set_xlabel('threshold probability'); ax.set_ylabel('net benefit (3-yr)')
    ax.set_title('decision-curve analysis', fontsize = 8)
    ax.legend(fontsize = 6.5, loc = 'upper right')
    style_ax(ax)

def draw_s1_nomogram(ax):
    d = load(os.path.join(mod, 'multivariable_cox.csv'))
    d = d[d['covariate'] != 'grade_z'] if False else d
    d = d.set_index('covariate')
    order = ['risk_z','stage_z','age_z','grade_z']
    lab = {'risk_z': 'risk score', 'stage_z': 'tumor stage', 'age_z': 'age', 'grade_z': 'grade'}
    d = d.reindex(order)
    pts = d['coef'] / d['coef'].max() * 100
    y = np.arange(len(order))[::-1]
    for yi, cov in zip(y, order):
        c = ATF3 if d.loc[cov, 'p'] < 0.05 else GRAY
        ax.plot([0, pts[cov]], [yi, yi], color = c, lw = 4, solid_capstyle = 'round')
        ax.text(pts[cov] + 2, yi, '%.0f pts' % pts[cov], va = 'center', fontsize = 6.5)
        sig_label(ax, pts[cov] + 26, yi, stars(d.loc[cov, 'p']), ha = 'left')
    ax.set_yticks(y); ax.set_yticklabels([lab[o] for o in order])
    ax.set_xlabel('points (per SD, multivariable Cox)')
    ax.set_xlim(0, 152)
    ax.set_title('nomogram contributions', fontsize = 8)
    style_ax(ax)

def figureS1():
    fig = plt.figure(figsize = (W2, mm(142)))
    gs = fig.add_gridspec(2, 2, hspace = 0.42, wspace = 0.38,
                          left = 0.1, right = 0.95, top = 0.93, bottom = 0.16)
    a = fig.add_subplot(gs[0, 0]); draw_s1_lasso(a); panel_label(a, 'A')
    b = fig.add_subplot(gs[0, 1]); draw_s1_calib(b); panel_label(b, 'B')
    c = fig.add_subplot(gs[1, 0]); draw_s1_dca(c); panel_label(c, 'C')
    d = fig.add_subplot(gs[1, 1]); draw_s1_nomogram(d); panel_label(d, 'D')
    star_key(fig)
    savefig(fig, 'FigureS1', W2)
    for nm, fn, w, h in [('S1A_lasso_cv', draw_s1_lasso, W1, mm(52)), ('S1B_calibration', draw_s1_calib, W1, mm(52)),
                         ('S1C_dca', draw_s1_dca, W1, mm(52)), ('S1D_nomogram', draw_s1_nomogram, W1, mm(50))]:
        f = plt.figure(figsize = (w, h)); ax = f.add_subplot(111); fn(ax); save_panel(f, nm)


def draw_s2_drug(ax):
    d = load(os.path.join(mod, 'riskgroup_drug.csv')).sort_values('delta')
    y = np.arange(len(d))
    cols = [RCC if p < 0.05 else GRAY for p in d['fdr']]
    ax.barh(y, d['delta'], color = cols, height = 0.6)
    ax.axvline(0, color = 'k', lw = 0.6)
    ax.set_yticks(y); ax.set_yticklabels(d['agent'])
    hstars(ax, y, d['delta'].values, d['fdr'].values)
    ax.set_xlabel('Δ predicted log-IC50 (high − low)')
    ax.set_title('drug sensitivity (oncoPredict)', fontsize = 8)
    style_ax(ax)

def draw_s2_tide(ax):
    d = load(os.path.join(mod, 'riskgroup_tide.csv')).sort_values('delta')
    y = np.arange(len(d))
    cols = [RCC if p < 0.05 else GRAY for p in d['fdr']]
    ax.barh(y, d['delta'], color = cols, height = 0.62)
    ax.axvline(0, color = 'k', lw = 0.6)
    ax.set_yticks(y); ax.set_yticklabels(d['metric'])
    hstars(ax, y, d['delta'].values, d['fdr'].values)
    ax.set_xlabel('Δ score (high − low risk)')
    ax.set_title('TIDE immune evasion', fontsize = 8)
    style_ax(ax)

def draw_s2_enrich(ax):
    sig = load(os.path.join(mod, 'enrichment_signature.csv')).head(6)
    deg = load(os.path.join(mod, 'enrichment_degs.csv')).head(6)
    def clean(t):
        import re
        return re.sub(r' \(GO:\d+\)', '', t)[:38]
    y1 = np.arange(len(deg))
    v = -np.log10(deg['Adjusted P-value'].values)
    ax.barh(y1, v, color = C1Q, height = 0.72)
    for yi, t in zip(y1, deg['Term']):
        ax.text(0.4, yi, clean(t), va = 'center', ha = 'left', fontsize = 6, color = 'white')
    ax.set_yticks([])
    ax.axvline(-np.log10(0.05), color = 'k', ls = ':', lw = 0.7)
    ax.set_xlabel('−log10 adjusted p')
    ax.set_title('enrichment: TAM-vs-mono DEGs', fontsize = 8)
    ax.set_xlim(0, v.max() * 1.05)
    ax.invert_yaxis()
    style_ax(ax)

def draw_s2_immune(ax):
    d = load(os.path.join(mod, 'riskgroup_immune.csv')).sort_values('delta')
    y = np.arange(len(d))
    cols = [RCC if p < 0.05 else GRAY for p in d['fdr']]
    ax.barh(y, d['delta'], color = cols, height = 0.62)
    ax.axvline(0, color = 'k', lw = 0.6)
    ax.set_yticks(y); ax.set_yticklabels([s.replace('_',' ') for s in d['signature']], fontsize = 6.5)
    hstars(ax, y, d['delta'].values, d['fdr'].values)
    ax.set_xlabel('Δ score (high − low risk)')
    ax.set_title('risk-group immune infiltration', fontsize = 8)
    style_ax(ax)

def figureS2():
    fig = plt.figure(figsize = (W2, mm(140)))
    gs = fig.add_gridspec(2, 2, hspace = 0.72, wspace = 0.5,
                          left = 0.13, right = 0.96, top = 0.94, bottom = 0.14)
    a = fig.add_subplot(gs[0, 0]); draw_s2_drug(a); panel_label(a, 'A')
    b = fig.add_subplot(gs[0, 1]); draw_s2_tide(b); panel_label(b, 'B')
    c = fig.add_subplot(gs[1, 0]); draw_s2_enrich(c); panel_label(c, 'C', x = -0.06)
    d = fig.add_subplot(gs[1, 1]); draw_s2_immune(d); panel_label(d, 'D', x = -0.4)
    star_key(fig)
    savefig(fig, 'FigureS2', W2)
    for nm, fn, w, h in [('S2A_drug', draw_s2_drug, W1, mm(50)), ('S2B_tide', draw_s2_tide, W1, mm(50)),
                         ('S2C_enrichment', draw_s2_enrich, W1, mm(52)), ('S2D_immune', draw_s2_immune, W1, mm(52))]:
        f = plt.figure(figsize = (w, h)); ax = f.add_subplot(111); fn(ax); save_panel(f, nm)


CNV = os.path.join(out, 'cnv', 'infercnv')
CNV_ORDER = ['Tumor','Pericyte','Endothelial','MSC','Myeloid_TAM','T_NK']
CNV_LAB = {'Tumor':'malignant','Pericyte':'pericyte','Endothelial':'endothelial',
           'MSC':'MSC','Myeloid_TAM':'myeloid/TAM','T_NK':'T/NK\n(reference)'}

def _cnv():
    d = load(os.path.join(CNV, 'cnv_score_per_cell.csv'))
    d['sample'] = d['cell'].str.replace(r'_[^_]+$', '', regex = True)
    return d

def draw_s3a(ax):
    d = _cnv()
    ref = d.loc[d.compartment == 'T_NK', 'cnv_score'].median()
    data = [d.loc[d.compartment == c, 'cnv_score'].values * 1e3 for c in CNV_ORDER]
    cols = [PCA if c == 'Tumor' else (C1Q if c == 'Myeloid_TAM' else (RCC if c == 'T_NK' else GRAY))
            for c in CNV_ORDER]
    bp = ax.boxplot(data, positions = np.arange(len(CNV_ORDER)), widths = 0.6, showfliers = False,
                    patch_artist = True, medianprops = dict(color = 'k', lw = 1.0),
                    whiskerprops = dict(lw = 0.7), capprops = dict(lw = 0.7),
                    boxprops = dict(lw = 0.7))
    for b, c in zip(bp['boxes'], cols):
        b.set_facecolor(c); b.set_alpha(0.85)
    ax.axhline(ref * 1e3, color = RCC, lw = 0.8, ls = '--')
    ax.set_xticks(np.arange(len(CNV_ORDER)))
    ax.set_xticklabels([CNV_LAB[c] for c in CNV_ORDER], rotation = 30, ha = 'right', fontsize = 6.5)
    ax.set_ylabel(r'CNV burden ($\times 10^{-3}$)')
    ax.set_title('per-cell inferred CNV burden', fontsize = 8)
    style_ax(ax)

def draw_s3b(ax):
    d = _cnv()
    ref = d.loc[d.compartment == 'T_NK', 'cnv_score'].median()
    med = d.groupby('compartment')['cnv_score'].median()
    order = [c for c in CNV_ORDER if c != 'T_NK']
    y = np.arange(len(order))[::-1]
    vals = [med[c] / ref for c in order]
    cols = [PCA if c == 'Tumor' else (C1Q if c == 'Myeloid_TAM' else GRAY) for c in order]
    ax.barh(y, vals, color = cols, height = 0.62)
    ax.axvline(1, color = RCC, lw = 0.8, ls = '--')
    for yi, v in zip(y, vals):
        ax.text(v + 0.05, yi, '%.2f' % v, va = 'center', fontsize = 6.5)
    ax.set_yticks(y); ax.set_yticklabels([CNV_LAB[c] for c in order], fontsize = 6.5)
    ax.set_xlabel('median CNV burden / T-NK reference')
    ax.set_xlim(0, max(vals) * 1.25)
    ax.set_title('fold over diploid reference', fontsize = 8)
    style_ax(ax)

def draw_s3c(ax):
    d = _cnv()
    rows = []
    for s_, g in d.groupby('sample'):
        t = g.loc[g.compartment == 'Tumor', 'cnv_score']
        m = g.loc[g.compartment == 'Myeloid_TAM', 'cnv_score']
        if len(t) >= 10 and len(m) >= 10:
            rows.append((s_.replace('RCC-', '').replace('-Tumor', ''), m.median() * 1e3, t.median() * 1e3))
    lim = max(max(r[1] for r in rows), max(r[2] for r in rows)) * 1.18
    ax.plot([0, lim], [0, lim], color = GRAY, lw = 0.8, ls = '--')
    for nm_, x_, y_ in rows:
        ax.scatter(x_, y_, s = 26, color = PCA, edgecolors = 'k', linewidths = 0.4, zorder = 3)
        ax.annotate(nm_, (x_, y_), textcoords = 'offset points', xytext = (4, 3), fontsize = 5.5)
    ax.set_xlim(0, lim); ax.set_ylim(0, lim)
    ax.set_xlabel(r'myeloid/TAM ($\times 10^{-3}$)')
    ax.set_ylabel(r'malignant ($\times 10^{-3}$)')
    ax.set_title('median CNV burden, per sample', fontsize = 8)
    style_ax(ax)

def figureS3():
    fig = plt.figure(figsize = (W2, mm(62)))
    gs = fig.add_gridspec(1, 3, wspace = 0.55, left = 0.09, right = 0.98, top = 0.86, bottom = 0.28)
    a = fig.add_subplot(gs[0, 0]); draw_s3a(a); panel_label(a, 'A')
    b = fig.add_subplot(gs[0, 1]); draw_s3b(b); panel_label(b, 'B')
    c = fig.add_subplot(gs[0, 2]); draw_s3c(c); panel_label(c, 'C')
    savefig(fig, 'FigureS3', W2)
    for nm_, fn, w, h in [('S3A_cnv_burden', draw_s3a, W1, mm(52)),
                          ('S3B_cnv_fold', draw_s3b, W1, mm(48)),
                          ('S3C_cnv_persample', draw_s3c, W1, mm(52))]:
        f = plt.figure(figsize = (w, h)); ax = f.add_subplot(111); fn(ax); save_panel(f, nm_)


if __name__ == '__main__':
    which = sys.argv[1:] or ['4','5','6','7','S1','S2','S3']
    if any(x in which for x in ['6']):
        prep_km()
    reg = {'4': figure4, '5': figure5, '6': figure6, '7': figure7, 'S1': figureS1,
           'S2': figureS2, 'S3': figureS3}
    for w in which:
        reg[w]()
        print('done figure', w)
