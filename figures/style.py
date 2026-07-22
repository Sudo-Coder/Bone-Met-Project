import glob
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

for pat in ['/autofs/projects-t3/hussain/scProj/envs/rcc_reinterp_venv/**/*iberation*.ttf',
            '/usr/share/fonts/**/*iberationSans*.ttf',
            '/usr/share/fonts/urw-base35/NimbusSans-*.otf']:
    for f in glob.glob(pat, recursive = True):
        try:
            fm.fontManager.addfont(f)
        except Exception:
            pass

_have = {f.name for f in fm.fontManager.ttflist}
_pref = [n for n in ['Arial','Helvetica','Liberation Sans','Nimbus Sans','DejaVu Sans'] if n in _have]
FONT = _pref[0]

matplotlib.rcParams.update({
    'font.family': _pref,
    'pdf.fonttype': 42,
    'ps.fonttype': 42,
    'svg.fonttype': 'none',
    'axes.titlesize': 8,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'legend.fontsize': 7,
    'font.size': 7,
    'axes.linewidth': 0.75,
    'xtick.major.width': 0.75,
    'ytick.major.width': 0.75,
    'xtick.direction': 'out',
    'ytick.direction': 'out',
    'axes.spines.top': False,
    'axes.spines.right': False,
    'lines.linewidth': 1.0,
    'legend.frameon': False,
    'figure.dpi': 200,
    'savefig.dpi': 600,
    'savefig.bbox': 'standard',
})

PCA = '#E69F00'
RCC = '#0072B2'
ATF3 = '#7B3294'
C1Q = '#1B9E77'
CLEC = '#0072B2'
GRAY = '#999999'
GRAD = ['#FDE0C8', '#FDBE85', '#E6550D', '#A63603', '#7F2704']

from matplotlib.colors import LinearSegmentedColormap
GRAD_CMAP = LinearSegmentedColormap.from_list('benign_tumor', GRAD)

def mm(x):
    return x / 25.4

W1 = mm(89)
W15 = mm(120)
W2 = mm(183)
HMAX = mm(247)

def panel_label(ax, letter, x = -0.16, y = 1.06):
    ax.text(x, y, letter, transform = ax.transAxes, fontsize = 15, fontweight = 'bold',
            va = 'bottom', ha = 'right')

def style_ax(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.tick_params(length = 2.5, width = 0.75)
    return ax

def save_panel(fig, name):
    import os
    d = os.path.join(os.path.dirname(__file__), 'panels')
    os.makedirs(d, exist_ok = True)
    fig.savefig(os.path.join(d, name + '.pdf'))
    fig.savefig(os.path.join(d, name + '.svg'))
    plt.close(fig)

def stars(p):
    if p < 1e-4: return '****'
    if p < 1e-3: return '***'
    if p < 1e-2: return '**'
    if p < 0.05: return '*'
    return 'ns'


# Analysis keys are kept stable in the scripts; these are the names that appear in figures.
MODULE_LABEL = {
    'complement_C1Q':  'Complement C1q',
    'CLEC_LAM8':       'CLEC-LAM8',
    'APOE_TREM2':      'Lipid handling',
    'MERTK_GPNMB':     'MERTK/GPNMB',
    'SPP1_TAM':        'SPP1$^+$ TAM',
    'panTAM':          'panTAM',
    'ATF3_NFkB':       'ATF3/NF-\u03baB',
    'Obradovic_TREM2': 'Obradovic-TREM2 Sig',
}

# Retained in the analysis as a sensitivity module but not reported: it was defined after
# seeing which genes were RCC-specific, and correlates r = 0.94 with complement C1q.
MODULE_HIDE = {'RCC_skew_CORE', 'complement_C1Q_C3', 'inflammatory_mono', 'MHC_I_APM'}

def mlabel(m):
    return MODULE_LABEL.get(m, m.replace('_', '/'))
