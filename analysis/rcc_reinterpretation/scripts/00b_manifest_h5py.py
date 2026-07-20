import h5py, numpy as np

paths = [
    "kidney-cancer/Cleaned_Data/integrated-with-stromal.h5ad",
    "prostate-cancer/Cleaned_Data/integrated.h5ad",
    "prostate-cancer/Cleaned_Data/integrated_with_kfoury_labels.h5ad",
    "prostate-cancer/integrated.h5ad",
]
CHECK=["ATF3","C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","MERTK","FOLR2","SELENOP","ABCA1","C3","SERPING1","CA9","PAX8","CD68","CD163"]

def cat_summary(grp, key):

    try:
        node = grp[key]
        if isinstance(node, h5py.Group) and 'categories' in node and 'codes' in node:
            cats = node['categories'][:]
            cats = [c.decode() if isinstance(c,bytes) else c for c in cats]
            codes = node['codes'][:]
            if len(cats) <= 40:
                vc = {}
                for c in codes:
                    if c<0: continue
                    vc[cats[c]] = vc.get(cats[c],0)+1
                return f"({len(cats)} uniq) " + ", ".join(f"{k}={v}" for k,v in sorted(vc.items(), key=lambda x:-x[1]))
            else:
                return f"({len(cats)} uniq categorical)"
        elif isinstance(node, h5py.Dataset):
            arr = node[:]
            uq = np.unique(arr)
            if len(uq)<=40:
                return f"({len(uq)} uniq) numeric/str"
            return f"({len(uq)} uniq)"
    except Exception as e:
        return f"err {e}"
    return "?"

for p in paths:
    print("="*90); print("FILE:",p)
    try:
        f = h5py.File(p,"r")
    except Exception as e:
        print("  ERR", e); continue

    var = f['var']
    idxkey = var.attrs.get('_index', b'_index')
    if isinstance(idxkey,bytes): idxkey=idxkey.decode()
    vn = var[idxkey][:]
    vn = [x.decode() if isinstance(x,bytes) else x for x in vn]
    obs = f['obs']
    oidx = obs.attrs.get('_index', b'_index')
    if isinstance(oidx,bytes): oidx=oidx.decode()
    n_obs = obs[oidx].shape[0]
    print(f"  n_obs={n_obs}  n_vars={len(vn)}")
    print("  obs keys:", [k for k in obs.keys()])
    for k in obs.keys():
        if k==oidx: continue
        s = cat_summary(obs,k)
        if s and 'uniq' in s and not s.startswith('(') is False:

            if 'categorical)' not in s and s.count('=')>0:
                print(f"    [{k}] {s}")
    vset=set(vn)
    print("  present:", [g for g in CHECK if g in vset])
    print("  ABSENT:", [g for g in CHECK if g not in vset])
    if 'raw' in f: print("  has raw, raw n_vars:", len(f['raw']['var'][ f['raw']['var'].attrs.get('_index','_index') if isinstance(f['raw']['var'].attrs.get('_index','_index'),str) else f['raw']['var'].attrs.get('_index','_index').decode() ][:]))
    print("  obsm:", list(f['obsm'].keys()) if 'obsm' in f else [])
    print("  layers:", list(f['layers'].keys()) if 'layers' in f else [])
    f.close()
