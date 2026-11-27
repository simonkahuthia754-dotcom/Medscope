#!/usr/bin/env python3
"""screen.py - MILLION-SCALE virtual screening funnel over the whole ChEMBL small-molecule collection (~2.4 M structures), for MANY diseases in one pass.
  python screen.py <disease key|all> [--limit 300000] [--top 15] [--keys lung-cancer,hiv,...]               (MedScope: screen ...)
Funnel:  stream every structure -> property window + no PAINS (parallel, all CPU cores) -> Morgan-count features -> per-disease random-forest potency model
         trained on MEASURED ChEMBL data for that target -> keep the best 500 per disease -> ADMET-AI organ-risk scoring of the survivors -> ranked table + CSV.
Honest scope: potency comes from a model trained on public measurements (error bars from scaffold-grouped CV), so it finds analogues of what is already known
to hit that target; targets with no measured data need docking/structure methods (much slower). Output = candidates to TEST, not drugs."""
import csv, gzip, heapq, os, re, sys, time, urllib.request
import state

FTP = "https://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/latest/"
OUT = os.path.join(os.getcwd(), "medscope_out")

def chemreps_url():
    html = urllib.request.urlopen(FTP, timeout=60).read().decode(); m = re.findall(r'href="(chembl_\d+_chemreps\.txt\.gz)"', html)
    if not m: raise RuntimeError("ChEMBL chemreps file not found on the EBI FTP")
    return FTP + sorted(set(m))[-1]

def stream(limit=None, url=None):
    """yield (chembl_id, smiles) for every structure, streaming (no 250 MB download to disk)"""
    r = urllib.request.urlopen(url or chemreps_url(), timeout=120); n = 0
    with gzip.GzipFile(fileobj=r) as f:
        next(f)
        for line in f:
            p = line.decode("utf-8", "ignore").rstrip("\n").split("\t")
            if len(p) >= 2 and p[1]:
                yield p[0], p[1]; n += 1
                if limit and n >= limit: return

_W = {}
def _init():
    from rdkit import RDLogger
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    RDLogger.DisableLog("rdApp.*"); p = FilterCatalogParams(); p.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS); _W["pains"] = FilterCatalog(p)
    from rdkit.Chem import rdFingerprintGenerator; _W["gen"] = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)

def featurize(chunk):
    """worker: [(id, smiles)] -> (ids, smiles, uint8 count matrix, n_total) for molecules inside the property window and without PAINS"""
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    if not _W: _init()
    ids, sm, X = [], [], []
    for cid, s in chunk:
        try:
            if len(s) > 120: continue
            m = Chem.MolFromSmiles(s)
            if m is None: continue
            frags = Chem.GetMolFrags(m)
            if len(frags) > 1: continue                           # salts / mixtures are handled by their parent entries
            mw = Descriptors.MolWt(m)
            if not 150 <= mw <= 650: continue
            if not -1.5 <= Crippen.MolLogP(m) <= 6 or rdMolDescriptors.CalcTPSA(m) > 150 or rdMolDescriptors.CalcNumHBD(m) > 5: continue
            if _W["pains"].HasMatch(m): continue
            ids.append(cid); sm.append(s); X.append(np.minimum(_W["gen"].GetCountFingerprintAsNumPy(m), 255).astype(np.uint8))
        except Exception: continue
    return ids, sm, (np.array(X, dtype=np.uint8) if X else np.zeros((0, 1024), np.uint8)), len(chunk)

def vec(smi):
    import numpy as np
    from rdkit import Chem
    from rdkit.Chem import rdFingerprintGenerator
    m = Chem.MolFromSmiles(smi)
    return None if m is None else np.minimum(rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024).GetCountFingerprintAsNumPy(m), 255).astype(np.uint8)

def train(key, n_trees=100, seed=7):
    """random forest on MEASURED ChEMBL pIC50 for this disease's target -> (model, training dict, strict-CV rmse) or None"""
    import numpy as np, diseases, lab, potency
    from sklearn.ensemble import RandomForestRegressor
    D = diseases.get(key); cid, nm = diseases.resolve_chembl(D)
    if not cid: return None
    data = potency.fetch(cid, overlay=lab.training_overlay(D["key"])); rows = [(s, y, vec(s)) for s, y in data.items()]; rows = [r for r in rows if r[2] is not None]
    if len(rows) < 80: return None
    X = np.array([r[2] for r in rows], dtype=np.float32); y = np.array([r[1] for r in rows])
    rf = RandomForestRegressor(n_estimators=n_trees, min_samples_leaf=1, max_features=.3, n_jobs=-1, random_state=seed).fit(X, y)
    cv = None
    try:
        M = potency.Model(data); cv = (M.rmse_c, M.auc_c)
    except Exception: pass
    return dict(key=D["key"], name=D["name"], target=nm, rf=rf, data=data, n=len(rows), cv=cv)

def scan(models, limit=None, top=500, block=60000, log=print):
    """one pass over the library for ALL diseases; returns ({key: [(score, id, smiles)]}, n_streamed, n_scored)"""
    import numpy as np, multiprocessing as mp
    heaps = {m["key"]: [] for m in models}; n_in = n_ok = 0; t0 = time.time(); buf_i, buf_s, buf_x = [], [], []
    def flush():
        nonlocal n_ok
        if not buf_x: return
        X = np.vstack(buf_x).astype(np.float32); n_ok += len(X)
        for m in models:
            p = m["rf"].predict(X); h = heaps[m["key"]]
            for j in np.argpartition(-p, min(top, len(p) - 1))[:top]:
                it = (float(p[j]), buf_i[j], buf_s[j])
                if len(h) < top: heapq.heappush(h, it)
                elif it[0] > h[0][0]: heapq.heapreplace(h, it)
        buf_i.clear(); buf_s.clear(); buf_x.clear()
    def chunks(it, k=4000):
        c = []
        for x in it:
            c.append(x)
            if len(c) == k: yield c; c = []
        if c: yield c
    cur = 0
    with mp.Pool(max(1, (os.cpu_count() or 2)), initializer=_init) as pool:
        for ids, sm, X, n in pool.imap_unordered(featurize, chunks(stream(limit)), chunksize=1):
            n_in += n; buf_i += ids; buf_s += sm
            if len(X): buf_x.append(X); cur += len(X)
            if cur >= block: flush(); cur = 0; log(f"  ... {n_in:,} structures streamed, {n_ok:,} passed filters ({time.time() - t0:.0f}s)")
        flush()
    log(f"  scan done: {n_in:,} structures, {n_ok:,} scored per disease, {time.time() - t0:.0f}s")
    return {k: sorted(h, reverse=True) for k, h in heaps.items()}, n_in, n_ok

def rescore(model, hits, keep=200):
    """ADMET-AI organ-risk scoring of the best potency hits -> list of dict, ranked by potency-and-safety"""
    import admet, stabilizer
    sm = [s for _, _, s in hits[:keep]]
    if not sm: return []
    try: ad = admet.predict(sm)
    except Exception as e: print(f"  (organ-risk scoring unavailable: {e}; ranking by potency only)"); ad = None
    out = []
    for sc, cid, s in hits[:keep]:
        rk, tox = stabilizer.risks(ad[s]) if ad else ({}, None); known = None
        try:
            import potency; st = potency.standardize(s); known = model["data"].get(st)
        except Exception: pass
        out.append(dict(id=cid, smiles=s, pic50=sc, measured=known, tox=tox, flags=admet.flags(ad[s]) if ad else [], rank=(min(1, max(0, (sc - 5) / 4)) * .6 + (1 - (tox if tox is not None else .5)) * .4)))
    return sorted(out, key=lambda r: -r["rank"])

def run(keys=None, limit=None, top=15, show=True):
    import diseases
    keys = keys or list(diseases.D); models = []
    for k in keys:
        try:
            m = train(k)
            if m: models.append(m); print(f"  trained {m['key']:<14} on {m['n']:>5} measured compounds ({m['target']}); strict CV rmse/AUC {m['cv']}")
            else: print(f"  skipped {k}: no measured target data (use docking/onboard)")
        except Exception as e: print(f"  skipped {k}: {type(e).__name__}: {e}")
    if not models: return {}
    hits, n_in, n_ok = scan(models, limit); res = {}; os.makedirs(OUT, exist_ok=True)
    for m in models:
        rows = rescore(m, hits[m["key"]]); res[m["key"]] = rows
        with open(os.path.join(OUT, f"screen_{m['key']}.csv"), "w", newline="") as f:
            w = csv.writer(f); w.writerow(["chembl_id", "smiles", "pred_pIC50", "measured_pIC50", "tox_index", "flags"])
            for r in rows: w.writerow([r["id"], r["smiles"], round(r["pic50"], 2), r["measured"], r["tox"] and round(r["tox"], 3), "; ".join(r["flags"])])
        if show:
            print(f"\n=== {m['key']}: {m['name']} - top {top} of {n_ok:,} screened structures (potency model + organ-risk) ===")
            for r in rows[:top]:
                print(f"  {r['id']:<14} pIC50 {r['pic50']:.1f}" + (f" (measured {r['measured']:.1f})" if r["measured"] else " (NOVEL for this target)") + f" | tox {r['tox'] if r['tox'] is None else round(r['tox'], 2)} | {r['smiles'][:70]}")
    state.record("screened_n", n_in); state.record("screen_diseases", len(models)); return res

def cli(a):
    fl = {a[i]: a[i + 1] for i in range(len(a) - 1) if a[i].startswith("--")}; rest = [x for i, x in enumerate(a) if not x.startswith("--") and not (i and a[i - 1] in fl)]
    keys = fl["--keys"].split(",") if "--keys" in fl else (None if not rest or rest[0] == "all" else [rest[0]])
    return run(keys, int(fl["--limit"]) if "--limit" in fl else None, int(fl.get("--top", 15)))

if __name__ == "__main__": cli(sys.argv[1:])
