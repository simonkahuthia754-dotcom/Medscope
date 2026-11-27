"""potency.py - MEASURED potency data + QSAR for a drug target (replaces assumed potency).
Training data = every ChEMBL IC50 (nM, '=' relation) measured on the target, cleaned (largest fragment = salts removed,
canonical SMILES, duplicates -> median pIC50), optionally overlaid with YOUR lab measurements (lab.py).
Model = similarity-weighted k-NN on Morgan fingerprints (+ optional RandomForest when scikit-learn exists), leave-one-out CV,
applicability domain, uncertainty, ROC AUC (potent pIC50>=7 vs weak <=5). A molecule outside the domain gets 'no reliable
prediction', never a guess. Typical QSAR error is ~0.6-0.9 log units (4-8x): never 'exact'."""
import math, random, statistics
import sources

def standardize(smi):
    """largest fragment (strip salts/counter-ions) -> canonical SMILES; falls back to the input if RDKit is missing"""
    try:
        from rdkit import Chem
        from rdkit.Chem.MolStandardize import rdMolStandardize
        m = Chem.MolFromSmiles(smi)
        return Chem.MolToSmiles(rdMolStandardize.LargestFragmentChooser().choose(m)) if m else None
    except Exception: return smi

def fetch(chembl_target, max_n=4000, overlay=None, types=("IC50",)):
    """-> {smiles: pIC50}: ChEMBL measurements, cleaned + deduplicated, then overlaid with lab data (lab values win)"""
    rows, off = {}, 0
    while off < max_n:
        a = sources.J(f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={chembl_target}&{'standard_type=IC50' if types == ('IC50',) else 'standard_type__in=' + ','.join(types)}&standard_units=nM"
                      f"&standard_relation=%3D&limit=1000&offset={off}")["activities"]
        for x in a:
            try:
                v = float(x["standard_value"]); s = x.get("canonical_smiles")
                if s and 0 < v < 1e9: rows.setdefault(s, []).append(9 - math.log10(v))
            except Exception: pass
        if len(a) < 1000: break
        off += 1000
    out = {}
    for s, v in rows.items():
        c = standardize(s)
        if c: out.setdefault(c, []).extend(v)
    out = {s: statistics.median(v) for s, v in out.items()}
    for s, y in (overlay or {}).items():
        c = standardize(s)
        if c: out[c] = y
    return out

def _rd_fp(smi):
    import chemistry
    m = chemistry.mol(smi); return chemistry.fp(m) if m else None

def _rd_sim(a, b):
    from rdkit import DataStructs
    return DataStructs.TanimotoSimilarity(a, b)

def roc_auc(pos, neg):
    """Mann-Whitney ROC AUC of predicted scores for truly-potent vs truly-weak compounds"""
    if not pos or not neg: return None
    return sum((p > n) + .5 * (p == n) for p in pos for n in neg) / (len(pos) * len(neg))

class Model:
    def __init__(self, data, fp=_rd_fp, sim=_rd_sim, k=5, domain=.40, seed=3):
        self.fp, self.sim, self.k, self.domain = fp, sim, k, domain
        self.items = [(s, y, f) for s, y in data.items() for f in [fp(s)] if f is not None]; self.n = len(self.items)
        self.loo = []; self.loo_c = []; self.rmse, self.r2 = self._cv(seed); self.auc = self._auc(); self.rmse_c = self.auc_c = None; self.rf = None; self.which = 'kNN'; self._cv_cluster(seed); self._rf(seed)

    def _nb(self, f, skip=None, ex=.999, cut=None):
        s = sorted(((self.sim(f, g), y, i) for i, (_, y, g) in enumerate(self.items) if i != skip), reverse=True)
        if ex: s = [t for t in s if t[0] < ex]          # drop the molecule itself (and exact duplicates): no leakage
        if cut: s = [t for t in s if t[0] < cut]        # cluster-aware: also drop the same chemical series
        return s[:self.k]

    def _est(self, sims):
        w = [max(s, 1e-6) ** 3 for s, _, _ in sims]; m = sum(wi * y for wi, (_, y, _) in zip(w, sims)) / sum(w)
        sd = math.sqrt(sum(wi * (y - m) ** 2 for wi, (_, y, _) in zip(w, sims)) / sum(w)); return m, sd, sims[0][0]

    def _cv(self, seed, cap=400):
        if self.n < 20: return None, None
        idx = random.Random(seed).sample(range(self.n), min(cap, self.n)); err, ys = [], []
        for i in idx:
            m, _, _ = self._est(self._nb(self.items[i][2], skip=i)); err.append(m - self.items[i][1]); ys.append(self.items[i][1]); self.loo.append((self.items[i][1], m))
        var = statistics.pvariance(ys) or 1e-9; mse = sum(e * e for e in err) / len(err); return math.sqrt(mse), 1 - mse / var

    def _cv_cluster(self, seed, cut=.55, cap=400):
        """STRICT validation: predict each compound with its whole chemical series (Tanimoto >= cut) removed - tests generalisation to NEW chemotypes, not interpolation inside a series"""
        if self.n < 20: return
        idx = random.Random(seed + 1).sample(range(self.n), min(cap, self.n)); err = []
        for i in idx:
            nb = self._nb(self.items[i][2], skip=i, cut=cut)
            if len(nb) < 3: continue
            m, _, _ = self._est(nb); err.append(m - self.items[i][1]); self.loo_c.append((self.items[i][1], m))
        if len(err) >= 20: self.rmse_c = math.sqrt(sum(e * e for e in err) / len(err))
        pos = [p for y, p in self.loo_c if y >= 7]; neg = [p for y, p in self.loo_c if y <= 5]; self.auc_c = roc_auc(pos, neg) if len(pos) >= 3 and len(neg) >= 3 else None

    def _vec(self, smi):
        from rdkit import Chem
        from rdkit.Chem import rdFingerprintGenerator
        m = Chem.MolFromSmiles(smi); return None if m is None else rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024).GetCountFingerprintAsNumPy(m).astype(float)

    def _rf(self, seed, cap=4000):
        """random-forest on Morgan counts, validated with SCAFFOLD-GROUPED 5-fold CV (whole scaffolds held out). Used instead of kNN only if it generalises better to new chemotypes."""
        try:
            import numpy as np
            from sklearn.ensemble import RandomForestRegressor
            from sklearn.model_selection import GroupKFold, cross_val_predict
            from rdkit import Chem
            from rdkit.Chem.Scaffolds import MurckoScaffold
        except Exception: return
        if self.n < 60: return
        rows = [(s, y) for s, y, _ in self.items][:cap]; X, y, g = [], [], []
        for s, v in rows:
            x = self._vec(s)
            if x is None: continue
            try: sc = MurckoScaffold.MurckoScaffoldSmiles(smiles=s)
            except Exception: sc = s
            X.append(x); y.append(v); g.append(hash(sc) % 100003)
        X, y, g = np.array(X), np.array(y), np.array(g)
        if len(set(g)) < 6: return
        rf = RandomForestRegressor(n_estimators=150, min_samples_leaf=1, n_jobs=-1, random_state=seed)
        pr = cross_val_predict(rf, X, y, groups=g, cv=GroupKFold(5)); rm = float(np.sqrt(np.mean((pr - y) ** 2)))
        pos = [a for a, b in zip(pr, y) if b >= 7]; neg = [a for a, b in zip(pr, y) if b <= 5]; auc = roc_auc(pos, neg) if len(pos) >= 3 and len(neg) >= 3 else None
        self.rmse_rf, self.auc_rf = rm, auc
        if self.rmse_c is None or rm < self.rmse_c:
            self.rf = rf.fit(X, y); self.which = "random-forest"; self.rmse_c, self.auc_c = rm, auc

    def _auc(self):
        pos = [p for y, p in self.loo if y >= 7]; neg = [p for y, p in self.loo if y <= 5]
        return roc_auc(pos, neg) if len(pos) >= 3 and len(neg) >= 3 else None

    def predict(self, smi, exclude_self=False):
        """-> dict(pic50, sd, ic50_nM, ic50_lo, ic50_hi, maxsim, in_domain, confidence). exclude_self=True for fair validation of known drugs"""
        f = self.fp(smi)
        if f is None or not self.items: return None
        sims = self._nb(f, ex=.999 if exclude_self else None)
        if not sims: return None
        m, sd, ms = self._est(sims); base = self.rmse or .8
        if self.rf is not None:
            x = self._vec(smi)
            if x is not None: m = float(self.rf.predict(x.reshape(1, -1))[0]); base = self.rmse_c or base   # strict-CV model chosen: its error bar applies
        sd_tot = max(base, sd) * (1 + 2 * max(0, .6 - ms)); ind = ms >= self.domain
        conf = "HIGH" if ind and ms >= .7 and sd_tot < 1.0 else "MEDIUM" if ind and ms >= .5 else "LOW"
        return {"pic50": m, "sd": sd_tot, "ic50_nM": 10 ** (9 - m), "ic50_lo": 10 ** (9 - m - sd_tot), "ic50_hi": 10 ** (9 - m + sd_tot),
                "maxsim": ms, "in_domain": ind, "confidence": conf if ind else "OUT-OF-DOMAIN"}

    def summary(self):
        a = f" | ROC AUC {self.auc:.2f} (potent>=7 vs weak<=5; target >0.8)" if self.auc else ""
        r = f"RMSE {self.rmse:.2f} log units (~{10 ** self.rmse:.1f}-fold), R2 {self.r2:.2f}{a}" if self.rmse else "too little data for cross-validation"
        c = f" | NEW-CHEMOTYPE CV [{self.which}] (series/scaffolds removed): RMSE {self.rmse_c:.2f}" + (f", ROC AUC {self.auc_c:.2f}" if self.auc_c else "") if self.rmse_c else ""
        return f"{self.n} measured compounds | leave-one-out CV: {r}{c}"
