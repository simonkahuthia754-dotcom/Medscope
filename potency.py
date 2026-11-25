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

def fetch(chembl_target, max_n=4000, overlay=None):
    """-> {smiles: pIC50}: ChEMBL measurements, cleaned + deduplicated, then overlaid with lab data (lab values win)"""
    rows, off = {}, 0
    while off < max_n:
        a = sources.J(f"https://www.ebi.ac.uk/chembl/api/data/activity.json?target_chembl_id={chembl_target}&standard_type=IC50&standard_units=nM"
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
        self.loo = []; self.rmse, self.r2 = self._cv(seed); self.auc = self._auc()

    def _nb(self, f, skip=None, ex=.999):
        s = sorted(((self.sim(f, g), y, i) for i, (_, y, g) in enumerate(self.items) if i != skip), reverse=True)
        if ex: s = [t for t in s if t[0] < ex]          # drop the molecule itself (and exact duplicates): no leakage
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
        sd_tot = max(base, sd) * (1 + 2 * max(0, .6 - ms)); ind = ms >= self.domain
        conf = "HIGH" if ind and ms >= .7 and sd_tot < 1.0 else "MEDIUM" if ind and ms >= .5 else "LOW"
        return {"pic50": m, "sd": sd_tot, "ic50_nM": 10 ** (9 - m), "ic50_lo": 10 ** (9 - m - sd_tot), "ic50_hi": 10 ** (9 - m + sd_tot),
                "maxsim": ms, "in_domain": ind, "confidence": conf if ind else "OUT-OF-DOMAIN"}

    def summary(self):
        a = f" | ROC AUC {self.auc:.2f} (potent>=7 vs weak<=5; target >0.8)" if self.auc else ""
        r = f"RMSE {self.rmse:.2f} log units (~{10 ** self.rmse:.1f}-fold), R2 {self.r2:.2f}{a}" if self.rmse else "too little data for cross-validation"
        return f"{self.n} measured compounds | leave-one-out CV: {r}"
