"""chemistry.py - RDKit cheminformatics: REAL structures, descriptors, filters, analog design, drawings.
Needs RDKit (pip install rdkit  /  conda install -c conda-forge rdkit). Not available in QPython/Android -> run on a PC/server/Colab.
Analogs come from (1) medicinal-chemistry edits via reaction SMARTS and (2) BRICS fragment recombination of real leads, then are
filtered by validity, QED, Lipinski/Veber, PAINS, similarity to a real lead. Potency is NOT predicted here (that needs docking/assays)."""
import itertools, os, random

try:
    from rdkit import Chem, DataStructs, RDLogger
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors, QED, AllChem, BRICS, rdFingerprintGenerator
    from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
    from rdkit.Chem.Draw import rdMolDraw2D
    RDLogger.DisableLog("rdApp.*"); HAVE = True
except Exception:
    HAVE = False

# well-known drug SMILES (validated at runtime; PubChem value is preferred when online)
LEAD_SMILES = {
    "zidovudine": "CC1=CN(C(=O)NC1=O)C2CC(C(O2)CO)N=[N+]=[N-]",
    "tenofovir": "CC(CN1C=NC2=C(N=CN=C21)N)OCP(=O)(O)O",
    "imatinib": "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1",
    "erlotinib": "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC",
    "chloroquine": "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12",
    "ciprofloxacin": "OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O",
    "linezolid": "CC(=O)NC[C@H]1CN(c2ccc(N3CCOCC3)c(F)c2)C(=O)O1",
}
EDITS = ["[c;H1:1]>>[c:1]F", "[c;H1:1]>>[c:1]Cl", "[c;H1:1]>>[c:1]C", "[c;H1:1]>>[c:1]OC", "[F:1]>>[Cl:1]", "[Cl:1]>>[F:1]",
         "[CH3:1][c:2]>>[CH2:1](F)[c:2]", "[OH:1][CH2:2]>>[F:1][CH2:2]", "[c;H1:1]>>[c:1]C#N", "[NH2:1][c:2]>>[OH:1][c:2]",
         # medicinal-chemistry transformations that typically cut hERG / liver / genotox liabilities: add polarity, remove lipophilic or reactive groups, scan ring nitrogens
         "[c;H1:1]>>[n:1]", "[c;H1:1]>>[c:1]O", "[c;H1:1]>>[c:1]C(N)=O", "[c;H1:1]>>[c:1]C(=O)O", "[c;H1:1]>>[c:1]C(F)(F)F", "[c;H1:1]>>[c:1]S(N)(=O)=O", "[c;H1:1]>>[c:1]CO",
         "[CH3:1][c:2]>>[CH2:1](O)[c:2]", "[CH3:1][c:2]>>[H][c:2]", "[c:1][Cl]>>[c:1]C(F)(F)F", "[c:1][Cl]>>[c:1]C#N", "[c:1][Br]>>[c:1]F", "[c:1][Br]>>[c:1]Cl", "[c:1][O:2][CH3]>>[c:1][O:2]",
         "[O;X2:1][CH3:2]>>[O:1]", "[N;X3;!$(N-C=O);!$(N-c):1][CH3]>>[N:1]", "[N;X3;!$(N-C=O);!$(N-c):1][CH3]>>[N:1]C(C)=O", "[CH2;!R:1]>>[C;!R:1](F)F",
         "[C:1]#[CH:2]>>[C:1]=[CH2:2]", "[CH2;R:1]>>[O;R:1]", "[CH2;R:1]>>[NH;R:1]", "[CH3;!$(C-[!#6;!#1]):1]>>[CH2:1]O", "[N;X3;H1;!$(N-C=O):1]>>[N:1]C(C)=O"]

def need():
    if not HAVE: raise RuntimeError("RDKit not installed. On a PC/Colab: pip install rdkit  (not available on QPython/Android)")

def mol(smi):
    need(); m = Chem.MolFromSmiles(smi) if smi else None; return m

_FC = None
def _catalog(kind):
    global _FC
    if _FC is None:
        _FC = {}
        for k, c in (("pains", FilterCatalogParams.FilterCatalogs.PAINS), ("brenk", FilterCatalogParams.FilterCatalogs.BRENK)):
            p = FilterCatalogParams(); p.AddCatalog(c); _FC[k] = FilterCatalog(p)
    return _FC[kind]

def alerts(m):
    """structural alerts: PAINS (assay interference) and Brenk (reactive / toxicophore / unstable groups)"""
    return {"pains": [e.GetDescription() for e in _catalog("pains").GetMatches(m)], "brenk": [e.GetDescription() for e in _catalog("brenk").GetMatches(m)]}

BASIC_AMINE = None
def describe(m):
    """real descriptors computed from the structure"""
    global BASIC_AMINE
    if BASIC_AMINE is None: BASIC_AMINE = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(NC=[O,S,N]);!$(N-a);!$(N-S(=O)=O);!$(N#*);!$(N=*)]")
    d = {"mw": Descriptors.MolWt(m), "logp": Crippen.MolLogP(m), "tpsa": rdMolDescriptors.CalcTPSA(m), "hbd": rdMolDescriptors.CalcNumHBD(m),
         "hba": rdMolDescriptors.CalcNumHBA(m), "rotb": rdMolDescriptors.CalcNumRotatableBonds(m), "qed": QED.qed(m),
         "basic_amine": m.HasSubstructMatch(BASIC_AMINE), "rings": rdMolDescriptors.CalcNumRings(m), "heavy": m.GetNumHeavyAtoms()}
    d["lipinski_viol"] = sum([d["mw"] > 500, d["logp"] > 5, d["hbd"] > 5, d["hba"] > 10]); d["veber_ok"] = d["tpsa"] <= 140 and d["rotb"] <= 10
    d["herg_flag"] = bool(d["basic_amine"] and d["logp"] > 3.7)   # pharmacophore flag only, NOT a hERG prediction
    return d

_GEN = None
def fp(m):
    global _GEN
    if _GEN is None: _GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    return _GEN.GetFingerprint(m)

def edits(m, rng, per=6):
    out = set()
    for sm in EDITS:
        try: rx = AllChem.ReactionFromSmarts(sm)
        except Exception: continue
        prods = set()
        for ps in rx.RunReactants((m,)):
            try: p = ps[0]; Chem.SanitizeMol(p); prods.add(Chem.MolToSmiles(p))
            except Exception: pass
        out |= set(rng.sample(sorted(prods), min(per, len(prods))))
    return out

def brics(leads, rng, tries=1500, nfrag=30):
    fr = set()
    for m in leads: fr |= BRICS.BRICSDecompose(m)
    fm = [x for x in (Chem.MolFromSmiles(f) for f in sorted(fr)) if x]; rng.shuffle(fm); out = set()
    for p in itertools.islice(BRICS.BRICSBuild(fm[:nfrag], maxDepth=2, scrambleReagents=True), tries):
        try: Chem.SanitizeMol(p); out.add(Chem.MolToSmiles(p))
        except Exception: pass
    return out

def candidates(leads, n=60, seed=7, mw=(200, 600), min_qed=.35, sim=(.25, .95)):
    """leads: {name: {"smiles":..., "ic50_nM":..., "data":...}} -> ranked real-structure candidates (leads included)"""
    need(); rng = random.Random(seed); L = {}
    for nm, v in leads.items():
        m = mol(v.get("smiles") or LEAD_SMILES.get(nm))
        if m: L[nm] = (m, v)
    if not L: return []
    lfp = {nm: fp(m) for nm, (m, _) in L.items()}; out, seen = [], set()
    def add(smi, name, parent, known):
        m = mol(smi); can = Chem.MolToSmiles(m) if m else None
        if not m or can in seen: return
        d = describe(m); a = alerts(m)
        sims = {nm: DataStructs.TanimotoSimilarity(fp(m), f) for nm, f in lfp.items()}; best = max(sims, key=sims.get); s = sims[best]
        if not known:
            if not (mw[0] <= d["mw"] <= mw[1]) or d["qed"] < min_qed or d["lipinski_viol"] > 1 or not d["veber_ok"] or a["pains"] or not (sim[0] <= s <= sim[1]): return
        seen.add(can); par = parent or best; pv = L[par][1]
        pen = min(1.0, .34 * len(a["brenk"]) + 1.0 * len(a["pains"]))
        score = round(min(98, 100 * (.45 * d["qed"] + .30 * min(1, s / .6) + .25 * (1 - pen))), 1)
        out.append({"name": name, "smiles": can, "known": known, "parent": par, "sim": round(s, 2), "alerts": a, "score": score, "ic50_nM": pv.get("ic50_nM") if known else None,
                    "data": pv.get("data", "RDKit"), **d})
    for nm, (m, _) in L.items(): add(Chem.MolToSmiles(m), nm, nm, True)
    pool = []
    for nm, (m, _) in L.items(): pool += [(s, nm) for s in edits(m, rng)]
    for s in brics([m for m, _ in L.values()], rng): pool.append((s, None))
    rng.shuffle(pool); i = 0
    for s, par in pool:
        if sum(1 for c in out if not c["known"]) >= n: break
        i += 1; add(s, f"cand-{i:03d}", par, False)
    return sorted(out, key=lambda c: -c["score"])

def svg(smi, path, size=(420, 320)):
    m = mol(smi); d = rdMolDraw2D.MolDraw2DSVG(*size); rdMolDraw2D.PrepareAndDrawMolecule(d, m); d.FinishDrawing()
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True); open(path, "w").write(d.GetDrawingText()); return path
