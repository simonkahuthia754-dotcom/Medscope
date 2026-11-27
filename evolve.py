#!/usr/bin/env python3
"""evolve.py - GENERATIVE design: invent NEW molecules (not in any database) for a disease target by evolutionary search.
  python evolve.py <disease key> [--gens 8] [--pop 200] [--seeds 20] [--top 10]                          (MedScope: evolve ...)
Start from the most potent MEASURED compounds, mutate them with medicinal-chemistry edits (chemistry.edits), score every child by the target's potency model
(strict-CV random forest) minus penalties for structural alerts / property violations, keep the fittest, repeat. Children found in ChEMBL are dropped, so what
remains is NOVEL structures. The score is a model prediction (error bar = new-chemotype CV): novel compounds need synthesis and a real assay."""
import random, sys
import state

def fitness(m, rf, vecf, smi, alerts_fn, props_fn):
    x = vecf(smi)
    if x is None: return None
    import numpy as np
    p = float(rf.predict(x.reshape(1, -1).astype(np.float32))[0])
    pen = .3 * len(alerts_fn(smi)); mw, lp, tp = props_fn(smi)
    pen += (.5 if mw > 600 else 0) + (.5 if lp > 5.5 else 0) + (.5 if tp > 140 else 0)
    return p - pen, p

def evolve(key, gens=8, pop=200, seeds=20, top=10, seed=11, show=True):
    import screen, chemistry, potency
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    M = screen.train(key)
    if not M: raise RuntimeError(f"no measured data for {key}: evolve needs a trained potency model")
    rng = random.Random(seed); known = set(M["data"]); rf = M["rf"]
    def props(s):
        m = Chem.MolFromSmiles(s); return Descriptors.MolWt(m), Crippen.MolLogP(m), rdMolDescriptors.CalcTPSA(m)
    def alerts(s):
        try: return list(chemistry.alerts(chemistry.mol(s)))
        except Exception: return []
    best_known = sorted(M["data"].items(), key=lambda kv: -kv[1])[:seeds]; cur = {}
    for s, y in best_known:
        f = fitness(M, rf, screen.vec, s, alerts, props)
        if f: cur[s] = f
    gen_hist = []; novel = {}
    for g in range(gens):
        kids = {}
        for s in list(cur)[:max(seeds, pop // 4)]:
            m = chemistry.mol(s)
            if m is None: continue
            for c in chemistry.edits(m, rng, per=3):
                try: c = potency.standardize(c)
                except Exception: continue
                if c in cur or c in kids or len(c) > 110: continue
                f = fitness(M, rf, screen.vec, c, alerts, props)
                if f: kids[c] = f
        cur.update(kids); cur = dict(sorted(cur.items(), key=lambda kv: -kv[1][0])[:pop])
        for s, f in cur.items():
            if s not in known: novel[s] = f
        gen_hist.append((g + 1, len(kids), max(f[1] for f in cur.values()), sum(f[1] for f in cur.values()) / len(cur)))
        if show: print(f"  generation {g + 1}: {len(kids)} new children, best predicted pIC50 {gen_hist[-1][2]:.2f}, population mean {gen_hist[-1][3]:.2f}")
    ranked = sorted(novel.items(), key=lambda kv: -kv[1][0])[:top]
    out = [dict(smiles=s, fit=f[0], pic50=f[1]) for s, f in ranked]
    state.record("evolve_novel", len(novel)); state.record("evolve_ok", True)
    if show:
        print(f"\n=== {key}: {len(novel)} NOVEL molecules generated (absent from ChEMBL training data for this target); top {len(out)} ===")
        print(f"   error bar: new-chemotype CV {M['cv']} (rmse, AUC) - predictions for unfamiliar structures are LESS reliable than for known series")
        for r in out: print(f"  pIC50~{r['pic50']:.1f}  {r['smiles']}")
    return dict(novel=len(novel), top=out, hist=gen_hist, model=M)

def cli(a):
    fl = {a[i]: a[i + 1] for i in range(len(a) - 1) if a[i].startswith("--")}; rest = [x for i, x in enumerate(a) if not x.startswith("--") and not (i and a[i - 1] in fl)]
    if not rest: print(__doc__); return
    return evolve(rest[0], int(fl.get("--gens", 8)), int(fl.get("--pop", 200)), int(fl.get("--seeds", 20)), int(fl.get("--top", 10)))

if __name__ == "__main__": cli(sys.argv[1:])
