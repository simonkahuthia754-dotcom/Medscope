#!/usr/bin/env python3
"""stabilizer.py - make a candidate SAFER before anyone touches a lab: lower predicted organ toxicity while keeping potency.
  python stabilizer.py "<SMILES|drug name>" [--ic50 nM] [--target lung-cancer] [--n 40] [--ancestry african]      (MedScope: stabilize ...)
Three levers, each reported with before/after numbers:
  1. STRUCTURE  - RDKit analogs (small edits) re-scored by ADMET-AI (hERG heart, DILI liver, Ames genotoxicity, ClinTox, carcinogenicity, skin), structural alerts,
                  and potency (measured-data QSAR when --target is given); keeps only analogs that do not lose >0.5 log potency.
  2. REGIMEN    - same free trough (efficacy) delivered with lower peaks: split dosing/extended release, PBPK-based; dose cut; organ peak exposure.
  3. POPULATION - mass_population: % of 10,000 virtual people with high exposure (liver/kidney disease, poor metabolisers, ...) before vs after.
SAFETY SCORE (0-100, transparent heuristic, NOT a probability of being safe): 40 x (1 - organ toxicity index) + 30 x (1 - %high-exposure population)
  + 15 x (therapeutic-margin proxy) + 15 x (no structural alerts). Predicted organ toxicity is a model output: animal studies and trials decide real safety."""
import math, sys
import pbpk, state

ENDP = {"heart": ("hERG", "hERG/QT"), "liver": ("DILI", "drug-induced liver injury"), "genome": ("AMES", "Ames mutagenicity"), "systemic": ("ClinTox", "clinical-trial toxicity"),
        "carcinogenic": ("Carcinogens", "carcinogenicity"), "skin": ("Skin_Reaction", "skin sensitisation")}
ORGAN_EXP = {"heart": "heart", "liver": "liver", "skin": "skin", "genome": None, "systemic": None, "carcinogenic": None}

def risks(ad, exp=None):
    """ADMET probabilities (+ organ exposure ratios from PBPK) -> per-organ risk and an overall toxicity index = 1 - prod(1 - risk)"""
    import admet
    out = {}
    for k, (key, label) in ENDP.items():
        p = admet.get(ad, key)
        if p is None: continue
        o = ORGAN_EXP[k]; w = 1.0 if not (exp and o and o in exp) else .5 + .5 * min(1.0, exp[o][1])   # a liability only matters where the drug actually goes
        out[k] = dict(p=p, risk=p * w, label=label)
    idx = 1 - math.prod(1 - v["risk"] for v in out.values()) if out else None; return out, idx

def regimen(p, ovr, ic50_nM, mw, fup=None, cover=1.0):
    """equal free trough (= efficacy) at tau = 24/12/8/6 h: dose, daily dose, peak-to-trough, peak free conc vs IC50"""
    rows = []
    for tau in (24.0, 12.0, 8.0, 6.0):
        d, r = pbpk.dose_for_target(p, "human", ic50_nM, mw, ovr=ovr, target_cover=cover, tau=tau)
        if not d: continue
        css = pbpk.steady_state(r, tau); sc = d / 100; pk, tr = max(css) * sc, min(css) * sc; fup = r["params"]["fup"]
        rows.append(dict(tau=tau, dose=d, daily=d * 24 / tau, peak_mgL=pk, trough_mgL=tr, ptr=pk / max(tr, 1e-12), peak_free_x_ic50=pk * fup / mw * 1e6 / ic50_nM))
    return rows

def analogs(smi, n=40, sim_min=.6, seed=7):
    import chemistry, potency
    from rdkit import DataStructs
    m = chemistry.mol(smi); fp0 = chemistry.fp(m); import random; out = {}
    for s in chemistry.edits(m, random.Random(seed), per=max(3, n // 8)):
        mm = chemistry.mol(s)
        if mm is None: continue
        sim = DataStructs.TanimotoSimilarity(fp0, chemistry.fp(mm))
        if sim >= sim_min and s != potency.standardize(smi): out[s] = sim
    return sorted(out.items(), key=lambda kv: -kv[1])[:n]

def alerts_of(smi):
    try:
        import chemistry; return list(chemistry.alerts(chemistry.mol(smi)))
    except Exception: return []

def profile(smi, ad, ic50_nM=None, model=None):
    import admet, ingredients
    from rdkit import Chem
    from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
    m = Chem.MolFromSmiles(smi); p = dict(logp=Crippen.MolLogP(m), mw=Descriptors.MolWt(m), tpsa=rdMolDescriptors.CalcTPSA(m)); ovr = admet.overrides(ad)
    po, iv, F = pbpk.bioavailability(p, "human", 100.0, 48.0, ovr=ovr); exp = ingredients.exposure(po); rk, idx = risks(ad, exp)
    pot = None
    if model is not None:
        try: pot = model.predict(smi)["pic50"]
        except Exception: pass
    return dict(smiles=smi, p=p, ovr=ovr, F=F, exp=exp, risks=rk, tox=idx, alerts=alerts_of(smi), pic50=pot, ic50_nM=ic50_nM)

def pop_high(prof, ic50, ancestry="african", n=3000):
    import mass_population as mp
    r = mp.run(prof["p"], prof["ovr"], {}, prof["p"]["mw"], 100.0, ic50, 24.0, n, ancestry, show=False); return r["overall"]["high"], r["overall"]["low"]

def safety_score(prof, high_pct, margin=None):
    t = 1 - (prof["tox"] if prof["tox"] is not None else .5); parts = dict(organ=40 * t, population=30 * (1 - high_pct / 100), margin=15 * (margin if margin is not None else .5), alerts=15 * (0 if prof["alerts"] else 1))
    return sum(parts.values()), parts

def stabilize(smi, ic50_nM=None, model=None, n=40, ancestry="african", show=True, ad_fn=None, popn=3000):
    import admet, potency
    smi = potency.standardize(smi); ad_fn = ad_fn or (lambda xs: admet.predict(xs))
    cand = analogs(smi, n); allsmi = [smi] + [s for s, _ in cand]; ad = ad_fn(allsmi); base = profile(smi, ad[smi], ic50_nM, model); sims = dict(cand)
    rows = []
    for s in allsmi[1:]:
        try: pr = profile(s, ad[s], ic50_nM, model)
        except Exception: continue
        if model is not None and base["pic50"] is not None and pr["pic50"] is not None and pr["pic50"] < base["pic50"] - .5: continue   # keep potency
        pr["sim"] = sims[s]; rows.append(pr)
    rows.sort(key=lambda r: (r["tox"] if r["tox"] is not None else 1, -len(r["alerts"]) * 0, -(r["pic50"] or 0)))
    best = rows[0] if rows and rows[0]["tox"] is not None and base["tox"] is not None and rows[0]["tox"] < base["tox"] else None
    ic = ic50_nM or 10.0; res = {"base": base, "best": best, "analogs_tested": len(allsmi) - 1, "kept": len(rows)}
    hb, lb = pop_high(base, ic, ancestry, popn); sb, pb = safety_score(base, hb); res.update(pop_before=(hb, lb), score_before=sb, parts_before=pb)
    tgt = best or base; ha, la = pop_high(tgt, ic, ancestry, popn) if best else (hb, lb); reg = regimen(tgt["p"], tgt["ovr"], ic, tgt["p"]["mw"])
    marg = None
    if reg: marg = max(0.0, min(1.0, 1 - min(r["peak_free_x_ic50"] for r in reg) / 50))   # lower peak free conc relative to IC50 = wider margin (proxy)
    sa, pa = safety_score(tgt, ha, marg); res.update(pop_after=(ha, la), score_after=sa, parts_after=pa, regimen=reg)
    red = (base["tox"] - best["tox"]) / base["tox"] if best and base["tox"] else 0.0; res["tox_reduction"] = red
    state.record("stab_ok", True); state.record("stab_reduction", round(red, 2)); state.record("stab_score", round(sa))
    if show: report(res, ic50_nM)
    return res

def report(r, ic50):
    b, x = r["base"], r["best"]; print(f"\nSTABILIZER: {r['analogs_tested']} analogs generated, {r['kept']} kept after the potency filter")
    f = lambda v: "n/a" if v is None else f"{v:.2f}"
    print(f"  PARENT  {b['smiles']}\n    organ toxicity index {f(b['tox'])} | " + ", ".join(f"{k} {v['risk']:.2f}" for k, v in b["risks"].items()) + f" | alerts: {', '.join(b['alerts']) or 'none'}" + (f" | pIC50 {b['pic50']:.1f}" if b["pic50"] else ""))
    if x: print(f"  BEST    {x['smiles']}  (similarity {x['sim']:.2f})\n    organ toxicity index {f(x['tox'])} | " + ", ".join(f"{k} {v['risk']:.2f}" for k, v in x["risks"].items()) + f" | alerts: {', '.join(x['alerts']) or 'none'}" + (f" | pIC50 {x['pic50']:.1f}" if x["pic50"] else ""))
    else: print("  no analog lowered the predicted toxicity without losing potency - the parent stays; use the regimen/population levers below")
    print(f"  predicted toxicity reduction: {100 * r['tox_reduction']:.0f}%")
    print("  REGIMEN (same free trough = same efficacy):  tau  dose   daily  peak/trough  peak free conc / IC50")
    for g in r["regimen"][:4]: print(f"                                              {g['tau']:>4g}h {g['dose']:>5.0f} mg {g['daily']:>5.0f} mg {g['ptr']:>6.1f}      {g['peak_free_x_ic50']:.1f}x")
    hb, lb = r["pop_before"]; ha, la = r["pop_after"]; print(f"  POPULATION high-exposure (AUC>3x or Cmax>2.5x reference): before {hb:.1f}% -> after {ha:.1f}%" + (f" | below effective trough: {lb:.1f}% -> {la:.1f}%" if lb is not None else ""))
    print(f"  SAFETY SCORE: before {r['score_before']:.0f}/100 -> after {r['score_after']:.0f}/100   (organ {r['parts_after']['organ']:.0f}/40, population {r['parts_after']['population']:.0f}/30, margin {r['parts_after']['margin']:.0f}/15, alerts {r['parts_after']['alerts']:.0f}/15)")
    print("  This is a heuristic on model predictions. Real organ safety needs in-vitro tox panels, animal studies (with ethics approval) and trials.\n")

def cli(a):
    fl = {a[i]: a[i + 1] for i in range(len(a) - 1) if a[i].startswith("--")}; rest = [x for i, x in enumerate(a) if not x.startswith("--") and not (i and a[i - 1] in fl)]
    if not rest: print(__doc__); return
    arg = " ".join(rest); smi = arg
    if not any(c in arg for c in "()=#[]") or arg.isalpha():
        import sources; smi = sources.pubchem(arg)["smiles"]
    model = None
    if "--target" in fl:
        import diseases, potency; D = diseases.get(fl["--target"]); model = potency.Model(potency.fetch(D["chembl"]))
    return stabilize(smi, float(fl["--ic50"]) if "--ic50" in fl else None, model, int(fl.get("--n", 40)), fl.get("--ancestry", "african"))

if __name__ == "__main__": cli(sys.argv[1:])
