#!/usr/bin/env python3
"""chemicals.py - candidate screening pipeline:  structures -> filters -> rat PBPK -> human PBPK -> dose & safety flags.
Structures: REAL (RDKit, chemistry.py) when RDKit is installed, else a property-only fallback.
Simulators: pbpk.py whole-body models (rat 250 g, human 70 kg). Scores are heuristic; potency of analogs is ASSUMED equal to the
parent's measured ChEMBL IC50 (unverified). Passing = worth real lab work. NOT medical advice; not a cure finder."""
import math, os, random, sys
import pbpk

NOTE = "SCREENING MODEL - heuristic scores, NOT cure probabilities. Unvalidated; not medical advice."
# fallback property table (approx): name: (MW, logP, HBD, HBA, class)
LIB = {"zidovudine": (267, 0.1, 2, 6, "antiviral"), "dolutegravir": (419, 1.0, 2, 7, "antiviral"), "tenofovir": (287, -1.6, 2, 8, "antiviral"),
       "imatinib": (494, 3.5, 2, 7, "kinase"), "erlotinib": (393, 3.3, 1, 7, "kinase"), "artemether": (298, 3.5, 0, 5, "antimalarial"),
       "chloroquine": (320, 4.6, 1, 3, "antimalarial"), "ciprofloxacin": (331, 0.3, 2, 6, "antibiotic"), "linezolid": (337, 0.6, 1, 6, "antibiotic")}
CLS = {"hiv": "antiviral", "virus": "antiviral", "hepatitis": "antiviral", "covid": "antiviral", "cancer": "kinase", "leukemia": "kinase",
       "tumor": "kinase", "malaria": "antimalarial", "tuberculosis": "antibiotic", "bacter": "antibiotic"}

def cls_for(d):
    d = d.lower(); return next((v for k, v in CLS.items() if k in d), None)

def traits(c):
    """property-only fallback score (used when RDKit is unavailable)"""
    fit = math.exp(-((c["logp"] - 2.5) / 2.2) ** 2) * math.exp(-((c["mw"] - 350) / 160) ** 2); v = sum([c["mw"] > 500, c["logp"] > 5, c["hbd"] > 5, c["hba"] > 10])
    c["lipinski_viol"] = v; c["alerts"] = {"pains": [], "brenk": []}; c["tpsa"] = c.get("tpsa") or min(200, 25 * c["hbd"] + 12 * c["hba"])
    c["score"] = round(min(98, 100 * (.5 * (.1 + .8 * fit) + .3 * (1 - .25 * v) + .2 * .85)), 1); return c

def fallback(cls, real, n=30, seed=7):
    r = random.Random(seed); out = []
    for nm, (mw, lp, hd, ha, k) in LIB.items():
        if cls and k != cls: continue
        x = real.get(nm) or {}; src = "PubChem/ChEMBL" if x else "built-in approx"
        if x: mw, lp, hd, ha = x["mw"], x["logp"], x["hbd"], x["hba"]
        base = dict(mw=mw, logp=lp, hbd=hd, hba=ha, tpsa=x.get("tpsa"), ic50_nM=x.get("ic50_nM"), data=src, smiles=x.get("smiles"), parent=nm)
        out.append(traits({**base, "name": nm, "known": True}))
        for i in range(n):
            out.append(traits({**base, "name": f"{nm}-v{i+1}", "known": False, "tpsa": None, "mw": round(mw * r.uniform(.8, 1.25)), "logp": round(lp + r.uniform(-1, 1), 1),
                               "hbd": max(0, hd + r.choice([-1, 0, 1])), "hba": max(1, ha + r.choice([-2, -1, 0, 1, 2]))}))
    return sorted(out, key=lambda c: -c["score"])

def band(cs, lo=70, hi=98): return [c for c in cs if lo <= c["score"] <= hi]

def props(c): return {"logp": c["logp"], "mw": c["mw"], "tpsa": c.get("tpsa") or min(200, 25 * c["hbd"] + 12 * c["hba"])}

ORGANS = ["liver", "kidney", "heart", "brain", "lung", "gut", "muscle", "fat", "skin"]

def exposure(r):
    """organ exposure: peak tissue conc (mg/L) and tissue:plasma AUC ratio from the PBPK run"""
    dt = r["t"][1] - r["t"][0]; ap = sum(r["cp"]) * dt or 1e-12
    return {o: (max(r["tissue"][o]), sum(r["tissue"][o]) * dt / ap) for o in ORGANS}

def rat_study(c, mgkg=5.0):
    p = props(c); kg = pbpk.SPECIES["rat"]["bw"]; po, iv, F = pbpk.bioavailability(p, "rat", mgkg * kg, 24.0)
    ok = F >= .1 and (po["t_half"] or 0) >= .5
    return {"mgkg": mgkg, "F": F, "cmax": po["cmax"], "t_half": po["t_half"], "auc": po["auc"], "exp": exposure(po), "pass": ok,
            "why": "" if ok else f"rat F={F:.2f}, t1/2={po['t_half'] and round(po['t_half'], 1)}h"}

def human_study(c, max_daily_mg=1000.0):
    p = props(c); po, iv, F = pbpk.bioavailability(p, "human", 100.0, 48.0); ic = c.get("ic50_nM"); dose = None
    if ic: dose, _ = pbpk.dose_for_target(p, "human", ic, c["mw"])
    th = po["t_half"] or 0; fl = []; ex = exposure(po); a = c["alerts"]
    if c.get("herg_flag"): fl.append("heart: basic amine + lipophilic (hERG/QT pharmacophore)")
    if a["brenk"]: fl.append("liver/general: reactive-group alerts: " + ", ".join(a["brenk"][:3]))
    if ex["brain"][1] > .5: fl.append("brain: substantial CNS exposure (brain:plasma %.1f)" % ex["brain"][1])
    if po["frac_renal"] > .5: fl.append("kidney: mainly renally cleared (dose-adjust in renal impairment)")
    if po["params"]["fa"] < .7: fl.append("gut: incomplete absorption predicted")
    why = [] if F >= .1 else [f"human F={F:.2f}"]
    if th and not 2 <= th <= 72: why.append(f"t1/2={th:.1f}h outside 2-72h")
    if dose and dose > max_daily_mg: why.append(f"needs {dose:.0f} mg/day (> {max_daily_mg:.0f})")
    return {"F": F, "t_half": th, "cmax": po["cmax"], "dose_mg_day": dose, "exp": ex, "flags": fl, "pass": not why, "why": "; ".join(why), "params": po["params"], "psrc": po["param_source"]}

def pipeline(d, top=5):
    cls = cls_for(d); real = {}; names = [n for n, v in LIB.items() if not cls or v[4] == cls]
    try:
        import sources; real = sources.enrich(names)
    except Exception: pass
    mode = "fallback (property-only; pip install rdkit for real structures)"; cs = None
    try:
        import chemistry
        if chemistry.HAVE:
            leads = {n: {"smiles": (real.get(n) or {}).get("smiles"), "ic50_nM": (real.get(n) or {}).get("ic50_nM"), "data": "PubChem/ChEMBL" if n in real else "built-in"} for n in names}
            cs = chemistry.candidates(leads); mode = "REAL STRUCTURES (RDKit)" if cs else mode
    except Exception as e: print("! chemistry error:", e)
    if not cs: cs = fallback(cls, real)
    for c in cs: c.setdefault("herg_flag", False)
    b = band(cs); print(f"\n{NOTE}\nDisease: {d} -> class: {cls or 'all'}\nStructures: {mode}")
    print(f"Data: {len(real)} leads from PubChem/ChEMBL" if real else "Data: offline - approximate built-in values")
    print(f"Candidates: {len(cs)} | in 70-98 band: {len(b)}\n\nRANKED (heuristic score):")
    for c in b[:top]: print(f"  {c['score']:>5}  {c['name']:<18} MW{c['mw']:.0f} logP{c['logp']:.1f} {'(known drug)' if c['known'] else '(analog of ' + c['parent'] + ')'}")
    surv = []
    for c in b[:top * 3]:
        if c["alerts"]["pains"]: c["stage"] = "failed structure filter: PAINS " + c["alerts"]["pains"][0]; continue
        r = rat_study(c)
        if not r["pass"]: c["stage"] = "failed rat PBPK: " + r["why"]; continue
        h = human_study(c)
        if not h["pass"]: c["stage"] = "failed human PBPK: " + h["why"]; continue
        c["stage"] = "PASSED simulated stages"; surv.append((c, r, h))
    print("\nATTRITION:")
    for c in b[:top * 3]: print(f"  {c['name']:<18} {c['stage']}")
    for c, r, h in surv[:top]:
        t = lambda x: "n/a" if x is None else f"{x:.1f}"
        print(f"\n>> {c['name']}  {c.get('smiles') or ''}")
        print(f"   rat {r['mgkg']} mg/kg po: Cmax {r['cmax']:.2f} mg/L, t1/2 {t(r['t_half'])} h, F {r['F']:.2f}")
        print(f"   human 100 mg po: Cmax {h['cmax']:.2f} mg/L, t1/2 {t(h['t_half'])} h, F {h['F']:.2f}; dose to hold free trough >= IC50: " +
              (f"{h['dose_mg_day']:.0f} mg/day (IC50 {c['ic50_nM']:.3g} nM from {c['data']}, assumed same for analogs)" if h["dose_mg_day"] else "unknown (no measured potency)"))
        print("   organ:plasma AUC ratio: " + ", ".join(f"{o} {v[1]:.1f}" for o, v in h["exp"].items()))
        for f in h["flags"]: print("   ! " + f)
    if surv:
        try:
            import chemistry
            if chemistry.HAVE:
                for c, _, _ in surv[:top]: chemistry.svg(c["smiles"], os.path.join("medscope_out", c["name"] + ".svg"))
                print("\nStructure drawings saved to ./medscope_out/*.svg")
        except Exception: pass
    print("\nParameters not measured are PREDICTED (low confidence). A 'pass' means: worth REAL lab testing. Not evidence of a cure.\n"); return surv

if __name__ == "__main__":
    pipeline(" ".join(sys.argv[1:]) or "hiv")
