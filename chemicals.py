#!/usr/bin/env python3
"""chemicals.py - TOY in-silico drug-candidate pipeline: generate -> rank -> lab -> rat -> human.
Stdlib only. EDUCATIONAL MODEL: scores are heuristic 'computational plausibility', NOT probabilities of
efficacy or safety. Candidates are property-space analogs of known drug-like leads, not real structures.
Real discovery needs docking/ADMET software, wet-lab assays, animal studies and clinical trials."""
import math, random, sys

NOTE = "SCREENING MODEL - heuristic scores, NOT cure probabilities. Unvalidated; not medical advice."
# approximate textbook-level properties: name: (MW, logP, HBD, HBA, class)
LIB = {"zidovudine": (267, 0.1, 2, 6, "antiviral"), "dolutegravir": (419, 1.0, 2, 7, "antiviral"),
       "tenofovir": (287, -1.6, 2, 8, "antiviral"), "imatinib": (494, 3.5, 2, 7, "kinase"),
       "erlotinib": (393, 3.3, 1, 7, "kinase"), "artemether": (298, 3.5, 0, 5, "antimalarial"),
       "chloroquine": (320, 4.6, 1, 3, "antimalarial"), "ciprofloxacin": (331, 0.3, 2, 6, "antibiotic"),
       "linezolid": (337, 0.6, 1, 6, "antibiotic")}
CLS = {"hiv": "antiviral", "virus": "antiviral", "hepatitis": "antiviral", "covid": "antiviral", "cancer": "kinase",
       "leukemia": "kinase", "tumor": "kinase", "malaria": "antimalarial", "tuberculosis": "antibiotic", "bacter": "antibiotic"}
ORG = ["liver", "kidney", "heart", "brain", "lungs", "gut", "blood", "immune", "skin"]

def cls_for(d):
    d = d.lower(); return next((v for k, v in CLS.items() if k in d), None)

def lip(c): return sum([c["mw"] > 500, c["logp"] > 5, c["hbd"] > 5, c["hba"] > 10])

def traits(c):
    """toy potency / toxicity / drug-likeness from properties"""
    fit = math.exp(-((c["logp"] - 2.5) / 2.2) ** 2) * math.exp(-((c["mw"] - 350) / 160) ** 2)
    c["pot"] = round(0.1 + 0.8 * fit, 3)
    if c.get("ic50_nM"): c["pot"] = round(max(.05, min(1, (9 - math.log10(c["ic50_nM"]) - 4) / 3.5)), 3)  # REAL ChEMBL potency
    elif c.get("anchor") is not None: c["pot"] = round(.5 * c["pot"] + .5 * c["anchor"], 3)  # analog inherits half of parent's real potency
    c["dl"] = round(1 - 0.25 * lip(c), 2)
    c["tox"] = round(min(1, 0.12 + 0.1 * max(0, c["logp"] - 3) + 0.2 * (c["mw"] > 500) + 0.05 * max(0, c["hbd"] - 3)), 3)
    c["score"] = round(min(98, 100 * (0.5 * c["pot"] + 0.3 * c["dl"] + 0.2 * (1 - c["tox"]))), 1); return c

def generate(cls=None, n=30, seed=7, real=None):
    """real: {name: {mw,logp,hbd,hba,ic50_nM}} from sources.enrich() overrides built-in approximations"""
    r = random.Random(seed); out = []; real = real or {}
    for nm, (mw, lp, hd, ha, k) in LIB.items():
        if cls and k != cls: continue
        x = real.get(nm); src = "PubChem/ChEMBL" if x else "built-in approx"
        if x: mw, lp, hd, ha = x["mw"], x["logp"], x["hbd"], x["hba"]
        p = traits({"name": nm, "mw": mw, "logp": lp, "hbd": hd, "hba": ha, "cls": k, "known": True, "ic50_nM": (x or {}).get("ic50_nM"), "data": src}); out.append(p)
        for i in range(n):
            out.append(traits({"name": f"{nm}-v{i+1}", "mw": round(mw * r.uniform(.8, 1.25)), "logp": round(lp + r.uniform(-1, 1), 1),
                               "hbd": max(0, hd + r.choice([-1, 0, 1])), "hba": max(1, ha + r.choice([-2, -1, 0, 1, 2])), "cls": k, "known": False,
                               "anchor": p["pot"] if x and x.get("ic50_nM") else None, "data": src}))
    return sorted(out, key=lambda c: -c["score"])

def band(cs, lo=70, hi=98): return [c for c in cs if lo <= c["score"] <= hi]

def lab(c, conc=(.001, .01, .1, 1, 10, 100)):
    """in-vitro sim: Hill dose-response on target cells + cytotoxicity on healthy cells (uM)"""
    ic50 = 10 ** (-(4 + 3.5 * c["pot"]) + 6); cc50 = 10 ** (-(4.2 - 1.2 * c["tox"] * 2) + 6)
    hill = lambda x, m: x ** 1.2 / (x ** 1.2 + m ** 1.2)
    return {"ic50_uM": ic50, "cc50_uM": cc50, "SI": cc50 / ic50, "curve": [(x, 100 * hill(x, ic50), 100 * (1 - hill(x, cc50))) for x in conc],
            "pass": cc50 / ic50 >= 10}

def pk(c, dose_mg, kg, vd_lkg, th, f=None, seed=1, h=24, dt=.25):
    """one-compartment oral PK -> (Cmax mg/L, AUC, [concentration series])"""
    f = f or max(.1, min(1, .9 * c["dl"] * (1.1 if c["logp"] > 0 else .7))); k = math.log(2) / th; ka = 1.5
    v = vd_lkg * kg; a = dose_mg * f; cs = []; t = 0
    while t <= h:
        cs.append(a * ka / (v * (ka - k)) * (math.exp(-k * t) - math.exp(-ka * t))); t += dt
    return max(cs), sum(cs) * dt, cs

def rat(c, mgkg=2, n=10, seed=3):
    """virtual rat study: n animals, inter-animal variability, toxicity threshold from tox trait"""
    r = random.Random(seed); lim = 2.5 * (1 - c["tox"]) ** 2 + .05; res = []; th0 = 1.5 + .8 * max(c["logp"], 0)
    for _ in range(n):
        kg = r.gauss(.25, .02); cm, auc, _ = pk(c, mgkg * kg, kg, .7, th0 * r.lognormvariate(0, .2))
        res.append((cm, auc, cm < lim * 2))
    return {"dose_mgkg": mgkg, "Cmax": sum(x[0] for x in res) / n, "AUC": sum(x[1] for x in res) / n, "survived": sum(x[2] for x in res),
            "n": n, "tox_limit": lim, "pass": sum(x[2] for x in res) >= .9 * n and sum(x[0] for x in res) / n > .05}

def human(c, mg=100, kg=70):
    """virtual human: PK + nine organ systems; each gets exposure-driven risk and optional reaction text"""
    th = (1.5 + .8 * max(c["logp"], 0)) * 3; cm, auc, _ = pk(c, mg, kg, 1.0, th); lim = 2.5 * (1 - c["tox"]) ** 2 + .05
    bbb = c["mw"] < 450 and c["logp"] > 1 and c["hbd"] <= 3; x = cm / lim
    risk = {"liver": x * (.5 + .15 * max(c["logp"], 0)), "kidney": x * (.4 + (.3 if c["logp"] < 0 else 0)), "heart": x * (.3 + .1 * max(c["logp"] - 3, 0)),
            "brain": x * (.5 if bbb else .05), "lungs": x * .15, "gut": x * (.6 if c["logp"] > 2 else .3), "blood": x * (.3 + .1 * c["tox"] * 5),
            "immune": x * .2, "skin": x * .15}
    msg = {"liver": "raised liver enzymes", "kidney": "reduced clearance / kidney strain", "heart": "QT-interval prolongation risk", "brain": "CNS effects (dizziness)",
           "lungs": "cough / airway irritation", "gut": "nausea, GI upset", "blood": "low blood-cell counts", "immune": "hypersensitivity", "skin": "rash"}
    rx = [f"{o}: {msg[o]}" for o in ORG if risk[o] > .6]; eff = min(90, 100 * c["pot"] * min(1, cm / .6))
    return {"Cmax": cm, "AUC": auc, "t_half_h": th, "crosses_BBB": bbb, "risk": {o: round(min(1, v), 2) for o, v in risk.items()},
            "reactions": rx, "efficacy_idx": round(eff, 1), "pass": not rx and eff > 40}

def pipeline(d, top=5):
    cls = cls_for(d); real = {}
    try:
        import sources; real = sources.enrich([n for n, v in LIB.items() if not cls or v[4] == cls])
    except Exception: pass
    cs = generate(cls, real=real); b = band(cs); print(f"\n{NOTE}\nDisease: {d} -> class: {cls or 'all'}")
    print(f"Data: {len(real)} leads use REAL PubChem/ChEMBL values; others use built-in approximations" if real else "Data: offline - built-in approximations (connect to internet for real PubChem/ChEMBL data)")
    print(f"Candidates: {len(cs)} | in 70-98 band: {len(b)}\n\nRANKED (plausibility score 0-100):")
    for c in b[:top]: print(f"  {c['score']:>5}  {c['name']:<18} MW{c['mw']} logP{c['logp']} {'(known drug)' if c['known'] else '(hypothetical analog)'}")
    surv = []
    for c in b[:top * 3]:
        l = lab(c); r = rat(c) if l["pass"] else None; h = human(c) if r and r["pass"] else None
        c["stage"] = "failed lab (poor selectivity)" if not l["pass"] else "failed rat tox/exposure" if not r["pass"] else "failed human sim: " + "; ".join(h["reactions"] or ["low efficacy"]) if not h["pass"] else "PASSED all simulated stages"
        if h and h["pass"]: surv.append((c, l, r, h))
    print("\nATTRITION (top candidates through virtual stages):")
    for c in b[:top * 3]: print(f"  {c['name']:<18} {c['stage']}")
    for c, l, r, h in surv[:top]:
        print(f"\n>> {c['name']}: IC50 {l['ic50_uM']:.3g}uM SI {l['SI']:.0f} | rat Cmax {r['Cmax']:.2f} survived {r['survived']}/{r['n']} | human t1/2 {h['t_half_h']:.1f}h efficacy-idx {h['efficacy_idx']} BBB {h['crosses_BBB']}")
        print("   organ risk:", ", ".join(f"{o} {v}" for o, v in h["risk"].items()))
    print("\nA 'pass' here means only: worth testing in a REAL lab. It is not evidence of a cure.\n"); return surv

if __name__ == "__main__":
    pipeline(" ".join(sys.argv[1:]) or "hiv")
