#!/usr/bin/env python3
"""mass_population.py - drug exposure across a MASS virtual population (10,000+) with age, sex, weight, ancestry-specific genes,
liver/kidney/heart disease, a main disease, concomitant enzyme inducers, food, and everything the PBPK can scale.
  python mass_population.py <drug name|SMILES> [--n 10000] [--dose 100] [--ic50 nM] [--ancestry african|european|east_asian] [--disease hiv|tb|cancer|inflammation|malaria]
                            [--peds] [--food fasted|light|heavy] [--rifampicin] [--tau 24]      (MedScope: population ...)
  python mass_population.py person age=70 sex=F weight=55 cyp2d6=PM ckd=3 liver=B food=heavy disease=tb <drug>
Method: PBPK is simulated in full for a stratified sample (<=300), then a log-linear surrogate (R2 reported, held-out) predicts exposure
for every subject. Scaling factors are APPROXIMATE literature values (verify before clinical use). Exposure thresholds are proxies, not tox data."""
import math, random, sys
import pbpk, state

# ancestry-specific genotype frequencies (approximate literature): CYP2D6/2C19/2C9 phenotype classes, CYP3A5 expresser, NAT2 slow, UGT1A1 reduced
GENO = {
 "african":    dict(CYP2D6=(.03, .30, .60, .07), CYP2C19=(.04, .30, .50, .16), CYP2C9=(.02, .20, .78), CYP3A5=.70, NAT2=.30, UGT1A1=.30),
 "european":   dict(CYP2D6=(.07, .10, .77, .06), CYP2C19=(.03, .26, .40, .31), CYP2C9=(.03, .30, .67), CYP3A5=.10, NAT2=.55, UGT1A1=.12),
 "east_asian": dict(CYP2D6=(.01, .50, .47, .02), CYP2C19=(.13, .46, .38, .03), CYP2C9=(.02, .10, .88), CYP3A5=.30, NAT2=.15, UGT1A1=.15)}
ACT = {"CYP2D6": (.05, .5, 1, 2.5), "CYP2C19": (.05, .5, 1, 1.6), "CYP2C9": (.2, .6, 1)}   # PM, IM, NM, UM/RM activity factors
LAB = {"CYP2D6": ("PM", "IM", "NM", "UM"), "CYP2C19": ("PM", "IM", "NM", "UM"), "CYP2C9": ("PM", "IM", "NM")}
# disease / condition effects on physiology: (qc, hepq, gfr, fup, ka, enzyme-mult) approximate
LIVER = {"A": dict(hepq=.9, fup=1.15, enz=.7), "B": dict(hepq=.7, fup=1.5, gfr=.85, enz=.4), "C": dict(hepq=.5, fup=2.0, gfr=.6, enz=.2)}
CKD = {"2": dict(gfr=.75), "3": dict(gfr=.45, fup=1.1), "4": dict(gfr=.2, fup=1.2, enz=.85)}
HF = dict(qc=.7, hepq=.6, gfr=.7, ka=.8)
MAIN = {"hiv": dict(bw=.93, fup=1.1, ka=.85), "tb": dict(bw=.88, fup=1.2, ka=.8), "cancer": dict(bw=.9, fup=1.3, enzall=.7),
        "inflammation": dict(fup=1.15, enzall=.65), "malaria": dict(hepq=.85, fup=1.2, ka=.7), "covid": dict(fup=1.2, enzall=.7), "diabetes": dict(ka=.85)}
INDUCER = {"CYP3A4": 3.0, "CYP2C19": 2.0, "CYP2C9": 1.8, "CYP2D6": 1.0}   # rifampicin
FOOD = {"fasted": (1, 1, 1), "light": (.7, 1, 1.07), "heavy": (.45, 1, 1.15)}   # ka, fa(base), hepatic flow

def pick(r, p, labels):
    a = 0.0
    for q, l in zip(p, labels):
        a += q
        if r <= a: return l
    return labels[-1]

def make(n=10000, seed=1, ancestry="african", disease=None, peds=False, food=None, rifampicin=False, fm=None, force=None):
    """virtual subjects; `force` = dict of fixed attributes (age, sex, weight, cyp2d6, ckd, liver, hf, food, disease)"""
    rng = random.Random(seed); G = GENO[ancestry]; fm = fm or {}; out = []; F = {k.lower(): v for k, v in (force or {}).items()}
    for _ in range(n):
        sex = F.get("sex") or rng.choice("MF"); age = float(F["age"]) if "age" in F else (rng.uniform(2, 85) if peds and rng.random() < .25 else min(90, max(18, rng.gauss(42, 18))))
        if age < 18: bw = 2 * age + 8 if age <= 10 else 28 + 3.3 * (age - 10); bw *= math.exp(rng.gauss(0, .12)); h = None
        else:
            h = rng.gauss(1.76, .07) if sex == "M" else rng.gauss(1.63, .06); bmi = math.exp(rng.gauss(math.log(24.5 if sex == "M" else 26.5), .2)); bw = bmi * h * h
        if "weight" in F: bw = float(F["weight"])
        bw = max(8, min(200, bw)); bmi = bw / (h * h) if h else None
        geno, enz = {}, {}
        for e in ("CYP2D6", "CYP2C19", "CYP2C9"):
            lab = F.get(e.lower()) or pick(rng.random(), G[e], LAB[e]); geno[e] = lab; enz[e] = ACT[e][LAB[e].index(lab)] * math.exp(rng.gauss(0, .15))
        geno["CYP3A5"] = "expresser" if rng.random() < G["CYP3A5"] else "non-expresser"; enz["CYP3A5"] = 1.0
        geno["NAT2"] = "slow" if rng.random() < G["NAT2"] else "rapid"; enz["NAT2"] = .4 if geno["NAT2"] == "slow" else 1.0
        geno["UGT1A1"] = "reduced" if rng.random() < G["UGT1A1"] else "normal"; enz["UGT1A1"] = .6 if geno["UGT1A1"] == "reduced" else 1.0
        enz["CYP3A4"] = math.exp(rng.gauss(0, .3)) * (1.5 if geno["CYP3A5"] == "expresser" else 1.0) ** .5
        mod = dict(qc=1.0, hepq=1.0, gfr=1.0, fup=1.0, ka=1.0, enz=1.0, bw=1.0, enzall=1.0); cond = []
        def apply(d):
            for k, v in d.items(): mod[k] = mod.get(k, 1.0) * v
        lv = F.get("liver") or (rng.choices("ABC", [.5, .35, .15])[0] if rng.random() < (.03 + .0004 * max(0, age - 40)) else None)
        ck = str(F["ckd"]) if "ckd" in F else (rng.choices("234", [.6, .32, .08])[0] if rng.random() < (.04 + .006 * max(0, age - 40)) else None)
        ck = None if ck in ("0", "none") else ck
        hf = F.get("hf") or rng.random() < (.01 + .004 * max(0, age - 50))
        if lv: apply(LIVER[lv]); cond.append(f"cirrhosis CP-{lv}")
        if ck: apply(CKD[ck]); cond.append(f"CKD{ck}")
        if hf and str(hf).lower() not in ("0", "false"): apply(HF); cond.append("heart failure")
        dz = F.get("disease") or disease
        if dz in MAIN: apply(MAIN[dz]); cond.append(dz)
        if bmi and bmi < 18.5: apply(dict(fup=1.25, enz=.85)); cond.append("underweight/malnourished")
        riff = F.get("rifampicin", rifampicin) or (dz == "tb" and False)
        if str(riff).lower() in ("1", "true", "yes"): cond.append("rifampicin")
        for e, m in INDUCER.items():
            if cond.count("rifampicin"): enz[e] *= m
        fd = F.get("food") or food or rng.choice(["fasted", "fasted", "light", "heavy"]); ka_s, _, hq = FOOD[fd]
        mat = age / (age + .5) if age < 18 else 1.0   # enzyme maturation (approximate)
        gfr_age = (1.0 if age < 18 else max(.25, 1 - .0085 * max(0, age - 40))) * math.exp(rng.gauss(0, .15))
        # hepatic clearance factor: weighted by metabolic fractions fm (enzyme -> fraction of hepatic CL); the rest is unaffected by genes
        tot = sum(fm.values()); g = 1 - tot
        for e, f in fm.items(): g += f * enz.get(e, 1.0) * mod["enzall"]
        g = (g if tot else mod["enzall"]) * mod["enz"] * mat * math.exp(rng.gauss(0, .25))
        out.append(dict(age=age, sex=sex, bw=bw * mod["bw"] if "weight" not in F else bw, bmi=bmi, geno=geno, cond=cond, food=fd, ancestry=ancestry,
                        qc_scale=mod["qc"] * math.exp(rng.gauss(0, .1)) * (1 - .004 * max(0, age - 30)), hepq_scale=mod["hepq"] * hq, gfr_scale=gfr_age * mod["gfr"],
                        fup_scale=mod["fup"] * math.exp(rng.gauss(0, .1)), ka_scale=mod["ka"] * ka_s, fa_scale=1.0, clint_scale=max(.02, g)))
    return out

def adjust_food(subs, p):
    """food effect depends on the drug: heavy meal raises absorption of lipophilic drugs, lowers that of very polar ones (approximate)"""
    for s in subs: s["fa_scale"] = (1.25 if p["logp"] > 3 else .9 if p["logp"] < 1 else 1.0) if s["food"] == "heavy" else (1.1 if s["food"] == "light" and p["logp"] > 3 else 1.0)
    return subs

FEAT = ["bw", "clint_scale", "gfr_scale", "fup_scale", "fa_scale", "ka_scale", "hepq_scale", "qc_scale"]
def X(s): return [1.0] + [math.log(max(1e-3, s[k])) for k in FEAT]

def solve(A, b):
    n = len(b); M = [r[:] + [b[i]] for i, r in enumerate(A)]
    for i in range(n):
        p = max(range(i, n), key=lambda r: abs(M[r][i])); M[i], M[p] = M[p], M[i]
        for r in range(n):
            if r != i and M[i][i]:
                f = M[r][i] / M[i][i]; M[r] = [a - f * c for a, c in zip(M[r], M[i])]
    return [M[i][n] / M[i][i] if M[i][i] else 0.0 for i in range(n)]

def fit(rows, ys, lam=1e-6):
    k = len(rows[0]); A = [[sum(r[i] * r[j] for r in rows) + (lam if i == j else 0) for j in range(k)] for i in range(k)]
    return solve(A, [sum(r[i] * y for r, y in zip(rows, ys)) for i in range(k)])

def r2(b, rows, ys):
    pr = [sum(x * w for x, w in zip(r, b)) for r in rows]; m = sum(ys) / len(ys)
    return 1 - sum((a - c) ** 2 for a, c in zip(ys, pr)) / max(1e-12, sum((a - m) ** 2 for a in ys))

def surrogate(p, subs, ovr=None, n_full=240, seed=5, hours=60.0):
    """full PBPK on a sample -> log-linear surrogate of AUC and Cmax; returns (coef_auc, coef_cmax, R2 held-out, ref-run)"""
    rng = random.Random(seed); samp = rng.sample(subs, min(n_full, len(subs))); rows, la, lc = [], [], []
    for s in samp:
        r = pbpk.simulate(p, "human", 100.0, hours, ovr=ovr, subj=s, every=.5); rows.append(X(s)); la.append(math.log(r["auc"])); lc.append(math.log(r["cmax"]))
    k = int(len(rows) * .8); ba, bc = fit(rows[:k], la[:k]), fit(rows[:k], lc[:k])
    ra = r2(ba, rows[k:], la[k:]) if len(rows) - k > 5 else None; ba, bc = fit(rows, la), fit(rows, lc)
    return ba, bc, ra, pbpk.simulate(p, "human", 100.0, hours, ovr=ovr, every=.5)

def predict(b, s): return math.exp(sum(x * w for x, w in zip(X(s), b)))

def groups(s):
    a = s["age"]; g = ["age " + ("<18" if a < 18 else "18-39" if a < 40 else "40-59" if a < 60 else "60-74" if a < 75 else "75+"), "sex " + s["sex"], "food " + s["food"]]
    if s["bmi"]: g.append("BMI " + ("<18.5" if s["bmi"] < 18.5 else "<25" if s["bmi"] < 25 else "<30" if s["bmi"] < 30 else "30+"))
    for e in ("CYP2D6", "CYP2C19", "CYP2C9"): g.append(f"{e} {s['geno'][e]}")
    g += ["CYP3A5 " + s["geno"]["CYP3A5"], "NAT2 " + s["geno"]["NAT2"], "UGT1A1 " + s["geno"]["UGT1A1"]] + ["condition: " + c for c in s["cond"]]
    return g or ["all"]

def run(p, ovr=None, fm=None, mw=300.0, dose=100.0, ic50_nM=None, tau=24.0, n=10000, ancestry="african", disease=None, peds=False, food=None, rifampicin=False, tox_auc=3.0, tox_cmax=2.5, show=True, record=True):
    subs = adjust_food(make(n, 1, ancestry, disease, peds, food, rifampicin, fm), p); b_a, b_c, R2, ref = surrogate(p, subs, ovr)
    fup0 = ref["params"]["fup"]; css = pbpk.steady_state(ref, tau); tr = min(css) / (sum(css) / len(css)); aref, cref = ref["auc"] * dose / 100, ref["cmax"] * dose / 100
    for s in subs:
        s["auc"] = predict(b_a, s) * dose / 100; s["cmax"] = predict(b_c, s) * dose / 100; s["fup"] = min(1, fup0 * s["fup_scale"])
        s["trough_nM"] = s["auc"] / tau * tr * s["fup"] / mw * 1e6; s["high"] = s["auc"] > tox_auc * aref or s["cmax"] > tox_cmax * cref
        s["low"] = bool(ic50_nM) and s["trough_nM"] < ic50_nM
    G = {}
    for s in subs:
        for g in groups(s): G.setdefault(g, []).append(s)
    pct = lambda L, f: 100 * sum(1 for s in L if s[f]) / len(L); med = lambda L, k: sorted(s[k] for s in L)[len(L) // 2]
    res = {"n": n, "r2": R2, "overall": dict(high=pct(subs, "high"), low=pct(subs, "low") if ic50_nM else None, auc_med=med(subs, "auc") / aref), "groups": {}, "subs": subs, "ref": ref}
    for g, L in G.items():
        if len(L) >= max(30, n // 200):
            adj = max(.25, min(2.0, aref / med(L, "auc"))); res["groups"][g] = dict(n=len(L), auc_ratio=med(L, "auc") / aref, high=pct(L, "high"), low=pct(L, "low") if ic50_nM else None, dose_adj=adj)
    if show: report(res, dose, ic50_nM, aref, ancestry)
    if record:
        state.record("pop_n", n)
        if R2 is not None: state.record("pop_r2", round(R2, 3))
    return res

def report(res, dose, ic50, aref, anc):
    o = res["overall"]; print(f"\nMASS POPULATION: {res['n']:,} virtual people ({anc} ancestry), {dose:g} mg | surrogate held-out R2 vs full PBPK: " + (f"{res['r2']:.3f}" if res["r2"] is not None else "n/a"))
    print(f"  median exposure {o['auc_med']:.2f}x healthy-adult reference | HIGH exposure (AUC>3x or Cmax>2.5x reference; a proxy for toxicity risk): {o['high']:.1f}%" + (f" | BELOW effective trough (free >= IC50 {ic50:g} nM): {o['low']:.1f}%" if o["low"] is not None else ""))
    print(f"  {'subgroup (most unusual first)':<34}{'n':>6}{'AUC/ref':>8}{'%high':>7}{'%low':>7}{'dose x':>8}")
    for g, v in sorted(res["groups"].items(), key=lambda kv: -abs(math.log(kv[1]["auc_ratio"] / o["auc_med"])))[:18]:
        print(f"  {g:<34}{v['n']:>6}{v['auc_ratio']:>8.2f}{v['high']:>7.1f}{'' if v['low'] is None else format(v['low'], '7.1f')}{v['dose_adj']:>8.2f}")
    print("  'dose x' = multiplier that would bring the group's median exposure to the healthy-adult reference (clamped 0.25-2). Scaling factors are approximate; confirm in clinical PK studies.")

def fm_from_admet(a):
    try:
        import admet; en = [e for e in admet.cyp_substrate(a) if "CYP" in e]
    except Exception: return {}
    out = {}
    for e in en:
        for k in ("CYP3A4", "CYP2D6", "CYP2C9", "CYP2C19"):
            if k in e: out[k] = .8 / len(en)
    return out

def drug(arg, dose=100.0):
    """name or SMILES -> (p, ovr, fm, mw, label)"""
    import ingredients
    if any(c in arg for c in "()=#[]") and not arg.isalpha():
        from rdkit import Chem
        from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
        m = Chem.MolFromSmiles(arg); p = dict(logp=Crippen.MolLogP(m), mw=Descriptors.MolWt(m), tpsa=rdMolDescriptors.CalcTPSA(m), ion=pbpk.ion_class(arg)); ovr = a = None
        try:
            import admet
            if admet.HAVE: a = admet.predict([arg])[arg]; ovr = admet.overrides(a)
        except Exception: pass
        return p, ovr, fm_from_admet(a) if a else {}, p["mw"], arg
    x = ingredients.profile(arg, dose); return x["props"], x["ovr"], fm_from_admet(x["admet"]) if x["admet"] else {}, x["props"]["mw"], arg

def selftest():
    """offline directional checks: known physiology must move exposure the right way; surrogate must track the full PBPK"""
    p = dict(logp=2.5, mw=350.0, tpsa=70.0); ok = []
    def auc(force=None, ovr=None, fm=None, food=None):
        s = adjust_food(make(1, 3, "european", force=force, fm=fm, food=food or "fasted"), p)[0]
        return pbpk.simulate(p, "human", 100.0, 72.0, ovr=ovr, subj=s, every=.5)
    ren = dict(clr=40.0, clint_hep=.3, _src={})
    a0 = auc(ovr=ren)["auc"]; a1 = auc({"ckd": 4}, ovr=ren)["auc"]; ok.append(("CKD4 raises exposure of a renally cleared drug", a1 > 1.5 * a0))
    fm = {"CYP2D6": .8}; s0 = adjust_food(make(1, 3, "european", force={"cyp2d6": "NM", "ckd": 0}, fm=fm, food="fasted"), p)[0]; s1 = adjust_food(make(1, 3, "european", force={"cyp2d6": "PM", "ckd": 0}, fm=fm, food="fasted"), p)[0]
    for s in (s0, s1): s.update(age=40.0, qc_scale=1, gfr_scale=1, fup_scale=1, hepq_scale=1, clint_scale=s["clint_scale"])
    ok.append(("CYP2D6 poor metabolisers have higher exposure to a CYP2D6 substrate", pbpk.simulate(p, "human", 100.0, 72.0, subj=s1, every=.5)["auc"] > 1.5 * pbpk.simulate(p, "human", 100.0, 72.0, subj=s0, every=.5)["auc"]))
    f0 = auc(food="fasted"); f1 = auc(food="heavy"); ok.append(("heavy meal delays/lowers Cmax", f1["cmax"] < f0["cmax"] and f1["tmax"] >= f0["tmax"]))
    c0 = auc(); c1 = auc({"liver": "C"}); ok.append(("severe cirrhosis raises exposure of a hepatically cleared drug", c1["auc"] > 1.3 * c0["auc"]))
    ok.append(("PBPK still conserves mass with all scalings", abs(f1["mass_balance"] - 1) < 1e-6 and abs(c1["mass_balance"] - 1) < 1e-6))
    res = run(p, mw=350.0, n=10000, ancestry="african", show=False); ok.append((f"10,000 subjects simulated, surrogate held-out R2 {res['r2']:.3f} >= 0.9", res["n"] >= 10000 and res["r2"] is not None and res["r2"] >= .9))
    ok.append(("ancestry changes genotype mix (CYP3A5 expresser: African > European)", sum(s["geno"]["CYP3A5"] == "expresser" for s in make(3000, 2, "african")) > 2 * sum(s["geno"]["CYP3A5"] == "expresser" for s in make(3000, 2, "european"))))
    for t, v in ok: print(f"  [{'PASS' if v else 'FAIL'}] {t}")
    state.record("pop_ok", all(v for _, v in ok)); return all(v for _, v in ok)

def cli(a):
    if not a or a[0] == "selftest": return selftest()
    kv = {x.split("=")[0]: x.split("=")[1] for x in a if "=" in x and not x.startswith("--")}; flags = {a[i]: a[i + 1] for i in range(len(a) - 1) if a[i].startswith("--") and not a[i + 1].startswith("--")}
    rest = [x for i, x in enumerate(a) if not x.startswith("--") and "=" not in x and not (i and a[i - 1].startswith("--") and a[i - 1] in flags)]
    person = bool(rest) and rest[0] == "person"; rest = rest[1:] if person else rest
    if not rest: print(__doc__); return
    dose = float(flags.get("--dose", 100)); p, ovr, fm, mw, nm = drug(" ".join(rest), dose); ic = float(flags["--ic50"]) if "--ic50" in flags else None
    if person:
        s = adjust_food(make(1, 1, kv.get("ancestry", "african"), kv.get("disease"), force=kv, fm=fm), p)[0]; r = pbpk.simulate(p, "human", dose, 72.0, ovr=ovr, subj=s, every=.5); ref = pbpk.simulate(p, "human", dose, 72.0, ovr=ovr, every=.5)
        print(f"{nm}: {dose:g} mg | person: age {s['age']:.0f}, {s['sex']}, {s['bw']:.0f} kg, {s['food']}, genes {s['geno']}, conditions {s['cond'] or 'none'}")
        print(f"  Cmax {r['cmax']:.3f} mg/L ({r['cmax'] / ref['cmax']:.2f}x reference adult) | AUC {r['auc']:.2f} ({r['auc'] / ref['auc']:.2f}x) | t1/2 {r['t_half'] and round(r['t_half'], 1)} h (reference {ref['t_half'] and round(ref['t_half'], 1)} h)")
        print(f"  suggested dose multiplier for equal exposure: {max(.25, min(2, ref['auc'] / r['auc'])):.2f}x"); return
    run(p, ovr, fm, mw, dose, ic, float(flags.get("--tau", 24)), int(flags.get("--n", 10000)), flags.get("--ancestry", "african"), flags.get("--disease"), "--peds" in a, flags.get("--food"), "--rifampicin" in a)

if __name__ == "__main__": cli(sys.argv[1:])
