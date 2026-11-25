"""humanbio.py - human biology layer: (1) virtual population (age, sex, body size, kidney decline, enzyme variability,
CYP genetic polymorphism phenotypes, albumin), (2) where the drug's human target is expressed in the body (Open Targets) vs where the
drug actually goes (PBPK organ exposure) -> on-target toxicity map. Frequencies/scalars are approximate literature values (European)."""
import json, math, random
import pbpk, sources

# (population frequency, metabolic activity factor, label) - approximate
PHENO = {"CYP2D6": [(.07, .1, "PM"), (.10, .5, "IM"), (.77, 1.0, "NM"), (.06, 2.5, "UM")],
         "CYP2C19": [(.03, .1, "PM"), (.26, .5, "IM"), (.40, 1.0, "NM"), (.31, 1.6, "RM/UM")],
         "CYP2C9": [(.03, .2, "PM"), (.30, .6, "IM"), (.67, 1.0, "NM")]}

def population(n=200, seed=1, cyp=None, fm=.7):
    """virtual adults. fm = fraction of clearance via `cyp` (polymorphic enzymes only matter in proportion to this)"""
    rng = random.Random(seed); out = []
    for _ in range(n):
        sex = rng.choice("MF"); age = rng.uniform(18, 80); h = rng.gauss(1.76, .07) if sex == "M" else rng.gauss(1.63, .06)
        bw = max(40, min(150, rng.lognormvariate(math.log(26), .17) * h * h)); ph = "NM"; cl = math.exp(rng.gauss(0, .3))
        if cyp in PHENO:
            r, acc = rng.random(), 0.0
            for f, fac, lab in PHENO[cyp]:
                acc += f
                if r <= acc: cl *= fm * fac + (1 - fm); ph = lab; break
        out.append({"sex": sex, "age": age, "bw": bw, "qc_scale": math.exp(rng.gauss(0, .1)) * (1 - .004 * max(0, age - 30)),
                    "gfr_scale": max(.25, 1 - .0085 * max(0, age - 40)) * math.exp(rng.gauss(0, .15)), "clint_scale": cl,
                    "fup_scale": math.exp(rng.gauss(0, .1)), "phenotype": ph})
    return out

def _pct(x, q): x = sorted(x); return x[min(len(x) - 1, int(q * (len(x) - 1) + .5))]

def population_pk(p, dose_mg, ovr=None, subs=None, ic50_nM=None, mw=None, hours=72.0):
    """simulate every virtual subject; returns percentiles + % of subjects whose steady-state free trough covers the IC50"""
    subs = subs or population(); rows = []
    for s in subs:
        r = pbpk.simulate(p, "human", dose_mg, hours, ovr=ovr, subj=s); row = {"cmax": r["cmax"], "auc": r["auc"], "t_half": r["t_half"] or float("nan"), "ph": s["phenotype"]}
        if ic50_nM and mw: row["trough_nM"] = min(pbpk.steady_state(r, 24.0)) * r["params"]["fup"] / mw * 1e6
        rows.append(row)
    S = lambda k: {q: _pct([r[k] for r in rows if r[k] == r[k]], v) for q, v in (("p5", .05), ("p50", .5), ("p95", .95))}
    out = {"n": len(rows), "cmax": S("cmax"), "auc": S("auc"), "t_half": S("t_half")}
    if ic50_nM and mw:
        out["trough_nM"] = S("trough_nM"); out["pct_covered"] = 100 * sum(r["trough_nM"] >= ic50_nM for r in rows) / len(rows)
        out["by_phenotype"] = {}
        for ph in sorted({r["ph"] for r in rows}):
            g = [r for r in rows if r["ph"] == ph]; out["by_phenotype"][ph] = (len(g), 100 * sum(r["trough_nM"] >= ic50_nM for r in g) / len(g))
    return out

ORGAN_MAP = {"liver": ["liver"], "kidney": ["kidney"], "heart": ["heart", "cardiac"], "brain": ["brain", "cortex", "cerebell"], "lung": ["lung"],
             "gut": ["intestin", "colon", "stomach", "duoden", "ileum"], "muscle": ["muscle"], "fat": ["adipose"], "skin": ["skin"]}
OT = "https://api.platform.opentargets.org/api/v4/graphql"

def _gql(q, v): return sources.J(OT, json.dumps({"query": q, "variables": v}), {"Content-Type": "application/json"})

def target_expression(symbol):
    """human target expression per organ (Open Targets RNA level) + recorded safety liabilities -> ({organ: level}, [events])"""
    h = _gql('query($q:String!){search(queryString:$q,entityNames:["target"],page:{index:0,size:1}){hits{id}}}', {"q": symbol})["data"]["search"]["hits"]
    if not h: return {}, []
    t = _gql("query($id:String!){target(ensemblId:$id){expressions{tissue{label} rna{level}} safetyLiabilities{event}}}", {"id": h[0]["id"]})["data"]["target"]
    lv = {}
    for e in t.get("expressions") or []:
        lab = e["tissue"]["label"].lower(); l = (e.get("rna") or {}).get("level") or 0
        for o, keys in ORGAN_MAP.items():
            if any(k in lab for k in keys): lv[o] = max(lv.get(o, 0), l)
    return lv, sorted({s["event"] for s in t.get("safetyLiabilities") or [] if s.get("event")})

def exposure_flags(exposure, levels, symbol, ratio=.5, level=3):
    """exposure: {organ: (peak mg/L, tissue:plasma AUC ratio)} from PBPK. Flag organs that BOTH express the target AND receive drug"""
    return [f"{o}: {symbol} expressed (level {lv}) and drug reaches it (tissue:plasma {exposure[o][1]:.1f}) -> on-target effects likely"
            for o, lv in sorted(levels.items()) if lv >= level and o in exposure and exposure[o][1] >= ratio]
