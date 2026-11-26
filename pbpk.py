"""pbpk.py - physiologically-based PK (PBPK) simulator: human (70 kg) and rat (250 g). Pure Python.
Perfusion-limited organs (gut, liver, kidney, heart, brain, muscle, fat, skin, rest, lung) + arterial/venous blood,
oral absorption, well-stirred hepatic metabolism, renal filtration, mass-balance checked.
Physiology: approximate literature values (Brown et al. 1997-style allometry). Drug parameters are PREDICTED from
structure descriptors (low confidence) unless you pass measured values via `ovr` (fup, clint, fa, ka, clr)."""
import math

SPECIES = {
    "human": dict(bw=70.0, qc=16.5, gfr=0.102, liver_g=25.7, mppgl=40, plasma=(0.96, 0.0035, 0.0023),
                  V=dict(gut=.017, liver=.026, kidney=.0044, heart=.0047, brain=.020, lung=.0076, muscle=.40, fat=.214, skin=.037, blood=.079),
                  Q=dict(gut=.151, liver_art=.076, kidney=.175, heart=.04, brain=.12, muscle=.17, fat=.05, skin=.058)),
    "rat": dict(bw=0.25, qc=14.1, gfr=0.312, liver_g=40.0, mppgl=45, plasma=(0.96, 0.0017, 0.0008),
                V=dict(gut=.027, liver=.034, kidney=.0073, heart=.0033, brain=.0057, lung=.005, muscle=.404, fat=.07, skin=.19, blood=.074),
                Q=dict(gut=.14, liver_art=.034, kidney=.141, heart=.049, brain=.02, muscle=.278, fat=.07, skin=.058)),
}
# approximate tissue composition (water, neutral lipid, phospholipid fractions) for Poulin-Theil-style partitioning
TC = {"gut": (.75, .04, .013), "liver": (.75, .014, .024), "kidney": (.78, .012, .029), "heart": (.78, .012, .016), "brain": (.79, .039, .005),
      "muscle": (.76, .010, .009), "fat": (.18, .79, .002), "skin": (.65, .06, .004), "rest": (.70, .04, .01), "lung": (.81, .02, .012)}
TIS = ["gut", "liver", "kidney", "heart", "brain", "muscle", "fat", "skin", "rest"]
ST = ["lumen", "gut", "liver", "kidney", "heart", "brain", "muscle", "fat", "skin", "rest", "lung", "art", "ven", "cl_liver", "cl_kidney", "fec"]
IX = {n: i for i, n in enumerate(ST)}
clamp = lambda x, a, b: max(a, min(b, x))

HEPC = {"human": 120.0, "rat": 117.0}   # hepatocytes, million cells per g liver (approx)

def predict(p, ovr=None):
    """drug parameters from descriptors p={logp,mw,tpsa}; ovr = measured/ADMET values replacing predictions:
    fup, clint (uL/min/mg microsomal), clint_hep (uL/min/1e6 hepatocytes), vss (L/kg), fa, ka, clr (L/h), efflux_brain"""
    o = ovr or {}; lp, mw, tp = p["logp"], p["mw"], p.get("tpsa", 80)
    fa = 1.0 if tp <= 100 and mw <= 500 else clamp(1 - max(tp - 100, 0) / 120 - (.3 if mw > 600 else 0), .05, 1)
    if lp < -1: fa *= .6
    d = {"fup": clamp(1 / (1 + 10 ** (.5 * lp - .2)), .01, 1), "clint": clamp(3 * 10 ** (.35 * lp), .5, 300), "fa": fa,
         "ka": 1.2 if tp <= 90 else .5, "efflux_brain": .3 if tp > 90 else 1.0, "clint_hep": None, "vss": None, "clr": None}
    src = {k: "rule-based" for k in d}
    for k, v in o.items():
        if k in d and v is not None: d[k] = v; src[k] = o.get("_src", {}).get(k, "measured/model")
    return d, src

RT = {"liver": .45, "kidney": .4, "heart": .3, "brain": .1, "muscle": .15, "fat": .05, "skin": .3, "rest": .3, "gut": .3, "lung": .3}  # tissue protein binding vs plasma (approx)

def kp(t, P, fup, plasma=None):
    """tissue:plasma partition. Unbound drug equilibrates with tissue water + lipids; plus tissue-protein binding:
    Kp = fup*(water + neutral-lipid*P + phospholipid*(0.3P+0.7)) + r_t*(1-fup)"""
    fw, fnl, fph = TC[t]
    return clamp(fup * (fw + fnl * P + fph * (.3 * P + .7)) + RT[t] * (1 - fup), .2, 60)

def build(p, species, ovr=None, subj=None):
    s = SPECIES[species]; sj = subj or {}; bw = sj.get("bw", s["bw"]); d, src = predict(p, ovr); d = dict(d)
    d["fup"] = clamp(d["fup"] * sj.get("fup_scale", 1), .003, 1); d["fa"] = clamp(d["fa"] * sj.get("fa_scale", 1), .03, 1); d["ka"] = d["ka"] * sj.get("ka_scale", 1); P = 10 ** clamp(p["logp"], -3, 6)
    CO = s["qc"] * bw ** .75 * sj.get("qc_scale", 1); V = {k: v * bw for k, v in s["V"].items()}
    V["art"], V["ven"] = V["blood"] / 3, V["blood"] * 2 / 3; V["lumen"] = 1.0
    V["rest"] = bw - sum(V[k] for k in ("gut", "liver", "kidney", "heart", "brain", "lung", "muscle", "fat", "skin", "blood"))
    Q = {k: v * CO for k, v in s["Q"].items()}
    for k in ("liver_art", "gut"): Q[k] *= sj.get("hepq_scale", 1)   # splanchnic/hepatic flow (cirrhosis, heart failure)
    Q["rest"] = CO - sum(Q.values())
    K = {t: kp(t, P, d["fup"]) for t in TIS + ["lung"]}
    if d.get("vss"):   # calibrate tissue partitioning so the model's Vss matches the supplied (measured / ADMET) Vss
        vt = sum(V[t] * K[t] for t in K); f = clamp((d["vss"] * bw - V["blood"]) / vt, .1, 10); K = {t: k * f for t, k in K.items()}
    if d["efflux_brain"] < 1: K["brain"] *= d["efflux_brain"]
    clint = (d["clint_hep"] * HEPC[species] if d.get("clint_hep") else d["clint"] * s["mppgl"]) * s["liver_g"] * bw * 1e-6 * 60 * sj.get("clint_scale", 1)
    clr = (d["clr"] or d["fup"] * s["gfr"] * bw) * sj.get("gfr_scale", 1)
    return dict(s=s, d=d, src=src, CO=CO, V=V, Q=Q, K=K, clint=clint, clr=clr, bw=bw)

def _deriv(y, m):
    V, Q, K, d, CO = m["V"], m["Q"], m["K"], m["d"], m["CO"]; ca = y[IX["art"]] / V["art"]; cven = y[IX["ven"]] / V["ven"]
    cv = {t: y[IX[t]] / V[t] / K[t] for t in TIS}; cvl = y[IX["lung"]] / V["lung"] / K["lung"]; out = [0.0] * len(ST)
    lum = y[0]; out[0] = -d["ka"] * lum; out[IX["fec"]] = (1 - d["fa"]) * d["ka"] * lum
    out[IX["gut"]] = d["fa"] * d["ka"] * lum + Q["gut"] * (ca - cv["gut"])
    qh = Q["liver_art"] + Q["gut"]; ml = m["clint"] * d["fup"] * cv["liver"]
    out[IX["liver"]] = Q["liver_art"] * ca + Q["gut"] * cv["gut"] - qh * cv["liver"] - ml; out[IX["cl_liver"]] = ml
    mk = m["clr"] * cv["kidney"]; out[IX["kidney"]] = Q["kidney"] * (ca - cv["kidney"]) - mk; out[IX["cl_kidney"]] = mk
    for t in ("heart", "brain", "muscle", "fat", "skin", "rest"): out[IX[t]] = Q[t] * (ca - cv[t])
    out[IX["lung"]] = CO * (cven - cvl); out[IX["art"]] = CO * (cvl - ca)
    out[IX["ven"]] = qh * cv["liver"] + sum(Q[t] * cv[t] for t in ("kidney", "heart", "brain", "muscle", "fat", "skin", "rest")) - CO * cven
    return out

def _matmul(A, B):
    Bt = list(zip(*B)); return [[sum(x * y for x, y in zip(r, c)) for c in Bt] for r in A]

def _expm(A, h):
    """matrix exponential exp(A*h) by scaling-and-squaring + Taylor series (exact for the linear PBPK system, no stiffness issues)"""
    n = len(A); B = [[a * h for a in r] for r in A]; nrm = max(sum(abs(a) for a in r) for r in B)
    sq = max(0, int(math.ceil(math.log2(nrm))) + 1) if nrm > .5 else 0; B = [[a / 2 ** sq for a in r] for r in B]
    I = [[1.0 if i == j else 0.0 for j in range(n)] for i in range(n)]; R = [r[:] for r in I]; T = [r[:] for r in I]
    for k in range(1, 20):
        T = _matmul(T, B); T = [[a / k for a in r] for r in T]; R = [[a + b for a, b in zip(r, t)] for r, t in zip(R, T)]
    for _ in range(sq): R = _matmul(R, R)
    return R

def simulate(p, species="human", dose_mg=100.0, hours=48.0, route="po", every=.1, ovr=None, subj=None, **_):
    """returns dict: t, cp (plasma mg/L), tissue conc (mg/L), PK metrics, mass balance, parameters + their source"""
    m = build(p, species, ovr, subj); n_st = len(ST); A = [[0.0] * n_st for _ in range(n_st)]
    for j in range(n_st):
        e = [0.0] * n_st; e[j] = 1.0; col = _deriv(e, m)
        for i in range(n_st): A[i][j] = col[i]
    Phi = _expm(A, every); y = [0.0] * n_st; y[IX["lumen" if route == "po" else "ven"]] = dose_mg; V = m["V"]
    steps = int(round(hours / every)); ts, cp, tis = [], [], {t: [] for t in TIS + ["lung"]}
    for i in range(steps + 1):
        ts.append(i * every); cp.append(max(0.0, y[IX["ven"]] / V["ven"]))
        for t in tis: tis[t].append(max(0.0, y[IX[t]] / V[t]))
        if i < steps: y = [sum(a * b for a, b in zip(r, y)) for r in Phi]
    ex = every; auc = sum((cp[i] + cp[i + 1]) / 2 * ex for i in range(len(cp) - 1)); cmax = max(cp); tmax = ts[cp.index(cmax)]
    late = [(t, math.log(c)) for t, c in zip(ts, cp) if t >= hours / 2 and c > 1e-12]; kz = None
    if len(late) > 5:
        mt = sum(a for a, _ in late) / len(late); mc = sum(b for _, b in late) / len(late)
        sl = sum((a - mt) * (b - mc) for a, b in late) / max(1e-12, sum((a - mt) ** 2 for a, _ in late)); kz = -sl if sl < 0 else None
    body = sum(y[IX[t]] for t in ST[:13]); mb = body + y[IX["cl_liver"]] + y[IX["cl_kidney"]] + y[IX["fec"]]
    return {"species": species, "route": route, "dose_mg": dose_mg, "t": ts, "cp": cp, "tissue": tis, "cmax": cmax, "tmax": tmax, "auc": auc,
            "kz": kz, "t_half": math.log(2) / kz if kz else None, "mass_balance": mb / dose_mg, "frac_liver_metab": y[IX["cl_liver"]] / dose_mg,
            "frac_renal": y[IX["cl_kidney"]] / dose_mg, "frac_unabsorbed": y[IX["fec"]] / dose_mg, "params": m["d"], "param_source": m["src"], "vss_L": m["V"]["blood"] + sum(m["V"][t] * m["K"][t] for t in m["K"]),
            "clint_Lh": m["clint"], "clr_Lh": m["clr"], "CO_Lh": m["CO"], "Kp": m["K"]}

def bioavailability(p, species="human", dose_mg=100.0, hours=48.0, ovr=None, subj=None):
    po = simulate(p, species, dose_mg, hours, "po", ovr=ovr, subj=subj); iv = simulate(p, species, dose_mg, hours, "iv", ovr=ovr, subj=subj)
    return po, iv, po["auc"] / iv["auc"] if iv["auc"] else 0.0

def steady_state(res, tau=24.0):
    """repeat-dose steady-state plasma profile over one interval by superposition (PK is linear)"""
    t, c, kz = res["t"], res["cp"], res["kz"] or 0.0; dtk = t[1] - t[0]; T = t[-1]; out = []
    for i in range(int(tau / dtk)):
        s, n = i * dtk, 0; tot = 0.0
        while True:
            x = s + n * tau
            if x <= T: j = int(round(x / dtk)); tot += c[min(j, len(c) - 1)]
            elif kz: tot += c[-1] * math.exp(-kz * (x - T))
            else: break
            if x > T + 8 / max(kz, 1e-9) or n > 400: break
            n += 1
        out.append(tot)
    return out

def dose_for_target(p, species, ic50_nM, mw, ovr=None, target_cover=1.0, tau=24.0, ref_mg=100.0, subj=None):
    """dose (mg per interval) so the trough FREE plasma conc >= target_cover x IC50 at steady state (linear PK scaling)"""
    r = simulate(p, species, ref_mg, 48.0, "po", ovr=ovr, subj=subj); fup = r["params"]["fup"]
    css = steady_state(r, tau); free_trough_nM = min(css) * fup / mw * 1e6
    if free_trough_nM <= 0: return None, r
    return ref_mg * (ic50_nM * target_cover) / free_trough_nM, r
