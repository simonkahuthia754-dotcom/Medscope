#!/usr/bin/env python3
"""recreation_bioengineering.py - rebuild a disease as a mechanistic within-host MODEL and test it against PUBLISHED clinical behaviour.
  python recreation_bioengineering.py [hiv|tb|lung-cancer|all] [--calibrate] [--regimen ic50_nM=.. conc_nM=..]      (MedScope: recreate <disease>)
What 'fidelity' means here (no inflated claims):
  * CALIBRATED features: published values the parameters were fitted to (fit quality only - passing is expected).
  * HELD-OUT features: published values the model was NOT fitted to - the real test of whether the mechanism is understood.
  * PRIOR features: parameters fixed from the literature (reported, not counted as evidence).
Passing every published summary statistic does NOT make the model identical to the real disease: it shows the mechanism reproduces those
behaviours. Reference ranges are approximate literature values (verify before citing). Defensive/therapeutic modelling of public
epidemiology and clinical data only. NOTE: the model structure was revised once after the first held-out results were seen (second-phase rate, reservoir-size and CD4-decline constraints added), so held-out features are no longer fully blind; the strictest test is a NEW disease/cohort. Calibrated parameters are stored in recreation_params.json (re-fit with --calibrate)."""
import json, math, os, random, sys
import state

HERE = os.path.dirname(os.path.abspath(__file__)); PF = os.path.join(HERE, "recreation_params.json")

def rk4(f, y, t, dt):
    k1 = f(t, y); k2 = f(t + dt / 2, [a + dt / 2 * b for a, b in zip(y, k1)]); k3 = f(t + dt / 2, [a + dt / 2 * b for a, b in zip(y, k2)]); k4 = f(t + dt, [a + dt * b for a, b in zip(y, k3)])
    return [max(0.0, a + dt / 6 * (b + 2 * c + 2 * d + e)) for a, b, c, d, e in zip(y, k1, k2, k3, k4)]

# ---------------------------------------------------------------- HIV: target cells, latent reservoir, infected, long-lived infected, CTL, immune damage
HIV0 = dict(lam=10.0, dT=0.01, beta=2.7e-8, p=1.5e6, c=23.0, delta=0.7, fL=2e-3, fM=.03, deltaM=.05, aL=3e-5, dL=4.9e-4, kE=.01, sE=.5, bE=.35, KE=.5, dE=.05, Emax=100.0, kD=6e-4, rD=4e-5, rA=100.0, Vs=2.0, Vk=1e4, Ik=1e-3)
HIV_FREE = {"fL": (1e-5, 5e-3), "fM": (.005, .1), "beta": (5e-9, 2e-7), "p": (1e4, 1e8), "kE": (.001, .2), "bE": (.05, 2.0), "kD": (1e-4, 1e-2)}
BODY = 2.5e8     # blood-volume-equivalent: 5e6 uL x ~50 for lymphoid tissue (order-of-magnitude, for reservoir counts)

def hiv_f(P, eps):
    def f(t, y):
        T, L, I, M, E, D = y; V = P["p"] * (I + M) / P["c"]; e = eps(t) if callable(eps) else eps; inf = (1 - e) * P["beta"] * T * V
        return [P["lam"] * (1 - D) - P["dT"] * T - inf, P["fL"] * inf - (P["aL"] + P["dL"]) * L,
                (1 - P["fL"] - P["fM"]) * inf + P["aL"] * L - P["delta"] * I - P["kE"] * E * I, P["fM"] * inf - P["deltaM"] * M,
                P["sE"] + P["bE"] * (I + M) / (I + M + P["KE"]) * E * (1 - E / (P["Emax"] * (1 - D) ** 3 + 1.0)) - P["dE"] * E, P["kD"] * (.5 * V / (V + P["Vk"]) + .5 * (I + M) / (I + M + P["Ik"])) * (1 - D) - P["rD"] * D * (1 + P["rA"] * P["Vs"] / (V + P["Vs"]))]
    return f

def hiv_run(P, days, eps=0.0, y0=None, dt=.1, I0=1e-5, every=1.0):
    y = y0 or [P["lam"] / P["dT"], 0.0, I0, 0.0, P["sE"] / P["dE"], 0.0]; f = hiv_f(P, eps); out = []; t = 0.0; nxt = 0.0
    while t <= days + 1e-9:
        if t >= nxt - 1e-9:
            V = P["p"] * (y[2] + y[3]) / P["c"]; out.append(dict(t=t, V=V, CD4=y[0] + y[1] + y[2] + y[3], L=y[1], I=y[2], D=y[5], y=y[:])); nxt += every
        y = rk4(f, y, t, dt); t += dt
    return out

def log10(x): return math.log10(max(x, 1e-9))

def hiv_features(P, long=True):
    """simulate stages: acute -> chronic -> ART -> interruption; return {name: (value, lo, hi, kind, unit, source)}"""
    inf = hiv_run(P, 1830, 0.0, every=1.0); pk = max(inf[:120], key=lambda r: r["V"]); sp = inf[365]; y1 = sp["y"]
    art = hiv_run(P, 730, .99, y0=y1, every=1.0); v = [r["V"] for r in art]
    r7 = art[7]["V"]; t_half1 = 7 * math.log(2) / math.log(max(v[0], 1e-9) / max(r7, 1e-9)) if r7 < v[0] else float("inf")
    sup = next((r["t"] for r in art if r["V"] < 50), float("inf")); gain = art[365]["CD4"] - art[0]["CD4"]
    L = [r["L"] for r in art[200:]]; kL = -math.log(L[-1] / L[0]) / (art[-1]["t"] - art[200]["t"]) if L[0] > 0 and L[-1] > 0 else 0; t_half_L = math.log(2) / kL if kL > 0 else float("inf")
    stop = hiv_run(P, 90, 0.0, y0=art[-1]["y"], every=.5); reb = next((r["t"] for r in stop if r["V"] > 50), float("inf"))
    F = {"acute peak viremia": (log10(pk["V"]), 5.5, 8.0, "calibrated", "log10 copies/mL", "Fiebig/acute HIV cohorts"), "time to acute peak": (pk["t"], 14, 35, "calibrated", "days", "acute HIV cohorts"),
         "set-point viremia (1 y)": (log10(sp["V"]), 3.5, 5.0, "calibrated", "log10 copies/mL", "MACS/CASCADE set-point"), "CD4 at 1 y (untreated)": (sp["CD4"], 450, 900, "calibrated", "cells/uL", "seroconverter cohorts"),
         "first-phase decay on ART (t1/2)": (t_half1, .8, 2.5, "calibrated", "days", "Perelson 1996/97"),
         "time to viral suppression <50 on ART": (sup, 40, 180, "held-out", "days", "ART trials (typically 8-24 weeks)"), "CD4 gain in 1st year of ART": (gain, 80, 320, "held-out", "cells/uL", "ART cohorts"),
         "viral rebound after stopping ART": (reb, 10, 42, "held-out", "days", "analytic treatment interruptions"), "latent reservoir decay (t1/2)": (t_half_L, 1000, 1800, "prior", "days", "Siliciano 2003 (44 months)")}
    if long:
        un = hiv_run(P, 365 * 15, 0.0, dt=.25, every=30.0); c200 = next((r["t"] / 365 for r in un if r["CD4"] < 200), float("inf"))
        F["time to CD4 <200 untreated"] = (c200, 6, 13, "held-out", "years", "natural-history cohorts (median ~10 y)")
    F["latent reservoir size on ART"] = (art[200]["L"] * BODY, 1e5, 1e7, "calibrated", "cells (body)", "Finzi/Chun/Siliciano reservoir estimates")
    F["CD4 decline rate, untreated years 2-5"] = ((inf[1825]["CD4"] - inf[730]["CD4"]) / 3, -100, -40, "calibrated", "cells/uL/year", "natural-history cohorts (~50-80/yr)"); return F

# ---------------------------------------------------------------- Bacterial infection with persisters (TB / MRSA-like): replicating, persister, immune clearance
TB0 = dict(g=.35, kI=.05, s=.03, w=.004, B0=1e4, Bmax=1e9, kill_r=.9, kill_p=.06, Hill=1.5)
def tb_f(P, c_over_mic):
    def f(t, y):
        Br, Bp = y; N = Br + Bp; g = P["g"] * max(0.0, 1 - N / P["Bmax"]); e = (lambda x: x ** P["Hill"] / (1 + x ** P["Hill"]))(c_over_mic(t) if callable(c_over_mic) else c_over_mic)
        return [g * Br - P["kI"] * Br - P["s"] * Br + P["w"] * Bp - P["kill_r"] * e * Br, P["s"] * Br - P["w"] * Bp - P["kill_p"] * e * Bp - P["kI"] * .3 * Bp]
    return f
def tb_run(P, days, c=0.0, y0=None, dt=.1, every=1.0):
    y = y0 or [P["B0"], 0.0]; f = tb_f(P, c); out = []; t = 0.0; nxt = 0.0
    while t <= days + 1e-9:
        if t >= nxt - 1e-9: out.append(dict(t=t, N=y[0] + y[1], Br=y[0], Bp=y[1], y=y[:])); nxt += every
        y = rk4(f, y, t, dt); t += dt
    return out
def tb_features(P):
    un = tb_run(P, 200); c = 4.0; tx = tb_run(P, 420, c, y0=un[200]["y"]); N0 = tx[0]["N"]; eba = (log10(N0) - log10(tx[2]["N"])) / 2
    late = (log10(tx[30]["N"]) - log10(tx[120]["N"])) / 90; cure = next((r["t"] for r in tx if r["N"] < 1.0), float("inf")); st = tb_run(P, 120, 0.0, y0=tx[120]["y"])
    return {"bacterial load at presentation": (log10(un[200]["N"]), 6.5, 9.0, "calibrated", "log10 CFU", "cavitary TB sputum loads"), "early bactericidal activity (days 0-2)": (eba, .25, .8, "calibrated", "log10 CFU/day", "Jindani/Mitchison EBA"),
            "late sterilizing phase decline": (late, .02, .12, "held-out", "log10 CFU/day", "biphasic kill, Mitchison"), "time to sterilization on therapy": (cure / 30.4, 3, 9, "held-out", "months", "standard 6-month regimens"),
            "relapse after stopping therapy at 2 months": (log10(st[-1]["N"]) - log10(max(tx[60]["N"], 1e-9)), -1, 12, "held-out", "log10 change (growth = relapse)", "short-course trials: high relapse")}

# ---------------------------------------------------------------- Solid tumour with a pre-existing resistant clone (Gompertz + targeted therapy)
CA0 = dict(a=.0007, K=1e12, N0=1e9, r0=1e-7, kill=.03, fitR=.85, pT=.4)
def ca_run(P, days, drug_from=0, dt=1.0, every=7.0, dose=1.0, on_until=1e9):
    S, R = P["N0"] * (1 - P["r0"]), P["N0"] * P["r0"]; out = []; t = 0.0
    while t <= days:
        N = S + R
        if t % every < dt: out.append(dict(t=t, N=N, R=R / max(N, 1)))
        g = P["a"] * math.log(P["K"] / max(N, 1)); on = drug_from <= t < on_until
        S = max(S + dt * (g * S - (P["kill"] * dose * S if on else 0)), 0.0); R = max(R + dt * (g * P["fitR"] * R), 0.0); t += dt
    return out
def ca_features(P):
    dbl = math.log(2) / (P["a"] * math.log(P["K"] / P["N0"])); tx = ca_run(P, 2400, drug_from=0); N0 = tx[0]["N"]; nad = min(tx, key=lambda r: r["N"])
    resp = 1 - nad["N"] / N0; pr = next((r for r in tx if r["t"] > nad["t"] and r["N"] >= 1.73 * nad["N"]), None)   # RECIST progression: +20% diameter = 1.73x volume over nadir
    prog = pr["t"] / 30.4 if pr else float("inf")
    return {"tumour volume doubling time": (dbl, 60, 300, "calibrated", "days", "NSCLC volumetric doubling (median ~150 d)"), "best objective response (shrinkage)": (resp, .3, .98, "calibrated", "fraction", "EGFR-TKI ORR, RECIST >=30%"),
            "time to progression on targeted therapy": (prog, 7, 22, "calibrated", "months", "EGFR-TKI PFS ~9-19 months (osimertinib ~19)"),
            "resistant fraction of tumour at progression": (pr["R"] if pr else 0.0, .4, 1.0, "held-out", "fraction", "T790M/bypass resistance in ~50-60% of progressions (Yu 2013, Oxnard 2011)"),
            "tumour-doubling time with sensitive clone only vs resistant (fitness cost)": (P["fitR"], .5, 1.0, "prior", "relative growth", "resistant clones carry a fitness cost")}

TB_FREE = {"g": (.1, .6), "kill_r": (.3, 3.0), "kI": (.005, .15)}; CA_FREE = {"kill": (.004, .2), "r0": (1e-8, 1e-2), "a": (3e-4, 3e-3)}   # physiological bounds (a resistant clone is not fitter than the sensitive one)
MODELS = {"hiv": ("HIV infection", HIV0, hiv_features, HIV_FREE), "tb": ("Tuberculosis", TB0, tb_features, TB_FREE), "lung-cancer": ("EGFR lung cancer", CA0, ca_features, CA_FREE)}

def load(key):
    P = dict(MODELS[key][1])
    try: P.update(json.load(open(PF)).get(key, {}))
    except Exception: pass
    return P

def score(F):
    d = {"calibrated": [0, 0], "held-out": [0, 0], "prior": [0, 0]}
    for k, (v, lo, hi, kind, *_r) in F.items():
        if kind in d: d[kind][1] += 1; d[kind][0] += int(lo <= v <= hi)
    return d

def penalty(F):
    s = 0.0
    for k, (v, lo, hi, kind, *_r) in F.items():
        if kind == "calibrated":
            if not (v == v) or v in (float("inf"),): s += 10; continue
            a, b = (log10(lo), log10(hi)) if lo > 0 and hi / lo > 50 else (lo, hi); x = log10(v) if lo > 0 and hi / lo > 50 else v; w = max(b - a, 1e-9); z = (x - (a + b) / 2) / (w / 2); s += z * z * .2 + 5 * (max(0, a - x) + max(0, x - b)) / w   # pull to the CENTRE of the published range, not its edge
    return s

def calibrate(key, iters=700, seed=11):
    return _calibrate(key, iters, seed) if key == "hiv" else _multi(key, iters)

def _multi(key, iters):
    best = None
    for sd in (11, 12, 13, 14):
        P, res = _calibrate(key, iters, sd, save=False)
        if best is None or res < best[1]: best = (P, res)
    _store(key, best[0]); return best[0]

def _store(key, P):
    free = MODELS[key][3]
    try: allp = json.load(open(PF))
    except Exception: allp = {}
    allp[key] = {k: P[k] for k in free}; json.dump(allp, open(PF, "w"), indent=1); print(f"  calibrated {MODELS[key][0]} (best of 4 seeds)")

def _calibrate(key, iters=700, seed=11, save=True):
    """random-restart coordinate search on CALIBRATED features only (held-out features are never used in the objective)"""
    name, base, ff, free = MODELS[key]; P = load(key); rng = random.Random(seed); best = penalty(ff(P, False) if key == "hiv" else ff(P))
    for i in range(iters):
        if best <= 1e-3: break
        Q = dict(P); k = rng.choice(list(free)); lo, hi = free[k]; Q[k] = min(hi, max(lo, P[k] * math.exp(rng.gauss(0, .5 if i < iters / 2 else .15))))
        s = penalty(ff(Q, False) if key == "hiv" else ff(Q))
        if s <= best: P, best = Q, s
    if save: _store(key, P); print(f"  residual {best:.3f} (lower = closer to the centre of every published range)"); return P
    return P, best

def show(key, calib=False):
    name, base, ff, free = MODELS[key]
    if calib: calibrate(key)
    P = load(key); F = ff(P); print(f"\n=== {name}: recreated disease vs published behaviour ===")
    for k, v in F.items():
        if k.startswith("_"): continue
        val, lo, hi, kind, unit, src = v; print(f"  [{'PASS' if lo <= val <= hi else 'FAIL'}] {k:<44}{val:>10.3g} {unit:<22} published {lo:g}-{hi:g} ({kind}; {src})")
    sc = score(F); cal, held = sc["calibrated"], sc["held-out"]
    print(f"  calibrated features reproduced: {cal[0]}/{cal[1]} ({100 * cal[0] / max(1, cal[1]):.0f}%)  |  HELD-OUT (predictive) features reproduced: {held[0]}/{held[1]} ({100 * held[0] / max(1, held[1]):.0f}%)")
    return sc

def virtual_patients(key="hiv", n=40, seed=3, sd=.25):
    """between-patient variability: perturb parameters (lognormal) and see what fraction of virtual patients stay inside every published range"""
    rng = random.Random(seed); base = load(key); ok = 0; ff = MODELS[key][2]; nf = 0
    for _ in range(n):
        P = {k: (v * math.exp(rng.gauss(0, sd)) if isinstance(v, float) and k not in ("fL", "aL", "dL") else v) for k, v in base.items()}
        F = ff(P, False) if key == "hiv" else ff(P); r = [lo <= v <= hi for k, (v, lo, hi, kind, *_a) in F.items() if kind in ("calibrated", "held-out")]; ok += all(r); nf += sum(r) / len(r)
    return ok / n, nf / n

def regimen(ic50_nM, conc_nM, hill=2.0, lra=1.0, years=3, stop_after=2.0):
    """candidate -> HIV model: ART efficacy from free concentration (Emax/Hill), simulate suppression, optional 'shock' (latency reversal x lra), what remains,
    and whether anything approaches sterilizing cure (< ~1 latent cell body-wide)"""
    P = load("hiv"); eps = conc_nM ** hill / (conc_nM ** hill + ic50_nM ** hill); print(f"  candidate: IC50 {ic50_nM:g} nM, free conc {conc_nM:g} nM -> efficacy {eps:.4f}")
    inf = hiv_run(P, 365, 0.0); y = inf[-1]["y"]; Pt = dict(P, aL=P["aL"] * lra, kE=P["kE"] * (3 if lra > 1 else 1)); run = hiv_run(Pt, 365 * years, eps, y0=y, dt=.25, every=30.0)
    Lc = run[-1]["L"] * BODY; print(f"  after {years} y: viral load {run[-1]['V']:.3g} copies/mL, CD4 {run[-1]['CD4']:.0f}, latent reservoir ~{Lc:.3g} cells (start ~{inf[-1]['L'] * BODY:.3g})")
    t_cure = None
    if Lc > 0 and run[-1]["L"] < run[0]["L"]:
        k = -math.log(run[-1]["L"] / run[0]["L"]) / (run[-1]["t"] - run[0]["t"]); t_cure = math.log(max(Lc, 1.0)) / k / 365 if k > 0 else None
    print("  extrapolated time to <1 latent cell body-wide: " + (f"{t_cure:.0f} years" if t_cure else "never at this rate") + (" (latency-reversal x%g applied)" % lra if lra > 1 else " (ART alone: reservoir decays too slowly -> lifelong therapy)"))
    return dict(eps=eps, reservoir=Lc, years_to_cure=t_cure)

def cli(a):
    flags = [x for x in a if x.startswith("--")]; args = [x for x in a if not x.startswith("--") and "=" not in x]; kv = dict(x.split("=") for x in a if "=" in x)
    if "ic50_nM" in kv: return regimen(float(kv["ic50_nM"]), float(kv.get("conc_nM", 100)), lra=float(kv.get("lra", 1)))
    keys = list(MODELS) if not args or args[0] == "all" else [args[0]]
    for k in keys:
        if k not in MODELS: print(f"no mechanistic model for '{k}' yet. Available: {', '.join(MODELS)}. (Honest limit: each disease needs its own validated model.)"); continue
        sc = show(k, "--calibrate" in flags)
        if k == "hiv": fr, mean = virtual_patients("hiv"); print(f"  virtual patients (parameters varied +-25%): {100 * fr:.0f}% reproduce ALL published ranges; mean {100 * mean:.0f}% of features per patient")
        state.record(f"recreate_{k}", [sc["calibrated"][0], sc["calibrated"][1], sc["held-out"][0], sc["held-out"][1]])
    state.record("recreate_ok", True)

if __name__ == "__main__": cli(sys.argv[1:])
