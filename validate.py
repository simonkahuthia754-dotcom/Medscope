#!/usr/bin/env python3
"""validate.py - benchmark MedScope's human PBPK against approximate literature human PK (pkbench.py: DEV = tuning set, TEST = blind set never used for tuning).
  python validate.py [--admet]     reports fold errors; the TEST numbers are the honest estimate for NEW molecules. Literature values are approximate."""
import math, sys
import pbpk, pkbench

fold = lambda pred, obs: max(pred / obs, obs / pred)
gm = lambda x: math.exp(sum(math.log(v) for v in x) / max(1, len(x)))
w = lambda x, k: 100 * sum(v <= k for v in x) / len(x)

def evaluate(S, ad=None, verbose=True):
    import admet
    ft, ff = [], []
    if verbose: print(f"{'drug':<14}{'t1/2 obs':>9}{'pred':>9}{'fold':>7}{'F obs':>8}{'pred':>7}")
    for nm, (mw, lp, tp, th, F, ion, smi) in S.items():
        ovr = admet.overrides(ad[smi]) if ad else None; p = {"logp": lp, "mw": mw, "tpsa": tp, "ion": ion}
        po, iv, f = pbpk.bioavailability(p, "human", 100.0, min(max(96.0, 6 * th), 2400.0), ovr=ovr); t = po["t_half"] or 1e4
        a = fold(t, th); b = fold(max(f, .01), F); ft.append(a); ff.append(b)
        if verbose: print(f"{nm:<14}{th:>9.1f}{t:>9.1f}{a:>7.1f}{F:>8.2f}{f:>7.2f}")
    return dict(t2=w(ft, 2), t3=w(ft, 3), tg=gm(ft), f2=w(ff, 2), fg=gm(ff))

def run(use_admet=False):
    import state
    ad = None
    if use_admet:
        import admet
        if not admet.HAVE: print("! admet-ai not installed (pip install admet-ai); using rule-based parameters"); use_admet = False
        else: ad = admet.predict([v[6] for S in (pkbench.DEV, pkbench.TEST) for v in S.values()])
    tag = "admet" if use_admet else "rule"; res = {}
    for nm, S in (("DEV (tuning set)", pkbench.DEV), ("TEST (blind, not used for tuning)", pkbench.TEST)):
        print(f"\n== {nm} [{tag}] =="); r = evaluate(S, ad, verbose="--v" in sys.argv or nm.startswith("TEST")); res[nm[:3]] = r
        print(f"  half-life within 2-fold {r['t2']:.0f}% | 3-fold {r['t3']:.0f}% | GMFE {r['tg']:.2f}   ||   bioavailability within 2-fold {r['f2']:.0f}% | GMFE {r['fg']:.2f}")
    state.record(f"pk_dev_{tag}", round(res["DEV"]["t2"])); state.record(f"pk_test_{tag}", round(res["TES"]["t2"])); state.record(f"pk_test_f_{tag}", round(res["TES"]["f2"]))
    state.record("pk_2fold_admet" if use_admet else "pk_2fold_rule", round(res["TES"]["t2"]))
    print("\nPublished PBPK 'good' benchmark: >=70% within 2-fold. Use TEST numbers as the expected accuracy on new molecules.")
    return res

if __name__ == "__main__": run("--admet" in sys.argv)
