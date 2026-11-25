#!/usr/bin/env python3
"""validate.py - benchmark MedScope's human PBPK against approximate literature human PK of well-known drugs.
Reports fold errors so you know the REAL accuracy of the simulator. Literature values are approximate (verify before citing).
Descriptors are approximate experimental values; with ADMET-AI installed, also pass `--admet` to use model-derived parameters."""
import math, sys
import pbpk

# name: (MW, logP, TPSA, t_half_h, F_oral, SMILES)  -- approximate literature values
REF = {
    "caffeine": (194.2, -0.07, 58.4, 5.0, 1.00, "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    "theophylline": (180.2, -0.02, 69.3, 8.0, 1.00, "Cn1c2nc[nH]c2c(=O)n(C)c1=O"),
    "propranolol": (259.3, 3.48, 41.5, 4.0, 0.26, "CC(C)NCC(O)COc1cccc2ccccc12"),
    "diazepam": (284.7, 2.82, 32.7, 43.0, 0.93, "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
    "warfarin": (308.3, 2.70, 63.6, 37.0, 0.99, "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
    "ibuprofen": (206.3, 3.97, 37.3, 2.0, 0.90, "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    "metformin": (129.2, -1.40, 91.5, 6.2, 0.55, "CN(C)C(=N)NC(=N)N"),
    "atenolol": (266.3, 0.16, 84.6, 6.5, 0.55, "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"),
    "zidovudine": (267.2, 0.05, 125.8, 1.1, 0.64, "CC1=CN(C(=O)NC1=O)C2CC(C(O2)CO)N=[N+]=[N-]"),
    "ciprofloxacin": (331.3, 0.28, 74.6, 4.0, 0.70, "OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O"),
    "linezolid": (337.4, 0.55, 71.1, 5.0, 1.00, "CC(=O)NC[C@H]1CN(c2ccc(N3CCOCC3)c(F)c2)C(=O)O1"),
    "imatinib": (493.6, 3.50, 86.3, 18.0, 0.98, "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
    "erlotinib": (393.4, 3.30, 74.7, 36.0, 0.60, "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"),
    "chloroquine": (319.9, 4.60, 28.2, 1000.0, 0.80, "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12"),
}
fold = lambda pred, obs: max(pred / obs, obs / pred)

def run(use_admet=False):
    ad = None
    if use_admet:
        import admet
        if not admet.HAVE: print("! admet-ai not installed (pip install admet-ai); using rule-based parameters"); use_admet = False
        else: ad = admet.predict([v[5] for v in REF.values()])
    rows = []
    for nm, (mw, lp, tp, th, F, smi) in REF.items():
        ovr = admet.overrides(ad[smi]) if ad else None
        po, iv, f = pbpk.bioavailability({"logp": lp, "mw": mw, "tpsa": tp}, "human", 100.0, max(96.0, 6 * th), ovr=ovr)
        rows.append((nm, th, po["t_half"] or float("nan"), F, f))
    print(f"{'drug':<14}{'t1/2 obs':>9}{'pred':>9}{'fold':>7}{'F obs':>8}{'pred':>7}")
    ft, ff = [], []
    for nm, th, pt, F, pf in rows:
        a = fold(pt, th) if pt == pt else float("inf"); b = fold(max(pf, .01), F); ft.append(a); ff.append(b)
        print(f"{nm:<14}{th:>9.1f}{pt:>9.1f}{a:>7.1f}{F:>8.2f}{pf:>7.2f}")
    gm = lambda x: math.exp(sum(math.log(v) for v in x if v != float('inf')) / max(1, len(x)))
    w2 = lambda x: 100 * sum(v <= 2 for v in x) / len(x); w3 = lambda x: 100 * sum(v <= 3 for v in x) / len(x)
    print(f"\nHalf-life : within 2-fold {w2(ft):.0f}% | within 3-fold {w3(ft):.0f}% | geometric-mean fold error {gm(ft):.2f}")
    print(f"Bioavail. : within 2-fold {w2(ff):.0f}% | within 3-fold {w3(ff):.0f}% | geometric-mean fold error {gm(ff):.2f}")
    print("Typical published PBPK 'good' benchmark: >=70% within 2-fold. Treat predictions for NEW molecules as less accurate than this.")
    return rows

if __name__ == "__main__":
    run("--admet" in sys.argv)
