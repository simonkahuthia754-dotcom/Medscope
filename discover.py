#!/usr/bin/env python3
"""discover.py - per-disease candidate discovery:  python discover.py <disease> [--nodock] [--n 60]
measured target data -> QSAR potency (with error bars) -> real structures -> ADMET-AI -> Vina docking (trust-gated) -> human PBPK +
virtual population + human target expression -> composite score (QED 25 + potency/docking 30 + bioavailability 15 + Lipinski 10 + tox 20).
Everything is labelled with its confidence. Output = candidates worth testing in a REAL lab (see lab.py to close the loop), not cures."""
import io, math, os, sys
import chemistry, diseases, lab, pbpk, potency, sources, humanbio
from chemicals import exposure
clamp = lambda x, a, b: max(a, min(b, x))
TOXINS = ["terfenadine", "astemizole", "cisapride", "troglitazone", "bromfenac", "trovafloxacin", "sertindole"]  # withdrawn/black-boxed drugs (hERG/liver)

class Tee:
    def __init__(self, *s): self.s = s
    def write(self, x): [f.write(x) for f in self.s]
    def flush(self): [f.flush() for f in self.s]

def composite(c):
    """transparent 0-100 ranking heuristic (NOT a probability of success)"""
    q = 25 * c["qed"]; pot = clamp((c["pic50"] - 5) / 4, 0, 1)
    if c.get("dock_pct") is not None and c.get("dock_trust"): pot = .6 * pot + .4 * c["dock_pct"] / 100
    p = 30 * pot; f = 15 * clamp(c["F"] / .8, 0, 1); notes = []
    if c["F"] < .2: f = 0; notes.append("F<0.2: dropped")
    lip = 10 if c["lipinski_viol"] == 0 else 5 if c["lipinski_viol"] == 1 else 0
    if c["logp"] > 5 or c["mw"] > 600 or c["tpsa"] > 140: lip = 0; notes.append("outside property window (logP>5 / MW>600 / TPSA>140)")
    t = 20.0; A = c.get("admet") or {}
    if (A.get("hERG") or 0) >= .5 or c.get("herg_flag"): t -= 15; notes.append("hERG risk -15")
    if (A.get("DILI") or 0) >= .5: t -= 15; notes.append("DILI risk -15")
    if (A.get("AMES") or 0) >= .5: t -= 10; notes.append("Ames risk -10")
    t -= 3 * min(3, len(c["alerts"]["brenk"])); t = max(0, t); tot = q + p + f + lip + t
    if c.get("toxin_sim", 0) > .8: tot = 0; notes.append(f"BLOCKED: {c['toxin_sim']:.2f} similar to a withdrawn toxic drug")
    if c.get("dose") and c["dose"] > 1000: tot *= .5; notes.append("dose >1000 mg/day: halved")
    if c.get("near_fail"): tot *= .8; notes.append("close analog of a molecule that FAILED in your lab: x0.8")
    return round(min(100, tot), 1), {"QED": round(q, 1), "potency": round(p, 1), "F": round(f, 1), "Lipinski": lip, "tox": round(t, 1)}, notes

def run(key, top=5, n=60, pop=100, dock_on=True):
    D = diseases.get(key); chemistry.need(); print(f"\n=== {D['key']}: {D['name']} | target: {D['target']} ===")
    if D.get("note"): print("NOTE:", D["note"])
    cid, nm = diseases.resolve_chembl(D); print(f"ChEMBL target {cid} = {nm}   <- check this is the protein you mean")
    ov = lab.training_overlay(D["key"]); data = potency.fetch(cid, overlay=ov); model = potency.Model(data); print("QSAR:", model.summary()); __import__("state").record("qsar_auc", model.auc) if model.auc else None; __import__("state").record("qsar_rmse", model.rmse) if model.rmse else None
    if ov: print(f"  includes {len(ov)} measurements from YOUR lab")
    if model.n < 30: print("Too little measured data for this target: potency cannot be predicted reliably. Measure leads in a lab (lab.py) or use docking only."); return []
    diseases.check(D["key"], model=model)
    leads = {}
    for s, y in sorted(data.items(), key=lambda kv: -kv[1]):
        m = chemistry.mol(s)
        if m and 200 <= chemistry.describe(m)["mw"] <= 650: leads[f"lead{len(leads) + 1}"] = {"smiles": s, "ic50_nM": 10 ** (9 - y), "data": "measured"}
        if len(leads) >= 12: break
    cs = chemistry.candidates(leads, n=n)
    for c in cs:
        if c["known"]: c["pic50"] = 9 - math.log10(c["ic50_nM"]); c["pot"] = {"confidence": "MEASURED", "sd": .3, "ic50_nM": c["ic50_nM"], "ic50_lo": c["ic50_nM"] / 2, "ic50_hi": c["ic50_nM"] * 2, "maxsim": 1.0, "in_domain": True}
        else:
            c["pot"] = model.predict(c["smiles"]); c["pic50"] = c["pot"]["pic50"] if c["pot"] else None
    cs = [c for c in cs if c.get("pot") and c["pot"]["in_domain"]]; __import__("state").record("gen_ok", bool(cs)); print(f"Candidates inside the QSAR applicability domain: {len(cs)}")
    sl = sorted(cs, key=lambda c: -(c["pic50"] - c["pot"]["sd"]))[:top * 3]
    if not sl: print("No candidates in domain."); return []
    try:
        import admet; pr = admet.predict([c["smiles"] for c in sl]) if admet.HAVE else {}
    except Exception as e: pr = {}; print("! ADMET-AI failed:", e)
    print("ADMET:", "ADMET-AI predictions used" if pr else "ADMET-AI not available -> rule-based parameters (LOW confidence; run validate.py to see the error)")
    tx = {}
    for t in TOXINS:
        try: tx[t] = chemistry.fp(chemistry.mol(sources.pubchem(t)["smiles"]))
        except Exception: pass
    fails = [chemistry.fp(chemistry.mol(s)) for s in lab.failed(D["key"]) if chemistry.mol(s)]
    from rdkit import DataStructs
    T = cal = None; rmsd = None
    if dock_on and D.get("pdb"):
        try:
            import docking; T = docking.prepare(D["pdb"]); r = docking.redock_rmsd(T); rmsd = r[1] if r else None
            cal = docking.calibrate(T, data); trust = rmsd is not None and rmsd <= 2.0 and cal["rho"] is not None and cal["rho"] >= .3
            print(f"Docking {D['pdb']} (ligand {T['lig']}): redock RMSD {rmsd and round(rmsd, 2)} A | rho(Vina, measured pIC50) {cal['rho'] and round(cal['rho'], 2)} -> "
                  + ("TRUSTED as a cross-check" if trust else "NOT trusted (reported, but excluded from the score)"))
            T["trust"] = trust; st = __import__('state'); st.record('redock_rmsd', round(rmsd, 2)) if rmsd is not None else None; st.record('dock_rho', round(cal['rho'], 2)) if cal['rho'] is not None else None
        except Exception as e: print("Docking unavailable:", e); T = None
    for c in sl:
        a = pr.get(c["smiles"], {}); c["admet"] = {k: admet.get(a, k) for k in ("hERG", "DILI", "AMES")} if a else {}
        c["ovr"] = admet.overrides(a) if a else None; c["flags"] = admet.flags(a) if a else []; fp = chemistry.fp(chemistry.mol(c["smiles"]))
        c["toxin_sim"] = max([DataStructs.TanimotoSimilarity(fp, t) for t in tx.values()] or [0]); c["near_fail"] = any(DataStructs.TanimotoSimilarity(fp, f) >= .7 for f in fails)
        p = {"logp": c["logp"], "mw": c["mw"], "tpsa": c["tpsa"]}; po, iv, c["F"] = pbpk.bioavailability(p, "human", 100.0, 48.0, ovr=c["ovr"]); c["po"] = po
        c["dose"], _ = pbpk.dose_for_target(p, "human", c["pot"]["ic50_nM"], c["mw"], ovr=c["ovr"]); c["dose_hi"], _ = pbpk.dose_for_target(p, "human", c["pot"]["ic50_hi"], c["mw"], ovr=c["ovr"])
        c["rat"] = pbpk.bioavailability(p, "rat", 5 * .25, 24.0, ovr=c["ovr"]); c["exp"] = exposure(po)
        if T:
            d = docking.dock_consistent(T, c["smiles"])
            if d: c["dock"] = d; c["dock_pct"] = docking.percentile(cal, d[0]); c["dock_trust"] = T["trust"] and d[2]
        c["score"], c["parts"], c["notes"] = composite(c)
    sl.sort(key=lambda c: -c["score"]); lv, safety = ({}, [])
    if D.get("symbol"):
        try: lv, safety = humanbio.target_expression(D["symbol"])
        except Exception as e: print("! target expression lookup failed:", e)
    print(f"\nTOP {top} (composite: QED 25 + potency/docking 30 + F 15 + Lipinski 10 + tox 20; a ranking aid, NOT a probability of success)")
    for c in sl[:top]:
        po = c["po"]; P = c["pot"]; kind = "KNOWN ACTIVE (measured)" if c["known"] else f"analog of {c['parent']}"
        print(f"\n>> {c['name']} [{kind}]  score {c['score']}  {c['parts']}\n   {c['smiles']}")
        print(f"   potency: pIC50 {c['pic50']:.1f} -> IC50 {P['ic50_nM']:.3g} nM (plausible range {P['ic50_lo']:.3g}-{P['ic50_hi']:.3g}) | confidence {P['confidence']} (nearest known active similarity {P['maxsim']:.2f})")
        if c.get("dock"): print(f"   docking: {c['dock'][0]:.1f} +/- {c['dock'][1]:.1f} kcal/mol ({'consistent' if c['dock'][2] else 'INCONSISTENT'}); better than {c['dock_pct']:.0f}% of known actives; {'used in score' if c['dock_trust'] else 'not used in score'}")
        th = po["t_half"]; print(f"   human PBPK 100 mg po: Cmax {po['cmax']:.2f} mg/L, t1/2 {'n/a' if not th else f'{th:.1f}'} h, F {c['F']:.2f} | params: " + ", ".join(f"{k}={v}" for k, v in po["param_source"].items() if k in ("fup", "clint_hep", "vss")))
        if c["dose"]:
            print(f"   dose to hold free trough >= IC50: {c['dose']:.0f} mg/day (pessimistic IC50: {c['dose_hi']:.0f} mg/day)")
            if c["dose"] <= 1000:
                cyp = (D.get("cyp") if D.get("cyp") in humanbio.PHENO else None); r = humanbio.population_pk({"logp": c["logp"], "mw": c["mw"], "tpsa": c["tpsa"]}, c["dose"], ovr=c["ovr"],
                    subs=humanbio.population(pop, cyp=cyp), ic50_nM=P["ic50_nM"], mw=c["mw"])
                print(f"   virtual population (n={r['n']}): effective trough in {r['pct_covered']:.0f}% of subjects; t1/2 {r['t_half']['p5']:.0f}-{r['t_half']['p95']:.0f} h (5th-95th pct)"
                      + (" | by " + cyp + " phenotype: " + ", ".join(f"{k} {v[1]:.0f}% (n={v[0]})" for k, v in r["by_phenotype"].items()) if cyp else ""))
        rp = c["rat"][0]; print(f"   rat 5 mg/kg po: Cmax {rp['cmax']:.2f} mg/L, F {c['rat'][2]:.2f} (human-derived parameters scaled to rat: approximate)")
        print("   organ:plasma AUC: " + ", ".join(f"{o} {v[1]:.1f}" for o, v in c["exp"].items()))
        for f in c["flags"] + c["notes"] + humanbio.exposure_flags(c["exp"], lv, D.get("symbol") or ""): print("   ! " + f)
        if c["exp"]["brain"][1] > .5: print("   ! CNS exposure not needed for a peripheral disease (brain:plasma %.1f)" % c["exp"]["brain"][1])
        if c["alerts"]["brenk"]: print("   ! structural alerts: " + ", ".join(c["alerts"]["brenk"][:3]))
    if safety: print("\nHuman target safety liabilities on record (Open Targets): " + ", ".join(safety[:8]))
    try:
        for c in sl[:top]: chemistry.svg(c["smiles"], os.path.join("medscope_out", f"{D['key']}_{c['name']}.svg"))
    except Exception: pass
    print("\nNEXT: test the top molecules in a lab, then log results:  lab add " + D["key"] + " <smiles> <ic50_nM> [cc50_nM] [predicted_pic50]\n"
          "These are hypotheses ranked by evidence quality, not cures. Animal/human work needs ethics approval and qualified partners.\n"); return sl[:top]

def main():
    a = [x for x in sys.argv[1:] if not x.startswith("--")]; os.makedirs("medscope_out", exist_ok=True)
    if not a: diseases.show(); return
    f = open(os.path.join("medscope_out", a[0] + "_report.txt"), "w"); old = sys.stdout; sys.stdout = Tee(old, f)
    try: run(a[0], dock_on="--nodock" not in sys.argv, n=int(sys.argv[sys.argv.index("--n") + 1]) if "--n" in sys.argv else 60)
    finally: sys.stdout = old; f.close()

if __name__ == "__main__": main()
