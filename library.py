"""library.py - EVERY existing approved drug, and drug repurposing.
  repurpose <disease>   score all ChEMBL approved small-molecule drugs (max phase 4) against that disease's target with the measured-data QSAR.
Why it matters: an approved drug already has human safety + PK data, so a real hit can reach patients years sooner than a new molecule.
Scores are PREDICTIONS with error bars, to be confirmed by an assay. Any drug by name/SMILES: see ingredients.py (`drug <name>`)."""
import json, os, sources, state, potency, diseases

def approved(limit=6000):
    """-> [{name, chembl, smiles}] approved small molecules (cached for 12 h by sources.J)"""
    out, off = [], 0
    while off < limit:
        r = sources.J(f"https://www.ebi.ac.uk/chembl/api/data/molecule.json?max_phase=4&molecule_type=Small%20molecule&limit=1000&offset={off}")["molecules"]
        for m in r:
            s = (m.get("molecule_structures") or {}).get("canonical_smiles")
            if s: out.append({"name": m.get("pref_name") or m["molecule_chembl_id"], "chembl": m["molecule_chembl_id"], "smiles": potency.standardize(s)})
        if len(r) < 1000: break
        off += 1000
    state.record("library_n", len(out)); return out

def repurpose(key, top=15):
    d = diseases.get(key); cid, nm = diseases.resolve_chembl(d); print(f"[{d['key']}] target {cid} = {nm}")
    import lab
    data = potency.fetch(cid, overlay=lab.training_overlay(d["key"])); m = potency.Model(data); print("QSAR:", m.summary())
    if m.n < 30: print("Too little measured data on this target for a reliable screen."); return []
    lib = approved(); print(f"Screening {len(lib)} approved drugs...")
    res = []
    for x in lib:
        p = m.predict(x["smiles"])
        if p and p["in_domain"]: res.append((p["pic50"] - p["sd"], p, x))
    res.sort(key=lambda t: -t[0]); print(f"{len(res)} drugs inside the model's domain. TOP {top} by conservative potency (pIC50 - sd):")
    for _, p, x in res[:top]:
        print(f"  {x['name']:<26} IC50 ~{p['ic50_nM']:.3g} nM (range {p['ic50_lo']:.3g}-{p['ic50_hi']:.3g}) | confidence {p['confidence']} | similarity {p['maxsim']:.2f} | {x['chembl']}")
    print("\nAlready approved = human safety/PK data exist. These are PREDICTIONS: confirm with an assay before any claim. Check each drug's current indication and dose/toxicity limits.")
    return res[:top]
