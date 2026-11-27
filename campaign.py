#!/usr/bin/env python3
"""campaign.py - run MedScope across EVERY disease at once and produce one master ranked table.
  python campaign.py [--limit 200000] [--evolve] [--keys a,b,c] [--top 5]                               (MedScope: campaign ...)
For each disease with measured target data: (1) screen the ChEMBL collection (screen.py), (2) optionally evolve novel molecules (evolve.py),
(3) mechanism analysis from disease genetics (reverse_bioengineering) so targets for diseases WITHOUT ChEMBL data are listed too.
Output: medscope_out/campaign.csv. Every row is a candidate TO TEST, not a drug."""
import csv, os, sys
import state

def run(keys=None, limit=None, top=5, do_evolve=False):
    import diseases, screen
    keys = keys or list(diseases.D); res = screen.run(keys, limit, top, show=True); rows = []
    for k, lst in res.items():
        for r in lst[:top]: rows.append(dict(disease=k, source="ChEMBL screen", id=r["id"], smiles=r["smiles"], pic50=round(r["pic50"], 2), tox=r["tox"] and round(r["tox"], 2), novel_for_target=r["measured"] is None))
    if do_evolve:
        import evolve
        for k in res:
            try:
                e = evolve.evolve(k, gens=5, pop=120, seeds=15, top=top, show=False)
                for r in e["top"]: rows.append(dict(disease=k, source="evolved (novel)", id="", smiles=r["smiles"], pic50=round(r["pic50"], 2), tox=None, novel_for_target=True))
            except Exception as ex: print(f"  evolve {k}: {ex}")
    os.makedirs(screen.OUT, exist_ok=True); p = os.path.join(screen.OUT, "campaign.csv")
    if rows:
        with open(p, "w", newline="") as f: w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    skipped = [k for k in keys if k not in res]
    print(f"\nCAMPAIGN: {len(res)} diseases screened, {len(rows)} candidates -> {p}")
    if skipped: print("  no measured target data (run `onboard`/docking path instead): " + ", ".join(skipped))
    state.record("campaign_diseases", len(res)); return rows

def cli(a):
    fl = {a[i]: a[i + 1] for i in range(len(a) - 1) if a[i].startswith("--")}
    return run(fl["--keys"].split(",") if "--keys" in fl else None, int(fl["--limit"]) if "--limit" in fl else None, int(fl.get("--top", 5)), "--evolve" in a)

if __name__ == "__main__": cli(sys.argv[1:])
