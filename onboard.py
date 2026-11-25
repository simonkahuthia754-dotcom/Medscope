"""onboard.py - add ANY disease (or a brand-new pathogen target) automatically; no hand-curation.
  onboard "alzheimer disease"              -> Open Targets: top disease-linked targets -> ChEMBL target ids -> PDB structures -> new keys
  onboard --uniprot P12345 --name "Novel virus protease"    -> for a NEW disease/pathogen once its protein sequence/UniProt entry exists
New diseases appear as ~/.medscope/diseases.json entries usable by check/discover/repurpose. Auto entries have an empty GOOD list
(add known drugs to enable per-disease validation). Data for brand-new targets will be sparse: the QSAR will say so (low n, out-of-domain)."""
import json, os, re, sys
import humanbio, sources, state

DF = os.path.join(os.path.expanduser("~"), ".medscope", "diseases.json")
slug = lambda s: re.sub(r"[^a-z0-9]+", "-", s.lower()).strip("-")

def targets_for(disease, n=15):
    h = humanbio._gql('query($q:String!){search(queryString:$q,entityNames:["disease"],page:{index:0,size:1}){hits{id name}}}', {"q": disease})["data"]["search"]["hits"]
    if not h: raise RuntimeError(f"'{disease}' not found in Open Targets; try a synonym")
    q = ("query($id:String!){disease(efoId:$id){name associatedTargets(page:{index:0,size:%d}){rows{score target{id approvedSymbol approvedName "
         "proteinIds{id source} tractability{modality value}}}}}}") % n
    rows = humanbio._gql(q, {"id": h[0]["id"]})["data"]["disease"]["associatedTargets"]["rows"]; out = []
    for r in rows:
        t = r["target"]; up = [p["id"] for p in t.get("proteinIds") or [] if p.get("source") == "uniprot_swissprot"]
        out.append({"symbol": t["approvedSymbol"], "name": t["approvedName"], "uniprot": up[0] if up else None, "score": round(r["score"], 2),
                    "sm": any(x.get("modality") == "SM" and x.get("value") for x in t.get("tractability") or [])})
    return h[0]["name"], out

def chembl_for(uniprot):
    for t in sources.J(f"https://www.ebi.ac.uk/chembl/api/data/target.json?target_components__accession={uniprot}&limit=5")["targets"]:
        if t.get("target_type") == "SINGLE PROTEIN": return t["target_chembl_id"], t["pref_name"]
    return None, None

def pdbs_for(uniprot, rows=5):
    q = {"query": {"type": "terminal", "service": "text", "parameters": {"attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
         "operator": "exact_match", "value": uniprot}}, "return_type": "entry", "request_options": {"paginate": {"start": 0, "rows": rows}, "sort": [{"sort_by": "rcsb_entry_info.resolution_combined", "direction": "asc"}]}}
    try: return [x["identifier"] for x in sources.J("https://search.rcsb.org/rcsbsearch/v2/query", json.dumps(q), {"Content-Type": "application/json"}).get("result_set", [])]
    except Exception: return []

def _save(entries):
    try: cur = json.load(open(DF))
    except Exception: cur = {}
    cur.update(entries); os.makedirs(os.path.dirname(DF), exist_ok=True); json.dump(cur, open(DF, "w"), indent=1)

def entry(name, target, symbol, uniprot, note=""):
    cid, pref = chembl_for(uniprot) if uniprot else (None, None); pd = pdbs_for(uniprot) if uniprot else []
    return {"name": name, "target": pref or target, "chembl": cid, "pdb": pd[0] if pd else None, "pdbs": pd, "pdb_ok": False, "symbol": symbol, "cyp": None, "good": [],
            "note": "AUTO-ONBOARDED: GOOD list empty (add known drugs to enable validation); structure ids unverified (redocking gate checks). " + note}

def onboard(disease, top=3):
    dn, ts = targets_for(disease); print(f"{dn}: top disease-linked targets (Open Targets association score; SM = small-molecule tractable)")
    for t in ts: print(f"  {t['symbol']:<10}{t['score']:<6}{'SM ' if t['sm'] else '-- '}{t['uniprot'] or '(no UniProt)':<12}{t['name']}")
    pick = [t for t in ts if t["sm"] and t["uniprot"]][:top] or [t for t in ts if t["uniprot"]][:top]; new = {}
    for t in pick:
        e = entry(dn, t["name"], t["symbol"], t["uniprot"]); k = f"{slug(disease)}-{t['symbol'].lower()}"; new[k] = e
        print(f"  -> {k}: ChEMBL {e['chembl'] or 'NONE (no measured data: docking/lab only)'} | PDB {', '.join(e['pdbs'][:3]) or 'none'}")
    _save(new); state.record("onboard_ok", True); print(f"Saved {len(new)} entries to {DF}. Next: discover {next(iter(new))}   (or repurpose {next(iter(new))})"); return new

def from_uniprot(name, uniprot):
    k = slug(name); e = entry(name, name, None, uniprot, "new pathogen/target"); _save({k: e}); state.record("onboard_ok", True)
    print(f"Saved '{k}': ChEMBL {e['chembl'] or 'NONE'} | PDB {', '.join(e['pdbs'][:3]) or 'none'} | until measured data exists, rely on docking + lab loop"); return {k: e}

def cli(a):
    if "--uniprot" in a:
        u = a[a.index("--uniprot") + 1]; n = a[a.index("--name") + 1] if "--name" in a else u; return from_uniprot(n, u)
    if not a: print(__doc__); return
    return onboard(" ".join(a))

if __name__ == "__main__": cli(sys.argv[1:])
