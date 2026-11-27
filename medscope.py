#!/usr/bin/env python3
"""MedScope - terminal research aggregator for any disease (HIV, cancer, malaria...).
Pulls live public data: PubMed (papers), ClinicalTrials.gov (trials), Open Targets (drugs/targets).
Stdlib only -> runs in QPython/Android. Research aid only, NOT medical advice."""
import json, os, sys, time, hashlib, textwrap, urllib.request, urllib.parse

H = os.path.join(os.path.expanduser("~"), ".medscope"); os.makedirs(H, exist_ok=True)
CF = os.path.join(H, "cache.json"); TTL = 6 * 3600
NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
CT = "https://clinicaltrials.gov/api/v2/studies"
OT = "https://api.platform.opentargets.org/api/v4/graphql"
DISC = "MedScope aggregates public research data. It is NOT medical advice and cannot confirm a cure."

def _cache():
    try: return json.load(open(CF))
    except Exception: return {}

def get(url, data=None, hdr=None):
    import sources
    try: return sources.J(url, data, hdr, ttl=TTL)
    except Exception as e: raise RuntimeError(f"network/API error: {e}")

def papers(q, n=10, recent=True):
    p = urllib.parse.urlencode({"db": "pubmed", "term": q, "retmode": "json", "retmax": n, "sort": "date" if recent else "relevance"})
    ids = get(NCBI + "esearch.fcgi?" + p)["esearchresult"]["idlist"]
    if not ids: return []
    s = get(NCBI + "esummary.fcgi?" + urllib.parse.urlencode({"db": "pubmed", "id": ",".join(ids), "retmode": "json"}))["result"]
    return [{"id": i, "title": s[i].get("title", ""), "journal": s[i].get("source", ""), "date": s[i].get("pubdate", ""),
             "authors": ", ".join(a["name"] for a in s[i].get("authors", [])[:3]), "url": f"https://pubmed.ncbi.nlm.nih.gov/{i}/"}
            for i in ids if i in s]

def trials(cond, n=10, status="RECRUITING"):
    p = {"query.cond": cond, "pageSize": n, "format": "json"}
    if status: p["filter.overallStatus"] = status
    out = []
    for s in get(CT + "?" + urllib.parse.urlencode(p)).get("studies", []):
        ps = s.get("protocolSection", {}); idm = ps.get("identificationModule", {})
        out.append({"id": idm.get("nctId", ""), "title": idm.get("briefTitle", ""),
                    "status": ps.get("statusModule", {}).get("overallStatus", ""),
                    "phase": "/".join(ps.get("designModule", {}).get("phases", [])) or "n/a",
                    "drugs": [i.get("name", "") for i in ps.get("armsInterventionsModule", {}).get("interventions", [])][:4],
                    "url": "https://clinicaltrials.gov/study/" + idm.get("nctId", "")})
    return out

GQ = """query($q:String!){search(queryString:$q,entityNames:["disease"],page:{index:0,size:1}){hits{id name}}}"""
GD = """query($id:String!){disease(efoId:$id){name description
 associatedTargets(page:{index:0,size:10}){rows{score target{approvedSymbol approvedName}}}
 drugAndClinicalCandidates{count rows{maxClinicalStage drug{name drugType}}}}}"""

def gql(q, v): return get(OT, json.dumps({"query": q, "variables": v}), {"Content-Type": "application/json"})

def biology(name):
    h = gql(GQ, {"q": name})["data"]["search"]["hits"]
    if not h: return None
    d = gql(GD, {"id": h[0]["id"]})
    if d.get("errors"): raise RuntimeError(d["errors"][0].get("message", "Open Targets query failed"))
    d = d["data"]["disease"]
    return {"id": h[0]["id"], "name": d["name"], "desc": d.get("description") or "",
            "targets": [(r["target"]["approvedSymbol"], r["target"]["approvedName"], round(r["score"], 2)) for r in d["associatedTargets"]["rows"]],
            "drugs": [(r["drug"]["name"], r["drug"].get("drugType", ""), r.get("maxClinicalStage", "")) for r in (d.get("drugAndClinicalCandidates") or {}).get("rows", [])]}

def w(t, i=2): return textwrap.fill(str(t), 78, initial_indent=" " * i, subsequent_indent=" " * (i + 2))

def sec(t): print("\n== " + t + " " + "=" * max(0, 74 - len(t)))

def safe(f, *a, **k):
    try: return f(*a, **k)
    except RuntimeError as e: print(w(f"! {e}")); return None

def show_bio(d):
    sec("BIOLOGY & DRUGS (Open Targets)")
    b = safe(biology, d)
    if b is None: print(w("no data")); return
    print(w(b["name"] + " [" + b["id"] + "]")); print(w(b["desc"][:500]))
    print("\n  Top disease-linked targets (gene, score 0-1):")
    for s, n, sc in b["targets"]: print(w(f"{s:<9}{sc}  {n}", 4))
    print("\n  Drugs / clinical candidates:")
    for n, t, ph in b["drugs"][:15]: print(w(f"{n} ({t}) stage: {ph}", 4))

def show_trials(d, n=8, st="RECRUITING"):
    sec(f"CLINICAL TRIALS ({st or 'any status'})")
    for t in safe(trials, d, n, st) or []:
        print(w(f"{t['id']} [{t['status']}, {t['phase']}] {t['title']}")); print(w("drugs: " + (", ".join(t["drugs"]) or "n/a") + "  " + t["url"], 4))

def show_papers(d, n=8, recent=True):
    sec("PAPERS (PubMed, " + ("newest" if recent else "most relevant") + ")")
    for p in safe(papers, d, n, recent) or []:
        print(w(f"{p['date']} | {p['title']}")); print(w(f"{p['authors']} - {p['journal']} {p['url']}", 4))

def overview(d):
    show_bio(d); show_trials(d); show_papers(d, 6); print("\n" + DISC)

def report(d, path=None):
    """save everything as a markdown report"""
    import io, contextlib
    b = io.StringIO()
    with contextlib.redirect_stdout(b): overview(d)
    path = path or os.path.join(os.getcwd(), d.replace(" ", "_") + "_report.md")
    open(path, "w").write(f"# MedScope report: {d}\n_{time.strftime('%Y-%m-%d')}_\n\n```\n{b.getvalue()}\n```\n"); return path

def watch():
    """emerging/new diseases: WHO outbreak news + newest PubMed on novel pathogens"""
    sec("WHO DISEASE OUTBREAK NEWS")
    r = safe(get, "https://www.who.int/api/news/diseaseoutbreaknews?sf_culture=en&$orderby=PublicationDateAndTime%20desc&$top=8")
    for i in (r or {}).get("value", []): print(w(f"{str(i.get('PublicationDateAndTime',''))[:10]} | {i.get('Title','')}"))
    show_papers("(novel OR emerging OR outbreak) AND (pathogen OR virus OR infection) AND (\"last 30 days\"[PDat])", 8)

LADDER = "Evidence ladder: lab/animal -> Ph1 safety -> Ph2 signal -> Ph3 proof -> approval. Only Ph3+/approval is real evidence of benefit."

def leads(d):
    """candidate leads ranked by evidence stage. These are HYPOTHESES, never confirmed cures."""
    sec("LEADS FOR " + d.upper()); print(w(LADDER))
    ts = []
    for ph in ("PHASE3", "PHASE2"):
        ts += [t for t in safe(trials, d + " " + ph, 10, None) or [] if ph in t["phase"]]
    for t in ts[:10]: print(w(f"[{t['phase']}|{t['status']}] {', '.join(t['drugs']) or 'n/a'} - {t['title']} {t['url']}"))
    show_papers(d + " drug repurposing", 6, False)
    print("\n  Treat everything above as leads to discuss with researchers/clinicians, not as cures.")

HELP = """MedScope commands (or run: python medscope.py <cmd> <disease>)
  overview <disease>   biology + drugs + trials + papers
  bio <disease>        targets & drugs (Open Targets)
  trials <disease> [status]   status: recruiting|completed|any (default recruiting)
  papers <disease> [best]     newest papers, or 'best' = most relevant
  report <disease>     save markdown report
  leads <disease>      Ph2/Ph3 trials + repurposing papers, ranked by evidence stage
  design <disease>     candidates -> rat PBPK -> human PBPK (chemicals.py; real structures with RDKit)
  search <query>       Europe PMC + OpenAlex in parallel (free, no key)
  brain <disease>      AI research brief (local Ollama = unlimited; or GEMINI_API_KEY / GROQ_API_KEY free tiers)
  diseases             list per-disease targets (hiv, covid, flu, tb, lung-cancer, ...)
  check <disease>      validate that disease's model: known GOOD vs BAD drugs (AUC)
  discover <disease>   full per-disease pipeline (needs RDKit; Vina/ADMET-AI optional)
  lab add|report|list  log lab results & track REAL prediction accuracy (closes the loop)
  score                evidence-based rating of MedScope (verified vs ceiling)
  onboard <disease>    add ANY disease automatically (targets, ChEMBL, PDB); or --uniprot ACC --name X for a new pathogen
  repurpose <disease>  screen ALL approved drugs against that disease's target
  screen <disease|all> [--limit N]  million-scale ChEMBL screen (potency model + organ risk), all diseases
  evolve <disease>     generate NOVEL molecules by evolutionary design
  campaign [--evolve]  every disease at once -> master table
  drug <name>          any drug/product: ingredients, structure, human PBPK, flags
  combo <a> <b> ...    combination/regimen checks (organ overlap, shared tox, CYP DDI)
  selftest             what works on this device + offline correctness checks
  watch                emerging diseases: WHO outbreak news + newest novel-pathogen papers
  help | quit   (menu mode remembers your last disease: just type a command alone)"""

def run(cmd, rest, last):
    if cmd in ("prefetch", "download", "audit"):
        import offline; offline.cli([cmd] + rest.split()); return last
    if cmd in ("screen", "evolve", "campaign"):
        try:
            import importlib; importlib.import_module(cmd).cli(rest.split())
        except Exception as e: print(w(f"! {e}"))
        return last
    if cmd in ("recreate", "stabilize", "population", "person"):
        try:
            if cmd == "recreate": import recreation_bioengineering as x; x.cli(rest.split())
            elif cmd == "stabilize": import stabilizer as x; x.cli(rest.split())
            else: import mass_population as x; x.cli((["person"] if cmd == "person" else []) + rest.split())
        except Exception as e: print(w(f"! {e}"))
        return last
    if cmd in ("mechanism", "how"):
        import reverse_bioengineering as rb
        try: rb.analyze(rest.replace("--save","").strip(), save="--save" in rest) if rest.replace("--save","").strip() else print("usage: mechanism <disease> [--save]")
        except Exception as e: print(w(f"! {e}"))
        return last
    if cmd == "score": import scorecard; scorecard.main(); return last
    if cmd == "onboard":
        import onboard
        try: onboard.cli(rest.split())
        except Exception as e: print(w(f"! {e}"))
        return last
    if cmd in ("drug", "combo"):
        import ingredients
        try: ingredients.cli(cmd, rest.split())
        except Exception as e: print(w(f"! {e}"))
        return last
    if cmd == "selftest": import selftest; selftest.main(); return last
    if cmd == "diseases": import diseases; diseases.show(); return last
    if cmd == "lab": import lab; lab.cli(rest.split()); return last
    a = rest.split(); st = {"recruiting": "RECRUITING", "completed": "COMPLETED", "any": None}
    if cmd == "trials" and a and a[-1].lower() in st: s = st[a.pop().lower()]
    else: s = "RECRUITING"
    best = cmd == "papers" and a and a[-1].lower() == "best"
    if best: a.pop()
    d = " ".join(a) or last
    if cmd == "watch": watch(); return last
    if cmd in ("overview", "bio", "trials", "papers", "report", "leads", "design", "search", "brain", "discover", "check", "repurpose") and not d: print(w("give a disease name, e.g. overview hiv")); return last
    if cmd == "overview": overview(d)
    elif cmd == "bio": show_bio(d)
    elif cmd == "trials": show_trials(d, 10, s)
    elif cmd == "papers": show_papers(d, 10, not best)
    elif cmd == "report": print("saved:", report(d))
    elif cmd == "leads": leads(d)
    elif cmd == "design": import chemicals; chemicals.pipeline(d)
    elif cmd == "discover":
        import discover
        try: discover.run(d.split()[0], dock_on="--nodock" not in rest)
        except Exception as e: print(w(f"! {e}"))
    elif cmd == "repurpose":
        import library
        try: library.repurpose(d.split()[0])
        except Exception as e: print(w(f"! {e}"))
    elif cmd == "check":
        import diseases
        try: diseases.check(d.split()[0])
        except Exception as e: print(w(f"! {e}"))
    elif cmd == "search":
        import sources; sec("ALL FREE ENGINES (Europe PMC + OpenAlex, parallel)")
        for x in sources.search_all(d, 8): print(w(f"{x['date']} [{x['src']}] {x['title']}")); print(w(f"{x['venue']} {x['url']}", 4))
    elif cmd == "brain":
        import brain; ev = [f"paper {p['date']}: {p['title']}" for p in safe(papers, d, 8) or []] + \
            [f"trial {t['phase']}/{t['status']}: {t['title']} drugs={','.join(t['drugs'])}" for t in safe(trials, d, 8, None) or []]
        try: n, t = brain.brief(d, ev); sec("AI BRIEF via " + n); print(t); print("\n" + DISC)
        except RuntimeError as e: print(w("! " + str(e)))
    elif cmd in ("help", "h", "?"): print(HELP); return last
    else: print("unknown command. type help"); return last
    return d

def main():
    if len(sys.argv) > 1: run(sys.argv[1].lower(), " ".join(sys.argv[2:]), ""); return
    print("MedScope - disease research explorer\n" + DISC + "\n" + HELP); last = ""
    while True:
        try: l = input("\nmedscope> ").strip()
        except (EOFError, KeyboardInterrupt): break
        if not l: continue
        if l.lower() in ("quit", "exit", "q"): break
        c, _, r = l.partition(" "); last = run(c.lower(), r, last)

if __name__ == "__main__": main()
