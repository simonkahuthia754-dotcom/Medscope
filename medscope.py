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
    k = hashlib.md5((url + (data or "")).encode()).hexdigest(); c = _cache()
    if k in c and time.time() - c[k][0] < TTL: return c[k][1]
    rq = urllib.request.Request(url, data=data.encode() if data else None, headers={"User-Agent": "MedScope/1.0", **(hdr or {})})
    try:
        with urllib.request.urlopen(rq, timeout=25) as r: j = json.loads(r.read().decode())
    except Exception as e:
        raise RuntimeError(f"network/API error: {e}")
    c[k] = [time.time(), j]
    try: json.dump(c, open(CF, "w"))
    except Exception: pass
    return j

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

HELP = """MedScope commands (or run: python medscope.py <cmd> <disease>)
  overview <disease>   biology + drugs + trials + papers
  bio <disease>        targets & drugs (Open Targets)
  trials <disease> [status]   status: recruiting|completed|any (default recruiting)
  papers <disease> [best]     newest papers, or 'best' = most relevant
  report <disease>     save markdown report
  help | quit   (menu mode remembers your last disease: just type a command alone)"""

def run(cmd, rest, last):
    a = rest.split(); st = {"recruiting": "RECRUITING", "completed": "COMPLETED", "any": None}
    if cmd == "trials" and a and a[-1].lower() in st: s = st[a.pop().lower()]
    else: s = "RECRUITING"
    best = cmd == "papers" and a and a[-1].lower() == "best"
    if best: a.pop()
    d = " ".join(a) or last
    if cmd in ("overview", "bio", "trials", "papers", "report") and not d: print(w("give a disease name, e.g. overview hiv")); return last
    if cmd == "overview": overview(d)
    elif cmd == "bio": show_bio(d)
    elif cmd == "trials": show_trials(d, 10, s)
    elif cmd == "papers": show_papers(d, 10, not best)
    elif cmd == "report": print("saved:", report(d))
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
