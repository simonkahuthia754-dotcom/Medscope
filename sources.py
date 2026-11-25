"""sources.py - free, keyless search engines + real compound data. Stdlib only, parallel, cached.
Europe PMC & OpenAlex (papers) | PubChem (real properties) | ChEMBL (measured IC50) | + PubMed/CT.gov in medscope.py"""
import json, os, time, hashlib, statistics, urllib.request, urllib.parse
from concurrent.futures import ThreadPoolExecutor

CF = os.path.join(os.path.expanduser("~"), ".medscope", "cache2.json"); TTL = 12 * 3600
Q = urllib.parse.quote

def J(url, data=None, hdr=None, ttl=TTL):
    k = hashlib.md5((url + (data or "")).encode()).hexdigest()
    try: c = json.load(open(CF))
    except Exception: c = {}
    if k in c and time.time() - c[k][0] < ttl: return c[k][1]
    rq = urllib.request.Request(url, data=data.encode() if data else None, headers={"User-Agent": "MedScope/2.0", **(hdr or {})})
    with urllib.request.urlopen(rq, timeout=25) as r: j = json.loads(r.read().decode())
    c[k] = [time.time(), j]
    try: os.makedirs(os.path.dirname(CF), exist_ok=True); json.dump(c, open(CF, "w"))
    except Exception: pass
    return j

def europepmc(q, n=10):
    r = J(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={Q(q)}&format=json&pageSize={n}&sort=P_PDATE_D%20desc")["resultList"]["result"]
    return [{"src": "EuropePMC", "title": x.get("title", ""), "date": x.get("firstPublicationDate", x.get("pubYear", "")), "venue": x.get("journalTitle", ""),
             "url": f"https://europepmc.org/article/{x.get('source','MED')}/{x.get('id','')}"} for x in r]

def openalex(q, n=10):
    r = J(f"https://api.openalex.org/works?search={Q(q)}&per-page={n}&sort=publication_date:desc")["results"]
    return [{"src": "OpenAlex", "title": x.get("display_name", ""), "date": x.get("publication_date", ""),
             "venue": ((x.get("primary_location") or {}).get("source") or {}).get("display_name", ""), "url": x.get("doi") or x.get("id", "")} for x in r]

def pubchem(name):
    """REAL molecular properties: (MW, XLogP, HBD, HBA)"""
    p = J(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{Q(name)}/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount/JSON")["PropertyTable"]["Properties"][0]
    return {"mw": float(p["MolecularWeight"]), "logp": float(p.get("XLogP", 0) or 0), "hbd": p["HBondDonorCount"], "hba": p["HBondAcceptorCount"], "cid": p["CID"]}

def chembl_ic50(name):
    """REAL measured potency: median IC50 (nM) across ChEMBL assays + number of measurements"""
    m = J(f"https://www.ebi.ac.uk/chembl/api/data/molecule/search.json?q={Q(name)}&limit=1")["molecules"]
    if not m: return None
    cid = m[0]["molecule_chembl_id"]
    a = J(f"https://www.ebi.ac.uk/chembl/api/data/activity.json?molecule_chembl_id={cid}&standard_type=IC50&standard_units=nM&limit=100")["activities"]
    v = [float(x["standard_value"]) for x in a if x.get("standard_value") and float(x["standard_value"]) > 0]
    return {"chembl": cid, "ic50_nM": statistics.median(v), "n": len(v)} if v else {"chembl": cid, "ic50_nM": None, "n": 0}

def search_all(q, n=8):
    """query every free paper index in parallel, dedupe by title"""
    with ThreadPoolExecutor(4) as ex: fs = [ex.submit(f, q, n) for f in (europepmc, openalex)]
    out, seen = [], set()
    for f in fs:
        try:
            for x in f.result():
                t = x["title"].lower()[:60]
                if t and t not in seen: seen.add(t); out.append(x)
        except Exception as e: out.append({"src": "error", "title": str(e), "date": "", "venue": "", "url": ""})
    return sorted(out, key=lambda x: x["date"], reverse=True)

def enrich(names):
    """fetch real props + measured IC50 for many compounds in parallel -> {name: {...}} (skips failures)"""
    def one(n):
        try:
            d = pubchem(n)
            try: d.update(chembl_ic50(n) or {})
            except Exception: pass
            return n, d
        except Exception: return n, None
    with ThreadPoolExecutor(6) as ex: return {n: d for n, d in ex.map(one, names) if d}
