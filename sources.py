"""sources.py - free, keyless search engines + real compound data. Stdlib only, parallel, cached.
Europe PMC & OpenAlex (papers) | PubChem (real properties) | ChEMBL (measured IC50) | + PubMed/CT.gov in medscope.py"""
import json, os, time, hashlib, statistics, urllib.request, urllib.parse, urllib.error
from concurrent.futures import ThreadPoolExecutor

import gzip, sqlite3, threading, urllib.parse as _up
HOME = os.environ.get("MEDSCOPE_HOME") or os.path.join(os.path.expanduser("~"), ".medscope")
CF = os.path.join(HOME, "cache.db"); TTL = 12 * 3600; STALE_OK = 365 * 86400; MAXB = 64 * 1024 * 1024
HOSTS = ("ebi.ac.uk", "ncbi.nlm.nih.gov", "openalex.org", "rcsb.org", "opentargets.org", "string-db.org", "uniprot.org", "nlm.nih.gov", "clinicaltrials.gov", "who.int", "wwwdev.ebi.ac.uk", "files.rcsb.org", "data.rcsb.org")
Q = urllib.parse.quote; _L = threading.Lock(); STATS = dict(hit=0, net=0, stale=0, fail=0)

def offline(): return os.environ.get("MEDSCOPE_OFFLINE", "") not in ("", "0")

def _secure_dir(p):
    os.makedirs(p, exist_ok=True)
    try: os.chmod(p, 0o700)
    except Exception: pass

def _db():
    _secure_dir(HOME); c = sqlite3.connect(CF, timeout=30); c.execute("PRAGMA journal_mode=WAL"); c.execute("CREATE TABLE IF NOT EXISTS c(k TEXT PRIMARY KEY, t REAL, v BLOB)"); return c

def _get(k):
    try:
        with _L, _db() as c: r = c.execute("SELECT t, v FROM c WHERE k=?", (k,)).fetchone()
        return (r[0], json.loads(gzip.decompress(r[1]).decode())) if r else None
    except Exception: return None

def _put(k, j):
    try:
        with _L, _db() as c: c.execute("REPLACE INTO c VALUES(?,?,?)", (k, time.time(), gzip.compress(json.dumps(j).encode(), 3)))
    except Exception: pass

def safe_url(url):
    u = _up.urlparse(url)
    if u.scheme != "https" or not u.hostname or not any(u.hostname == h or u.hostname.endswith("." + h) for h in HOSTS): raise ValueError(f"blocked URL (https + allow-listed hosts only): {url[:80]}")
    return url

def J(url, data=None, hdr=None, ttl=TTL):
    """cached JSON fetch: fresh cache -> network (retry, size cap, https allow-list) -> STALE cache if the network is down; MEDSCOPE_OFFLINE=1 never touches the network"""
    safe_url(url); k = hashlib.sha1((url + (data or "")).encode()).hexdigest(); c = _get(k)
    if c and (time.time() - c[0] < ttl or offline()): STATS["hit"] += 1; return c[1]
    if offline(): STATS["fail"] += 1; raise RuntimeError("offline and not cached: run `prefetch` once on a connection")
    rq = urllib.request.Request(url, data=data.encode() if data else None, headers={"User-Agent": "Mozilla/5.0 (compatible; MedScope)", "Accept-Encoding": "gzip", **(hdr or {})}); err = None; j = None
    for a in range(4):   # retry transient server/network errors (5xx, timeouts), never client errors
        try:
            with urllib.request.urlopen(rq, timeout=40) as r:
                raw = r.read(MAXB + 1)
                if len(raw) > MAXB: raise ValueError("response too large")
                if r.headers.get("Content-Encoding") == "gzip": raw = gzip.decompress(raw)
                j = json.loads(raw.decode()); break
        except urllib.error.HTTPError as e:
            err = e
            if e.code < 500 or a == 3: break
        except (TimeoutError, urllib.error.URLError, ConnectionError, OSError) as e:
            err = e
            if a == 3: break
        time.sleep(1.5 * 2 ** a)
    if j is not None:
        STATS["net"] += 1; _put(k, j); return j
    if c and time.time() - c[0] < STALE_OK: STATS["stale"] += 1; return c[1]     # weak/no connection: serve the last good answer
    STATS["fail"] += 1; raise err or RuntimeError("fetch failed")

def cache_info():
    try:
        with _db() as c: n, = c.execute("SELECT COUNT(*) FROM c").fetchone()
        return n, os.path.getsize(CF) / 1e6
    except Exception: return 0, 0.0


def europepmc(q, n=10):
    r = J(f"https://www.ebi.ac.uk/europepmc/webservices/rest/search?query={Q(q)}&format=json&pageSize={n}&sort=P_PDATE_D%20desc")["resultList"]["result"]
    return [{"src": "EuropePMC", "title": x.get("title", ""), "date": x.get("firstPublicationDate", x.get("pubYear", "")), "venue": x.get("journalTitle", ""),
             "url": f"https://europepmc.org/article/{x.get('source','MED')}/{x.get('id','')}"} for x in r]

def openalex(q, n=10):
    r = J(f"https://api.openalex.org/works?search={Q(q)}&per-page={n}&sort=publication_date:desc")["results"]
    return [{"src": "OpenAlex", "title": x.get("display_name", ""), "date": x.get("publication_date", ""),
             "venue": ((x.get("primary_location") or {}).get("source") or {}).get("display_name", ""), "url": x.get("doi") or x.get("id", "")} for x in r]

def pubchem(name):
    """REAL molecular properties + structure: mw, XLogP, HBD, HBA, TPSA, SMILES"""
    for smi in ("ConnectivitySMILES", "CanonicalSMILES"):
        try: p = J(f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{Q(name)}/property/MolecularWeight,XLogP,HBondDonorCount,HBondAcceptorCount,TPSA,{smi}/JSON")["PropertyTable"]["Properties"][0]; break
        except Exception as e: err = e
    else: raise err
    return {"mw": float(p["MolecularWeight"]), "logp": float(p.get("XLogP", 0) or 0), "hbd": p["HBondDonorCount"], "hba": p["HBondAcceptorCount"],
            "tpsa": float(p["TPSA"]) if p.get("TPSA") is not None else None, "smiles": next((v for k, v in p.items() if "SMILES" in k), None), "cid": p["CID"]}

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
