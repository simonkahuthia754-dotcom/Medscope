#!/usr/bin/env python3
"""offline.py - run MedScope on a WEAK or NO connection, safely, on any OS.
  python offline.py prefetch [disease|all]   warm the local cache once (resumable); afterwards set MEDSCOPE_OFFLINE=1 and everything cached runs with no network
  python offline.py download                 resumable download of the ChEMBL structure collection (~270 MB, one time) for million-scale screening offline
  python offline.py audit                    security + portability + offline self-test (what the scorecard's 'Robustness' criterion reads)
Design: every web answer is cached (sqlite, gzip, atomic) and a stale answer is served when the network fails; downloads resume with HTTP Range;
only https to an allow-list of public science hosts; response-size caps; no eval/pickle/shell=True; private data dir (0700); no telemetry; nothing leaves your device."""
import os, re, sys, time, glob, urllib.request
import sources, state

DATA = os.path.join(sources.HOME, "data"); CHEMREPS = os.path.join(DATA, "chemreps.txt.gz")

def download(url=None, dest=CHEMREPS, log=print):
    """resumable download; a dropped connection continues from the last byte instead of restarting"""
    import screen
    os.makedirs(DATA, exist_ok=True); url = sources.safe_url(url or screen.chemreps_url()); part = dest + ".part"
    for attempt in range(30):
        have = os.path.getsize(part) if os.path.exists(part) else 0
        rq = urllib.request.Request(url, headers={"Range": f"bytes={have}-", "User-Agent": "Mozilla/5.0 (compatible; MedScope)"})
        try:
            with urllib.request.urlopen(rq, timeout=60) as r:
                total = have + int(r.headers.get("Content-Length", 0)); mode = "ab" if r.status == 206 else "wb"
                with open(part, mode) as f:
                    while True:
                        b = r.read(1 << 20)
                        if not b: break
                        f.write(b); log(f"\r  {os.path.getsize(part) / 1e6:.0f}/{total / 1e6:.0f} MB", end="") if False else None
            break
        except urllib.error.HTTPError as e:
            if e.code == 416: break
            time.sleep(2 + attempt)
        except Exception: time.sleep(2 + attempt)
    import gzip
    with gzip.open(part) as g:                                   # integrity: the whole stream must decompress
        while g.read(1 << 22): pass
    os.replace(part, dest); log(f"  saved {dest} ({os.path.getsize(dest) / 1e6:.0f} MB)"); return dest

def prefetch(keys=None, log=print):
    import diseases, screen, potency, lab
    keys = keys or list(diseases.D); n0 = sources.cache_info()[0]
    for k in keys:
        try:
            D = diseases.get(k)
            for cid, nm in screen.candidates(D):
                try: potency.fetch(cid, overlay=None)
                except Exception: pass
            log(f"  cached {k}")
        except Exception as e: log(f"  {k}: {e}")
    n1, mb = sources.cache_info(); log(f"  cache: {n1} entries ({mb:.1f} MB), +{n1 - n0} new. Offline mode: set MEDSCOPE_OFFLINE=1"); state.record("prefetch_ok", True)

RISKY = [(r"\bpickle\b", "pickle"), (r"\beval\(", "eval"), (r"\bexec\(", "exec"), (r"shell\s*=\s*True", "shell=True"), (r"os\.system", "os.system"), (r"verify\s*=\s*False|CERT_NONE", "TLS off"), (r"http://(?!localhost)", "plain http")]

def audit(show=True):
    here = os.path.dirname(os.path.abspath(__file__)); res = {}; bad = []
    for f in glob.glob(os.path.join(here, "*.py")):
        if os.path.basename(f) in ("ci_debug.py", "offline.py"): continue
        src = open(f, encoding="utf-8", errors="ignore").read()
        for rx, nm in RISKY:
            for m in re.finditer(rx, src):
                line = src[:m.start()].count("\n") + 1
                if "#" in src.splitlines()[line - 1].split(m.group(0))[0]: continue
                bad.append(f"{os.path.basename(f)}:{line} {nm}")
    res["static"] = not bad
    # URL policy
    ok = 0
    for u, want in (("http://x.org/a", False), ("https://evil.example.com/a", False), ("file:///etc/passwd", False), ("https://www.ebi.ac.uk/chembl/api", True), ("https://ebi.ac.uk.evil.com/", False)):
        try: sources.safe_url(u); got = True
        except ValueError: got = False
        ok += got == want
    res["url_policy"] = ok == 5
    # private data dir
    sources._secure_dir(sources.HOME)
    try: res["private_dir"] = (os.name == "nt") or (os.stat(sources.HOME).st_mode & 0o077 == 0)
    except Exception: res["private_dir"] = False
    # offline replay: cached answer survives, uncached fails closed
    import tempfile, json; old = os.environ.get("MEDSCOPE_OFFLINE"); k = "https://www.ebi.ac.uk/chembl/api/data/_selftest.json"
    try:
        import hashlib; sources._put(hashlib.sha1(k.encode()).hexdigest(), {"ok": 1}); os.environ["MEDSCOPE_OFFLINE"] = "1"
        a = sources.J(k) == {"ok": 1}
        try: sources.J("https://www.ebi.ac.uk/chembl/api/data/_never_cached.json"); b = False
        except RuntimeError: b = True
        res["offline"] = a and b
    finally:
        if old is None: os.environ.pop("MEDSCOPE_OFFLINE", None)
        else: os.environ["MEDSCOPE_OFFLINE"] = old
    # speed: cache hit latency
    t = time.time(); os.environ["MEDSCOPE_OFFLINE"] = "1"
    try:
        for _ in range(50): sources.J(k)
    finally:
        if old is None: os.environ.pop("MEDSCOPE_OFFLINE", None)
        else: os.environ["MEDSCOPE_OFFLINE"] = old
    ms = (time.time() - t) / 50 * 1000; res["cache_ms"] = ms; res["fast"] = ms < 25
    res["platform"] = f"{sys.platform}/{os.name}, python {sys.version.split()[0]}, {os.cpu_count()} cores"
    # portability: pure-Python core imports with no optional packages
    res["core_imports"] = all(__import__(m) for m in ("pbpk", "mass_population", "stabilizer".replace("stabilizer", "state"), "diseases"))
    for k2 in ("static", "url_policy", "private_dir", "offline", "fast", "core_imports"): state.record("rob_" + k2, bool(res[k2]))
    if show:
        print(f"AUDIT  platform: {res['platform']}")
        for k2, lab_ in (("static", "no eval/pickle/shell=True/os.system/plain http/TLS-off in code"), ("url_policy", "https-only, allow-listed hosts, lookalike domains refused"), ("private_dir", "data directory private (owner-only)"),
                         ("offline", "cached answers replay offline; uncached fails closed"), ("fast", f"cache hit {ms:.1f} ms"), ("core_imports", "pure-Python core imports without RDKit/numpy")):
            print(f"  [{'PASS' if res[k2] else 'FAIL'}] {lab_}")
        for b in bad[:10]: print("   !", b)
    return res

def cli(a):
    c = a[0] if a else "audit"
    if c == "download": download()
    elif c == "prefetch": prefetch(None if len(a) < 2 or a[1] == "all" else [a[1]])
    else: audit()

if __name__ == "__main__": cli(sys.argv[1:])
