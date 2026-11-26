#!/usr/bin/env python3
"""verify.py - run every live check that earns scorecard points, on a machine WITH internet + RDKit/Vina/ADMET-AI (PC, Colab, GitHub Actions).
  python verify.py        -> runs each step, tolerates failures, prints scorecard, saves medscope_out/verify.json + scorecard.txt"""
import contextlib, io, json, os, subprocess, sys, time, traceback
os.makedirs("medscope_out", exist_ok=True); R = {}
def step(name, f):
    t = time.time(); buf = io.StringIO()
    try:
        with contextlib.redirect_stdout(buf): f()
        R[name] = {"ok": True, "s": round(time.time() - t)}
    except BaseException as e:
        R[name] = {"ok": False, "err": f"{type(e).__name__}: {e}", "tb": traceback.format_exc()[-900:]}
    out = buf.getvalue(); open(f"medscope_out/{name}.log", "w").write(out); print("\n".join("    | "+l[:200] for l in out.splitlines() if l.strip() and "it/s" not in l and "GPU" not in l and "TPU" not in l and "Tip:" not in l and "Warning" not in l)[-3500:], flush=True); print(f"[{'OK ' if R[name]['ok'] else 'FAIL'}] {name}" + ("" if R[name]["ok"] else "  " + R[name]["err"]), flush=True)
def env():
    for m in ("vina", "meeko", "rdkit", "admet_ai"):
        try: __import__(m); print("import ok", m)
        except BaseException as e: print("IMPORT FAIL", m, type(e).__name__, e)
    try: from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate; print("meeko names ok")
    except BaseException as e: print("MEEKO NAMES FAIL", e)
    import humanbio, json
    def t(q, lab):
        try: print(lab, json.dumps(humanbio._gql(q, {})["data"])[:300])
        except BaseException as e: print(lab, "FAIL", getattr(e, "read", lambda: b"")()[:400] if hasattr(e, "read") else e)
    try:
        f = humanbio._gql('{__type(name:"Target"){fields{name type{name kind ofType{name}}}}}', {})["data"]["__type"]["fields"]
        print("Target expression-like fields:", [(x["name"], x["type"]["name"] or x["type"]["ofType"] and x["type"]["ofType"]["name"]) for x in f if "xpress" in x["name"].lower()])
    except BaseException as e: print("introspection FAIL", e)
def sel(): import selftest; selftest.main()
def val(): import validate, admet; validate.run(False); admet.HAVE and validate.run(True)
def onb(): import onboard; onboard.onboard("tuberculosis", 2)
def lib(): import library; assert len(library.approved(1500)) >= 1000
def drg(): import ingredients; ingredients.cli("drug", ["metformin"])
def cmb(): import ingredients; ingredients.cli("combo", ["abacavir", "dolutegravir", "lamivudine"])
def dsc(): import discover; discover.run("lung-cancer", 5, 30, 50, True)
def rep(): import library; library.repurpose("blood-cancer", 10)
def mech(): import reverse_bioengineering as rb; rb.analyze("HIV infection", key="hiv")
def disc_chk(): import diseases; diseases.check("lung-cancer")
for n, f in [("env", env), ("selftest", sel), ("validate", val), ("onboard", onb), ("library", lib), ("drug", drg), ("combo", cmb), ("discover", dsc), ("repurpose", rep), ("mechanism", mech), ("check", disc_chk)]: step(n, f)
import scorecard; buf = io.StringIO()
with contextlib.redirect_stdout(buf): scorecard.main(refresh=False)
print(buf.getvalue()); open("medscope_out/scorecard.txt", "w").write(buf.getvalue()); json.dump(R, open("medscope_out/verify.json", "w"), indent=1)
