#!/usr/bin/env python3
"""verify.py - run every live check that earns scorecard points, on a machine WITH internet + RDKit/Vina/ADMET-AI (PC, Colab, GitHub Actions).
  python verify.py        -> runs each step, tolerates failures, prints scorecard, saves medscope_out/verify.json + scorecard.txt"""
import os, contextlib, io, json, os, subprocess, sys, time, traceback
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
    try: lv, ev = humanbio.target_expression("EGFR"); print("expression EGFR:", lv, "| safety events:", ev[:5])
    except BaseException as e: print("EXPRESSION FAIL", type(e).__name__, e)
    try:
        f = humanbio._gql('{__type(name:"Target"){fields{name type{name kind ofType{name}}}}}', {})["data"]["__type"]["fields"]
        print("Target expression-like fields:", [(x["name"], x["type"]["name"] or x["type"]["ofType"] and x["type"]["ofType"]["name"]) for x in f if "xpress" in x["name"].lower()])
    except BaseException as e: print("introspection FAIL", e)
def dck():
    import diseases, docking, potency, state, traceback
    D = diseases.get("lung-cancer"); print("pdb", D.get("pdb"), "chembl", D.get("chembl"))
    try:
        T = docking.prepare(D["pdb"]); print("prepared", T["lig"], T["center"], T["size"])
        r = docking.redock_rmsd(T); print("redock (score, RMSD A, cluster size, runs):", r)
        if r: state.record("redock_rmsd", round(r[1], 2))
        rate, n, rows = docking.panel_redock(); print("PANEL redock success (<=2 A):", rate, "of", n)
        if rate is not None and n >= 5: state.record("redock_rate", round(rate, 2)); state.record("redock_n", n)
        data = potency.fetch(D["chembl"]); cal = docking.calibrate(T, data); print("calibration", {k: v for k, v in cal.items() if k != "pts"}); 
        if cal.get("rho") is not None: state.record("dock_rho", round(cal["rho"], 2))
    except BaseException: print(traceback.format_exc()[-1800:]); raise
def popn(): import mass_population as m; assert m.selftest(); m.cli(["metformin", "--n", "10000", "--ic50", "1000"])
def rec(): import recreation_bioengineering as r; r.cli(["all"])
def stab(): import stabilizer; stabilizer.cli(["erlotinib", "--ic50", "20", "--target", "lung-cancer", "--n", "40"])
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
def scl():
    import screen, evolve, diseases, state; r = screen.run(None, limit=int(os.environ.get("SCREEN_LIMIT", "300000")), top=5, show=False); tot = 0
    for k in list(r)[:6]:
        try: e = evolve.evolve(k, gens=8, pop=300, seeds=25, top=3, show=False); tot += e["novel"]; print(f"  evolved {k}: {e['novel']} novel molecules")
        except Exception as ex: print(f"  evolve {k}: {ex}")
    state.record("evolve_novel", tot); print(f"  screened diseases {len(r)}, novel generated {tot}")
def rob():
    import offline; offline.prefetch(["lung-cancer"]); offline.audit()
for n, f in [("env", env), ("selftest", sel), ("validate", val), ("onboard", onb), ("library", lib), ("drug", drg), ("combo", cmb), ("population", popn), ("recreate", rec), ("stabilize", stab), ("docking", dck), ("discover", dsc), ("repurpose", rep), ("mechanism", mech), ("scale", scl), ("robust", rob), ("check", disc_chk)]: step(n, f)
import scorecard; buf = io.StringIO()
with contextlib.redirect_stdout(buf): scorecard.main(refresh=False)
print(buf.getvalue()); open("medscope_out/scorecard.txt", "w").write(buf.getvalue()); json.dump(R, open("medscope_out/verify.json", "w"), indent=1)
