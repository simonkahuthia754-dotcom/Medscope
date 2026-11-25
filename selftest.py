#!/usr/bin/env python3
"""selftest.py - tells you EXACTLY what works on this device (Termux, PC, Colab) and runs offline correctness checks.
python selftest.py"""
import importlib, math, os, shutil, sys, tempfile

def have(m):
    try: importlib.import_module(m); return True
    except Exception: return False

def net():
    try:
        import sources; sources.J("https://www.ebi.ac.uk/chembl/api/data/status.json", ttl=0); return True
    except Exception: return False

def main():
    termux = "com.termux" in os.environ.get("PREFIX", "") or "TERMUX_VERSION" in os.environ
    print(f"Python {sys.version.split()[0]} | {'Termux/Android' if termux else sys.platform}")
    cap = {"rdkit": have("rdkit"), "vina": have("vina"), "meeko": have("meeko"), "obabel": bool(shutil.which("obabel")), "admet-ai": have("admet_ai"),
           "scikit-learn": have("sklearn"), "internet (ChEMBL)": net()}
    for k, v in cap.items(): print(f"  {'OK ' if v else '-- '} {k}")
    t1 = True; t2 = cap["rdkit"] and cap["internet (ChEMBL)"]; t3 = t2 and cap["admet-ai"] and cap["vina"] and cap["meeko"] and cap["obabel"]
    print("\nCapability tiers:\n  Tier 1 (pure Python): whole-body PBPK, virtual population, lab loop, literature search ........ AVAILABLE" if t1 else "")
    print(f"  Tier 2 (+RDKit +internet): real structures, QSAR potency from ChEMBL, per-disease validation ... {'AVAILABLE' if t2 else 'missing: ' + ', '.join(k for k in ('rdkit', 'internet (ChEMBL)') if not cap[k])}")
    print(f"  Tier 3 (+ADMET-AI +Vina): trained ADMET models, docking with trust gates ........................ {'AVAILABLE' if t3 else 'missing: ' + ', '.join(k for k in ('admet-ai', 'vina', 'meeko', 'obabel') if not cap[k])}")
    if termux and not cap["rdkit"]: print("\nTermux tip: RDKit/Vina/torch are hard to build on Android. Run Tier 2/3 on a PC, a free Google Colab notebook, or a cheap cloud VM; keep Termux for Tier 1 and for logging lab results.")
    ok = []
    def chk(name, f):
        try: r = bool(f())
        except Exception as e: r = False; print("   error:", e)
        ok.append(r); print(f"  [{'PASS' if r else 'FAIL'}] {name}")
    print("\nOffline checks:")
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    import pbpk, humanbio, diseases, lab
    P = {"logp": 2.0, "mw": 350, "tpsa": 80}
    chk("PBPK conserves mass (human)", lambda: abs(pbpk.simulate(P, "human", 100, 72)["mass_balance"] - 1) < 1e-3)
    chk("PBPK conserves mass (rat)", lambda: abs(pbpk.simulate(P, "rat", .5, 72)["mass_balance"] - 1) < 1e-3)
    chk("rat clears faster than human", lambda: pbpk.simulate(P, "rat", .5, 72)["t_half"] < pbpk.simulate(P, "human", 100, 72)["t_half"])
    chk("Vss override honoured", lambda: abs(pbpk.simulate(P, "human", 100, 48, ovr={"vss": 1.0})["vss_L"] - 70) < 1)
    chk("virtual population runs", lambda: humanbio.population_pk(P, 100, subs=humanbio.population(20), hours=48)["n"] == 20)
    chk("every disease has its own target + good list, no good/bad overlap", lambda: all(v["target"] and v["good"] and not (set(v["good"]) & set(diseases.bad_for(k))) for k, v in diseases.D.items() if k in diseases.BUILTIN))
    def labloop():
        old_home, old_f = os.environ.get("HOME"), lab.F; d = tempfile.mkdtemp(); lab.F = os.path.join(d, "lab.jsonl")
        try:
            lab.add("x", "CCO", 500, 8000, 7.0); lab.add("x", "CCN", 20000, None, 8.0); return bool(lab.training_overlay("x")) and lab.failed("x") == ["CCN"] and len(lab.load("x")) == 2
        finally:
            lab.F = old_f
            if old_home is not None: os.environ["HOME"] = old_home
    chk("lab loop logs, overlays measured data, tracks failures", labloop)
    import state; state.record("offline_pass", all(ok)); state.record("tier2", bool(t2)); state.record("tier3", bool(t3))
    print(f"\n{sum(ok)}/{len(ok)} offline checks passed")
    return all(ok)

if __name__ == "__main__": sys.exit(0 if main() else 1)
