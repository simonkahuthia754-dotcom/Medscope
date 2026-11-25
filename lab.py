"""lab.py - the REAL feedback loop: predict -> lab measures -> log -> retrain -> predict better.
Every run is logged (who / when / molecule / prediction / measurement / cost / days) to ~/.medscope/lab_history.jsonl.
Measured lab values REPLACE ChEMBL values in that disease's training set; FAILED molecules down-weight their close analogs.
HONEST METRIC: `lab report` shows the PROSPECTIVE error of your predictions (predicted vs later-measured). That number - not a
promised '98%' - is your real accuracy, and it only improves if the model actually learns from your data. Ethics/IRB approval and a
qualified lab partner are required before any animal or human work."""
import getpass, json, math, os, time

F = os.path.join(os.path.expanduser("~"), ".medscope", "lab_history.jsonl")

def _rows():
    try: return [json.loads(l) for l in open(F) if l.strip()]
    except Exception: return []

def add(disease, smiles, ic50_nM=None, cc50_nM=None, predicted_pic50=None, cost_usd=None, days=None, note="", who=None, outcome=None):
    """log one tested molecule. outcome auto = PASS if IC50 <= 1 uM and selectivity (CC50/IC50) >= 10, else FAILED"""
    if outcome is None:
        ok = ic50_nM is not None and ic50_nM <= 1000 and (cc50_nM is None or cc50_nM / ic50_nM >= 10); outcome = "PASS" if ok else "FAILED"
    r = {"ts": time.strftime("%Y-%m-%d %H:%M:%S"), "who": who or getpass.getuser(), "disease": disease, "smiles": smiles, "ic50_nM": ic50_nM, "cc50_nM": cc50_nM,
         "predicted_pic50": predicted_pic50, "outcome": outcome, "cost_usd": cost_usd, "days": days, "note": note}
    os.makedirs(os.path.dirname(F), exist_ok=True)
    with open(F, "a") as f: f.write(json.dumps(r) + "\n")
    return r

def load(disease=None): return [r for r in _rows() if disease in (None, r["disease"])]

def training_overlay(disease):
    """{smiles: pIC50} measured in YOUR lab (replaces ChEMBL values for the same structure)"""
    return {r["smiles"]: 9 - math.log10(r["ic50_nM"]) for r in load(disease) if r.get("ic50_nM")}

def failed(disease): return [r["smiles"] for r in load(disease) if r["outcome"] == "FAILED"]

def report(disease=None):
    R = load(disease); print(f"\nLAB HISTORY{' for ' + disease if disease else ''}: {len(R)} tests logged in {F}")
    if not R: return
    P = [r for r in R if r["outcome"] == "PASS"]; print(f"  hit rate (PASS = IC50<=1uM & selective): {len(P)}/{len(R)} = {100 * len(P) / len(R):.0f}%")
    pr = [(r["predicted_pic50"], 9 - math.log10(r["ic50_nM"])) for r in R if r.get("predicted_pic50") is not None and r.get("ic50_nM")]
    if len(pr) >= 3:
        e = [p - m for p, m in pr]; rmse = math.sqrt(sum(x * x for x in e) / len(e)); w3 = 100 * sum(abs(x) <= math.log10(3) for x in e) / len(e)
        print(f"  PROSPECTIVE prediction error: RMSE {rmse:.2f} log units (~{10 ** rmse:.1f}-fold) | within 3-fold {w3:.0f}% (n={len(pr)})")
        if len(pr) >= 8:
            h = len(pr) // 2; r1 = math.sqrt(sum(x * x for x in e[:h]) / h); r2 = math.sqrt(sum(x * x for x in e[h:]) / (len(e) - h))
            print(f"  learning trend: first half RMSE {r1:.2f} -> second half {r2:.2f} ({'improving' if r2 < r1 else 'NOT improving'})")
    else: print("  (log predicted_pic50 together with each measurement to track true prediction accuracy; need >= 3)")
    c = [r["cost_usd"] for r in R if r.get("cost_usd")]; d = [r["days"] for r in R if r.get("days")]
    if c: print(f"  cost: total ${sum(c):,.0f} | ${sum(c) / len(c):,.0f} per test")
    if d: print(f"  time per loop: {sum(d) / len(d):.1f} days average")

def cli(a):
    """lab add <disease> <smiles> <ic50_nM> [cc50_nM|-] [predicted_pic50|-] [cost_usd|-] [days|-] [note...] | lab report [disease] | lab list [disease]"""
    if not a or a[0] == "report": return report(a[1] if len(a) > 1 else None)
    if a[0] == "list":
        for r in load(a[1] if len(a) > 1 else None): print(f"{r['ts']} {r['who']} {r['disease']} {r['outcome']:<6} IC50={r['ic50_nM']} {r['smiles']}")
        return
    if a[0] == "add" and len(a) >= 4:
        v = lambda i: float(a[i]) if len(a) > i and a[i] != "-" else None
        r = add(a[1], a[2], v(3), v(4), v(5), v(6), v(7), " ".join(a[8:])); print("logged:", r["outcome"], r["smiles"]); return
    print(cli.__doc__)
