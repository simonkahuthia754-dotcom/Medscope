import sys, traceback
def sect(n): print("\n##", n, flush=True)
sect("expression")
try:
    import humanbio
    for s in ("EGFR", "CCR5", "ABL1"): print(s, humanbio.target_expression(s))
except BaseException: print(traceback.format_exc()[-900:])
sect("pk")
try:
    import validate; validate.run(False); validate.run(True)
except BaseException: print(traceback.format_exc()[-1500:])
sect("docking")
try:
    import docking
    T = docking.prepare("1M17"); print("consensus redock (score, RMSD, cluster, runs):", docking.redock_rmsd(T))
except BaseException: print(traceback.format_exc()[-1500:])
