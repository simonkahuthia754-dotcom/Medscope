"""state.py - remembers what has ACTUALLY been verified on this machine (feeds scorecard.py). Only real runs write here."""
import json, os, time
F = os.path.join(os.path.expanduser("~"), ".medscope", "state.json")

def load():
    try: return json.load(open(F))
    except Exception: return {}

def record(key, value):
    s = load(); s[key] = value; s["_t_" + key] = time.strftime("%Y-%m-%d %H:%M")
    try: os.makedirs(os.path.dirname(F), exist_ok=True); json.dump(s, open(F, "w"))
    except Exception: pass
