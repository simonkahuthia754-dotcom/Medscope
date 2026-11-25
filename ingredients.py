"""ingredients.py - drugs are made of INGREDIENTS. Brand/product -> active ingredients (RxNorm), each ingredient -> structure (PubChem)
-> human PBPK + ADMET flags; combinations -> overlap checks.
  drug <name>              any existing drug or product: ingredients, structures, properties, human PBPK, flags
  combo <name> <name> ...  fixed-dose combinations / regimens: ingredients, organ-exposure overlap, shared toxicity flags, CYP interaction risk
NOT modelled: excipients (binders, coatings, solubilisers) - they change dissolution/absorption and need formulation studies;
SYNERGY / antagonism - cannot be predicted here, needs combination assays (e.g. checkerboard)."""
import sys
import pbpk, sources, state
from chemicals import exposure

RX = "https://rxnav.nlm.nih.gov/REST"

def ingredients(product):
    """product/brand/generic name -> list of active-ingredient names (RxNorm). Falls back to the name itself."""
    try:
        r = sources.J(f"{RX}/rxcui.json?name={sources.Q(product)}&search=2")["idGroup"].get("rxnormId") or []
        if not r: return [product]
        names = []
        for g in sources.J(f"{RX}/rxcui/{r[0]}/related.json?tty=IN+MIN+PIN")["relatedGroup"].get("conceptGroup", []):
            names += [c["name"] for c in g.get("conceptProperties", [])]
        state.record("rxnorm_ok", True); return sorted(set(n.lower() for n in names)) or [product]
    except Exception: return [product]

def profile(name, dose_mg=100.0):
    """structure + properties + human PBPK + ADMET flags for ONE ingredient"""
    pc = sources.pubchem(name); state.record("pubchem_ok", True); p = {"logp": pc["logp"], "mw": pc["mw"], "tpsa": pc.get("tpsa") or 80}
    out = {"name": name, "smiles": pc["smiles"], "props": p, "flags": [], "cyp_sub": [], "ovr": None, "admet": {}}
    try:
        import admet
        if admet.HAVE and pc["smiles"]:
            a = admet.predict([pc["smiles"]])[pc["smiles"]]; out.update(ovr=admet.overrides(a), flags=admet.flags(a), cyp_sub=admet.cyp_substrate(a), admet=a)
    except Exception: pass
    po, iv, F = pbpk.bioavailability(p, "human", dose_mg, 48.0, ovr=out["ovr"]); out.update(po=po, F=F, exp=exposure(po)); return out

def show_profile(x):
    th = x["po"]["t_half"]; print(f"\n- {x['name']}  {x['smiles']}\n  MW {x['props']['mw']:.0f} logP {x['props']['logp']:.1f} TPSA {x['props']['tpsa']:.0f} | human 100 mg po: Cmax {x['po']['cmax']:.2f} mg/L, "
          f"t1/2 {'n/a' if not th else f'{th:.1f}'} h, F {x['F']:.2f} ({'ADMET-AI' if x['ovr'] else 'rule-based, LOW confidence'})")
    print("  organ:plasma AUC: " + ", ".join(f"{o} {v[1]:.1f}" for o, v in x["exp"].items()))
    for f in x["flags"]: print("  ! " + f)

def drug(name):
    ing = ingredients(name); print(f"'{name}' -> active ingredient(s): {', '.join(ing)}")
    for i in ing:
        try: show_profile(profile(i))
        except Exception as e: print(f"\n- {i}: lookup failed ({e})")
    print("\nExcipients are not modelled. Predictions: see validate.py for the simulator's measured accuracy.")

def combo(names):
    ing = []
    for n in names: ing += [i for i in ingredients(n) if i not in ing]
    print("Combination ingredients: " + ", ".join(ing)); P = []
    for i in ing:
        try: x = profile(i); show_profile(x); P.append(x)
        except Exception as e: print(f"\n- {i}: lookup failed ({e})")
    if len(P) < 2: print("Need at least two resolvable ingredients."); return
    print("\nCOMBINATION CHECKS")
    for o in ("liver", "kidney", "heart", "brain", "gut"):
        hi = [x["name"] for x in P if x["exp"][o][1] >= 1.0]
        if len(hi) > 1: print(f"  ! organ overlap: {', '.join(hi)} all concentrate in the {o} (tissue:plasma >= 1): additive burden possible")
    for k in ("hERG", "liver injury", "mutagenic", "carcinogen"):
        hit = [x["name"] for x in P if any(k.lower() in f.lower() for f in x["flags"])]
        if len(hit) > 1: print(f"  ! shared liability '{k}': {', '.join(hit)}")
    done = False
    for a in P:
        for e in [f.split("DDI: ")[1].split(" ")[0] for f in a["flags"] if "DDI: " in f]:
            for b in P:
                if b is not a and e in b["cyp_sub"]: print(f"  ! DDI: {a['name']} inhibits {e}; {b['name']} is a {e} substrate -> {b['name']} exposure may rise"); done = True
    if any(x["ovr"] for x in P): state.record("ddi_ok", True)
    else: print("  (ADMET-AI not available: CYP interaction and shared-toxicity checks skipped)")
    print("  Synergy/antagonism cannot be predicted here (needs combination assays). Excipients not modelled.")

def cli(cmd, a):
    if cmd == "drug": drug(" ".join(a))
    else: combo(a)

if __name__ == "__main__": cli(sys.argv[1], sys.argv[2:])
