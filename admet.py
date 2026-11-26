"""admet.py - ADMET-AI (Chemprop-RDKit models trained on ~41 TDC/ChEMBL measured ADMET datasets) -> PBPK parameters + safety flags.
Install: pip install admet-ai   (PC/Colab/server; needs torch; not for QPython). Predictions are model outputs with typical
2-3x error on clearance/Vss; use `validate.py --admet` to measure accuracy on known drugs, and replace with measured data when available."""
try:
    from admet_ai import ADMETModel
    HAVE = True
except Exception:
    HAVE = False
_M = None
clamp = lambda x, a, b: max(a, min(b, x))

def predict(smiles):
    """list of SMILES -> {smiles: {property: value}} (model loads once)"""
    global _M
    if not HAVE: raise RuntimeError("ADMET-AI not installed: pip install admet-ai")
    _M = _M or ADMETModel(); r = _M.predict(smiles=list(smiles))
    return {s: dict(v) for s, v in r.to_dict("index").items()} if hasattr(r, "to_dict") else {smiles[0]: dict(r)}

def get(pr, *subs):
    """find a property by substring(s), ignoring percentile columns"""
    for k, v in pr.items():
        kl = k.lower()
        if "percentile" not in kl and all(s.lower() in kl for s in subs):
            try: return float(v)
            except Exception: return None
    return None

def overrides(pr):
    """ADMET predictions -> PBPK parameter overrides (+ provenance)"""
    o = {"_src": {}}; ppb = get(pr, "PPBR"); clh = get(pr, "Clearance_Hepatocyte"); vss = get(pr, "VDss"); hia = get(pr, "HIA")
    if ppb is not None: o["fup"] = clamp(1 - ppb / 100, .005, 1); o["_src"]["fup"] = "ADMET-AI (PPBR)"
    if clh is not None and clh > 0: o["clint_hep"] = clh; o["_src"]["clint_hep"] = "ADMET-AI (hepatocyte CL)"
    if vss is not None and vss > 0: o["vss"] = vss; o["_src"]["vss"] = "ADMET-AI (VDss)"
    th = get(pr, "Half_Life")
    if th is not None and th > 0: o["t_half"] = th; o["_src"]["t_half"] = "ADMET-AI (Half_Life_Obach)"
    if hia is not None: o["fa"] = clamp(.2 + .8 * hia, .2, 1); o["_src"]["fa"] = "ADMET-AI (HIA)"
    return o

FLAGS = [("hERG", "heart: hERG/QT liability"), ("AMES", "genotoxicity: Ames mutagenic"), ("DILI", "liver: drug-induced liver injury"),
         ("Carcinogens", "carcinogenicity alert"), ("ClinTox", "clinical-trial toxicity"), ("Skin_Reaction", "skin: sensitisation"),
         ("CYP3A4_Veith", "DDI: CYP3A4 inhibitor"), ("CYP2D6_Veith", "DDI: CYP2D6 inhibitor"), ("CYP2C9_Veith", "DDI: CYP2C9 inhibitor"), ("Pgp", "P-gp interaction")]

def flags(pr, thr=.5):
    """probability-type endpoints above thr -> human-readable flags with the model probability"""
    out = []
    for key, text in FLAGS:
        v = get(pr, key)
        if v is not None and v >= thr: out.append(f"{text} (p={v:.2f})")
    return out

def cyp_substrate(pr, thr=.5):
    """which CYP enzymes likely metabolise it (drives genetic-polymorphism scenarios)"""
    return [e for e in ("CYP2D6", "CYP3A4", "CYP2C9", "CYP2C19") if (get(pr, e, "Substrate") or 0) >= thr]
