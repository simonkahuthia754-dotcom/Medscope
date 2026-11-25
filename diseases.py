"""diseases.py - PER-DISEASE targeting. No more lumped 'virus' / 'cancer' / 'bacteria'.
Each disease = its own molecular target, its own GOOD list (known effective drugs for THAT target), its own BAD list (drugs for
other diseases), its own validation (AUC + separation) and its own human PK context. ChEMBL target ids are resolved by name at
runtime unless fixed here - ALWAYS read the printed target name. `pdb_ok` = I am confident that structure matches the target;
others are best guesses and are verified at run time by the redocking RMSD gate (<= 2 A) before any docking score is trusted.
Add your own diseases in ~/.medscope/diseases.json (same fields)."""
import json, math, os
import sources

BAD_COMMON = ["dolutegravir", "nirmatrelvir", "oseltamivir", "vancomycin", "isoniazid", "imatinib", "celecoxib", "fluconazole", "ciprofloxacin", "erlotinib", "sofosbuvir"]

D = {
    "hiv": dict(name="HIV-1 infection", target="HIV-1 integrase", chembl=None, pdb="1QS4", pdb_ok=False, symbol=None, cyp="CYP3A4",
                good=["dolutegravir", "raltegravir", "elvitegravir", "bictegravir", "cabotegravir"]),
    "covid": dict(name="COVID-19", target="SARS-CoV-2 3C-like proteinase (Mpro)", chembl=None, pdb="6LU7", pdb_ok=True, symbol=None, cyp="CYP3A4",
                  good=["nirmatrelvir", "ensitrelvir"]),
    "flu": dict(name="Influenza", target="Influenza neuraminidase", chembl=None, pdb="3B7E", pdb_ok=False, symbol=None, cyp=None,
                good=["oseltamivir", "zanamivir", "peramivir", "laninamivir"]),
    "hepc": dict(name="Hepatitis C", target="HCV NS5B RNA-dependent RNA polymerase", chembl=None, pdb=None, pdb_ok=False, symbol=None, cyp=None,
                 good=["sofosbuvir", "dasabuvir"]),
    "mrsa": dict(name="MRSA (ribosome-targeting)", target="Bacterial 50S ribosome / 23S rRNA (oxazolidinones)", chembl=None, pdb=None, pdb_ok=False, symbol=None, cyp=None,
                 good=["linezolid", "tedizolid"], note="vancomycin acts on the cell wall (D-Ala-D-Ala), NOT the ribosome - it belongs to a separate cell-wall entry"),
    "tb": dict(name="Tuberculosis (rifamycin target)", target="DNA-directed RNA polymerase beta subunit (Mycobacterium tuberculosis)", chembl=None, pdb="5UHB", pdb_ok=False,
               symbol=None, cyp="CYP3A4", good=["rifampicin", "rifapentine", "rifabutin"], note="isoniazid targets InhA (mycolic-acid synthesis), not RNA polymerase"),
    "ecoli": dict(name="E. coli / Gram-negative infection", target="DNA gyrase", chembl=None, pdb="2XCT", pdb_ok=False, symbol=None, cyp=None,
                  good=["ciprofloxacin", "levofloxacin", "moxifloxacin", "norfloxacin"]),
    "lung-cancer": dict(name="EGFR-mutant lung cancer", target="Epidermal growth factor receptor", chembl="CHEMBL203", pdb="1M17", pdb_ok=True, symbol="EGFR", cyp="CYP3A4",
                        good=["erlotinib", "gefitinib", "osimertinib", "afatinib"], skip=["imatinib"]),
    "breast-cancer": dict(name="HER2+ breast cancer", target="Receptor tyrosine-protein kinase erbB-2", chembl="CHEMBL1824", pdb="3PP0", pdb_ok=False, symbol="ERBB2", cyp="CYP3A4",
                          good=["lapatinib", "neratinib", "tucatinib"], skip=["imatinib", "erlotinib"]),
    "blood-cancer": dict(name="CML (BCR-ABL)", target="Tyrosine-protein kinase ABL", chembl="CHEMBL1862", pdb="1IEP", pdb_ok=True, symbol="ABL1", cyp="CYP3A4",
                         good=["imatinib", "dasatinib", "nilotinib", "ponatinib", "bosutinib"], skip=["erlotinib"]),
    "malaria": dict(name="Malaria (P. falciparum)", target="Plasmodium falciparum P-type ATPase 4 (PfATP4)", chembl=None, pdb=None, pdb_ok=False, symbol=None, cyp=None,
                    good=["cipargamin"], note="measured PfATP4 data in ChEMBL is sparse: expect LOW confidence potency predictions"),
    "fungi": dict(name="Fungal infection", target="Lanosterol 14-alpha demethylase (CYP51)", chembl=None, pdb=None, pdb_ok=False, symbol=None, cyp="CYP3A4",
                  good=["fluconazole", "voriconazole", "itraconazole", "posaconazole"]),
    "inflammation": dict(name="Inflammation / pain", target="Cyclooxygenase-2 (PTGS2)", chembl="CHEMBL230", pdb="3LN1", pdb_ok=True, symbol="PTGS2", cyp="CYP2C9",
                         good=["celecoxib", "etoricoxib", "valdecoxib"], skip=[]),
}
BUILTIN = set(D)
try:
    D.update(json.load(open(os.path.join(os.path.expanduser("~"), ".medscope", "diseases.json"))))
except Exception: pass

def get(key):
    k = key.lower().replace("_", "-").replace(" ", "-")
    if k not in D: raise KeyError(f"unknown disease '{key}'. Available: {', '.join(sorted(D))}")
    return {**D[k], "key": k}

def bad_for(key):
    d = get(key); ex = set(d["good"]) | set(d.get("skip", []))
    return [b for b in BAD_COMMON if b not in ex]

def show():
    print(f"{'key':<14}{'target':<58}PDB")
    for k, v in D.items(): print(f"{k:<14}{v['target'][:56]:<58}{(v['pdb'] or '-') + ('' if v.get('pdb_ok') or not v['pdb'] else '?')}")
    print("\n'?' = structure id not verified by me; the redocking gate checks it at run time. Notes:")
    for k, v in D.items():
        if v.get("note"): print(f"  {k}: {v['note']}")

def resolve_chembl(d):
    """ChEMBL target id for a disease (fixed id, else name search). Returns (id, preferred name)"""
    if d.get("chembl"): return d["chembl"], d["target"]
    r = sources.J(f"https://www.ebi.ac.uk/chembl/api/data/target/search.json?q={sources.Q(d['target'])}&limit=10")["targets"]
    for t in r:
        if t.get("target_type") in ("SINGLE PROTEIN", "PROTEIN COMPLEX"): return t["target_chembl_id"], t["pref_name"]
    raise RuntimeError(f"no ChEMBL target found for '{d['target']}' - set 'chembl' for this disease in ~/.medscope/diseases.json")

def auc(good, bad):
    """ROC AUC = P(random good drug scores above random bad drug); 0.5 = no discrimination"""
    if not good or not bad: return None
    w = sum((g > b) + .5 * (g == b) for g in good for b in bad); return w / (len(good) * len(bad))

def separation(good, bad):
    """standardised mean difference (Cohen's d) between good and bad scores"""
    if len(good) < 2 or len(bad) < 2: return None
    mg, mb = sum(good) / len(good), sum(bad) / len(bad)
    vg = sum((x - mg) ** 2 for x in good) / (len(good) - 1); vb = sum((x - mb) ** 2 for x in bad) / (len(bad) - 1)
    sp = math.sqrt((vg + vb) / 2) or 1e-9; return (mg - mb) / sp

def check(key, model=None, T=None, verbose=True):
    """per-disease validation: score known GOOD vs BAD drugs for THIS target. Drugs are removed from the training data when scored (no leakage)."""
    d = get(key); cid = None
    if model is None:
        import potency
        cid, nm = resolve_chembl(d); print(f"[{d['key']}] ChEMBL target: {cid} = {nm}"); model = potency.Model(potency.fetch(cid)); print("  " + model.summary())
    res = {"good": {}, "bad": {}}
    for grp, names in (("good", d["good"]), ("bad", bad_for(key))):
        for n in names:
            try:
                smi = sources.pubchem(n)["smiles"]; p = model.predict(smi, exclude_self=True)
                if p: res[grp][n] = p["pic50"] if p["in_domain"] else None
            except Exception as e:
                if verbose: print(f"  ! {n}: {e}")
    g = [v for v in res["good"].values() if v is not None]; b = [v for v in res["bad"].values() if v is not None]
    res["auc"], res["sep"] = auc(g, b), separation(g, b)
    if verbose:
        f = lambda x: "n/a" if x is None else f"{x:.2f}"
        print(f"  GOOD pIC50: " + ", ".join(f"{n}={f(v)}" for n, v in res["good"].items()))
        print(f"  BAD  pIC50: " + ", ".join(f"{n}={f(v)}" for n, v in res["bad"].items()) + "   (n/a = outside model's domain)")
        print(f"  per-disease validation: AUC {f(res['auc'])} | separation {f(res['sep'])} (n good {len(g)}, bad {len(b)}) - AUC>0.8 is good; tiny n = noisy")
    return res
