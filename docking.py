"""docking.py - AutoDock Vina docking against REAL human/viral protein structures from the RCSB PDB.
Needs: pip install vina meeko rdkit  + Open Babel (`obabel` on PATH) for receptor prep. PC/Colab/server only.
Trust gates (reported, never skipped): (1) REDOCKING - re-dock the crystal ligand; RMSD <= 2 A means the box/protocol reproduces reality.
(2) CALIBRATION - dock known actives with measured pIC50; Spearman rho tells you how much this target's Vina scores mean anything.
Vina scores rank poses/compounds roughly; they are NOT binding affinities. Use them to cross-check the QSAR potency, not replace it."""
import math, os, shutil, subprocess, urllib.request
try:
    from vina import Vina
    from meeko import MoleculePreparation, PDBQTWriterLegacy, PDBQTMolecule, RDKitMolCreate
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    HAVE = True
except Exception:
    HAVE = False
HOME = os.path.join(os.path.expanduser("~"), ".medscope", "pdb")

def need():
    if not HAVE: raise RuntimeError("docking needs: pip install vina meeko rdkit (and Open Babel: conda install -c conda-forge openbabel)")
    if not shutil.which("obabel"): raise RuntimeError("Open Babel `obabel` not found on PATH (needed to prepare the receptor)")

def fetch_pdb(pid):
    os.makedirs(HOME, exist_ok=True); p = os.path.join(HOME, pid.upper() + ".pdb")
    if not os.path.exists(p): urllib.request.urlretrieve(f"https://files.rcsb.org/download/{pid.upper()}.pdb", p)
    return p

SOLV = {"HOH", "WAT", "SO4", "PO4", "GOL", "EDO", "PEG", "PGE", "MG", "ZN", "CA", "NA", "CL", "K", "MN", "FE", "NI", "CO", "CU", "CD", "ACT", "DMS", "FMT", "TRS", "MES", "EPE", "BME", "IOD", "BR", "NO3", "ACY", "IPA", "MPD", "PG4", "1PE", "P6G"}

def _pq(s):
    try: return PDBQTMolecule(s, is_dock=True)   # meeko <=0.5
    except TypeError: return PDBQTMolecule(s)         # meeko >=0.6 auto-detects Vina output

def prepare(pdb_id, lig_id=None):
    """-> dict(rec, center, size, crystal, lig): protein only, hydrogens at pH 7.4, box centred on the co-crystal ligand
    (lig_id None = auto-pick the largest non-solvent HETATM group)"""
    need(); src = fetch_pdb(pdb_id); prot, groups = [], {}
    for l in open(src):
        if l.startswith("ATOM"):
            if l[16] in " A": prot.append(l)    # keep a single alternate location: duplicated overlapping atoms would clash and block the pocket
        elif l.startswith("HETATM") and l[16] in " A" and l[17:20].strip() not in SOLV and l[76:78].strip() != "H": groups.setdefault((l[17:20].strip(), l[21]), []).append(l)
    cand = {k: v for k, v in groups.items() if lig_id is None or k[0] == lig_id}
    if not cand: raise RuntimeError(f"no ligand{' ' + lig_id if lig_id else ''} found in {pdb_id}")
    (lig_id, chain), lig = max(cand.items(), key=lambda kv: len(kv[1]))
    xyz = [(float(l[30:38]), float(l[38:46]), float(l[46:54])) for l in lig]
    c = tuple(sum(p[i] for p in xyz) / len(xyz) for i in range(3)); ext = max(max(p[i] for p in xyz) - min(p[i] for p in xyz) for i in range(3))
    d = os.path.join(HOME, f"{pdb_id.upper()}_{lig_id}"); os.makedirs(d, exist_ok=True)
    open(d + "/protein.pdb", "w").writelines(prot); open(d + "/ligand.pdb", "w").writelines(lig + ["END\n"])
    subprocess.run(["obabel", d + "/protein.pdb", "-O", d + "/protein.pdbqt", "-xr", "-p", "7.4"], check=True, capture_output=True)
    return {"rec": d + "/protein.pdbqt", "center": c, "size": max(20.0, ext + 12.0), "crystal": d + "/ligand.pdb", "lig": lig_id}

def native_smiles(code):
    """SMILES of a PDB ligand code from RCSB (template for bond orders in the redocking check)"""
    import sources
    r = sources.J(f"https://data.rcsb.org/rest/v1/core/chemcomp/{code}").get("rcsb_chem_comp_descriptor", {})
    L = {k.lower().replace("_", ""): v for k, v in r.items()}; return L.get("smilesstereo") or L.get("smiles")

def _lig(smi, seed=7):
    m = Chem.AddHs(Chem.MolFromSmiles(smi))
    if AllChem.EmbedMolecule(m, randomSeed=seed) != 0: return None
    try: AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    except Exception: pass
    s, ok, err = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m)[0]); return s if ok else None

def dock(T, smi, exh=16, seed=7):
    """-> (best Vina score kcal/mol, pose pdbqt string) or None"""
    s = _lig(smi, seed)
    if not s: return None
    v = Vina(sf_name="vina", seed=seed); v.set_receptor(T["rec"]); v.set_ligand_from_string(s)
    v.compute_vina_maps(center=list(T["center"]), box_size=[T["size"]] * 3); v.dock(exhaustiveness=exh, n_poses=5)
    return float(v.energies(n_poses=1)[0][0]), v.poses(n_poses=1)

def redock_rmsd(T, native=None):
    """redock the crystal ligand; returns (score, heavy-atom RMSD to crystal pose in A). <=2 A = protocol trustworthy"""
    native = native or native_smiles(T["lig"]); r = dock(T, native)
    if not r: return None
    crystal = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(native), Chem.MolFromPDBFile(T["crystal"]))
    pose = RDKitMolCreate.from_pdbqt_mol(_pq(r[1]))[0]
    return r[0], rdMolAlign.CalcRMS(Chem.RemoveHs(pose), Chem.RemoveHs(crystal))

def _rank(x):
    o = sorted(range(len(x)), key=lambda i: x[i]); r = [0.0] * len(x)
    for k, i in enumerate(o): r[i] = k
    return r

def calibrate(T, data, n=24):
    """dock n known actives spread over the measured pIC50 range; Spearman rho(-score, pIC50) = how informative Vina is for this target"""
    items = sorted(data.items(), key=lambda kv: kv[1]); step = max(1, len(items) // n); sc, py = [], []
    for s, y in items[::step][:n]:
        r = dock(T, s, exh=8)
        if r: sc.append(-r[0]); py.append(y)
    if len(sc) < 8: return {"rho": None, "scores": sc, "n": len(sc)}
    a, b = _rank(sc), _rank(py); ma, mb = sum(a) / len(a), sum(b) / len(b)
    rho = sum((x - ma) * (y - mb) for x, y in zip(a, b)) / math.sqrt(sum((x - ma) ** 2 for x in a) * sum((y - mb) ** 2 for y in b))
    return {"rho": rho, "scores": sorted(sc), "n": len(sc)}

def percentile(cal, e):
    """share of calibration actives that this candidate scores better than (0-100)"""
    s = cal["scores"]; return 100 * sum(x <= -e for x in s) / len(s) if s else None

def dock_consistent(T, smi, runs=3):
    """repeat docking with different seeds: (mean score, sd, reliable). Reliable = sd <= 0.5 kcal/mol (a stable answer, not a lucky run)"""
    sc = [r[0] for r in (dock(T, smi, seed=sd) for sd in (7, 11, 13, 17, 19)[:runs]) if r]
    if not sc: return None
    m = sum(sc) / len(sc); sd = math.sqrt(sum((x - m) ** 2 for x in sc) / max(1, len(sc) - 1)); return m, sd, len(sc) >= 2 and sd <= .5
