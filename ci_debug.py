import sys, json, traceback
def sect(n): print("\n##", n, flush=True)
sect("OT BaselineExpression")
try:
    import humanbio
    q = 'query($n:String!){__type(name:$n){fields{name type{name kind ofType{name kind ofType{name kind ofType{name}}}}}}}'
    print("rows type:", [f["type"] for f in humanbio._gql(q, {"n": "BaselineExpression"})["data"]["__type"]["fields"] if f["name"] == "rows"])
    for tn in ("ExpressionRow", "BaselineExpressionRow", "TissueExpression", "BaselineTissueExpression", "Expression2"):
        r = humanbio._gql(q, {"n": tn})["data"]["__type"]
        if r: print(tn, [(f["name"], f["type"]["name"] or (f["type"]["ofType"] or {}).get("name")) for f in r["fields"]])
except BaseException: print(traceback.format_exc()[-900:])
sect("docking")
try:
    try: import vina, meeko; print("vina/meeko import ok")
    except BaseException: print(traceback.format_exc()[-600:])
    import docking
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    T = docking.prepare("1M17"); native = docking.native_smiles(T["lig"]); print("native smiles", native, "center", T["center"], "size", T["size"])
    crystal = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(native), Chem.MolFromPDBFile(T["crystal"])); print("crystal atoms", crystal.GetNumHeavyAtoms(), "centroid", crystal.GetConformer().GetPositions().mean(0))
    for exh, seed in ((16, 7), (32, 7), (32, 11), (64, 3)):
        s = docking._lig(native, seed); v = docking.Vina(sf_name="vina", seed=seed); v.set_receptor(T["rec"]); v.set_ligand_from_string(s)
        v.compute_vina_maps(center=list(T["center"]), box_size=[T["size"]] * 3); v.dock(exhaustiveness=exh, n_poses=10); en = v.energies(n_poses=10)
        mols = docking.RDKitMolCreate.from_pdbqt_mol(docking._pq(v.poses(n_poses=10))); out = []
        for i, m in enumerate(mols):
            try: out.append((round(float(en[i][0]), 2), round(rdMolAlign.CalcRMS(Chem.RemoveHs(m), Chem.RemoveHs(crystal)), 2)))
            except Exception as e: out.append(("err", str(e)[:80]))
        print(f"exh {exh} seed {seed}: (score, RMSD A) =", out, flush=True)
except BaseException: print(traceback.format_exc()[-1500:])
