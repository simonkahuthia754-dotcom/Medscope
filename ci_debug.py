import sys, json, traceback, math
def sect(n): print("\n##", n, flush=True)
sect("OT BaselineExpression")
try:
    import humanbio
    for tn in ("BaselineExpression", "Expression", "TissueSpecificity"):
        r = humanbio._gql('query($n:String!){__type(name:$n){fields{name type{name kind ofType{name kind ofType{name}}}}}}', {"n": tn})["data"]["__type"]
        print(tn, [(f["name"], f["type"]["name"] or (f["type"]["ofType"] or {}).get("name")) for f in (r or {}).get("fields", [])])
except BaseException: print(traceback.format_exc()[-900:])
sect("docking")
try:
    import docking
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdMolAlign
    T = docking.prepare("1M17"); native = docking.native_smiles(T["lig"]); print("native smiles", native, "center", T["center"], "size", T["size"])
    crystal = AllChem.AssignBondOrdersFromTemplate(Chem.MolFromSmiles(native), Chem.MolFromPDBFile(T["crystal"])); print("crystal atoms", crystal.GetNumHeavyAtoms())
    cc = crystal.GetConformer().GetPositions(); print("crystal centroid", cc.mean(0))
    for exh, seed in ((16, 7), (32, 7), (32, 11), (64, 3)):
        s = docking._lig(native, seed); v = docking.Vina(sf_name="vina", seed=seed); v.set_receptor(T["rec"]); v.set_ligand_from_string(s)
        v.compute_vina_maps(center=list(T["center"]), box_size=[T["size"]] * 3); v.dock(exhaustiveness=exh, n_poses=10); en = v.energies(n_poses=10)
        out = []
        for i in range(min(10, len(en))):
            pq = v.poses(n_poses=10)
            break
        pm = docking._pq(v.poses(n_poses=10)); mols = docking.RDKitMolCreate.from_pdbqt_mol(pm)
        for i, m in enumerate(mols):
            try: out.append((round(float(en[i][0]), 2), round(rdMolAlign.CalcRMS(Chem.RemoveHs(m), Chem.RemoveHs(crystal)), 2), round(float(m.GetConformer().GetPositions().mean(0)[0]), 1)))
            except Exception as e: out.append(("err", str(e)[:60]))
        print(f"exh {exh} seed {seed}: (score, RMSD A, pose x-centroid) =", out, flush=True)
except BaseException: print(traceback.format_exc()[-1500:])
