"""pkbench.py - human PK benchmark sets (approximate literature: MW, logP, TPSA, t1/2 h, F oral, ionisation class at pH 7.4, SMILES).
DEV = the 14 drugs used to diagnose and tune the rule-based model. TEST = 20 further drugs written down BEFORE the ionisation-class rules were added;
the model has never been tuned on them, so TEST is the honest estimate of accuracy on new molecules. Values approximate - verify before citing."""
# name: (MW, logP, TPSA, t_half_h, F, ion, SMILES)
DEV = {
 "caffeine": (194.2, -0.07, 58.4, 5.0, 1.00, "neutral", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"), "theophylline": (180.2, -0.02, 69.3, 8.0, 1.00, "neutral", "Cn1c2nc[nH]c2c(=O)n(C)c1=O"),
 "propranolol": (259.3, 3.48, 41.5, 4.0, 0.26, "base", "CC(C)NCC(O)COc1cccc2ccccc12"), "diazepam": (284.7, 2.82, 32.7, 43.0, 0.93, "neutral", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
 "warfarin": (308.3, 2.70, 63.6, 37.0, 0.99, "acid", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"), "ibuprofen": (206.3, 3.97, 37.3, 2.0, 0.90, "acid", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
 "metformin": (129.2, -1.40, 91.5, 6.2, 0.55, "base", "CN(C)C(=N)NC(=N)N"), "atenolol": (266.3, 0.16, 84.6, 6.5, 0.55, "base", "CC(C)NCC(O)COc1ccc(CC(N)=O)cc1"),
 "zidovudine": (267.2, 0.05, 125.8, 1.1, 0.64, "neutral", "CC1=CN(C(=O)NC1=O)C2CC(C(O2)CO)N=[N+]=[N-]"), "ciprofloxacin": (331.3, 0.28, 74.6, 4.0, 0.70, "zwitterion", "OC(=O)C1=CN(C2CC2)c2cc(N3CCNCC3)c(F)cc2C1=O"),
 "linezolid": (337.4, 0.55, 71.1, 5.0, 1.00, "neutral", "CC(=O)NC[C@H]1CN(c2ccc(N3CCOCC3)c(F)c2)C(=O)O1"), "imatinib": (493.6, 3.50, 86.3, 18.0, 0.98, "base", "Cc1ccc(NC(=O)c2ccc(CN3CCN(C)CC3)cc2)cc1Nc1nccc(-c2cccnc2)n1"),
 "erlotinib": (393.4, 3.30, 74.7, 36.0, 0.60, "base", "COCCOc1cc2ncnc(Nc3cccc(C#C)c3)c2cc1OCCOC"), "chloroquine": (319.9, 4.60, 28.2, 1000.0, 0.80, "base", "CCN(CC)CCCC(C)Nc1ccnc2cc(Cl)ccc12")}
TEST = {
 "acetaminophen": (151.2, 0.46, 49.3, 2.5, 0.88, "neutral", "CC(=O)Nc1ccc(O)cc1"), "midazolam": (325.8, 3.30, 30.2, 2.5, 0.44, "neutral", "Cc1ncc2n1-c1ccc(Cl)cc1C(c1ccccc1F)=NC2"),
 "verapamil": (454.6, 3.80, 63.9, 5.0, 0.22, "base", "COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC"), "phenytoin": (252.3, 2.50, 58.2, 22.0, 0.90, "neutral", "O=C1NC(=O)C(c2ccccc2)(c2ccccc2)N1"),
 "lidocaine": (234.3, 2.40, 32.3, 1.8, 0.35, "base", "CCN(CC)CC(=O)Nc1c(C)cccc1C"), "naproxen": (230.3, 3.20, 46.5, 14.0, 0.95, "acid", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
 "omeprazole": (345.4, 2.20, 77.1, 1.0, 0.60, "neutral", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1"), "diclofenac": (296.1, 4.40, 49.3, 2.0, 0.55, "acid", "OC(=O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
 "amlodipine": (408.9, 3.00, 99.6, 40.0, 0.64, "base", "CCOC(=O)C1=C(COCCN)NC(C)=C(C(=O)OC)C1c1ccccc1Cl"), "sildenafil": (474.6, 2.70, 113.0, 4.0, 0.40, "base", "CCCc1nn(C)c2c1nc([nH]c2=O)-c1cc(ccc1OCC)S(=O)(=O)N1CCN(C)CC1"),
 "losartan": (422.9, 4.30, 92.5, 2.0, 0.33, "acid", "CCCCc1nc(Cl)c(CO)n1Cc1ccc(cc1)-c1ccccc1-c1nn[nH]n1"), "fluconazole": (306.3, 0.50, 81.7, 30.0, 0.90, "neutral", "OC(Cn1cncn1)(Cn1cncn1)c1ccc(F)cc1F"),
 "lamotrigine": (256.1, 1.90, 90.7, 29.0, 0.98, "neutral", "Nc1nnc(c(N)n1)-c1cccc(Cl)c1Cl"), "carbamazepine": (236.3, 2.45, 46.3, 15.0, 0.85, "neutral", "NC(=O)N1c2ccccc2C=Cc2ccccc12"),
 "acyclovir": (225.2, -1.60, 114.0, 3.0, 0.20, "neutral", "Nc1nc2n(COCCO)cnc2c(=O)[nH]1"), "ranitidine": (314.4, 1.10, 86.0, 2.5, 0.50, "base", "CNC(=C[N+](=O)[O-])NCCSCc1ccc(CN(C)C)o1"),
 "citalopram": (324.4, 3.50, 36.3, 35.0, 0.80, "base", "CN(C)CCCC1(OCc2cc(C#N)ccc12)c1ccc(F)cc1"), "nifedipine": (346.3, 2.20, 110.0, 2.0, 0.50, "neutral", "COC(=O)C1=C(C)NC(C)=C(C(=O)OC)C1c1ccccc1[N+](=O)[O-]"),
 "metoprolol": (267.4, 1.90, 50.7, 3.5, 0.50, "base", "COCCc1ccc(OCC(O)CNC(C)C)cc1"), "levofloxacin": (361.4, -0.40, 73.3, 7.0, 0.99, "zwitterion", "CC1COc2c(N3CCN(C)CC3)c(F)cc3c(=O)c(C(O)=O)cn1c23")}
