# MedScope
Terminal research explorer for any disease. Python stdlib only (works in QPython/Android).

**Not medical advice. Cannot discover or confirm a cure.** Scores are heuristic screening outputs.

## Run
`python medscope.py` (menu) or `python medscope.py <command> <disease>`

| command | what it does |
|---|---|
| overview / bio / trials / papers | Open Targets, ClinicalTrials.gov, PubMed |
| search | Europe PMC + OpenAlex in parallel (free, no key) |
| leads | Ph2/Ph3 trials + repurposing papers by evidence stage |
| watch | WHO outbreak news + newest novel-pathogen papers |
| design | `chemicals.py`: real structures (RDKit) -> filters -> rat & human whole-body PBPK -> dose + safety flags |
| brain | AI brief of the evidence |
| diseases / check / discover | PER-DISEASE targets, per-disease GOOD-vs-BAD validation (AUC), full discovery pipeline |
| lab add/report/list | log lab results, retrain on them, track REAL prospective accuracy |
| selftest | what works on this device (Termux/PC/Colab) + offline checks |
| report | save markdown report |

## Data & "brains"
- Real compound data when online: PubChem (properties), ChEMBL (measured IC50). Offline falls back to approximations.
- Brain backends: local **Ollama** (free, unlimited; set `OLLAMA_HOST`, `MEDSCOPE_MODEL`), or free tiers `GEMINI_API_KEY` / `GROQ_API_KEY` (rate-limited).
- Responses are cached in `~/.medscope`.

## Chemistry & simulators
- `chemistry.py` (needs **RDKit**: `pip install rdkit`; PC/Colab/server, not QPython): real SMILES, descriptors, QED, Lipinski/Veber, PAINS/Brenk alerts,
  analogs by medicinal-chemistry edits + BRICS recombination, SVG drawings. Without RDKit, `design` falls back to a property-only model.
- `pbpk.py` (pure Python, runs anywhere): rat (250 g) and human (70 kg) PBPK with 9 organs + blood, oral absorption, liver metabolism, kidney clearance,
  exact linear solver, mass-balance checked. Drug parameters are predicted from structure unless you pass measured ones (fup, clint, fa, ka, clr).

## Any disease, every drug, ingredients, rating
- `onboard "<disease>"` adds ANY disease automatically (Open Targets -> ChEMBL -> PDB); `onboard --uniprot ACC --name X` for a new pathogen once its protein is known.
- `repurpose <disease>` screens ALL approved drugs against that target (approved drugs already have human safety/PK data).
- `drug <name>` / `combo <a> <b>`: product -> active ingredients (RxNorm) -> structure, human PBPK, flags, organ overlap, CYP interaction. Excipients and synergy are NOT modelled.
- `score` (scorecard.py): evidence-based rating. Only things actually run and measured on your machine earn points; ceiling is 73/100 for software alone.

## Per-disease pipeline (`discover <disease>`)
`diseases.py` gives every disease its own target, GOOD list, BAD list and validation (13 diseases: hiv, covid, flu, hepc, mrsa, tb, ecoli, lung-cancer, breast-cancer, blood-cancer, malaria, fungi, inflammation; add yours in `~/.medscope/diseases.json`).
Flow: measured ChEMBL data (cleaned, no leakage) -> QSAR potency with error bars + ROC AUC -> RDKit analogs -> ADMET-AI -> Vina (trust-gated by redocking RMSD <= 2 A and Spearman rho vs measured pIC50) -> human PBPK + virtual population (CYP phenotypes, age, kidney function) + where the human target is expressed -> composite score (QED 25, potency/docking 30, F 15, Lipinski 10, tox 20).
Tiers: **1** pure Python (PBPK, population, lab loop, search) | **2** + RDKit + internet (structures, QSAR) | **3** + ADMET-AI + Vina. Run `python selftest.py` to see what your device supports.
**Lab loop:** predict -> test in a lab -> `lab add` -> retrain -> `lab report` shows your REAL prospective error. No accuracy number is promised; measure it.

## Limits
Not validated against clinical data. Predicted parameters are low confidence; no docking/potency prediction for new analogs (they inherit the parent's measured IC50 as an assumption). A "pass" only means "worth real lab testing".
