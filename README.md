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

## Limits
Not validated against clinical data. Predicted parameters are low confidence; no docking/potency prediction for new analogs (they inherit the parent's measured IC50 as an assumption). A "pass" only means "worth real lab testing".
