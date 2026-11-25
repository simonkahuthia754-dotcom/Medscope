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
| design | `chemicals.py`: candidates -> lab -> rat -> human simulators |
| brain | AI brief of the evidence |
| report | save markdown report |

## Data & "brains"
- Real compound data when online: PubChem (properties), ChEMBL (measured IC50). Offline falls back to approximations.
- Brain backends: local **Ollama** (free, unlimited; set `OLLAMA_HOST`, `MEDSCOPE_MODEL`), or free tiers `GEMINI_API_KEY` / `GROQ_API_KEY` (rate-limited).
- Responses are cached in `~/.medscope`.

## Limits
Simulators use simple property rules, not docking or validated ADMET models. A "pass" only means "worth real lab testing".
