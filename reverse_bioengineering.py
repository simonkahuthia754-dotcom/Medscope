#!/usr/bin/env python3
"""reverse_bioengineering.py - reconstruct HOW a disease forms, from public human-biology data, to find where it can be broken.
  python reverse_bioengineering.py "<disease name>" [--save]      (MedScope: mechanism <disease name>)
Pipeline: Open Targets causal evidence (genetics vs downstream) -> pathway modules (Reactome) -> protein-interaction network
(STRING) with hub/bottleneck analysis -> therapeutic window (where the target is also expressed in vital organs) -> modality hint
(small molecule / antibody / PROTAC / other) -> non-redundant target COMBINATIONS (block escape routes) -> CURE BLOCKERS (why most
diseases are hard to eliminate completely) -> ranked shortlist handed to onboard/discover/repurpose.
Scope = defensive therapeutic target discovery from public databases. Associations are NOT proof of causation; genetics is the strongest
evidence. A target shortlist is a hypothesis list for experiments, not a cure."""
import json, math, sys
import humanbio, sources, state

GEN = ("genetic_association", "somatic_mutation")         # causal-leaning evidence types in Open Targets
VITAL = ("heart", "brain", "liver", "kidney", "lung")

# textbook-level knowledge (verify against current literature): lifecycle / formation and why elimination is hard
KNOW = {
 "hiv": {"forms": "virus enters CD4 T cells (CD4 + CCR5/CXCR4, gp41 fusion) -> reverse transcription -> integration of provirus into host DNA -> assembly/maturation",
         "blockers": ["latent reservoir: integrated provirus in resting memory CD4 T cells is invisible to drugs and immunity", "provirus is part of the host genome",
                      "high viral diversity -> escape/resistance", "lifelong therapy needed because drugs block steps, they do not remove the provirus"],
         "strategy": "combine >=2 independent drug classes (entry/RT/integrase/protease/capsid) for suppression; elimination needs reservoir-targeting approaches (latency reversal + clearance, or excision/gene therapy) - unsolved"},
 "covid": {"forms": "spike binds ACE2 -> entry -> polyprotein cleaved by Mpro/PLpro -> RNA polymerase copies genome -> assembly",
           "blockers": ["rapid variant evolution", "reinfection / waning immunity", "needs early dosing (replication peaks early)"],
           "strategy": "target conserved enzymes (Mpro, RdRp) in combination; host-directed adjuncts for hyper-inflammation"},
 "flu": {"forms": "hemagglutinin binds sialic acid -> entry -> polymerase (PA/PB1/PB2) replicates RNA -> neuraminidase releases virions",
         "blockers": ["antigenic drift/shift", "resistance mutations", "short treatment window (<48 h)"], "strategy": "conserved targets (polymerase PA, neuraminidase) in combination; early dosing"},
 "hepc": {"forms": "virus replicates in hepatocytes via NS3/4A protease, NS5A, NS5B polymerase", "blockers": ["largely SOLVED: direct-acting antivirals cure >95%; remaining issues are access, cost, reinfection"],
          "strategy": "already curable: focus on access"},
 "mrsa": {"forms": "Staphylococcus aureus carries mecA (PBP2a) beta-lactam resistance; protein synthesis at the 50S ribosome", "blockers": ["biofilms", "efflux and target mutation", "persister cells"],
          "strategy": "combine different mechanisms; anti-biofilm adjuncts; novel classes"},
 "tb": {"forms": "M. tuberculosis survives inside macrophages and granulomas, goes dormant (low metabolism)", "blockers": ["dormant persisters tolerate drugs", "thick lipid cell wall/poor drug penetration into granulomas", "6+ month therapy -> non-adherence -> resistance"],
        "strategy": "multi-drug regimens hitting different pathways (RNA polymerase, cell wall, ATP synthase); drugs active on dormant cells to shorten therapy"},
 "ecoli": {"forms": "Gram-negative: outer membrane + efflux limit drug entry; gyrase/topoisomerase IV drive DNA replication", "blockers": ["efflux pumps", "plasmid-borne resistance (ESBL, carbapenemase)"], "strategy": "novel classes with outer-membrane penetration; combinations with efflux/beta-lactamase inhibitors"},
 "malaria": {"forms": "sporozoites infect liver -> merozoites invade red cells (blood stage, symptoms) -> gametocytes transmit via mosquito; P. vivax keeps dormant liver hypnozoites",
             "blockers": ["hypnozoites cause relapse", "gametocytes keep transmission going", "artemisinin partial resistance in Africa and Asia"], "strategy": "combination therapy across stages (blood schizonticide + hypnozoitocide + transmission blocker)"},
 "fungi": {"forms": "ergosterol synthesis (CYP51) maintains the fungal membrane", "blockers": ["azole resistance (target mutation, efflux)", "biofilms", "host immunocompromise"], "strategy": "novel targets; combination with immune support"},
 "lung-cancer": {"forms": "driver mutations (e.g. EGFR) keep growth signalling switched on; clones evolve", "blockers": ["tumour heterogeneity and clonal evolution", "on-treatment resistance (e.g. EGFR C797S/T790M, bypass MET/HER2)", "metastasis and immune evasion"],
                 "strategy": "match drug to driver; combination/sequencing to pre-empt bypass; add immunotherapy"},
 "breast-cancer": {"forms": "HER2 amplification drives growth signalling", "blockers": ["resistance via downstream PI3K activation", "brain metastasis", "heterogeneity"], "strategy": "dual HER2 blockade + downstream pathway inhibitors"},
 "blood-cancer": {"forms": "BCR-ABL fusion kinase drives leukaemic cell growth", "blockers": ["quiescent leukaemic stem cells survive kinase inhibitors", "ABL kinase-domain resistance mutations (e.g. T315I)"], "strategy": "next-generation inhibitors for mutants; stem-cell-directed combinations"},
 "inflammation": {"forms": "COX-2 makes prostaglandins that drive pain and swelling", "blockers": ["symptom control only: the upstream cause is separate", "redundant inflammatory pathways"], "strategy": "treat cause where known; combine mechanisms"},
}
GENERIC_BLOCKERS = ["redundant/compensating pathways", "target also needed in healthy vital organs (narrow therapeutic window)", "late diagnosis after irreversible damage",
                    "chronic disease with no single 'off' switch", "no small-molecule-tractable target (needs antibody, gene therapy, cell therapy)"]

def _g(q, v): return humanbio._gql(q, v)["data"]

def causal_targets(disease, n=25):
    h = _g('query($q:String!){search(queryString:$q,entityNames:["disease"],page:{index:0,size:1}){hits{id name}}}', {"q": disease})["search"]["hits"]
    if not h: raise RuntimeError(f"'{disease}' not found in Open Targets; try a synonym")
    q = ("query($id:String!){disease(efoId:$id){name description associatedTargets(page:{index:0,size:%d}){rows{score datatypeScores{id score} "
         "target{id approvedSymbol approvedName proteinIds{id source} tractability{modality value} pathways{pathway pathwayId}}}}}}") % n
    d = _g(q, {"id": h[0]["id"]})["disease"]; rows = []
    for r in d["associatedTargets"]["rows"]:
        t = r["target"]; ds = {x["id"]: x["score"] for x in r["datatypeScores"]}
        mods = sorted({x["modality"] for x in t.get("tractability") or [] if x.get("value")})
        rows.append({"symbol": t["approvedSymbol"], "name": t["approvedName"], "assoc": r["score"], "genetic": max([ds.get(k, 0) for k in GEN] or [0]), "drugs": ds.get("known_drug", 0),
                     "uniprot": next((p["id"] for p in t.get("proteinIds") or [] if p.get("source") == "uniprot_swissprot"), None), "mods": mods,
                     "pathways": sorted({p["pathway"] for p in t.get("pathways") or []})})
    return h[0]["name"], (d.get("description") or ""), rows

def network(symbols, min_score=.7):
    """STRING protein-protein interactions among the targets -> (edges, degree, betweenness)"""
    try:
        r = sources.J("https://string-db.org/api/json/network?species=9606&required_score=%d&identifiers=%s" % (int(min_score * 1000), "%0d".join(symbols)))
    except Exception: return [], {}, {}
    E = {(x["preferredName_A"], x["preferredName_B"]) for x in r}; adj = {s: set() for s in symbols}
    for a, b in E:
        if a in adj and b in adj and a != b: adj[a].add(b); adj[b].add(a)
    return sorted(E), {s: len(v) for s, v in adj.items()}, betweenness(adj)

def betweenness(adj):
    """Brandes algorithm (unweighted, undirected); high value = chokepoint on many shortest paths"""
    bc = {v: 0.0 for v in adj}
    for s in adj:
        S, P, sig, dist, Q = [], {v: [] for v in adj}, dict.fromkeys(adj, 0), dict.fromkeys(adj, -1), [s]; sig[s] = 1; dist[s] = 0
        while Q:
            v = Q.pop(0); S.append(v)
            for w in adj[v]:
                if dist[w] < 0: dist[w] = dist[v] + 1; Q.append(w)
                if dist[w] == dist[v] + 1: sig[w] += sig[v]; P[w].append(v)
        dl = dict.fromkeys(adj, 0.0)
        while S:
            w = S.pop()
            for v in P[w]: dl[v] += sig[v] / sig[w] * (1 + dl[w])
            if w != s: bc[w] += dl[w]
    m = max(bc.values() or [0]) or 1; return {k: v / m for k, v in bc.items()}

def window(symbol):
    """therapeutic window: vital organs where the human target is also strongly expressed (= where on-target toxicity could occur)"""
    try: lv, _ = humanbio.target_expression(symbol)
    except Exception: return None, []
    bad = sorted(o for o in VITAL if lv.get(o, 0) >= 3); return 1 - len(bad) / len(VITAL), bad

def shortlist(rows, deg, bet, top=8):
    """transparent ranking: genetics .30 + overall association .20 + tractability .15 + network chokepoint .10 + selectivity window .15 + existing-drug evidence .10"""
    out = []
    for r in rows[:15]:
        sel, bad = window(r["symbol"]); mx = max(deg.values() or [1]) or 1
        s = (.30 * r["genetic"] + .20 * r["assoc"] + .15 * (1 if "SM" in r["mods"] else .6 if r["mods"] else 0) + .10 * (bet.get(r["symbol"], 0) * .5 + deg.get(r["symbol"], 0) / mx * .5)
             + .15 * (sel if sel is not None else .5) + .10 * r["drugs"])
        out.append({**r, "score": round(100 * s, 1), "vital_expr": bad, "hub": bet.get(r["symbol"], 0)})
    return sorted(out, key=lambda x: -x["score"])[:top]

def combinations(sl, k=3):
    """non-redundant target sets: highest-ranked targets that sit in DIFFERENT pathways, so one escape route cannot defeat the regimen"""
    chosen, used = [], set()
    for t in sl:
        p = set(t["pathways"])
        if not chosen or len(p & used) / (len(p) or 1) < .5: chosen.append(t); used |= p
        if len(chosen) == k: break
    return chosen

def analyze(disease, top=8, key=None, save=False):
    name, desc, rows = causal_targets(disease); print(f"\n=== HOW '{name}' FORMS: reverse-engineered from public human data ===")
    if desc: print(desc[:420] + ("..." if len(desc) > 420 else ""))
    k = key if key in KNOW else None
    if k: print(f"\nFormation (textbook-level, verify): {KNOW[k]['forms']}")
    print("\n1) CAUSAL EVIDENCE: genetics-first (genetic/somatic evidence is the closest public proxy for causation)")
    for r in sorted(rows, key=lambda r: -r["genetic"])[:8]: print(f"   {r['symbol']:<9} genetic {r['genetic']:.2f} | overall {r['assoc']:.2f} | modalities {','.join(r['mods']) or '-'} | {r['name'][:46]}")
    pw = {}
    for r in rows:
        for p in r["pathways"]: pw[p] = pw.get(p, 0) + 1
    print("\n2) PATHWAY MODULES the disease runs through (targets per pathway):")
    for p, c in sorted(pw.items(), key=lambda x: -x[1])[:6]: print(f"   {c:>2}  {p}")
    E, deg, bet = network([r["symbol"] for r in rows]); print(f"\n3) INTERACTION NETWORK (STRING, high confidence): {len(E)} links among {len(rows)} targets" if E else "\n3) INTERACTION NETWORK: unavailable (offline?)")
    for s, b in sorted(bet.items(), key=lambda x: -x[1])[:4]:
        if b > 0: print(f"   chokepoint {s}: centrality {b:.2f}, {deg[s]} partners")
    sl = shortlist(rows, deg, bet, top); print(f"\n4) SHORTLIST (genetics .30 + association .20 + tractability .15 + chokepoint .10 + window .15 + existing drugs .10; a ranking aid, not proof)")
    for t in sl:
        print(f"   {t['score']:>5}  {t['symbol']:<9} {','.join(t['mods']) or 'no tractable modality':<14}" + (f" ! also expressed in vital organ(s): {', '.join(t['vital_expr'])} (narrow window)" if t["vital_expr"] else ""))
    cb = combinations(sl); print("\n5) NON-REDUNDANT COMBINATION (different pathways, so one escape route cannot defeat it): " + ", ".join(t["symbol"] for t in cb))
    print("\n6) CURE BLOCKERS - why complete elimination is hard:")
    for b in (KNOW[k]["blockers"] if k else GENERIC_BLOCKERS): print("   - " + b)
    if k: print(f"\nStrategy (textbook-level): {KNOW[k]['strategy']}")
    if any(t["genetic"] >= .6 and not t["mods"] for t in sl): print("\nNote: strong genetic evidence with no small-molecule target suggests antibody, gene-therapy or cell-therapy routes.")
    if save:
        import onboard
        for t in sl[:3]:
            if t["uniprot"]: print(f"   saving {t['symbol']} ..."); onboard._save({f"{onboard.slug(disease)}-{t['symbol'].lower()}": onboard.entry(name, t["name"], t["symbol"], t["uniprot"], "from reverse_bioengineering")})
    print("\nNEXT: onboard \"" + disease + "\" -> discover/repurpose <key> to find molecules for the top targets. These are hypotheses to test in experiments, not cures.\n")
    state.record("mechanism_ok", True); return sl

def main():
    a = [x for x in sys.argv[1:] if not x.startswith("--")]
    if not a: print(__doc__); return
    analyze(" ".join(a), save="--save" in sys.argv)

if __name__ == "__main__": main()
