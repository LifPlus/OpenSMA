"""Quick standalone runner for CRISPR gRNA design + 5-year simulation.

References:
  - Lorson CE et al., PNAS 1999; 96(11):6307-6311
  - Finkel RS et al., NEJM 2017; 377(18):1723-1732
  - Baranello G et al., NEJM 2021; 384(10):915-923
"""

import math
import json

# ── CRISPR gRNA Design ──────────────────────────────────────────────
# Genomic context: intron6-tail | Exon 7 (SMN2, T at +6) | intron7-start
SEQ = ("TAGTTCAGCTTCAGCCTTATATTGTGATCTTAATTTCGATTTGGTTTTAG"
       "ACAAAATCAAAAAGAAGGAAGGTGCTCACATTCCTTAAATTAAGGAG"
       "TAAGTCTGCCAGCATTATGAAAGT")

# Find the critical T at position +6 of Exon 7 (SMN2-specific)
idx = SEQ.find("GGTTTTAG")
edit_site = idx + 5  # 0-indexed position of the critical T

print(f"Edit site index: {edit_site}, nucleotide: {SEQ[edit_site]}")
print(f"Context around edit: ...{SEQ[edit_site-5:edit_site+6]}...")

def score_grna(g):
    gc = (g.count('G') + g.count('C')) / len(g) * 100
    s = 50
    if 40 <= gc <= 70: s += 15
    if g[-1] == 'G':   s += 5
    if 'TTTT' in g:    s -= 25
    if g[0] == 'G':    s += 5
    if g[-2:] == 'GG': s += 8
    return int(min(100, s))

def seed_gc_risk(g):
    seed = g[-12:]
    gc = (seed.count('G') + seed.count('C')) / 12 * 100
    return ("HIGH" if gc > 75 else "MODERATE" if gc > 55 else "LOW"), round(gc, 1)

candidates = []
for i in range(len(SEQ) - 23):
    pam = SEQ[i+20:i+23]
    if pam[1] == 'G' and pam[2] == 'G':
        g = SEQ[i:i+20]
        if len(g) < 20: continue
        ep = edit_site - i
        be = (3 <= ep <= 9)
        risk, sgc = seed_gc_risk(g)
        candidates.append({
            'gRNA': g, 'PAM': pam, 'Score': score_grna(g),
            'EditPos_in_gRNA': ep, 'BE_Suitable': be,
            'OffTarget_Risk': risk, 'SeedGC': sgc
        })

candidates.sort(key=lambda x: (-x['Score'], not x['BE_Suitable']))

print("\n=== TOP SpCas9 gRNA CANDIDATES (NGG PAM) ===")
print(f"{'gRNA (5->3)'.ljust(22)} {'PAM'} {'Score':>5} {'EditPos':>8} {'BE?':>4} {'OffTgt':>8}")
print('─'*60)
for c in candidates[:10]:
    print(f"{c['gRNA'].ljust(22)} {c['PAM']}  {c['Score']:>4}  {str(c['EditPos_in_gRNA']):>7}  "
          f"{'YES' if c['BE_Suitable'] else 'no':>4}  {c['OffTarget_Risk']:>8}")

be_candidates = [c for c in candidates if c['BE_Suitable']]
print(f"\n✓ {len(be_candidates)} Base-Editing-suitable gRNAs identified")
if be_candidates:
    print("\n=== BASE EDITING CANDIDATES (CBE window pos 4-9 from PAM-proximal) ===")
    for c in be_candidates[:3]:
        print(f"  gRNA:     {c['gRNA']}")
        print(f"  PAM:      {c['PAM']}  |  Score: {c['Score']}  |  EditPos: {c['EditPos_in_gRNA']}  |  OffTarget: {c['OffTarget_Risk']}")
        print(f"  EFFECT:   CBE corrects T→C at position {c['EditPos_in_gRNA']} of gRNA = SMN2 treated as SMN1")
        print()
else:
    print("\n  NOTE: On + strand no BE-window hit. Top HDR/knockout candidate:")
    c = candidates[0]
    print(f"  gRNA: {c['gRNA']} | PAM: {c['PAM']} | Score: {c['Score']}")

print("\n" + "="*70)
print("  5-YEAR PATIENT OUTCOME SIMULATION")
print("="*70)

def smn_traj(scenario, t_months):
    if scenario == "Untreated":
        return 15.0
    elif scenario == "Nusinersen":
        return 15 + 25 * (1 - math.exp(-0.8 * t_months))
    elif scenario == "OpenSMA_Phase1":
        return 15 + 30 * (1 - math.exp(-0.9 * t_months))
    elif scenario == "Full_Cure":
        if t_months <= 3:
            return 15 + 30 * (1 - math.exp(-0.9 * t_months))
        else:
            bridge = 15 + 30 * (1 - math.exp(-2.7))
            # Realistic base editing: ~50% efficiency → ~60% SMN target
            return bridge + (60 - bridge) * (1 - math.exp(-0.6 * (t_months - 3)))
    return 15.0

def mn_at(scenario, months, steps=100):
    """Simulate motor neuron survival."""
    decay = {"Untreated": 0.05, "Nusinersen": 0.018,
              "OpenSMA_Phase1": 0.014, "Full_Cure": 0.003}
    mn = 60.0
    dt = months / steps
    for step in range(steps):
        t = step * dt
        smn = smn_traj(scenario, t)
        prot = smn / 100.0
        rate = decay[scenario] * (1 - prot * 0.92)
        mn = max(mn * math.exp(-rate * dt), 2.0)
    return round(mn, 1)

def prob_milestone(smn, threshold):
    return round(100 / (1 + math.exp(-0.12 * (smn - threshold))), 1)

scenarios = ["Untreated", "Nusinersen", "OpenSMA_Phase1", "Full_Cure"]
checkpoints = [12, 24, 60]

print(f"\n{'Scenario'.ljust(22)} {'Month':>6} {'SMN%':>6} {'MN%':>6} {'Breath':>7} {'Sit':>6} {'Stand':>7} {'Walk':>6}")
print('─'*75)

summary = {}
for sc in scenarios:
    summary[sc] = {}
    for cp in checkpoints:
        smn = round(min(100, smn_traj(sc, cp)), 1)
        mn  = mn_at(sc, cp)
        breath = prob_milestone(smn, 38)
        sit    = prob_milestone(smn, 50)
        stand  = prob_milestone(smn, 62)
        walk   = prob_milestone(smn, 75)
        summary[sc][f"month_{cp}"] = {
            "SMN_%": smn, "MotorNeurons_%": mn,
            "Breathing_%": breath, "Sit_Unaided_%": sit,
            "Stand_w_Support_%": stand, "Walk_Indep_%": walk
        }
        print(f"{sc.ljust(22)} {str(cp)+'m':>6} {smn:>6} {mn:>6} {breath:>7} {sit:>6} {stand:>7} {walk:>6}")
    print()

with open("full_cure_functional_outcomes.json","w") as f:
    json.dump(summary, f, indent=2)

print("Functional outcomes saved to 'full_cure_functional_outcomes.json'")
print("\n=== WINNER ===")
smn_fc = round(min(100, smn_traj("Full_Cure", 60)), 1)
print(f"Full Cure Protocol at 5 years: SMN = {smn_fc}%")
print(f"  → Walking probability:    {prob_milestone(smn_fc, 75)}%")
print(f"  → Breathing unaided:      {prob_milestone(smn_fc, 38)}%")
print(f"  → Sitting unaided:        {prob_milestone(smn_fc, 50)}%")
print(f"  → Standing with support:  {prob_milestone(smn_fc, 62)}%")
