"""
CRISPR / Base Editing Guide RNA Designer (Phase 6 — Module E)
===============================================================
This module designs and evaluates guide RNAs (gRNAs) targeting:
  1. SpCas9 (NGG PAM) — Classic Cas9 for reference
  2. SaCas9 (NNGRRT PAM) — Smaller, fits AAV better
  3. Base Editing (CBE: C→T correction at SMN2 Exon 7 +6 position)
     Using cytosine base editor (BE4max or equivalent) to correct
     the SMN2-specific C→T back to C (making it functionally SMN1-like)

KEY TARGET:
  SMN2 Exon 7 position +6: The single C→T nucleotide change that
  distinguishes SMN2 from SMN1, disrupts the ESE (Exonic Splicing
  Enhancer), and causes Exon 7 skipping → disease.

Goal: Design a gRNA that positions a base editor EXACTLY over
position +6 of Exon 7 to perform the T→C correction
(in genomic/coding strand terms: editing the T back to C at the
SMN2-specific locus without cutting the DNA).
"""

import pandas as pd
import re
from Bio import pairwise2
from Bio.Seq import Seq

# ══════════════════════════════════════════════════════════════════════
# Reference sequences (from GenBank NG_008728.1 / NC_000005)
# SMN2 genomic context around Exon 7 (critical editing window)
#
# References:
#   - Lorson CE et al., PNAS 1999; 96(11):6307-6311
#   - Cartegni L & Krainer AR, Nat Genet 2002; 30:377-384
#   - GenBank NG_008728.1 (SMN2 genomic region)
# ══════════════════════════════════════════════════════════════════════

"""
SMN1 Exon 7 (reference, normal):   ...GGTTT[C]AGACAAAAT...
SMN2 Exon 7 (disease, to correct): ...GGTTT[T]AGACAAAAT...

Position convention:
  +1 = first nucleotide of Exon 7
  +6 = the critical C (SMN1) / T (SMN2) position (Lorson et al. 1999)

Genomic context: ~30bp intron 6 tail + Exon 7 (54 nt) + intron 7 start
"""

# Genomic context centred on the SMN2 Exon 7 editing window
# This is a ~115bp window around the edit site, coding strand (5'→3')
# Intron6-tail | Exon7 (54 nt, with T at +6) | Intron7-start
# The T at the +6 position within Exon 7 is the SMN2-specific mutation
SMN2_GENOMIC_CONTEXT = (
    "TAGTTCAGCTTCAGCCTTATATTGTGATCTTAATTTCGATTTGGTTTTAG"
    "ACAAAATCAAAAAGAAGGAAGGTGCTCACATTCCTTAAATTAAGGAG"
    "TAAGTCTGCCAGCATTATGAAAGT"
)

# Mark the target T (SMN2 position +6) within the context
# In our string, let's find the "GGTTTTTGAC" pattern and locate position +6
def find_edit_site(genomic_seq):
    """
    Find the SMN2-specific T at position +6 of Exon 7 within the genomic context.
    SMN2 has GGTTTT (T at +6) while SMN1 has GGTTTC (C at +6).
    The pattern TTTTAG marks the end of exon 7 positions +3 through +7.
    """
    # Find "GGTTTTAG" which contains the SMN2-specific T at +6
    # In SMN1 this would be "GGTTTCAG" (C at position +6)
    pattern = "GGTTTTAG"
    idx = genomic_seq.find(pattern)
    if idx != -1:
        return idx + 5  # 0-indexed position of the critical T (pos +6 of exon 7)
    # Fallback: search for the broader exon 7 start context
    pattern = "TTTTAG"
    idx = genomic_seq.find(pattern)
    if idx != -1:
        return idx + 3  # position of the critical T
    return None

EDIT_SITE_INDEX = find_edit_site(SMN2_GENOMIC_CONTEXT)

# ══════════════════════════════════════════════════════════════════════
# PAM Site Finder
# ══════════════════════════════════════════════════════════════════════
def find_pam_sites(seq, pam_pattern, strand='+', window=25):
    """
    Find all PAM sites and extract the adjacent 20nt protospacer (guide RNA).
    For SpCas9: NGG (look for NGG downstream of protospacer)
    For SaCas9: NNGRRT
    For base editing: edit site must fall in editing window (pos 4-8 from PAM-proximal end)
    """
    results = []
    # Compile regex for PAM
    pam_re = re.compile(pam_pattern)

    for m in pam_re.finditer(seq):
        pam_start = m.start()
        pam_end   = m.end()
        pam_len   = pam_end - pam_start

        # gRNA = the 20nt UPSTREAM of the PAM (on the same strand)
        grna_start = pam_start - 20
        grna_end   = pam_start

        if grna_start < 0:
            continue  # not enough sequence upstream

        grna = seq[grna_start:grna_end]
        pam  = seq[pam_start:pam_end]

        # Distance from edit site to gRNA (for base editing window check)
        if EDIT_SITE_INDEX is not None:
            # Position of edit site within the gRNA (1-indexed from 5' end)
            # Edit site must be at position 4-8 from PAM-proximal end for CBE
            edit_in_grna = EDIT_SITE_INDEX - grna_start  # 0-indexed
            be_window_pos = 20 - edit_in_grna  # distance from PAM end
            be_suitable = (3 <= edit_in_grna <= 9)  # CBE editing window: nt 4-9
        else:
            be_window_pos = -1
            be_suitable = False

        results.append({
            'strand':        strand,
            'gRNA_seq':      grna,
            'PAM':           pam,
            'gRNA_start':    grna_start,
            'gRNA_end':      grna_end,
            'PAM_start':     pam_start,
            'edit_pos_in_gRNA': edit_in_grna if EDIT_SITE_INDEX else 'N/A',
            'BE_window_suitable': be_suitable
        })

    return results

def reverse_complement(seq):
    complement = str.maketrans('ACGT', 'TGCA')
    return seq.translate(complement)[::-1]

# ══════════════════════════════════════════════════════════════════════
# Efficiency Score (Doench 2016 Rule Set 2 — simplified heuristics)
# ══════════════════════════════════════════════════════════════════════
def score_grna_efficiency(grna):
    """
    Simplified Doench 2016 Rule Set 2 heuristics.
    Full model requires a trained ML model; here we apply rule-based surrogates.
    """
    score = 50.0  # baseline

    # GC content 40-70% is optimal
    gc = (grna.count('G') + grna.count('C')) / len(grna) * 100
    if 40 <= gc <= 70:
        score += 15
    elif gc < 30 or gc > 80:
        score -= 20

    # G at position 20 (PAM-proximal, improves Cas9 binding)
    if grna[-1] == 'G':
        score += 5

    # No poly-T run (>4T causes RNA Pol III termination in gRNA expression)
    if 'TTTT' in grna:
        score -= 25

    # Avoid runs of >4 same base
    for base in ['A','C','G','T']:
        if base * 4 in grna:
            score -= 10

    # G at position 1 (5' G helps U6 promoter transcription initiation)
    if grna[0] == 'G':
        score += 5

    # GG dinucleotide at PAM-proximal end (positions 19-20)
    if grna[-2:] == 'GG':
        score += 8

    return round(min(100, max(0, score)), 1)

# ══════════════════════════════════════════════════════════════════════
# Off-Target Risk (Genomic Alignment Heuristic)
# ══════════════════════════════════════════════════════════════════════
def off_target_risk(grna):
    """
    Phase 2 Upgrade: Uses pairwise sequence alignment against a mock "off-target"
    genome library to estimate homology-based off-target risks, combined with
    seed-region GC heuristics.
    """
    seed = grna[-12:]  # last 12 nt (PAM-proximal = seed region)
    gc_seed = (seed.count('G') + seed.count('C')) / 12 * 100

    # Mock critical off-target regions (e.g. tumor suppressor genes prone to editing)
    mock_off_targets = [
        "GGTTTCAGACAAAATCAAAA",  # Near-perfect match (SMN1 itself!)
        "GCTTTCAGACAAAATCAAAA",  # 1 mismatch
        "TAGTTCAGCTTCAGCCTTAT",  # Random genome noise
        "CCGTTCAGACAAAATCAAAA"   # 2 mismatches
    ]
    
    max_homology_score = 0
    for ot in mock_off_targets:
        # Local alignment: match=2, mismatch=-1, gap_open=-2, gap_extend=-1
        alignments = pairwise2.align.localms(grna, ot, 2, -1, -2, -1)
        if alignments:
            best_score = alignments[0].score # Max score is 2 * len(grna) = 40
            pct_homology = (best_score / (2 * len(grna))) * 100
            if pct_homology > max_homology_score:
                max_homology_score = pct_homology

    # High homology to off-targets overrides seed heuristics
    if max_homology_score > 85:
        return "HIGH (Homology Alert)", gc_seed, max_homology_score
    elif max_homology_score > 70:
        return "MODERATE (Homology Alert)", gc_seed, max_homology_score

    # Very high GC in seed → more potential off-targets (non-specific hybridization)
    if gc_seed > 75:
        return "HIGH (GC-Rich Seed)", gc_seed, max_homology_score
    elif gc_seed > 55:
        return "MODERATE (GC-Rich Seed)", gc_seed, max_homology_score
    else:
        return "LOW", gc_seed, max_homology_score

# ══════════════════════════════════════════════════════════════════════
# Main: Run full gRNA design
# ══════════════════════════════════════════════════════════════════════
def design_crispr_grnas():
    seq  = SMN2_GENOMIC_CONTEXT
    rcseq = reverse_complement(seq)

    print(f"Genomic context length: {len(seq)} bp")
    print(f"Edit site index (SMN2-specific T): position {EDIT_SITE_INDEX}")
    print(f"Context at edit site: ...{seq[max(0,EDIT_SITE_INDEX-5):EDIT_SITE_INDEX+6]}...")
    print()

    all_grnas = []

    # SpCas9: NGG PAM on + strand
    sp_plus  = find_pam_sites(seq,  r'[ACGT]GG', strand='+')
    # SpCas9: NGG PAM on - strand (RC of seq)
    sp_minus = find_pam_sites(rcseq, r'[ACGT]GG', strand='-')

    for g in sp_plus + sp_minus:
        g['Cas_variant'] = 'SpCas9'
        g['PAM_type'] = 'NGG'
        all_grnas.append(g)

    # SaCas9: NNGRRT PAM (more compact, better AAV)
    sa_plus  = find_pam_sites(seq,  r'[ACGT]{2}G[AG][AG]T', strand='+')
    for g in sa_plus:
        g['Cas_variant'] = 'SaCas9'
        g['PAM_type'] = 'NNGRRT'
        all_grnas.append(g)

    # Score all candidates
    results = []
    for g in all_grnas:
        grna = g['gRNA_seq']
        if len(grna) < 20:
            continue
        efficiency = score_grna_efficiency(grna)
        ot_risk, seed_gc, seq_homology = off_target_risk(grna)
        be_ok = g.get('BE_window_suitable', False)

        results.append({
            'Cas_Variant':          g['Cas_variant'],
            'PAM_Type':             g['PAM_type'],
            'Strand':               g['strand'],
            'gRNA_Sequence_5_3':    grna,
            'PAM':                  g['PAM'],
            'Edit_Pos_in_gRNA':     g.get('edit_pos_in_gRNA', 'N/A'),
            'BE_Window_Suitable':   be_ok,
            'Efficiency_Score':     efficiency,
            'OffTarget_Risk':       ot_risk,
            'Seed_GC_%':            round(seed_gc, 1),
            'Max_OffTarget_Homology_%': round(seq_homology, 1)
        })

    df = pd.DataFrame(results)

    # Priority: BE-suitable, sorted by efficiency desc, off-target asc
    df_priority = df[df['BE_Window_Suitable'] == True].sort_values(
        ['Efficiency_Score', 'OffTarget_Risk'],
        ascending=[False, True]
    )

    print("="*78)
    print("  TOP BASE-EDITING-SUITABLE gRNAs (Edit site in CBE window)")
    print("  These can correct SMN2 T→C at Exon 7 +6 WITHOUT double-strand breaks")
    print("="*78)
    if not df_priority.empty:
        print(df_priority[['Cas_Variant','gRNA_Sequence_5_3','PAM',
                            'Edit_Pos_in_gRNA','Efficiency_Score',
                            'OffTarget_Risk']].head(5).to_string(index=False))
    else:
        print("No BE-window-suitable gRNAs found in this context window.")
        print("All SpCas9 / SaCas9 gRNAs (for reference knockout / HDR):")
        print(df.sort_values('Efficiency_Score', ascending=False).head(5)[
            ['Cas_Variant','gRNA_Sequence_5_3','PAM','Efficiency_Score','OffTarget_Risk']
        ].to_string(index=False))

    print("\n--- All gRNA Candidates ---")
    print(df.sort_values('Efficiency_Score', ascending=False).head(10).to_string(index=False))

    df.to_csv('crispr_grna_candidates.csv', index=False)
    print(f"\nAll {len(df)} gRNA candidates saved to 'crispr_grna_candidates.csv'")
    return df, df_priority

if __name__ == "__main__":
    df_all, df_be = design_crispr_grnas()
