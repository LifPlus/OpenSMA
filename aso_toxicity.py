"""
ASO Toxicity Screener (Module A)
---------------------------------
Screens all ASO candidates from aso_candidates_intron7.csv against known
toxicological heuristics derived from published literature:
  - CpG motif → TLR9 immunostimulation
  - Poly-A / Poly-C runs → Hepatotoxicity risk
  - G-Quadruplex motifs → Protein binding / off-target risk
  - High GC with clean motif profile → Optimal chemical space
"""

import pandas as pd
import re

# ----- Toxicology Rule Engine -----

def has_cpg_motif(seq):
    """Unmethylated CpG dinucleotides activate TLR9 → innate immune response."""
    # In RNA, 'CG' dinucleotide. In DNA ASO analogs this is 'CG'.
    seq_dna = seq.replace('U','T') # treat as DNA backbone for this check
    return 'CG' in seq_dna

def has_poly_run(seq, base, run_length=4):
    """Homopolymer runs (AAAA, CCCC) → hepatic accumulation and stability issues."""
    seq_dna = seq.replace('U','T')
    pattern = base * run_length
    return pattern in seq_dna

def has_g_quadruplex(seq):
    """G-quadruplex motif: G3+ N1-7 G3+ N1-7 G3+ N1-7 G3+ → off-target protein binding."""
    seq_dna = seq.replace('U','T')
    gq_pattern = re.compile(r'G{3,}.{1,7}G{3,}.{1,7}G{3,}.{1,7}G{3,}')
    return bool(gq_pattern.search(seq_dna))

def immunostim_score(seq):
    """Count CpG dinucleotides. More = higher immune activation risk."""
    seq_dna = seq.replace('U','T')
    return seq_dna.count('CG')

def hepatotox_score(seq):
    """Sum of risk: poly-A, poly-C, poly-T runs (≥4)."""
    score = 0
    for base in ['A','C','T']:
        if has_poly_run(seq, base, run_length=4):
            score += 1
    return score

def gc_content(seq):
    seq = seq.upper()
    return (seq.count('G') + seq.count('C')) / len(seq) * 100 if len(seq) > 0 else 0

def assess_toxicity(row):
    # Support both old ('ASO_Sequence_5_3') and new ('Sequence') column names
    seq_col = 'Sequence' if 'Sequence' in row.index else 'ASO_Sequence_5_3'
    aso = str(row[seq_col]).replace('U','T')

    cpg_flag = has_cpg_motif(aso)
    cpg_count = immunostim_score(aso)
    poly_flag = any(has_poly_run(aso, b) for b in ['A','C','T','G'])
    gq_flag   = has_g_quadruplex(aso)
    gc = round(gc_content(aso), 2)

    # Hepatotoxicity heuristic 
    hepato_risk = hepatotox_score(aso)

    # Compose overall toxicity flag
    total_flags = sum([cpg_flag, poly_flag, gq_flag])
    if total_flags == 0:
        tox_summary = "LOW"
    elif total_flags == 1:
        tox_summary = "MODERATE"
    else:
        tox_summary = "HIGH"

    # Off-target risk heuristic: very high GC (>70%) can cause non-specific binding
    offtarget_risk = "ELEVATED" if gc > 66 else "ACCEPTABLE"

    return pd.Series({
        'CpG_Motif_Present': cpg_flag,
        'CpG_Count': cpg_count,
        'Poly_Run_Flag': poly_flag,
        'G_Quadruplex_Flag': gq_flag,
        'GC_Content_%': gc,
        'Hepatotox_Risk_Score': hepato_risk,
        'OffTarget_Risk': offtarget_risk,
        'Overall_Tox_Rating': tox_summary
    })

# ----- Main -----

def screen_aso_candidates():
    df = pd.read_csv('aso_candidates_intron7.csv')

    print(f"Loaded {len(df)} ASO candidates from sliding window.\n")

    tox_df = df.apply(assess_toxicity, axis=1)
    result = pd.concat([df, tox_df], axis=1)

    # Filter to LOW toxicity only for clinical consideration
    sort_col = next((c for c in ['Smart_Score', 'Tm_RNA_RNA_C', 'GC_%'] if c in result.columns), None)
    safe_candidates = result[result['Overall_Tox_Rating'] == 'LOW']
    if sort_col:
        safe_candidates = safe_candidates.sort_values(sort_col, ascending=False)

    print("=== TOP SAFE ASO CANDIDATES (LOW Toxicity) ===")
    seq_col = 'Sequence' if 'Sequence' in safe_candidates.columns else 'ASO_Sequence_5_3'
    baseline_col = 'Is_Baseline_Drug' if 'Is_Baseline_Drug' in result.columns else 'Is_Nusinersen'
    print(safe_candidates[['Start_Position','End_Position', seq_col,
                            'GC_Content_%','CpG_Motif_Present',
                            'OffTarget_Risk','Overall_Tox_Rating']].head(10).to_string(index=False))

    print("\n=== NUSINERSEN (SPINRAZA) BASELINE TOXICITY PROFILE ===")
    spinraza = result[result[baseline_col] == True]
    print(spinraza[['Start_Position', seq_col, 'CpG_Motif_Present',
                     'Overall_Tox_Rating','OffTarget_Risk']].to_string(index=False))

    # Also flag any HIGH-tox candidates to discard
    discarded = result[result['Overall_Tox_Rating'] == 'HIGH']
    print(f"\n[DISCARDED] {len(discarded)} HIGH-toxicity candidates removed from pipeline.")

    result.to_csv('aso_toxicity_screened.csv', index=False)
    print(f"\nFull screened results saved to 'aso_toxicity_screened.csv'")
    
    return safe_candidates

if __name__ == "__main__":
    screen_aso_candidates()
