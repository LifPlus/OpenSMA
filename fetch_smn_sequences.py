import os
from Bio import Entrez, SeqIO

# Set your NCBI email (required by Entrez)
Entrez.email = "mrsmn.research@example.com" 

def fetch_sequence(accession_id):
    print(f"Fetching sequence for {accession_id} from NCBI...")
    try:
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        return str(record.seq)
    except Exception as e:
        print(f"Error fetching {accession_id}: {e}")
        return None

def analyze_smn_sequences():
    """
    Analyze the critical SMN1/SMN2 Exon 7 and Intron 7 sequences.
    
    References:
        - Lorson CE et al., PNAS 1999; 96(11):6307-6311 (C→T transition at Exon 7 +6)
        - Cartegni L & Krainer AR, Nat Genet 2002; 30:377-384 (ESE disruption mechanism)
        - Singh NN et al., Mol Cell Biol 2006; 26(4):1333-1346 (ISS-N1 characterization)
        - Hua Y et al., Am J Hum Genet 2008; 82(4):834-848 (Nusinersen target mapping)
    
    Gene accessions:
        - SMN1 mRNA: NM_000344.4
        - SMN2 mRNA: NM_017411.4  
        - SMN2 genomic: NG_008728.1
    """
    print("Initializing SMN Exon 7 and Intron 7 consensus sequences based on established literature.")
    print("  References: Lorson et al. PNAS 1999, Cartegni & Krainer Nat Genet 2002, Singh et al. MCB 2006")
    
    # SMN Exon 7 = 54 nucleotides (coding strand, 5'→3')
    # Critical C→T (C6U in mRNA) transition at position +6 (1-indexed)
    # SMN1 +6 = C → creates ESE recognized by SF2/ASF → exon 7 included
    # SMN2 +6 = T → disrupts ESE, creates ESS → exon 7 skipped → disease
    smn1_exon7 = "GGTTTCAGACAAAATCAAAAAGAAGGAAGGTGCTCACATTCCTTAAATTAAGGA"
    smn2_exon7 = "GGTTTTAG" + "ACAAAATCAAAAAGAAGGAAGGTGCTCACATTCCTTAAATTAAGGA"
    
    # Intron 7 begins right after Exon 7 (splice donor GT at +1,+2).
    # The first 50 nucleotides of Intron 7 play a massive role.
    # ISS-N1 (Intronic Splicing Silencer N1) at positions +10 to +24 (Singh et al. MCB 2006)
    # SMN1 and SMN2 intron 7 sequences are identical in this region.
    # Source: NG_008728.1 (GenBank), Singh et al. 2006 Fig. 1
    intron7_start = "GTAAGTCTGCCAGCATTATGAAAGTGAATCTTACTTTTGTAAAACTTTAT"
    
    print("\n--- Sequence Analysis ---")
    print(f"SMN1 Exon 7: {smn1_exon7}")
    print(f"SMN2 Exon 7: {smn2_exon7}")
    
    # C→T transition check at position +6 (1-indexed) = index [5] (0-indexed)
    assert len(smn1_exon7) == 54, f"SMN1 Exon 7 must be 54 nt, got {len(smn1_exon7)}"
    assert len(smn2_exon7) == 54, f"SMN2 Exon 7 must be 54 nt, got {len(smn2_exon7)}"
    if smn1_exon7[5] == 'C' and smn2_exon7[5] == 'T':
        print("✓ Critical Exon 7 +6 C→T transition verified (Lorson et al. 1999).")
        print("  SMN1[+6]=C → SF2/ASF ESE intact → exon 7 included")
        print("  SMN2[+6]=T → ESE disrupted → exon 7 skipping → SMA")
    else:
        print("⚠ WARNING: C→T transition verification FAILED at position +6!")
        print(f"  SMN1[5]={smn1_exon7[5]}, SMN2[5]={smn2_exon7[5]}")
    
    # ISS-N1 mapping
    print("\n--- Intron 7 Regulatory Elements ---")
    print(f"Intron 7 (first 50bp): {intron7_start}")
    
    # ISS-N1 is located at positions +10 to +24 of Intron 7 (15 nt)
    # 0-indexed in Python: [9:24] → positions 9 through 23 inclusive
    # Reference: Singh et al. Mol Cell Biol 2006; 26(4):1333-1346
    iss_n1 = intron7_start[9:24]
    print(f"ISS-N1 Region (+10 to +24): {iss_n1}")
    assert iss_n1 == "CCAGCATTATGAAAG", f"ISS-N1 validation failed: got {iss_n1}"
    print("✓ ISS-N1 sequence validated against Singh et al. 2006")
    
    # Spinraza (Nusinersen) targets ISS-N1 region.
    # Nusinersen is an 18-mer 2'-O-methoxyethyl phosphorothioate ASO
    # targeting Intron 7 from +10 to +27 (Hua et al. Am J Hum Genet 2008)
    spinraza_target = intron7_start[9:27]
    print(f"Spinraza (Nusinersen) Target Region (+10 to +27): {spinraza_target}")
    assert len(spinraza_target) == 18, f"Nusinersen target must be 18-mer, got {len(spinraza_target)}"
    
    # Generating the antisense sequence for Spinraza
    # RNA transcription (Intron 7 is copied to pre-mRNA as U instead of T)
    pre_mrna_spinraza_target = spinraza_target.replace('T', 'U')
    
    # ASO is DNA/RNA analog, complement of the pre-mRNA target.
    # A->U, U->A, C->G, G->C (but ASO uses T often in MOE/PMO, let's just use RNA complement A-U G-C for binding math)
    complement_map = str.maketrans('ACGU', 'UGCA')
    spinraza_aso_rna = pre_mrna_spinraza_target.translate(complement_map)[::-1] # 5' to 3' antisense
    print(f"Nusinersen Analog RNA Sequence (5'->3'): {spinraza_aso_rna}")
    
    # Save to a fasta file for structural analysis
    with open("smn2_target_region.fasta", "w") as f:
        f.write(">SMN2_Exon7_Intron7_Flank\n")
        f.write(smn2_exon7 + intron7_start + "\n")
        f.write(">ISS_N1_Region\n")
        f.write(iss_n1 + "\n")
        
    print("\nMapping complete. Seed sequences saved to 'smn2_target_region.fasta'.")

if __name__ == "__main__":
    analyze_smn_sequences()
