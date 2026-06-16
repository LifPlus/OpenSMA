import json

def parse_mock_vcf(vcf_path):
    # In a real scenario, we would use pysam.VariantFile or similar to parse raw VCFs
    # Here we simulate reading a patient VCF and returning the key genotypes.
    print(f"Parsing VCF data from {vcf_path}...")
    
    # Mocking patient genetic profile extraction
    # A real tool would look at locus 5q13.2 (chr5:70,220,768-70,249,769 for SMN1 and similar for SMN2)
    patient_data = {
        "Patient_ID": "SMA_Pt_001",
        "SMN1_Copy_Number": 0,    # 0 = SMA diagnosed
        "SMN2_Copy_Number": 2,    # Determines SMA severity type (Type 1/2 usually have 2-3 copies)
        "SMN2_c.859G>C_Modifier": False, # Rare variant creating an ESE, reducing severity
        "PLS3_Expression_Level": "Normal", # Plastin 3 - protective modifier if overexpressed
        "CORO1C_Variant": "WildType"       # Coronin 1C - protective modifier
    }
    return patient_data

def generate_intervention_strategy(patient_data):
    print("\n--- Generating Personalized Intervention Strategy ---")
    smn1_cn = patient_data["SMN1_Copy_Number"]
    smn2_cn = patient_data["SMN2_Copy_Number"]
    
    if smn1_cn > 0:
        return "Patient has SMN1 copies. Unlikely to have classical SMA. No intervention recommended."
        
    print(f"Diagnosis: SMA Patient. (SMN1 Copies: {smn1_cn}, SMN2 Copies: {smn2_cn})")
    print(f"Modifiers: PLS3 ({patient_data['PLS3_Expression_Level']}), SMN2 859G>C ({patient_data['SMN2_c.859G>C_Modifier']})")
    
    strategy = {
        "Primary_Recommendation": "",
        "Rationale": "",
        "ASO_Candidate": "N/A",
        "Small_Molecule_Candidate": "N/A"
    }
    
    # Logic matrix based on current standards + our novel findings
    if smn2_cn <= 2:
        # Severe (Type 1 or severe Type 2) - Needs immediate high-efficacy splicing correction
        # BBB permeability is critical for early CNS rescue.
        strategy["Primary_Recommendation"] = "Combination Therapy / High Affinity ASO"
        strategy["Rationale"] = "Low SMN2 copy number indicates severe phenotype. ASO designed for maximum RNA affinity and a small molecule helper is recommended."
        strategy["ASO_Candidate"] = "Novel ASO (+4 to +21 window) with optimized Tm: 61.31°C"
        strategy["Small_Molecule_Candidate"] = "Base_Risdiplam or Fluoro_Analog (for rapid systemic distribution)"
    elif smn2_cn == 3:
        # Type 2 or 3
        if patient_data['PLS3_Expression_Level'] == 'High' or patient_data['SMN2_c.859G>C_Modifier']:
            strategy["Primary_Recommendation"] = "Monotherapy (Small Molecule)"
            strategy["Rationale"] = "Protective modifiers present. Oral small molecule splicing modifier should suffice to maintain systemic SMN levels without invasive intrathecal ASO injections."
            strategy["ASO_Candidate"] = "None recommended (invasive)"
            strategy["Small_Molecule_Candidate"] = "Base_Risdiplam"
        else:
            strategy["Primary_Recommendation"] = "Standard Splicing Correction (ASO or Small Molecule)"
            strategy["Rationale"] = "Moderate phenotype. Either Spinraza analog or Risdiplam analog."
            strategy["ASO_Candidate"] = "Nusinersen Analog (+10 to +27)"
            strategy["Small_Molecule_Candidate"] = "Base_Risdiplam"
    else:
        # Type 3 or 4 (4+ copies)
        strategy["Primary_Recommendation"] = "Watchful Waiting / Low-dose systemic therapy"
        strategy["Rationale"] = "High SMN2 copy numbers generally result in a milder phenotype. Aggressive early intervention might not outweigh risks if asymptomatic."
        strategy["Small_Molecule_Candidate"] = "Low-dose Methoxy_Analog (simulated safe systemic profile)"

    print("\n[AI DECISION MATRIX OUTPUT]")
    print(json.dumps(strategy, indent=4))
    
    return strategy

if __name__ == "__main__":
    # Simulate a VCF parse and decision creation
    mock_data = parse_mock_vcf("dummy_patient.vcf")
    generate_intervention_strategy(mock_data)
