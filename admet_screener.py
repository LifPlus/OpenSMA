"""
Advanced ADMET & Toxicity Screener for Small Molecules (Module B)
------------------------------------------------------------------
Evaluates Risdiplam analogs against:
  - PAINS filtering (pan-assay interference patterns)
  - CNS MPO score (Wager et al. 2010) for BBB permeability
  - hERG cardiac toxicity heuristic
  - Reactive group / hepatotoxicity substructure screening
  - Ro5 (Lipinski's Rule of 5) compliance
"""
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

def ro5_check(mol):
    """Lipinski's Rule of Five: drug-like oral bioavailability check."""
    mw   = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    hbd  = Descriptors.NumHDonors(mol)
    hba  = Descriptors.NumHAcceptors(mol)
    violations = sum([mw > 500, logp > 5, hbd > 5, hba > 10])
    return violations, mw, logp, hbd, hba

def cns_mpo_score(mol):
    """
    CNS MPO scoring (Wager et al. 2010, ACS Chem. Neurosci.)
    Parameters: LogP, LogD(≈LogP), MW, TPSA, HBD, pKa-most-basic
    Score 0-6; ≥4 = CNS penetrant.
    """
    logp = Crippen.MolLogP(mol)
    mw   = Descriptors.MolWt(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd  = Descriptors.NumHDonors(mol)
    rings = rdMolDescriptors.CalcNumRings(mol)

    score = 0
    # LogP (optimal 2-4 for CNS)
    if 2.0 <= logp <= 4.0: score += 1
    elif logp < 2.0: score += 0.5
    # MW < 400 best, < 500 ok
    if mw < 400: score += 1
    elif mw < 500: score += 0.5
    # TPSA < 76 best, < 90 acceptable
    if tpsa < 76: score += 1
    elif tpsa < 90: score += 0.5
    # HBD ≤ 3
    if hbd <= 3: score += 1
    # Aromatic rings ≤ 3 (reduces metabolic burden)
    if rings <= 3: score += 1
    # Rotatable bonds ≤ 8 (Veber rule for CNS)
    rot_bonds = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if rot_bonds <= 8: score += 1

    return round(score, 2)

def pains_filter(mol):
    """Flag molecules with pan-assay interference patterns."""
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    return catalog.HasMatch(mol)

def herg_alert(mol):
    """
    Simple hERG cardiac toxicity heuristic.
    Alert triggered if molecule has:
    - Basic nitrogen (common in hERG blockers)
    - AND LogP > 3
    - AND aromatic ring count >= 2
    This is an approximation - formal hERG prediction needs ML models.
    """
    logp = Crippen.MolLogP(mol)
    rings = rdMolDescriptors.CalcNumAromaticRings(mol)
    # Check for any nitrogen (basic N is an hERG alert)
    has_basic_n = any(atom.GetAtomicNum() == 7 
                      and atom.GetFormalCharge() == 0 
                      and atom.GetTotalDegree() < 4 
                      for atom in mol.GetAtoms())
    
    if has_basic_n and logp > 3 and rings >= 2:
        return "ALERT"
    return "CLEAR"

def reactive_group_flag(mol):
    """
    Flag reactive substructures that can cause idiosyncratic hepatotoxicity.
    Checks for: epoxides, Michael acceptors, aldehydes, quinones.
    """
    alert_smarts = {
        "Aldehyde": "[CH]=O",
        "Epoxide": "[C]1O[C]1",
        "Michael_Acceptor": "C=CC=O",
        "Quinone": "O=C1C=CC(=O)C=C1"
    }
    flags = []
    for name, smarts in alert_smarts.items():
        pattern = Chem.MolFromSmarts(smarts)
        if pattern and mol.HasSubstructMatch(pattern):
            flags.append(name)
    return flags if flags else ["None"]

def analyze_molecules():
    # Try pipeline-generated molecules first, fall back to manual analogs
    input_file = None
    for fname in ['generated_molecules.csv', 'small_molecule_analogs.csv']:
        if os.path.exists(fname):
            input_file = fname
            break
    if input_file is None:
        print("Error: No molecule input file found.")
        return pd.DataFrame()
    df = pd.read_csv(input_file)
    print(f"Loaded {len(df)} small molecule candidates from {input_file}.\n")

    results = []
    for idx, row in df.iterrows():
        comp_name = row['Compound'] if 'Compound' in row else f"DeNovo_Cand_{idx+1}"
        mol = Chem.MolFromSmiles(row['SMILES'])
        if mol is None:
            print(f"Warning: Could not parse SMILES for {comp_name}")
            continue

        ro5_violations, mw, logp, hbd, hba = ro5_check(mol)
        cns_score = cns_mpo_score(mol)
        pains_flag = pains_filter(mol)
        herg = herg_alert(mol)
        react = reactive_group_flag(mol)
        tpsa = Descriptors.TPSA(mol)
        rot = rdMolDescriptors.CalcNumRotatableBonds(mol)

        # Overall CNS/BBB suitability
        bbb_ok = (cns_score >= 4 and not pains_flag and ro5_violations == 0)
        safety_flag = "SAFE" if (herg == "CLEAR" and react == ["None"] and not pains_flag) else "REVIEW"

        results.append({
            "Compound": comp_name,
            "SMILES": row['SMILES'],
            "MW": round(mw, 2),
            "LogP": round(logp, 2),
            "TPSA": round(tpsa, 2),
            "RotBonds": rot,
            "HBD": hbd,
            "HBA": hba,
            "Ro5_Violations": ro5_violations,
            "CNS_MPO_Score": cns_score,
            "PAINS_Flag": pains_flag,
            "hERG_Alert": herg,
            "Reactive_Groups": ", ".join(react),
            "BBB_Suitable": bbb_ok,
            "Safety_Rating": safety_flag
        })

    result_df = pd.DataFrame(results)

    print("=== SMALL MOLECULE ADMET ANALYSIS ===")
    print(result_df.to_string(index=False))

    result_df.to_csv('admet_screened_molecules.csv', index=False)
    print("\nFull ADMET results saved to 'admet_screened_molecules.csv'")
    return result_df

if __name__ == "__main__":
    analyze_molecules()
