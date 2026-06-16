import os
import random
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem import Crippen
from rdkit.Chem import QED

import pandas as pd
import numpy as np

# Force field for 3D conformers
from rdkit.Chem import rdForceFieldHelpers

def score_molecule(mol):
    """
    Calculates a fitness score for an SMA drug candidate.
    Prioritizes Blood-Brain Barrier (BBB) penetration, QED (Drug-likeness), and small MW.
    """
    if mol is None:
        return 0.0
        
    try:
        Chem.SanitizeMol(mol)
    except:
        return 0.0

    mw = Descriptors.MolWt(mol)
    logp = Crippen.MolLogP(mol)
    tpsa = Descriptors.TPSA(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    qed = QED.qed(mol)

    score = qed * 100 # Base score from QED (0-100)
    
    # BBB Penetrability Heuristics
    if 2.0 <= logp <= 5.0:
        score += 20
    if tpsa < 90:
        score += 20
    if mw < 500:
        score += 10
        
    # Penalties
    if hbd > 3: score -= 10
    if hba > 7: score -= 10
    if logp > 5.5: score -= 20
    if mw > 600: score -= 30
    
    # --- PHASE 2: 3D Conformer & Docking Heuristic ---
    mol_3d = Chem.AddHs(mol)
    try:
        # Generate 3D coordinates
        AllChem.EmbedMolecule(mol_3d, randomSeed=42)
        # Minimize energy using MMFF94 (Proxy for internal steric strain / docking prep)
        ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_3d, rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_3d))
        if ff:
            ff.Minimize()
            energy = ff.CalcEnergy()
            # Favorable (lower) energy gives a slight boost, high strain penalizes
            if energy < 150:
                score += 15
            elif energy > 300:
                score -= 10
    except:
        # If 3D embedding fails, penalize structurally (likely too rigid/strained)
        score -= 15
    # --------------------------------------------------
        
    return max(0, score)

def mutate_molecule(mol):
    """
    Applies a random SMARTS reaction to mutate the molecule (Generative step).
    """
    # Define common medicinal chemistry transformations
    rxn_smarts = [
        '[c:1][H]>>[c:1]F',           # Add Fluorine to aromatic Carbon
        '[c:1][H]>>[c:1]C',           # Add Methyl
        '[c:1][H]>>[c:1]OC',          # Add Methoxy
        '[c:1][H]>>[c:1]Cl',          # Add Chlorine
        '[c:1][H]>>[n:1]',            # C to N in aromatic ring
        '[CX3:1](=[OX1:2])[H]>>[CX3:1](=[OX1:2])N' # Aldehyde to Amide (if any)
    ]
    
    rxn = AllChem.ReactionFromSmarts(random.choice(rxn_smarts))
    try:
        products = rxn.RunReactants((mol,))
        valid_mols = []
        for prod in products:
            p_mol = prod[0]
            try:
                Chem.SanitizeMol(p_mol)
                valid_mols.append(p_mol)
            except:
                pass
        
        if valid_mols:
            return random.choice(valid_mols)
    except:
        pass
        
    return None

def run_generative_algorithm(base_smiles, generations=5, population_size=50):
    print(f"Initializing Generative AI (Genetic Algorithm) for De Novo Drug Design...")
    base_mol = Chem.MolFromSmiles(base_smiles)
    
    if base_mol is None:
        print("Error: Invalid Base SMILES.")
        return []

    base_score = score_molecule(base_mol)
    print(f"Base Molecule QED+BBB Fitness Score: {base_score:.2f}")

    # Initial population: start with copies of the base molecule
    population = [Chem.Mol(base_mol) for _ in range(population_size)]
    
    best_candidates = {}
    best_candidates[base_smiles] = base_mol

    for gen in range(generations):
        print(f"--- Generation {gen + 1}/{generations} ---")
        new_population = []
        
        # Mutate
        for mol in population:
            # Apply 1 to 3 mutations
            mutant = Chem.Mol(mol)
            num_mutations = random.randint(1, 3)
            for _ in range(num_mutations):
                new_mutant = mutate_molecule(mutant)
                if new_mutant is not None:
                    mutant = new_mutant
            
            if mutant is not None:
                new_population.append(mutant)
            
        # Add some base molecules back to preserve core traits (Elitism)
        new_population.extend([Chem.Mol(base_mol) for _ in range(int(population_size * 0.2))])
        
        # Evaluate and Select
        scored_pop = []
        for mol in new_population:
            if mol is None:
                continue
            smiles = Chem.MolToSmiles(mol)
            if smiles not in best_candidates: # Unique only
                sc = score_molecule(mol)
                if sc > 0:
                    scored_pop.append((sc, smiles, mol))
                    
        # Sort by score descending
        scored_pop.sort(key=lambda x: x[0], reverse=True)
        
        # Keep top N for next generation
        population = [x[2] for x in scored_pop[:population_size]]
        
        # Save to hall of fame
        for sc, smiles, mol in scored_pop[:10]:
            best_candidates[smiles] = mol
            
    # Final evaluation of unique hall of fame
    final_results = []
    for smiles, mol in best_candidates.items():
        score = score_molecule(mol)
        mw = Descriptors.MolWt(mol)
        logp = Crippen.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        qed = QED.qed(mol)
        
        bbb_score = "HIGH" if (2.0 <= logp <= 5.0 and tpsa < 90 and mw < 500) else "MODERATE/LOW"
        
        # Calculate 3D Energy for display
        mol_3d = Chem.AddHs(mol)
        energy_val = "N/A"
        try:
            if AllChem.EmbedMolecule(mol_3d, randomSeed=42) == 0:
                ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol_3d, rdForceFieldHelpers.MMFFGetMoleculeProperties(mol_3d))
                if ff:
                    ff.Minimize()
                    energy_val = round(ff.CalcEnergy(), 1)
        except:
            pass

        final_results.append({
            "SMILES": smiles,
            "Fitness_Score": round(score, 2),
            "QED": round(qed, 3),
            "MW": round(mw, 2),
            "LogP": round(logp, 2),
            "TPSA": round(tpsa, 2),
            "3D_Energy_kcal/mol": energy_val,
            "BBB": bbb_score
        })
        
    return pd.DataFrame(final_results)

def generate_risdiplam_analogs():
    # Known robust Risdiplam core representation
    risdiplam_smiles = "CC1=CN=C2N1N=CC(=C2C3=CC4=C(C=C3)C(=O)N=C(N4)N5CCC6(CC5)CC6)C"
    
    # Run generative algorithm
    df = run_generative_algorithm(risdiplam_smiles, generations=5, population_size=100)
    
    # Sort and take top 10 novel molecules
    df_sorted = df.sort_values(by="Fitness_Score", ascending=False).head(10)
    
    print("\n=== Top 10 De Novo Generated SMA Candidates ===")
    print(df_sorted.to_string(index=False))
    
    df_sorted.to_csv("generated_molecules.csv", index=False)
    print("\nSaved Top De Novo candidates to 'generated_molecules.csv'.")
    
    # Generate an image of the very best novel generator
    top_smiles = df_sorted.iloc[0]["SMILES"]
    top_mol = Chem.MolFromSmiles(top_smiles)
    if top_mol:
        img = Draw.MolToImage(top_mol, size=(400, 400))
        img.save("top_denovo_molecule_2d.png")
        print("Saved Top De Novo molecule 2D structure to 'top_denovo_molecule_2d.png'.")

if __name__ == "__main__":
    generate_risdiplam_analogs()
