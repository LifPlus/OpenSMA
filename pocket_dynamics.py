"""
OpenSMA Phase 7 — Metropolis Monte Carlo (MMC) Molecular Dynamics Binding Stability Simulator
=============================================================================================
Performs a simulated 310 K (body temperature) thermodynamic binding run of small molecules
inside the SMN2 pocket. Tracks ligand root-mean-square deviation (RMSD) to detect structural
stability or pocket escape.
"""

import os
import json
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem

def parse_pdb_coords(pdb_path):
    coords = []
    if not os.path.exists(pdb_path):
        return coords
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    element = line[76:78].strip()
                    if not element:
                        element = line[12:16].strip()[0]
                    coords.append((x, y, z, element))
                except Exception:
                    continue
    return coords

def get_aligned_ligand_coords(smiles):
    # 1. Get Risdiplam bioactive reference conformer
    ref_smiles = "CC1=CN=C2N1N=CC(=C2C3=CC4=C(C=C3)C(=O)N=C(N4)N5CCC6(CC5)CC6)C"
    ref_mol = Chem.MolFromSmiles(ref_smiles)
    ref_mol = Chem.AddHs(ref_mol)
    AllChem.EmbedMolecule(ref_mol, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(ref_mol)
    
    # 2. Prepare candidate molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None, None
    mol = Chem.AddHs(mol)
    
    # 3. Generate 3D conformer and align to Risdiplam using O3A
    conf_ids = AllChem.EmbedMultipleConfs(mol, numConfs=10, randomSeed=42)
    if not conf_ids:
        if AllChem.EmbedMolecule(mol, randomSeed=42) != 0:
            return None, None, None
        conf_ids = [0]
    AllChem.MMFFOptimizeMoleculeConfs(mol)
    
    best_rmsd = 999.0
    best_conf_id = 0
    from rdkit.Chem import rdMolAlign
    for conf_id in conf_ids:
        try:
            o3a = rdMolAlign.GetO3A(mol, ref_mol, prbCid=conf_id, refCid=0)
            rmsd = o3a.Align()
            if rmsd < best_rmsd:
                best_rmsd = rmsd
                best_conf_id = conf_id
        except Exception:
            continue
            
    # 4. Get coordinates of the best aligned conformer
    conformer = mol.GetConformer(best_conf_id)
    ligand_atoms = []
    for i, atom in enumerate(mol.GetAtoms()):
        pos = conformer.GetAtomPosition(i)
        ligand_atoms.append({
            "idx": i,
            "element": atom.GetSymbol(),
            "x": pos.x,
            "y": pos.y,
            "z": pos.z
        })
    return mol, best_conf_id, ligand_atoms

def calculate_interaction_energy(ligand_coords, pocket_coords):
    total_energy = 0.0
    for l in ligand_coords:
        lx, ly, lz = l["x"], l["y"], l["z"]
        for px, py, pz, pelem in pocket_coords:
            dx = lx - px
            dy = ly - py
            dz = lz - pz
            r2 = dx*dx + dy*dy + dz*dz
            r = np.sqrt(r2)
            if r < 1.2:
                total_energy += 150.0  # heavy clash penalty
                continue
            
            # Lennard-Jones 12-6
            sigma = 3.2
            eps = 0.15
            term6 = (sigma / r) ** 6
            term12 = term6 ** 2
            v_lj = 4.0 * eps * (term12 - term6)
            total_energy += v_lj
            
            # Simple hydrogen bond attraction
            if l["element"] in ["O", "N"] and pelem in ["O", "N"]:
                if 2.5 <= r <= 3.5:
                    total_energy -= 0.6  # HB stabilization
                    
    return total_energy

def apply_perturbation(ligand_coords, dx, dy, dz, rx, ry, rz):
    xs = [l["x"] for l in ligand_coords]
    ys = [l["y"] for l in ligand_coords]
    zs = [l["z"] for l in ligand_coords]
    cx = np.mean(xs)
    cy = np.mean(ys)
    cz = np.mean(zs)
    
    # Rotation matrices
    cx_rot = np.cos(rx); sx_rot = np.sin(rx)
    Rx = np.array([[1, 0, 0], [0, cx_rot, -sx_rot], [0, sx_rot, cx_rot]])
    
    cy_rot = np.cos(ry); sy_rot = np.sin(ry)
    Ry = np.array([[cy_rot, 0, sy_rot], [0, 1, 0], [-sy_rot, 0, cy_rot]])
    
    cz_rot = np.cos(rz); sz_rot = np.sin(rz)
    Rz = np.array([[cz_rot, -sz_rot, 0], [sz_rot, cz_rot, 0], [0, 0, 1]])
    
    R = Rz @ Ry @ Rx
    
    new_coords = []
    for l in ligand_coords:
        ox = l["x"] - cx
        oy = l["y"] - cy
        oz = l["z"] - cz
        
        v = np.array([ox, oy, oz])
        v_rot = R @ v
        
        nx = v_rot[0] + cx + dx
        ny = v_rot[1] + cy + dy
        nz = v_rot[2] + cz + dz
        
        new_coords.append({
            "idx": l["idx"],
            "element": l["element"],
            "x": nx,
            "y": ny,
            "z": nz
        })
    return new_coords

def calculate_rmsd(coords_a, coords_b):
    sq_dist_sum = 0.0
    for a, b in zip(coords_a, coords_b):
        dx = a["x"] - b["x"]
        dy = a["y"] - b["y"]
        dz = a["z"] - b["z"]
        sq_dist_sum += dx*dx + dy*dy + dz*dz
    return np.sqrt(sq_dist_sum / len(coords_a))

def run_metropolis_monte_carlo(smiles, receptor_pdb, num_steps=1000, temp=310.0):
    """
    Runs a Metropolis Monte Carlo molecular dynamics stability simulation at 310 K.
    """
    np.random.seed(42)
    # 1. Parse receptor coordinates
    receptor_coords = parse_pdb_coords(receptor_pdb)
    if not receptor_coords:
        print("Error: Receptor coords not found.")
        return None
        
    # 2. Align ligand and get 3D coordinates
    mol, best_conf, init_coords = get_aligned_ligand_coords(smiles)
    if not init_coords:
        print("Error: Conformer generation or alignment failed.")
        return None
        
    # 3. Filter receptor pocket atoms near the ligand center to speed up
    l_xs = [l["x"] for l in init_coords]
    l_ys = [l["y"] for l in init_coords]
    l_zs = [l["z"] for l in init_coords]
    lc_x, lc_y, lc_z = np.mean(l_xs), np.mean(l_ys), np.mean(l_zs)
    
    pocket_coords = []
    for rx, ry, rz, relem in receptor_coords:
        dist = np.sqrt((rx - lc_x)**2 + (ry - lc_y)**2 + (rz - lc_z)**2)
        if dist <= 14.0:
            pocket_coords.append((rx, ry, rz, relem))
            
    if not pocket_coords:
        print("Warning: No pocket atoms found within 14 Å. Using full receptor.")
        pocket_coords = receptor_coords
        
    # 4. MMC Parameters
    kb = 0.001987  # kcal/(mol*K)
    kt = kb * temp  # ~0.616 kcal/mol at 310 K
    
    current_coords = list(init_coords)
    current_energy = calculate_interaction_energy(current_coords, pocket_coords)
    
    rmsd_traj = []
    energy_traj = []
    
    accepted_steps = 0
    
    for step in range(num_steps):
        # Sample small translation and rotation
        dx, dy, dz = np.random.uniform(-0.06, 0.06, 3)
        rx, ry, rz = np.random.uniform(-0.025, 0.025, 3)
        
        trial_coords = apply_perturbation(current_coords, dx, dy, dz, rx, ry, rz)
        trial_energy = calculate_interaction_energy(trial_coords, pocket_coords)
        
        dE = trial_energy - current_energy
        
        if dE <= 0.0:
            accept = True
        else:
            p = np.exp(-dE / kt)
            accept = np.random.uniform(0.0, 1.0) < p
            
        if accept:
            current_coords = trial_coords
            current_energy = trial_energy
            accepted_steps += 1
            
        rmsd = calculate_rmsd(current_coords, init_coords)
        rmsd_traj.append(rmsd)
        energy_traj.append(current_energy)
        
    final_rmsd = rmsd_traj[-1]
    stability = "STABLE (Good binding)" if final_rmsd < 2.2 else "STRETCHED (Partial drift)" if final_rmsd < 4.0 else "UNSTABLE (Pocket escape)"
    
    return {
        "rmsd_traj": rmsd_traj,
        "energy_traj": energy_traj,
        "final_rmsd": round(final_rmsd, 2),
        "avg_rmsd": round(float(np.mean(rmsd_traj)), 2),
        "stability": stability,
        "acceptance_ratio": round(accepted_steps / num_steps, 2)
    }

def run_pocket_simulation_pipeline():
    print("=" * 80)
    print("  OPENSMA Phase 7 — Metropolis Monte Carlo 3D Binding Stability Simulator")
    print("=" * 80)
    
    # Load docked molecules
    if not os.path.exists("docked_candidates_ranked.csv"):
        print("Error: docked_candidates_ranked.csv not found. Run docking first.")
        return
        
    df = pd.read_csv("docked_candidates_ranked.csv")
    df_active = df[df.get("Fitness_Score", "") != "REF"].copy()
    if df_active.empty:
        df_active = df
        
    receptor_pdb = "receptor_clean.pdb"
    if not os.path.exists(receptor_pdb):
        print("Error: receptor_clean.pdb not found.")
        return
        
    results = []
    
    # Run MMC on each candidate
    for idx, row in df_active.head(5).iterrows():
        smiles = str(row["SMILES"])
        comp_name = row.get("Compound")
        if pd.isna(comp_name) or not comp_name:
            comp_name = f"Candidate_{idx+1}"
            
        print(f"\nRunning 310 K Metropolis Monte Carlo Simulation for {comp_name}...")
        res = run_metropolis_monte_carlo(smiles, receptor_pdb, num_steps=1000)
        if res:
            print(f"  → Final RMSD: {res['final_rmsd']} Å")
            print(f"  → Binding Stability: {res['stability']}")
            print(f"  → Acceptance Ratio: {res['acceptance_ratio']}")
            
            results.append({
                "Compound": comp_name,
                "SMILES": smiles,
                "Vina_Affinity_kcal/mol": row.get("Vina_Affinity_kcal/mol"),
                "MMC_Final_RMSD_A": res["final_rmsd"],
                "MMC_Avg_RMSD_A": res["avg_rmsd"],
                "MMC_Stability": res["stability"],
                "Acceptance_Ratio": res["acceptance_ratio"]
            })
            
            # Plot individual RMSD trajectory
            try:
                import matplotlib.pyplot as plt
                plt.figure(figsize=(7, 4))
                plt.plot(res["rmsd_traj"], color='#8E24AA', lw=1.5, label='MMC Trajectory')
                plt.axhline(y=2.2, color='green', linestyle='--', alpha=0.6, label='Stable threshold (<2.2 Å)')
                plt.axhline(y=4.0, color='red', linestyle='--', alpha=0.6, label='Escape threshold (>4.0 Å)')
                plt.xlabel("Monte Carlo Step")
                plt.ylabel("RMSD from Docked Pose (Å)")
                plt.title(f"310 K Binding Stability Run: {comp_name}")
                plt.legend(facecolor='#1E2732', labelcolor='white')
                plt.grid(True, alpha=0.15)
                # save plot
                fname = f"mmc_rmsd_{comp_name.replace(' ', '_')}.png"
                plt.savefig(fname, dpi=120, bbox_inches='tight')
                plt.close()
            except Exception as e:
                pass
                
    df_out = pd.DataFrame(results)
    df_out.to_csv("mmc_stability_results.csv", index=False)
    print("\nSaved MMC Stability results to 'mmc_stability_results.csv'")

if __name__ == "__main__":
    run_pocket_simulation_pipeline()
