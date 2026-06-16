import os
import numpy as np
import json

def prepare_receptor_and_box(pdb_path):
    print(f"Loading {pdb_path}...")
    
    # Read the PDB file and split ligand vs receptor, taking only the first MODEL
    receptor_lines = []
    ligand_lines = []
    in_model_1 = False
    model_count = 0
    
    with open(pdb_path, "r") as f:
        for line in f:
            if line.startswith("MODEL"):
                model_count += 1
                if model_count == 1:
                    in_model_1 = True
                else:
                    break # Stop after MODEL 1
            if model_count > 1:
                break
                
            if line.startswith("ATOM"):
                receptor_lines.append(line)
            elif line.startswith("HETATM") and "Y59" in line:
                ligand_lines.append(line)
            elif line.startswith("TER"):
                receptor_lines.append(line)
            elif line.startswith("ENDMDL"):
                break # Stop after MODEL 1
    
    if not ligand_lines:
        print("Error: Ligand Y59 not found in PDB file. Searching for any HETATM as fallback...")
        with open(pdb_path, "r") as f:
            for line in f:
                if line.startswith("HETATM"):
                    ligand_lines.append(line)
    
    # Save receptor PDB
    receptor_pdb = "receptor_clean.pdb"
    with open(receptor_pdb, "w") as f:
        f.writelines(receptor_lines)
    print(f"Saved cleaned receptor to {receptor_pdb}")
    
    # Extract ligand coordinates to find the search box center and size
    ligand_coords = []
    for line in ligand_lines:
        try:
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            ligand_coords.append([x, y, z])
        except ValueError:
            continue
    
    if not ligand_coords:
        print("Error: Could not extract ligand coordinates.")
        return None

    ligand_coords = np.array(ligand_coords)
    center = np.mean(ligand_coords, axis=0)
    # Define box size (ligand span + padding)
    size = np.max(ligand_coords, axis=0) - np.min(ligand_coords, axis=0) + 12.0 # 12A padding
    
    print(f"Binding Pocket Center: {center}")
    print(f"AutoDock Vina Search Box Size: {size}")
    
    box_data = {
        "center_x": float(center[0]),
        "center_y": float(center[1]),
        "center_z": float(center[2]),
        "size_x": float(size[0]),
        "size_y": float(size[1]),
        "size_z": float(size[2])
    }
    
    with open("docking_box_real.json", "w") as f:
        json.dump(box_data, f, indent=4)
    print("Saved search box parameters to docking_box_real.json")
    
    # Manual PDBQT converter for receptor (add zero charge and atom types from element)
    receptor_pdbqt = "receptor_clean.pdbqt"
    with open(receptor_pdbqt, "w") as f:
        for line in receptor_lines:
            if line.startswith("ATOM"):
                atom_name = line[12:16].strip()
                element = line[76:78].strip() if len(line) > 76 else atom_name[0]
                # PDBQT line: PDB line up to 66 + partial charge + AD atom type
                pdbqt_line = line[:66].ljust(66) + f"  0.00  0.00    {element:>2}\n"
                f.write(pdbqt_line)
        f.write("TER\nEND\n")
    print(f"Saved receptor PDBQT to {receptor_pdbqt}")
    
    return box_data

if __name__ == "__main__":
    prepare_receptor_and_box("smn2_receptor_real.pdb")
