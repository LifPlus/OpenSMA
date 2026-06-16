import requests
import os

def download_pdb(pdb_id, output_path):
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    print(f"Downloading {pdb_id} from RCSB...")
    response = requests.get(url)
    if response.status_code == 200:
        with open(output_path, "w") as f:
            f.write(response.text)
        print(f"Successfully downloaded {pdb_id} to {output_path}")
        return True
    else:
        print(f"Failed to download {pdb_id}. Status code: {response.status_code}")
        return False

if __name__ == "__main__":
    # 8R62 is the Solution structure of Risdiplam bound to the RNA duplex
    pdb_id = "8R62"
    output_file = "smn2_receptor_real.pdb"
    download_pdb(pdb_id, output_file)
