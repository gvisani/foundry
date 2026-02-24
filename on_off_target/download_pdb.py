


"""Download PDB files from the RCSB Protein Data Bank."""

import argparse
import urllib.request
from pathlib import Path


def download_pdb(pdb_id: str, output_dir: str) -> Path:
    """
    Download a PDB file from RCSB and save it to the specified directory.

    Args:
        pdb_id: 4-character PDB identifier (e.g., '1abc')
        output_dir: Directory where the PDB file will be saved

    Returns:
        Path to the downloaded file

    Raises:
        ValueError: If pdb_id is not 4 characters
        urllib.error.URLError: If download fails
    """
    pdb_id = pdb_id.lower()
    if len(pdb_id) != 4:
        raise ValueError(f"PDB ID must be 4 characters, got: {pdb_id}")

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
    filepath = output_path / f"{pdb_id}.pdb"

    urllib.request.urlretrieve(url, filepath)
    print(f"Downloaded {pdb_id}.pdb to {filepath}")

    return filepath


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Download a PDB file from RCSB")
    parser.add_argument("--pdbid", required=True, help="4-character PDB identifier (e.g., 1abc)")
    parser.add_argument("--out_dir", required=False, default="./input_pdbs", help="Directory to save the PDB file")

    args = parser.parse_args()
    download_pdb(args.pdbid, args.out_dir)