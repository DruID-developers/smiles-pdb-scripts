# DruID SMILES/PDB Scripts

## Dependencies

- Python 3.6+
- Biopython
- rdkit
- tqdm

## Installation

`git clone https://github.com/jeffkinnison/smiles-pdb-scripts`

## Downloading PDBs

1. Put PDB IDs into a file with one ID per line
2. `python pull_pdbs.py --input <id_file> --output <directory>`
3. The pulled PDBs will now be in `<directory>`

## Extracting Individual Chains

`python --input <pdb_file> --output <output_pdb_file> --chain <chain_id>`
