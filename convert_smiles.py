"""Convert smiles strings in a CSV file into files."""
import argparse
import os
import re
import sys

from rdkit import Chem


def parse_args(args=None):
    p = argparse.ArgumentParser(description=sys.modules[__name__].__doc__)

    p.add_argument('--input', type=str,
        help='path to input file with one SMILES string per line')
    p.add_argument('--output', type=str,
        help='output directory to write files to')
    p.add_argument('--format', choices=['smi', 'sdf'], default='smi',
        help='output file format')
    p.add_argument('--delimiter', type=str, default=' ',
        help='input file column delimiter')
    
    return p.parse_args(args)


def convert(infile, outdir, fformat='smi', delimiter=' '):
    """Convert smiles strings in a CSV file into files.

    Parameters
    ----------
    infile : str
        Input CSV file.
    outdir : str
        Directory to which to write output files.
    fformat : {'smi','sdf'}
        File format in which to save strings. Default: ``'smi'``
    delimiter : str, optional
        Field delimiter in ``infile``. Default: ``','``
    """
    os.makedirs(outdir, exist_ok=True)
    
    with open(infile, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        # Assumes the name is the first field and string is the fifth
        name, _, _, _, smiles = line.split(delimiter)
        
        # Replace spaces with underscores for easier linux-ing
        name = re.sub(r'[\s]', '_', name)
        
        # If 'smi', save the string to file. If 'sdf', use rdkit tools to
        # convert the string to SDF format.
        if fformat == 'smi':
            with open(os.path.join(outdir, f'{name}', 'structure.smi')) as f:
                f.write(smiles)
        elif fformat == 'sdf':
            mol = Chem.MolFromSmiles(smiles)
            writer = Chem.SDWriter(os.path.join(outdir, f'{name}', 'structure.sdf'))
            writer.write(mol)


if __name__ == '__main__':
    args = parse_args()
    convert(args.input, args.output, fformat=args.format,
            delimiter=args.delimiter)
