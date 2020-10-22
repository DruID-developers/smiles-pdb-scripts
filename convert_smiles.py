import argparse
import os
import sys

from rdkit import Chem


def parse_args(args=None):
    p = argparse.ArgumentParser(description=sys.modules[__name__].__doc__)

    p.add_argument('--input', type=str,
        help='path to input file with one SMILES string per line')
    p.add_argument('--output', type=str,
        help='output directory to write SDFs to')
    p.add_argument('--delimiter', type=str, default=' ',
        help='input file column delimiter')
    
    return p.parse_args(args)


def convert(infile, outdir, delimiter=' '):
    os.makedirs(outdir, exist_ok=True)
    
    with open(infile, 'r') as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        smiles = line  # , name = line.split(delimiter)
        mol = Chem.MolFromSmiles(smiles)
        writer = Chem.SDWriter(os.path.join(outdir, f'{str(i).zfill(5)}.sdf'))
        writer.write(mol)


if __name__ == '__main__':
    args = parse_args()
    convert(args.input, args.output, delimiter=args.delimiter)
