"""Extract the atoms in a chain of a PDB file."""
import argparse
import logging
import os
import re

import tqdm


logging.basicConfig(format='%(level)s %(asctime)s: %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S%p')


def parse_args(args=None):
    """Get command line arguments.

    Parameters
    ----------
    args : str, optional
        The arguments to parse. If not provided, parses args from `sys.argv`.

    Returns
    -------
    args : `argparse.Namespace`
        An object with ``input``, ``output``, and ``chain`` attributes.
    """
    p = argparse.ArgumentParser(
        description='extract the atoms in a chain of a PDB file')

    p.add_argument('--input', type=str, help='path to input PDB file')
    p.add_argument('--output', type=str, help='path to output PDB file')
    p.add_argument('--chain', type=str, help='chain to extract (e.g., A, B)')

    return p.parse_args(args)


def main(infile, outfile, chain):
    """Extract data about a chain in a PDB file.

    Parameters
    ----------
    infile : str
      Input PDB file to parse.
    outfile : str
      Output PDB to write substructure to.
    chain : str
      The chain to extract from the PDB file.
    """
    lines = []
    
    logging.info(f'Reading {infile}')
    with open(infile, 'r') as f:
        lines = list(f.readlines())
    
    if len(lines) == 0:
      raise ValueError('Empty file provided as input')

    logging.info('Setting up line and data collection')
    outlines = []
    atoms = []
    field_names = {
      'REMARK': 0,
      'HET': 0,
      'HELIX': 0,
      'SHEET': 0,
      'SITE': 0,
      'ORIGX': 0,
      'SCALE': 0,
      'MTRIX': 0,
      'ATOM': 0,
      'HETATM': 0,
      'TER': 0,
      'CONECT': 0,
      'SEQRES': 0
    }
    last_link = -1
    terminus = False
    
    logging.info('Parsing for chain {chain}')
    for idx, line in enumerate(tqdm(lines)):
        fields = line.split()

        if fields[0] in ['DBREF', 'SEQRES', 'HET'] and fields[2] == chain:
            outlines.append(line)
            if fields[0] in field_names:
                field_names[fields[0]] += 1
        elif fields[0] in ['HETNAM', 'FORMUL']:
            outlines.append(line)
            if fields[0] in field_names:
                field_names[fields[0]] += 1
        elif fields[0] == 'HELIX' and fields[4] == chain and fields[7] == chain:
            outlines.append(line)
            if fields[0] in field_names:
                field_names[fields[0]] += 1
        elif fields[0] == 'LINK' and fields[3] == chain and fields[7] == chain and last_link == -1 or idx == last_link + 1:
            outlines.append(line)
            if fields[0] in field_names:
                field_names[fields[0]] += 1
            last_link = idx
        elif re.search(r'^(CRYST|ORIGX|SCALE|MTRIX)[\d]+$', fields[0]):
            outlines.append(line)
            if fields[0][:-1] in field_names:
                field_names[fields[0][:-1]] += 1
        elif fields[0] in ['ATOM', 'HETATM'] and fields[4] == chain and not terminus:
            outlines.append(line)
            if fields[0] in field_names:
                field_names[fields[0]] += 1
            atoms.append(fields[1])
        elif fields[0] == 'TER' and fields[3] == chain and not terminus:
            outlines.append(line)
            if fields[0] in field_names:
                field_names[fields[0]] += 1
            terminus = True
        elif fields[0] == 'CONECT' and fields[1] in atoms:
            outlines.append(line)
            if fields[0] in field_names:
                field_names[fields[0]] += 1

    logging.info('Creating MASTER metadata line.')
    masterline = ''.join([
        'MASTER    ',
        str(field_names['REMARK']).rjust(5, ' '),
        '0'.rjust(5, ' '),
        str(field_names['HET']).rjust(5, ' '),
        str(field_names['HELIX']).rjust(5, ' '),
        str(field_names['SHEET']).rjust(5, ' '),
        '0'.rjust(5, ' '),
        str(field_names['SITE']).rjust(5, ' '),
        str(field_names['ORIGX'] + field_names['SCALE'] + field_names['MTRIX']).rjust(5, ' '),
        str(field_names['ATOM'] + field_names['HETATM']).rjust(5, ' '),
        str(field_names['TER']).rjust(5, ' '),
        str(field_names['CONECT']).rjust(5, ' '),
        str(field_names['SEQRES']).rjust(5, ' '),
        '\n',
    ])

    outlines.append(masterline)
    outlines.append('END   ')

    logging.info(f'Writing chain {chain} to {outfile}.')
    with open(outfile, 'w') as f:
        f.write(''.join(outlines))

    logging.info('Done')
    

if __name__ == '__main__':
    args = parse_args()
    main(args.input, args.output, args.chain)
    
