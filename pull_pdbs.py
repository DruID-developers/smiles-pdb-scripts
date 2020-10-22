"""Download PDB files by ID."""
import argparse
import logging
import sys

import Bio
from Bio.PDB import PDBList
import tqdm


logging.basicConfig(format='%(level)s %(asctime)s: %(message)s',
                    datefmt='%Y-%m-%d %I:%M:%S%p')


def parse_args(args=None):
    p = argparse.ArgumentParser(description=sys.modules[__name__].__doc__)

    p.add_argument('--input', type=str,
                   help='text file with one PDB ID per line')
    p.add_argument('--output', type=str,
                   help='directory to write output files to')
    p.add_argument('--server', type=str, default='ftp://ftp.wwpdb.org',
                   help='database URL to pull PDBS from')
    p.add_argument('--file_format', type=str, default='mmCiF',
                   choices=['mmCiF', 'pdb', 'xml', 'mmtf', 'bundle'],
                   help='the file format to pull')
    
    return p.parse_args(args)


def download_pdbs(id_file, out_dir, server='ftp://ftp.wwpdb.org',
                  file_format=None):
    """Download molecular files by PDB ID.

    Parameters
    ----------
    id_file : str
        Path to a text file with one PDB ID per line.
    out_dir : str
        Path to the directory to save downloaded files to.
    server : str
        URL of the server to download files from.
    file_format : {'mmCiF', 'pdb', 'xml', 'mmtf', 'bundle}
        The file format to download.

    Notes
    -----
    This function will directly download and save the files to ``out_dir`` with
    no return value or other non-logging side effects.
    """
    logging.info(f'Connecting to {server}')
    pdbl = PDBList(server=server)

    PDBlist2 = []

    logging.info('Reading PDB IDs from ')
    with open(id_file, 'r') as f:
        PDBlist2 = f.readlines()

    logging.info('Downloading pdbs...')
    pdb_ids = tqdm(PDBList2)
    for pdb_id in pdb_ids:
        pdb_ids.set_postfix_str(f'Downloading {pdb_id}')
        pdbl.retrieve_pdb_file(pdb_id.strip(), pdir=out_dir)

    logging.info('Done')

if __name__ == '__main__':
    args = parse_args()
    download_pdbs(
        args.input,
        args.output,
        database=args.database,
        file_format=args.file_format)
