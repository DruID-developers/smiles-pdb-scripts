import argparse

import Bio
from Bio.PDB import PDBList

pdbl = PDBList()

PDBlist2=[]

with open('./druid_pdblist.txt', 'r') as f:
    PDBlist2 = f.readlines()

for i in PDBlist2:
    pdbl.retrieve_pdb_file(i.strip(), pdir='pdbxs')
