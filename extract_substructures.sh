#!/usr/bin/env bash

module load python/3.7.3

mkdir -p substructures

for pdb in pdbs/*; do
    if [ -f $pdb ]; then
        echo "Extracting from $pdb"
        fname="$(basename $pdb)"
        python extract_substructure.py \
            --input $pdb \
            --output substructures/$fname \
            --chain A
    fi
done
