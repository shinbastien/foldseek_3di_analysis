#!/usr/bin/env python3
import os
import sys
from glob import glob

# Add project root to path so we can import training_3di_gpu
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

from training_3di_gpu_new.extract_pdb_features import (
    get_coords_from_pdb,
    find_nearest_residues,
    calc_angles_forloop
)
from Bio.PDB import PDBParser

INPUT_DIR = os.path.dirname(__file__)

pdb_files = sorted(glob(os.path.join(INPUT_DIR, '*.pdb')))
if not pdb_files:
    print('No PDB files found in', INPUT_DIR)
    sys.exit(1)

parser = PDBParser(QUIET=True)

for pdb_path in pdb_files:
    base = os.path.splitext(os.path.basename(pdb_path))[0]
    out_path = os.path.join(INPUT_DIR, base + '.3di.txt')

    print('Processing', pdb_path, '->', out_path)
    # read coords using existing helper
    coords, valid_mask = get_coords_from_pdb(pdb_path, full_backbone=True)

    # find partner index for each residue (1st nearest neighbor)
    partner_idx = find_nearest_residues(coords, valid_mask, k=1, min_seq_dist=1)

    angles, new_valid_mask = calc_angles_forloop(coords, partner_idx, valid_mask)

    # write output: one line per residue (only residues with valid features)
    with open(out_path, 'w') as fo:
        header = ['idx', 'resname', 'resnum', 'partner_idx'] + [f'f{i+1}' for i in range(angles.shape[1])]
        fo.write('\t'.join(header) + '\n')

        # need residue metadata: parse PDB to get residues
        structure = parser.get_structure(base, pdb_path)
        model = structure[0]
        chain = list(model.get_chains())[0]
        residues = [r for r in chain.get_residues() if len(r.id[0].strip()) == 0]

        n_res = len(residues)
        for i in range(n_res):
            if not new_valid_mask[i]:
                continue
            res = residues[i]
            resname = res.resname
            resnum = res.id[1]
            p_idx = int(partner_idx[i]) if not (partner_idx is None) else -1
            feats = angles[i]
            line = [str(i), resname, str(resnum), str(p_idx)] + [f"{float(x):.6g}" for x in feats]
            fo.write('\t'.join(line) + '\n')

print('Done.')
