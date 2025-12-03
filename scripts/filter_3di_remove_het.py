#!/usr/bin/env python3
"""
Filter a 3Di FASTA produced from a PDB so that HETATM positions are removed.
Usage:
  python scripts/filter_3di_remove_het.py <pdb_file> <3di_fasta> [--out OUT]

This reads the PDB (first chain), builds a mask of residues with empty hetero flag
(standard residues), and filters the 3Di sequence accordingly. The filtered
sequence overwrites the input 3Di FASTA by default (i.e. the original
`*_3di.fasta` will be replaced). Use `--out` to write to a different path.

This is a small post-processing convenience to fix existing generated 3Di files
that include placeholders for HETATM positions ("X...").
"""
import sys
from pathlib import Path
from Bio.PDB import PDBParser


def load_fasta(path):
    with open(path, 'r') as f:
        lines = [l.rstrip('\n') for l in f]
    if not lines:
        raise SystemExit('Empty file: %s' % path)
    header = lines[0] if lines[0].startswith('>') else None
    seq = ''.join(lines[1:]) if header else ''.join(lines)
    return header, seq


def write_fasta(path, header, seq):
    with open(path, 'w') as f:
        if header:
            f.write(header + '\n')
        f.write(seq + '\n')


def main():
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument('pdb')
    p.add_argument('fasta')
    p.add_argument('--out', '-o', default=None,
                   help='Output path. If omitted, the input FASTA is overwritten')
    args = p.parse_args()

    pdb_path = Path(args.pdb)
    fasta_path = Path(args.fasta)
    if not pdb_path.exists():
        raise SystemExit('PDB not found: %s' % pdb_path)
    if not fasta_path.exists():
        raise SystemExit('3Di FASTA not found: %s' % fasta_path)

    header, seq = load_fasta(fasta_path)

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('m', str(pdb_path))
    model = structure[0]
    chain = list(model.get_chains())[0]
    residues = list(chain.get_residues())

    if len(residues) != len(seq):
        # It's common that the generated 3Di has placeholders; but to filter we
        # expect per-residue mapping by chain entries. If lengths differ, we
        # still proceed but will align by chain length vs seq length trimming/padding.
        print('Warning: chain residue count (%d) != 3Di seq length (%d)'
              % (len(residues), len(seq)))

    non_het_mask = [len(r.id[0].strip()) == 0 for r in residues]

    filtered = []
    for i, keep in enumerate(non_het_mask):
        if i >= len(seq):
            break
        if keep:
            filtered.append(seq[i])

    filtered_seq = ''.join(filtered)

    out = args.out
    if out is None:
        # Overwrite the original 3Di FASTA as requested
        out = str(fasta_path)

    write_fasta(out, header, filtered_seq)
    print('Wrote filtered 3Di to', out)


if __name__ == '__main__':
    main()
