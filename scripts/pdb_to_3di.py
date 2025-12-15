#!/usr/bin/env python3
"""
Convert PDB file to 3Di sequence.

Usage:
    python pdb_to_3di.py <pdb_file> <model_dir>
    
Example:
    python pdb_to_3di.py my_protein.pdb training_3di_gpu_new/tmp
    
The model_dir should contain:
    - encoder.pt (trained encoder model)
    - states.txt (cluster centroids)
"""

import numpy as np
import sys
import os
import argparse
from pathlib import Path

import torch

# Add project root to path
ROOT = os.path.abspath(os.path.dirname(__file__))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

import importlib

# Try to import extract_pdb_features from one of the training_3di_gpu variants that exist
extract_pdb_features = None
_candidates = [
    'encoders_and_tools.training_3di_gpu_9f',
    'encoders_and_tools.training_3di_gpu_8f',
    'encoders_and_tools.training_3di_gpu_10f',
]
for _cand in _candidates:
    try:
        mod = importlib.import_module(_cand + '.extract_pdb_features')
        extract_pdb_features = mod
        # attach a note for debugging
        # print(f"[DEBUG] using extract_pdb_features from {_cand}")
        break
    except Exception:
        continue

if extract_pdb_features is None:
    raise ModuleNotFoundError(
        "Could not find a package providing extract_pdb_features.\n"
        "Searched: %s\n" % (', '.join(_candidates)) +
        "Make sure one of those folders is on PYTHONPATH and contains extract_pdb_features.py"
    )

# 50 letters (X/x are missing)
LETTERS = 'ABCDEFGHIJKLMNOPQRSTUVWYZabcdefghijklmnopqrstuvwyz'


def validate_model_files(model_dir):
    """
    Validate that encoder.pt and states.txt exist in model_dir.
    Returns paths to encoder and states files, or raises error.
    """
    model_dir = Path(model_dir)
    
    if not model_dir.exists():
        raise FileNotFoundError(f"❌ Error: Model directory does not exist: {model_dir}")
    
    encoder_path = model_dir / "encoder.pt"
    states_path = model_dir / "states.txt"
    
    errors = []
    if not encoder_path.exists():
        errors.append(f"  - encoder.pt not found in {model_dir}")
    if not states_path.exists():
        errors.append(f"  - states.txt not found in {model_dir}")
    
    if errors:
        error_msg = "❌ Error: Required model files are missing:\n" + "\n".join(errors)
        raise FileNotFoundError(error_msg)
    
    print(f"✓ Found encoder model: {encoder_path}")
    print(f"✓ Found centroids: {states_path}")
    
    return encoder_path, states_path


def validate_pdb_file(pdb_path):
    """
    Validate that PDB file exists.
    """
    pdb_path = Path(pdb_path)
    
    if not pdb_path.exists():
        raise FileNotFoundError(f"❌ Error: PDB file does not exist: {pdb_path}")
    
    if pdb_path.suffix.lower() != '.pdb':
        raise ValueError(f"❌ Error: File must have .pdb extension: {pdb_path}")
    
    print(f"✓ Found PDB file: {pdb_path}")
    return pdb_path


def extract_features(pdb_path, virt_cb=(270, 0, 2)):
    """
    Extract 3D structural features from PDB file.
    HETATM records are automatically filtered out.
    """
    print(f"\n[Step 1] Reading PDB structure...")
    print(f"  - Input file: {pdb_path}")
    
    # Read coordinates (HETATM filtering happens in get_atom_coordinates)
    coords, valid_mask = extract_pdb_features.get_coords_from_pdb(
        str(pdb_path), 
        full_backbone=True
    )
    print(f"  - Total residues: {len(coords)}")
    print(f"  - Valid residues (non-HETATM): {sum(valid_mask)}")
    
    print(f"\n[Step 2] Adjusting virtual C-beta positions...")
    coords = extract_pdb_features.move_CB(coords, virt_cb=virt_cb)
    
    print(f"\n[Step 3] Finding nearest neighbor residues...")
    partner_idx = extract_pdb_features.find_nearest_residues(coords, valid_mask)
    
    print(f"\n[Step 4] Calculating angular features...")
    features, valid_mask2 = extract_pdb_features.calc_angles_forloop(
        coords, partner_idx, valid_mask
    )
    print(f"  - Residues with complete features: {sum(valid_mask2)}")
    
    # By default historically we used first 8 feature columns (angular descriptors).
    # If caller provided expected_dim it will be handled by caller; here we keep
    # the full features and let caller slice/pad as needed.
    vae_features = features
    
    return vae_features, valid_mask2


def predict_with_encoder(encoder, features):
    """
    Encode features using the trained encoder model.
    """
    print(f"\n[Step 5] Encoding features with neural network...")
    encoder.eval()
    with torch.no_grad():
        z = encoder(torch.tensor(features, dtype=torch.float32))
        encoded = z.detach().numpy()
    print(f"  - Encoded feature shape: {encoded.shape}")
    return encoded


def discretize_to_states(encoded_features, centroids):
    """
    Discretize encoded features to 3Di states by finding nearest centroids.
    """
    print(f"\n[Step 6] Discretizing to 3Di states...")
    
    # Calculate distances to all centroids
    distances = np.sqrt(
        np.sum((encoded_features[:, np.newaxis, :] - centroids[np.newaxis, :, :])**2, axis=-1)
    )
    
    # Find nearest centroid for each residue
    states = np.argmin(distances, axis=1)
    print(f"  - Number of unique states used: {len(np.unique(states))}")
    
    return states


def states_to_sequence(states, invalid_state='X'):
    """
    Convert state indices to 3Di sequence string.
    """
    sequence = ''.join([
        LETTERS[state] if state != -1 else invalid_state
        for state in states
    ])
    return sequence


def pdb_to_3di(pdb_path, model_dir, output_path=None, no_header: bool=False):
    """
    Main function to convert PDB to 3Di sequence.
    """
    # Validate inputs
    print("="*60)
    print("PDB to 3Di Converter")
    print("="*60)
    
    pdb_path = validate_pdb_file(pdb_path)
    encoder_path, states_path = validate_model_files(model_dir)
    
    # Load model files
    print(f"\n[Loading Models]")
    print(f"  - Loading encoder...")
    encoder = torch.load(encoder_path, map_location=torch.device('cpu'))
    print(f"  - Loading centroids...")
    centroids = np.loadtxt(states_path)
    print(f"  - Centroids shape: {centroids.shape}")
    # Infer encoder expected input dimension so we can adapt extracted features.
    def infer_encoder_input_dim(enc) -> int:
        """Try to infer encoder input feature dimension (number of input channels).

        Strategy:
         - If enc is an nn.Module, search for first nn.Linear and use its in_features.
         - Otherwise inspect named parameters/state_dict and take the second dimension of
           the first 2D parameter (assumed weight shape (out,in)).
        """
        import torch.nn as nn
        # nn.Module: search modules
        try:
            if isinstance(enc, nn.Module):
                for m in enc.modules():
                    if isinstance(m, nn.Linear):
                        return int(m.in_features)
        except Exception:
            pass

        # Fallback: inspect parameters
        try:
            # state-dict like mapping
            if isinstance(enc, dict):
                items = enc.items()
            else:
                items = getattr(enc, 'state_dict', lambda: {})().items()
            for k, v in items:
                arr = np.asarray(v)
                if arr.ndim == 2:
                    # weight shape (out, in)
                    return int(arr.shape[1])
        except Exception:
            pass

        # Last resort: try named_parameters if module
        try:
            for name, p in getattr(enc, 'named_parameters', lambda: [])():
                if p.ndim == 2:
                    return int(p.shape[1])
        except Exception:
            pass

        raise RuntimeError('Could not infer encoder input dimension')

    try:
        input_dim = infer_encoder_input_dim(encoder)
        print(f"  - Inferred encoder input dimension: {input_dim}")
    except Exception as e:
        print('[WARN] Could not infer encoder input dimension:', e)
        # default to 8 (historical)
        input_dim = 8
        print('  - Falling back to input_dim = 8')

    # Extract features from PDB (full per-residue features)
    features, valid_mask = extract_features(pdb_path)

    # Ensure features have the correct number of columns expected by encoder
    if features.shape[1] >= input_dim:
        used_features = features[:, :input_dim]
    else:
        # pad with zeros
        pad_width = input_dim - features.shape[1]
        pad = np.zeros((features.shape[0], pad_width), dtype=features.dtype)
        used_features = np.hstack([features, pad])

    # Encode only valid residues
    valid_features = used_features[valid_mask]
    print(f"\n  - Features to encode: {valid_features.shape}")

    encoded = predict_with_encoder(encoder, valid_features)
    valid_states = discretize_to_states(encoded, centroids)
    
    # Map back to full sequence
    print(f"\n[Step 7] Creating final 3Di sequence...")
    states = np.full(len(valid_mask), -1)
    states[valid_mask] = valid_states
    
    sequence_3di = states_to_sequence(states)
    
    # Determine output path
    # protein name (used for header) — always derive from pdb filename
    protein_name = pdb_path.stem  # filename without extension
    if output_path is None:
        output_path = pdb_path.parent / f"{protein_name}_3di.txt"
    else:
        output_path = Path(output_path)
    
    # Prepare output sequence: filter out HETATM / non-standard residues
    # so that the output 3Di sequence length matches the canonical
    # amino-acid sequence in the FASTA (common expectation).
    print(f"\n[Step 8] Preparing output sequence (filtering HETATM by default)...")
    try:
        from Bio.PDB import PDBParser
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure(protein_name, str(pdb_path))
        model = structure[0]
        chain = list(model.get_chains())[0]
        residues = list(chain.get_residues())
        non_het_mask = [len(r.id[0].strip()) == 0 for r in residues]
        filtered_sequence = ''.join([
            sequence_3di[i]
            for i, keep in enumerate(non_het_mask)
            if keep
        ])
    except Exception:
        # If PDB parsing or mask creation fails, fall back to full sequence
        filtered_sequence = sequence_3di

    # Save result
    print(f"\n[Step 9] Saving results...")
    with open(output_path, 'w') as f:
        if not no_header:
            f.write(f">{protein_name}\n")
        f.write(filtered_sequence + "\n")

    print(f"  - Output saved to: {output_path}")

    # Summary
    print(f"\n{'='*60}")
    print(f"✓ Conversion completed successfully!")
    print(f"{'='*60}")
    print(f"  Protein name: {protein_name}")
    print(f"  Total residues (chain entries): {len(sequence_3di)}")
    print(f"  Valid 3Di states: {sum(states != -1)}")
    print(f"  Invalid positions (X): {sum(states == -1)}")
    print(f"  Output 3Di length (filtered): {len(filtered_sequence)}")
    print(f"  Output file: {output_path}")
    print(f"{'='*60}\n")

    return filtered_sequence, output_path


def main():
    parser = argparse.ArgumentParser(
        description='Convert PDB file to 3Di sequence',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python pdb_to_3di.py protein.pdb training_3di_gpu_new/tmp
  python pdb_to_3di.py protein.pdb training_3di_gpu_original/tmp -o output.txt
        """
    )
    
    parser.add_argument('pdb_file', type=str,
                       help='Path to input PDB file')
    parser.add_argument('model_dir', type=str,
                       help='Directory containing encoder.pt and states.txt')
    parser.add_argument('-o', '--output', type=str, default=None,
                       help='Output file path (default: <protein_name>_3di.txt)')
    parser.add_argument('--no-header', action='store_true', help='Do not write a FASTA header (omit >protein_name line)')

    args = parser.parse_args()

    try:
        pdb_to_3di(args.pdb_file, args.model_dir, args.output, no_header=args.no_header)
        return 0
    except Exception as e:
        print(f"\n{'='*60}")
        print(f"❌ Error occurred:")
        print(f"{'='*60}")
        print(f"{e}")
        print(f"{'='*60}\n")
        return 1


if __name__ == '__main__':
    sys.exit(main())
