#!/usr/bin/env python3
"""
Flatten pdbstyle-2.08 directory and copy/rename files to scope_pdb_backup.

pdbstyle-2.08 structure:
  pdbstyle-2.08/
    XX/  (2-letter prefix folders)
      YY/  (additional folders)
        domain_id.ent

Output:
  scope_pdb_backup/
    domain_id.pdb (renamed from .ent)
"""

import os
import shutil
from pathlib import Path
import logging

logging.basicConfig(level=logging.INFO, format='%(message)s')

def flatten_pdbstyle(source_dir, target_dir):
    """Flatten pdbstyle-2.08 and copy files to scope_pdb_backup."""
    
    source_path = Path(source_dir)
    target_path = Path(target_dir)
    
    target_path.mkdir(exist_ok=True)
    
    logging.info("="*80)
    logging.info(f"Flattening {source_dir} to {target_dir}")
    logging.info("="*80)
    
    # Find all .ent files recursively
    ent_files = list(source_path.rglob("*.ent"))
    logging.info(f"Found {len(ent_files)} .ent files")
    
    copied = 0
    skipped = 0
    
    for ent_file in ent_files:
        # Get domain ID from filename (remove .ent extension)
        domain_id = ent_file.stem
        
        # Create target path with .pdb extension
        target_file = target_path / f"{domain_id}.pdb"
        
        try:
            # Skip if already exists
            if target_file.exists():
                skipped += 1
                continue
            
            # Copy and rename
            shutil.copy2(ent_file, target_file)
            copied += 1
            
            if copied % 1000 == 0:
                logging.info(f"  Copied {copied} files...")
        
        except Exception as e:
            logging.error(f"Error copying {ent_file}: {e}")
    
    logging.info("="*80)
    logging.info(f"✓ Flattening completed!")
    logging.info(f"  Copied: {copied}")
    logging.info(f"  Skipped (already existed): {skipped}")
    logging.info(f"  Total in target dir: {len(list(target_path.glob('*.pdb')))}")
    logging.info("="*80)

if __name__ == "__main__":
    source = "/mnt/scratch/jugipalace/foldseek_new_3di/scope40/pdbstyle-2.08"
    target = "/mnt/scratch/jugipalace/foldseek_new_3di/scope_pdb_backup"
    
    flatten_pdbstyle(source, target)
