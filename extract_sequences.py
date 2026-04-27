#!/usr/bin/env python3
"""
Part 5: Extract sequences from PDB files with caching support

This module extracts amino acid sequences from PDB structure files.
It handles two types of PDB files:
    1. Single PDB files (.pdb) - standard PDB format
    2. PDB bundles (.tar.gz) - compressed archives containing multiple PDB files

Results are cached to avoid re-parsing PDB files in subsequent runs.
"""

from pathlib import Path
import tarfile
from collections import defaultdict
import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB import is_aa
from Bio.SeqUtils import seq1
from tqdm import tqdm
from utils import timer


def parse_chain_id_mapping(mapping_file):
    """
    Parse chain ID mapping file from a PDB bundle.
    
    Bundle PDBs contain multiple structure files with renamed chains.
    This mapping file tells us which original chain corresponds to each new chain.
    
    Format of mapping file:
        bundle_name.pdb:
            A B    # New chain A maps to original chain B
            C D    # New chain C maps to original chain D
    
    Args:
        mapping_file: Path to chain-id-mapping.txt file
    
    Returns:
        Dictionary: {pdb_filename: {new_chain: original_chain}}
    
    Example:
        {
            '182l-pdb-bundle1.pdb': {'A': 'B', 'C': 'D'},
            '182l-pdb-bundle2.pdb': {'A': 'E', 'F': 'G'}
        }
    """
    mapping = defaultdict(dict)
    current_bundle = None
    
    with open(mapping_file) as f:
        for line in f:
            line = line.rstrip()
            
            # Line ending with ".pdb:" indicates a new bundle file
            # Example: "182l-pdb-bundle1.pdb:"
            if line.endswith(".pdb:"):
                current_bundle = line[:-1]  # Remove trailing colon
                continue
            
            # Skip empty lines and header lines
            if not line.strip() or line.strip().startswith("New chain"):
                continue
            
            # Parse chain mapping: "new_chain original_chain"
            parts = line.split()
            if len(parts) == 2 and current_bundle:
                new_chain, original_chain = parts
                mapping[current_bundle][new_chain] = original_chain
    
    return mapping


@timer
def extract_sequences_from_single_pdbs(pdb_folder, assembly_ids, cache_file=None):
    """
    Extract sequences from single PDB files.
    
    Args:
        pdb_folder: Path to directory containing PDB files
        assembly_ids: List of PDB IDs to extract sequences for
        cache_file: Optional path to cache file (CSV)
    
    Returns:
        DataFrame with columns: ['Assembly ID', 'Chain', 'Sequence']
    
    Logic:
        1. Check if cached version exists (fast path)
        2. If not, parse each PDB file using Bio.PDB
        3. Extract one-letter amino acid sequences for each chain
        4. Save to cache for future runs
    """
    # Fast path: load from cache if available
    if cache_file and Path(cache_file).exists():
        print(f"  Loading cached single PDB sequences from: {cache_file}")
        df = pd.read_csv(cache_file)
        print(f"  Loaded {len(df)} sequences from cache")
        return df
    
    print("Extracting sequences from single PDB files...")
    parser = PDBParser(QUIET=True)  # QUIET suppresses unnecessary warnings
    
    pdb_sequences = []
    
    # Iterate through each PDB ID
    for pdb_id in tqdm(assembly_ids, desc="  Processing single PDBs", unit="files"):
        pdb_file = pdb_folder / f"{pdb_id}.pdb"
        
        if pdb_file.exists():
            try:
                # Parse the PDB file
                structure = parser.get_structure(pdb_id, pdb_file)
                
                # Iterate through models (usually just one)
                for model in structure:
                    # Iterate through chains
                    for chain in model.get_chains():
                        chain_id = chain.id
                        
                        # Extract amino acid residues only (skip water, ligands, etc.)
                        residues = [res.get_resname() for res in chain if is_aa(res)]
                        
                        if not residues:
                            continue
                        
                        # Convert three-letter codes to one-letter sequence
                        # Example: ['ALA', 'GLY', 'LYS'] → 'AGK'
                        oneletterseq = seq1("".join(residues))
                        
                        pdb_sequences.append({
                            "Assembly ID": pdb_id,
                            "Chain": chain_id,
                            "Sequence": oneletterseq
                        })
                        
            except Exception as e:
                tqdm.write(f"  Warning: Could not parse {pdb_id}: {e}")
    
    df = pd.DataFrame(pdb_sequences)
    print(f"  Extracted {len(df)} sequences from single PDB files")
    
    # Save to cache for faster future runs
    if cache_file and len(df) > 0:
        df.to_csv(cache_file, index=False)
        print(f"  Saved to cache: {cache_file}")
    
    return df


@timer
def extract_sequences_from_bundles(pdb_folder, assembly_ids, cache_file=None):
    """
    Extract sequences from PDB bundles (tar.gz archives).
    
    Bundles contain multiple PDB files (e.g., different conformations or assemblies).
    Each bundle has a mapping file to track chain IDs across files.
    
    Args:
        pdb_folder: Path to directory containing tar.gz bundle files
        assembly_ids: List of bundle IDs to extract sequences for
        cache_file: Optional path to cache file (CSV)
    
    Returns:
        DataFrame with columns: ['Assembly ID', 'Original Chain ID', 'Sequence']
    
    Logic:
        1. Check cache first
        2. For each bundle, extract tar.gz to temporary folder
        3. Parse chain mapping file to understand chain renaming
        4. Extract sequences from each PDB in the bundle
        5. Group by (Assembly ID, Original Chain ID) and concatenate sequences
    """
    # Fast path: load from cache if available
    if cache_file and Path(cache_file).exists():
        print(f"  Loading cached bundle sequences from: {cache_file}")
        df = pd.read_csv(cache_file)
        print(f"  Loaded {len(df)} sequences from cache")
        return df
    
    print("Extracting sequences from bundled PDB files...")
    parser = PDBParser(QUIET=True)
    
    all_sequences = []  # Collect all sequences before grouping
    
    for entry_id in tqdm(assembly_ids, desc="  Processing bundles", unit="bundles"):
        tar_file = pdb_folder / f"{entry_id}-pdb-bundle.tar.gz"
        
        if not tar_file.exists():
            tqdm.write(f"  Warning: Tar file not found: {tar_file}")
            continue
        
        # Create temporary extraction directory
        extract_folder = pdb_folder / f"{entry_id}-pdb-bundle"
        extract_folder.mkdir(exist_ok=True)
        
        try:
            # Extract the tar.gz archive
            with tarfile.open(tar_file) as tar:
                tar.extractall(path=extract_folder)
            
            # Find the chain mapping file
            mapping_file = list(extract_folder.rglob(f"{entry_id}-chain-id-mapping.txt"))
            if not mapping_file:
                tqdm.write(f"  Warning: Mapping file not found for {entry_id}")
                continue
            
            # Parse the mapping file
            chain_mapping = parse_chain_id_mapping(mapping_file[0])
            
            # Process each PDB file in the extracted bundle
            for pdb_file in extract_folder.rglob("*-pdb-bundle*.pdb"):
                pdb_name = pdb_file.name
                
                try:
                    structure = parser.get_structure(pdb_name, pdb_file)
                    
                    for model in structure:
                        for chain in model.get_chains():
                            new_chain_id = chain.id
                            
                            # Map back to original chain ID using the mapping file
                            original_chain_id = chain_mapping.get(pdb_name, {}).get(new_chain_id, new_chain_id)
                            
                            # Extract amino acid residues
                            residues = [res.get_resname() for res in chain if is_aa(res)]
                            if not residues:
                                continue
                            
                            # Convert to one-letter sequence
                            oneletterseq = seq1("".join(residues))
                            
                            all_sequences.append({
                                "Assembly ID": entry_id,
                                "New Chain ID": new_chain_id,
                                "Original Chain ID": original_chain_id,
                                "Sequence": oneletterseq,
                            })
                            
                except Exception as e:
                    tqdm.write(f"  Warning: Error parsing {pdb_name}: {e}")
                    
        except Exception as e:
            tqdm.write(f"  Warning: Error processing {entry_id}: {e}")
    
    # Convert to DataFrame
    df_sequences = pd.DataFrame(all_sequences)
    
    if len(df_sequences) > 0:
        # Group by Assembly ID and Original Chain ID
        # Multiple PDB files in the same bundle may contribute to the same chain
        # Concatenate sequences in order (across files)
        diff_seq = df_sequences.groupby(
            ["Assembly ID", "Original Chain ID"], as_index=False
        ).agg({"Sequence": lambda x: "".join(x)})
        
        print(f"  Extracted {len(diff_seq)} sequences from {len(assembly_ids)} bundles")
        
        # Save to cache
        if cache_file and len(diff_seq) > 0:
            diff_seq.to_csv(cache_file, index=False)
            print(f"  Saved to cache: {cache_file}")
        
        return diff_seq
    else:
        print("  No sequences extracted from bundles")
        return pd.DataFrame(columns=["Assembly ID", "Original Chain ID", "Sequence"])


@timer
def extract_all_sequences(pdb_folder, qbio_df, cache_dir=None):
    """
    Main function to extract all sequences with caching.
    
    Determines which PDB IDs are single files vs bundles by checking
    the file system. Extracts sequences from both types and returns
    them as separate DataFrames.
    
    Args:
        pdb_folder: Path to folder containing PDB files
        qbio_df: Q-BioLiP DataFrame with 'Assembly ID' column
        cache_dir: Directory to store cache files (if None, use pdb_folder)
    
    Returns:
        Tuple of (easy_seq_df, diff_seq_df)
        - easy_seq_df: Sequences from single PDB files
        - diff_seq_df: Sequences from bundles (already grouped by chain)
    
    Logic:
        1. Setup cache directory
        2. Check if both caches exist (fast path)
        3. Get all available PDB files in the folder
        4. Separate assembly IDs into single files vs bundles
        5. Extract sequences from each type
    """
    print(f"Extracting sequences from PDB files in: {pdb_folder}")
    
    # Setup cache directory
    if cache_dir is None:
        cache_dir = pdb_folder
    else:
        cache_dir = Path(cache_dir)
        cache_dir.mkdir(exist_ok=True)
    
    # Cache file paths
    single_cache = cache_dir / "single_pdb_sequences.csv"
    bundle_cache = cache_dir / "bundle_pdb_sequences.csv"
    
    # FAST PATH: If both caches exist, load them directly
    if single_cache.exists() and bundle_cache.exists():
        print("  Found cached sequence data!")
        print(f"  Loading single PDB sequences from: {single_cache}")
        easy_seq = pd.read_csv(single_cache)
        print(f"  Loading bundle sequences from: {bundle_cache}")
        diff_seq = pd.read_csv(bundle_cache)
        print(f"  Loaded {len(easy_seq)} single and {len(diff_seq)} bundle sequences from cache")
        return easy_seq, diff_seq
    
    # SLOW PATH: Need to extract sequences from PDB files
    
    # Step 1: Get all single PDB files in the folder
    all_pdbs_assembly = []
    for pdb_path in pdb_folder.glob("*.pdb"):
        all_pdbs_assembly.append(pdb_path.stem)
    all_pdbs_assembly_set = set(all_pdbs_assembly)
    
    # Step 2: Get unique assembly IDs from Q-BioLiP data
    assembly_ids = set(qbio_df["Assembly ID"].astype(str))
    
    # Step 3: Separate IDs based on file type
    single_file_ids = []   # IDs that exist as .pdb files
    bundle_ids = []        # IDs that don't (assumed to be bundles)
    
    for id_val in assembly_ids:
        if id_val in all_pdbs_assembly_set:
            single_file_ids.append(id_val)
        else:
            bundle_ids.append(id_val)
    
    print(f"  Found {len(single_file_ids):,} single PDB files and {len(bundle_ids):,} bundles")
    
    # Step 4: Extract sequences from single PDB files
    easy_seq = extract_sequences_from_single_pdbs(pdb_folder, single_file_ids, single_cache)
    
    # Step 5: Extract sequences from bundles
    diff_seq = extract_sequences_from_bundles(pdb_folder, bundle_ids, bundle_cache)
    
    return easy_seq, diff_seq