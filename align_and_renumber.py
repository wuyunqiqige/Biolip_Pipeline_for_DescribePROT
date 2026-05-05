#!/usr/bin/env python3
"""
Align sequences to describePROT and renumber binding sites

This module performs the critical task of mapping binding site positions from
original PDB numbering to DescribePROT reference sequence coordinates.

Key concepts:
    - Each binding site is a string like "L84" (amino acid L at position 84)
    - Positions are relative to the PDB sequence
    - We need to map them to DescribePROT coordinates via BLAST alignment
    - Uses pairwise BLAST (query vs single target) for maximum speed

Architecture:
    1. Group by (UniProt_ID, Sequence) to avoid redundant BLAST runs
    2. For each unique pair, run BLAST to align to DescribePROT
    3. Create position mapping dictionary from alignment
    4. Map each binding site using the position map
    5. Validate amino acid matches (only keep correct matches)
    6. Cache results for faster subsequent runs
"""

import pandas as pd
import numpy as np
from pathlib import Path
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import tempfile
import re
from tqdm import tqdm
import pickle
from typing import Dict, List, Optional, Tuple, Any


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def filter_failed_alignments(df: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    """
    Remove records that completely failed BLAST alignment.
    
    alignment_status values:
        - 'no_hit': BLAST found no alignment
        - 'no_sequence': No sequence available for alignment
        - Other values: Successful alignment (keep)
    
    Args:
        df: DataFrame with 'alignment_status' column
        dataset_name: Name for logging (QBioLiP or BioLiP2)
    
    Returns:
        Filtered DataFrame with failed alignments removed
    """
    if df is None or len(df) == 0:
        return df
    
    original_count = len(df)
    # Keep only rows that are NOT 'no_hit' or 'no_sequence'
    df_filtered = df[~df['alignment_status'].isin(['no_hit', 'no_sequence'])].copy()
    
    removed = original_count - len(df_filtered)
    print(f"\n  {dataset_name}: Removed {removed:,} records with NO BLAST hit, kept {len(df_filtered):,}")
    
    if len(df_filtered) > 0 and 'alignment_status' in df_filtered.columns:
        status_counts = df_filtered['alignment_status'].value_counts()
        status_counts_dict = {k: int(v) for k, v in status_counts.items()}
        print(f"    Kept breakdown: {status_counts_dict}")
    
    return df_filtered


def count_binding_site_mismatches(df, binding_col, seq_col='describePROT_sequence', label=""):
    """
    Validate binding sites against reference sequence.
    
    Checks if the amino acid at each binding site position matches
    the expected amino acid in the reference sequence.
    
    Args:
        df: DataFrame containing binding sites and sequences
        binding_col: Column name containing binding site strings (e.g., "L84 V87")
        seq_col: Column name containing reference sequence
        label: Description for printing (e.g., "Original validation")
    
    Returns:
        Dictionary with counts or None if column missing
    
    This is diagnostic only - no filtering is applied.
    """
    if binding_col not in df.columns:
        print(f"    Warning: {binding_col} not found, skipping mismatch count")
        return None
    
    has_binding = df[binding_col].notna() & (df[binding_col] != '')
    has_sequence = df[seq_col].notna() & (df[seq_col] != '')
    valid_rows = has_binding & has_sequence
    
    total_sites = 0
    match_count = 0
    mismatch_count = 0
    out_of_range_count = 0
    
    for idx, row in df[valid_rows].iterrows():
        binding_sites = str(row[binding_col])
        sequence = str(row[seq_col])
        
        # Parse binding sites: "L84 V87" -> [('L',84), ('V',87)]
        sites = []
        for part in binding_sites.strip().split():
            if len(part) >= 2:
                aa = part[0]                    # First char = amino acid
                try:
                    pos = int(part[1:])         # Rest = position number
                    sites.append((aa, pos))
                except ValueError:
                    continue
        
        for aa, pos in sites:
            total_sites += 1
            idx_pos = pos - 1                   # Convert to 0-based index
            
            if 0 <= idx_pos < len(sequence):
                actual_aa = sequence[idx_pos].upper()
                if actual_aa == aa.upper():
                    match_count += 1
                else:
                    mismatch_count += 1
            else:
                out_of_range_count += 1
    
    if total_sites > 0:
        print(f"    {label}: {match_count:,}/{total_sites:,} matches ({match_count/total_sites*100:.2f}%)")
        print(f"      Mismatches: {mismatch_count:,} ({mismatch_count/total_sites*100:.2f}%)")
        print(f"      Out of range: {out_of_range_count:,} ({out_of_range_count/total_sites*100:.2f}%)")
    else:
        print(f"    {label}: No binding sites found")
    
    return {
        'total_sites': total_sites,
        'match_count': match_count,
        'mismatch_count': mismatch_count,
        'out_of_range_count': out_of_range_count
    }


# =============================================================================
# BLAST ALIGNMENT FUNCTIONS
# =============================================================================

def run_pairwise_blast(query_seq: str, query_id: str, target_seq: str, target_id: str) -> List[Dict]:
    """
    Run BLAST between a query sequence and a single target sequence.
    
    Args:
        query_seq: Query amino acid sequence (from PDB)
        query_id: Identifier for query (UniProt ID)
        target_seq: Target sequence (from DescribePROT)
        target_id: Identifier for target (UniProt ID)
    
    Returns:
        List of hit dictionaries containing alignment information
    
    BLAST parameters:
        - evalue: 0.001 (stringent threshold)
        - max_hsps: 1 (only best hit)
        - num_threads: 8 (parallel processing)
    
    File handling:
        - Creates temporary FASTA files for query and target
        - Runs BLAST, outputs XML
        - Cleans up temporary files after completion
    """
    # Create temporary FASTA files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">{query_id}\n{query_seq}\n")
        query_file = f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as f:
        f.write(f">{target_id}\n{target_seq}\n")
        target_file = f.name
    
    output_file = tempfile.NamedTemporaryFile(suffix='.xml', delete=False).name
    
    try:
        # Run BLAST
        blast_cline = NcbiblastpCommandline(
            query=query_file,
            subject=target_file,
            outfmt=5,           # XML output format
            out=output_file,
            evalue=0.001,       # E-value threshold
            max_hsps=1,         # Only one HSP per alignment
            num_threads=8       # Use 8 CPU cores
        )
        stdout, stderr = blast_cline()
        
        # Parse XML output
        hits = []
        with open(output_file) as result_handle:
            blast_records = NCBIXML.parse(result_handle)
            for record in blast_records:
                for alignment in record.alignments:
                    for hsp in alignment.hsps:
                        hits.append({
                            'hit_def': alignment.hit_def,
                            'hit_id': target_id,
                            'evalue': hsp.expect,
                            'score': hsp.score,
                            'identity': (hsp.identities / hsp.align_length) * 100,
                            'align_length': hsp.align_length,
                            'query_start': hsp.query_start,
                            'query_end': hsp.query_end,
                            'sbjct_start': hsp.sbjct_start,
                            'sbjct_end': hsp.sbjct_end,
                            'query_seq': hsp.query,
                            'sbjct_seq': hsp.sbjct,
                            'match_seq': hsp.match
                        })
        return hits
    
    finally:
        # Clean up temporary files
        Path(query_file).unlink(missing_ok=True)
        Path(target_file).unlink(missing_ok=True)
        Path(output_file).unlink(missing_ok=True)


def parse_blast_alignment(hit: Dict, query_seq: str, subject_seq: str) -> Optional[Dict]:
    """
    Parse BLAST alignment to create position mapping.
    
    The mapping tells us: for each position in the query sequence,
    what position in the subject sequence it aligns to.
    
    Example:
        Query:   M A V T G
        Subject: M A - T G
        Mapping: 1→1, 2→2, 3→(gap), 4→3, 5→4
    
    Args:
        hit: BLAST hit dictionary from run_pairwise_blast()
        query_seq: Original query sequence
        subject_seq: Original subject sequence
    
    Returns:
        Dictionary with alignment metadata and position_map
    """
    query_start = hit['query_start']
    query_end = hit['query_end']
    sbjct_start = hit['sbjct_start']
    sbjct_end = hit['sbjct_end']
    identity = hit['identity']
    
    position_map = {}
    query_seq_aligned = hit['query_seq']
    sbjct_seq_aligned = hit['sbjct_seq']
    
    query_pos = query_start
    sbjct_pos = sbjct_start
    
    # Walk through aligned sequences, tracking positions
    for q_char, s_char in zip(query_seq_aligned, sbjct_seq_aligned):
        if q_char != '-':               # Query has a residue
            if s_char != '-':           # Subject also has a residue (match/mismatch)
                position_map[query_pos] = sbjct_pos
                query_pos += 1
                sbjct_pos += 1
            else:                       # Query residue, subject gap
                query_pos += 1
        else:                           # Query gap
            if s_char != '-':           # Subject residue, query gap
                sbjct_pos += 1
    
    return {
        'subject_id': hit['hit_id'],
        'identity': identity,
        'evalue': hit['evalue'],
        'score': hit['score'],
        'alignment_length': hit['align_length'],
        'query_start': query_start,
        'query_end': query_end,
        'subject_start': sbjct_start,
        'subject_end': sbjct_end,
        'position_map': position_map
    }

def find_best_alignment_for_sites(binding_sites_str, hits, query_seq, target_dp_seq, subject_sequence):
    """
    Try multiple BLAST hits to find one that maps the most binding sites
    """
    
    original_sites = parse_binding_sites(binding_sites_str)
    if not original_sites:
        return None, 0
    
    best_alignment = None
    best_mapped_count = 0
    best_mapped_sites = None
    
    for hit in hits[:5]:  # Try top 5 hits
        alignment = parse_blast_alignment(hit, query_seq, target_dp_seq)
        if not alignment:
            continue
        
        # Count how many binding sites map with this alignment
        mapped_count = 0
        mapped_sites = []
        
        for aa, position in original_sites:
            if position in alignment['position_map']:
                mapped_position = alignment['position_map'][position]
                idx = mapped_position - 1
                if 0 <= idx < len(subject_sequence):
                    actual_aa = subject_sequence[idx].upper()
                    if actual_aa == aa.upper():
                        mapped_count += 1
                        mapped_sites.append((aa, mapped_position))
        
        if mapped_count > best_mapped_count:
            best_mapped_count = mapped_count
            best_alignment = alignment
            best_mapped_sites = mapped_sites
    
    return best_alignment, best_mapped_sites

def create_perfect_match_result(uniprot_id: str, target_seq: str, query_seq: str) -> Dict:
    """
    Create result for when sequences are identical (no BLAST needed).
    
    This is an optimization: if query and target are exactly the same,
    we can create a 1:1 position map without running BLAST.
    
    Args:
        uniprot_id: UniProt identifier
        target_seq: Target sequence
        query_seq: Query sequence (must equal target_seq)
    
    Returns:
        Alignment dictionary with identity=100% and 1:1 position map
    """
    # Map position i in query to position i in subject
    position_map = {i+1: i+1 for i in range(len(query_seq))}
    
    return {
        'subject_id': uniprot_id,
        'identity': 100.0,
        'evalue': 0.0,
        'score': 999999,
        'alignment_length': len(query_seq),
        'query_start': 1,
        'query_end': len(query_seq),
        'subject_start': 1,
        'subject_end': len(query_seq),
        'position_map': position_map
    }


# =============================================================================
# BINDING SITE MAPPING FUNCTIONS
# =============================================================================

def parse_binding_sites(binding_sites_str: str) -> List[Tuple[str, int]]:
    """
    Parse binding sites string into list of (amino_acid, position).
    
    Args:
        binding_sites_str: Space-separated binding sites like "L84 V87 Y88"
    
    Returns:
        List of tuples: [('L',84), ('V',87), ('Y',88)]
    """
    if pd.isna(binding_sites_str) or binding_sites_str == '':
        return []
    
    sites = []
    for part in str(binding_sites_str).strip().split():
        part = part.strip()
        if len(part) >= 2:
            aa = part[0]                # First character = amino acid
            try:
                pos = int(part[1:])     # Rest = position number
                sites.append((aa, pos))
            except ValueError:
                continue
    return sites


def format_binding_sites(sites: List[Tuple[str, int]]) -> str:
    """
    Format list of (amino_acid, position) back to space-separated string.
    
    Args:
        sites: List of tuples like [('L',84), ('V',87)]
    
    Returns:
        String: "L84 V87" or None if empty
    """
    if not sites:
        return None
    return ' '.join([f"{aa}{pos}" for aa, pos in sites])


def map_binding_sites_with_alignment(
    binding_sites_str: str, 
    alignment: Dict,
    subject_sequence: str = None,
    debug: bool = False,
    uniprot_id: str = None,
    query_seq: str = None
) -> Dict[str, Any]:
    """
    Map binding sites using the alignment position map.
    
    This is the core function that renumbers binding sites:
        1. Parse original binding sites
        2. For each site, find mapped position using alignment
        3. Validate that amino acid matches subject sequence
        4. Only keep sites that pass validation
    
    Args:
        binding_sites_str: Original binding sites (e.g., "L84 V87")
        alignment: Alignment dictionary with 'position_map'
        subject_sequence: DescribePROT sequence for validation
        debug: Enable debug printing
        uniprot_id: For debug messages
        query_seq: For debug messages
    
    Returns:
        Dictionary with mapping results
    """
    # Handle empty input
    if pd.isna(binding_sites_str) or binding_sites_str == '':
        return {
            'mapped_sites': [],
            'mapped_sites_str': None,
            'all_mapped': False,
            'mapped_count': 0,
            'total_sites': 0,
            'removed_mismatches': 0,
            'removed_out_of_range': 0
        }
    
    original_sites = parse_binding_sites(binding_sites_str)
    
    if not original_sites:
        return {
            'mapped_sites': [],
            'mapped_sites_str': None,
            'all_mapped': False,
            'mapped_count': 0,
            'total_sites': 0,
            'removed_mismatches': 0,
            'removed_out_of_range': 0
        }
    
    position_map = alignment['position_map']
    mapped_sites = []
    removed_mismatches = 0
    removed_out_of_range = 0
    
    for aa, position in original_sites:
        if position in position_map:
            mapped_position = position_map[position]
            
            # Validate amino acid match
            if subject_sequence:
                idx = mapped_position - 1
                if 0 <= idx < len(subject_sequence):
                    actual_aa = subject_sequence[idx].upper()
                    if actual_aa == aa.upper():
                        mapped_sites.append((aa, mapped_position))
                    else:
                        removed_mismatches += 1
                else:
                    removed_out_of_range += 1
            else:
                mapped_sites.append((aa, mapped_position))
        # position not in map - silently dropped (gaps in alignment)
    
    mapped_sites_str = format_binding_sites(mapped_sites)
    
    if debug and (removed_mismatches > 0 or removed_out_of_range > 0):
        print(f"\n  DEBUG MAPPING for {uniprot_id}:")
        print(f"    Original: {binding_sites_str}")
        print(f"    Mapped: {mapped_sites_str}")
        print(f"    Removed mismatches: {removed_mismatches}, out of range: {removed_out_of_range}")
    
    return {
        'mapped_sites': mapped_sites,
        'mapped_sites_str': mapped_sites_str,
        'all_mapped': len(mapped_sites) == len(original_sites),
        'mapped_count': len(mapped_sites),
        'total_sites': len(original_sites),
        'removed_mismatches': removed_mismatches,
        'removed_out_of_range': removed_out_of_range
    }


# =============================================================================
# MAIN ALIGNMENT FUNCTION
# =============================================================================

def align_and_renumber_dataset(
    df: pd.DataFrame, 
    describeprot_df: pd.DataFrame, 
    dataset_name: str, 
    use_cache: bool = True
) -> pd.DataFrame:
    """
    Align all sequences to describePROT using pairwise BLAST.
    
    This is the main orchestration function. Steps:
        1. Setup: Load describePROT sequences, check cache
        2. Identify columns to preserve (affinity, metadata)
        3. Group by (UniProt_ID, Sequence) to avoid redundant BLAST
        4. For each unique pair, run BLAST and create position map
        5. Apply mapping to all rows in the group
        6. Restore preserved columns
        7. Filter out failed alignments
    
    Args:
        df: Input DataFrame with sequences and binding sites
        describeprot_df: DescribePROT DataFrame with 'ACC' and 'seq'
        dataset_name: "QBioLiP" or "BioLiP2" (affects column detection)
        use_cache: Whether to use/save alignment cache
    
    Returns:
        DataFrame with renumbered binding sites and describePROT sequences
    """
    
    print(f"\n{'='*60}")
    print(f"  ALIGNMENT FOR {dataset_name}")
    print(f"{'='*60}")
    
    print(f"\n  Initializing alignment for {dataset_name}...")
    
    # =========================================================================
    # STEP 1: Setup describePROT sequences
    # =========================================================================
    subject_sequences = dict(zip(describeprot_df['ACC'], describeprot_df['seq']))
    print(f"  Loaded {len(subject_sequences):,} describePROT sequences")
    
    # =========================================================================
    # STEP 2: Load or create alignment cache
    # =========================================================================
    cache_file = Path(f"{dataset_name}_alignment_cache.pkl")
    
    alignment_cache = {}
    if use_cache and cache_file.exists():
        try:
            with open(cache_file, 'rb') as f:
                alignment_cache = pickle.load(f)
            print(f"  Loaded {len(alignment_cache)} cached alignments")
        except Exception as e:
            print(f"  Error loading cache: {e}, starting fresh")
            alignment_cache = {}
    
    # =========================================================================
    # STEP 3: Initialize new columns in DataFrame
    # =========================================================================
    df = df.copy()
    df['describePROT_ID'] = None
    df['describePROT_sequence'] = None
    df['alignment_identity'] = None
    df['alignment_evalue'] = None
    df['alignment_score'] = None
    df['alignment_coverage'] = None
    df['renumbered_binding_sites'] = None
    df['binding_sites_mapped_count'] = None
    df['binding_sites_total_count'] = None
    df['binding_sites_all_mapped'] = False
    df['alignment_status'] = 'pending'
    
    # =========================================================================
    # STEP 4: Identify sequence column
    # =========================================================================
    seq_col = None
    possible_seq_cols = ['Sequence', 'Receptor_sequence', 'seq']
    
    for col in possible_seq_cols:
        if col in df.columns:
            seq_col = col
            print(f"  Using sequence column: '{seq_col}'")
            break
    
    if seq_col is None:
        raise ValueError(f"No sequence column found. Available: {df.columns.tolist()}")
    
    # =========================================================================
    # STEP 5: Identify binding sites column (dataset-specific)
    # =========================================================================
    sites_col = None
    
    if dataset_name == "QBioLiP":
        # Q-BioLiP uses Binding_site_pdb from flattening
        possible_cols = ['Binding_site_pdb', 'Binding_site_original', 'Binding_sites', 'Binding Site', 'QBioLiP_Renumbered']
        for col in possible_cols:
            if col in df.columns:
                sites_col = col
                print(f"  Using Q-BioLiP binding sites column: '{sites_col}'")
                break
        
        if sites_col is None:
            print(f"  WARNING: No binding sites column found for QBioLiP")
            print(f"    Available columns with 'binding' or 'site': {[col for col in df.columns if 'binding' in col.lower() or 'site' in col.lower()][:10]}")
    
    elif dataset_name == "BioLiP2":
        # BioLiP2 has pre-renumbered binding sites
        possible_cols = ['Binding_site_renumbered', 'Binding_site_pdb', 'Binding_sites', 'Binding_site_original']
        for col in possible_cols:
            if col in df.columns:
                sites_col = col
                print(f"  Using BioLiP2 binding sites column: '{sites_col}'")
                break
        
        if sites_col is None:
            print(f"  WARNING: No binding sites column found for BioLiP2")
            print(f"    Available columns with 'binding' or 'site': {[col for col in df.columns if 'binding' in col.lower() or 'site' in col.lower()][:10]}")
    
    # =========================================================================
    # STEP 6: Identify columns to preserve (affinity, metadata)
    # =========================================================================
    preserve_cols = []
    
    # Patterns for identifying columns to preserve
    affinity_patterns = ['affinity', 'BindingMOAD', 'PDBbind', 'BindingDB', 'MOAD', 'PDBbind-CN']
    metadata_patterns = ['DrugBank', 'Ligand', 'PDB_ID', 'Assembly', 'Resolution', 'Stoichiometry']
    
    for col in df.columns:
        col_lower = col.lower()
        # Skip columns used in alignment
        if col in ['UniProt_ID', seq_col, sites_col]:
            continue
        # Skip columns we're creating
        if col.startswith('alignment_') or col in ['describePROT_ID', 'describePROT_sequence', 'renumbered_binding_sites',
                                                    'binding_sites_mapped_count', 'binding_sites_total_count', 
                                                    'binding_sites_all_mapped', 'alignment_status']:
            continue
        
        is_affinity = any(pattern in col for pattern in affinity_patterns)
        is_metadata = any(pattern.lower() in col_lower for pattern in metadata_patterns)
        
        if is_affinity or is_metadata:
            preserve_cols.append(col)
    
    if preserve_cols:
        print(f"  Preserving {len(preserve_cols)} additional columns: {preserve_cols[:5]}..." if len(preserve_cols) > 5 else f"  Preserving columns: {preserve_cols}")
        # Store preserved data as dictionary for each row
        # This ensures affinity values stay attached to their rows through grouping
        df['_preserved_data'] = df[preserve_cols].apply(lambda x: x.to_dict(), axis=1)
    
    # =========================================================================
    # STEP 7: Pre-alignment validation statistics
    # =========================================================================
    print(f"\n  {'-'*40}")
    print(f"  BEFORE ALIGNMENT (original binding sites against original sequences)")
    print(f"  {'-'*40}")
    
    if sites_col and sites_col in df.columns:
        non_null_count = df[sites_col].notna().sum()
        print(f"    Binding sites column '{sites_col}' has {non_null_count:,} non-null values")
        if non_null_count > 0:
            sample = df[sites_col].dropna().iloc[0]
            print(f"    Sample binding site: {sample}")
    
    # =========================================================================
    # STEP 8: Group by (UniProt_ID, Sequence) to avoid redundant BLAST
    # =========================================================================
    print(f"\n  Finding unique (UniProt_ID, Sequence) pairs...")
    grouped = df.groupby(['UniProt_ID', seq_col])
    total_pairs = len(grouped)
    print(f"  Found {total_pairs:,} unique pairs out of {len(df):,} total rows")
    
    # =========================================================================
    # STEP 9: Process each unique pair
    # =========================================================================
    pair_results = {}
    cached_count = 0
    direct_match_count = 0
    blast_count = 0
    no_target_count = 0
    
    total_original_sites = 0
    total_mapped_sites = 0
    total_removed_mismatches = 0
    total_removed_out_of_range = 0
    
    for idx, ((uniprot_id, query_seq), group_indices) in enumerate(tqdm(grouped, desc=f"  Aligning unique pairs")):
        cache_key = (uniprot_id, query_seq)
        
        if cache_key in alignment_cache:
            mapping_result = alignment_cache[cache_key]
            cached_count += 1
        else:
            target_dp_seq = subject_sequences.get(uniprot_id)
            
            if target_dp_seq is None:
                mapping_result = None
                no_target_count += 1
            elif query_seq == target_dp_seq:
                # Sequences identical - create perfect match without BLAST
                best_alignment = create_perfect_match_result(uniprot_id, target_dp_seq, query_seq)
                mapping_result = {
                    'best_overall_alignment': best_alignment,
                    'has_alignment': True
                }
                direct_match_count += 1
                alignment_cache[cache_key] = mapping_result
            else:
                # Run BLAST and use the top hit
                hits = run_pairwise_blast(query_seq, uniprot_id, target_dp_seq, uniprot_id)
                
                if hits:
                    # Use only the top hit (first one)
                    best_alignment = parse_blast_alignment(hits[0], query_seq, target_dp_seq)
                    if best_alignment:
                        mapping_result = {
                            'best_overall_alignment': best_alignment,
                            'has_alignment': True
                        }
                    else:
                        mapping_result = None
                    blast_count += 1
                else:
                    mapping_result = None
                    blast_count += 1
            
            alignment_cache[cache_key] = mapping_result
        
        pair_results[(uniprot_id, query_seq)] = mapping_result
        
        # Save checkpoint
        if use_cache and (idx + 1) % 1000 == 0:
            with open(cache_file, 'wb') as f:
                pickle.dump(alignment_cache, f)
            print(f"\n  Checkpoint saved: {idx+1}/{total_pairs} pairs")
    
    # Final cache save
    if use_cache:
        with open(cache_file, 'wb') as f:
            pickle.dump(alignment_cache, f)
        print(f"\n  Final cache saved: {len(alignment_cache)} alignments to {cache_file}")
    
    print(f"\n  Processing stats:")
    print(f"    From cache: {cached_count:,}")
    print(f"    Direct matches (no BLAST): {direct_match_count:,}")
    print(f"    Pairwise BLAST runs: {blast_count:,}")
    print(f"    No target sequence found: {no_target_count:,}")
    
    # =========================================================================
    # STEP 10: Apply results to each row
    # =========================================================================
    print(f"\n  Applying results to {len(df):,} rows...")
    
    success_count = 0
    no_hits_count = 0
    
    for idx, row in tqdm(df.iterrows(), total=len(df), desc=f"  Applying results"):
        uniprot_id = row['UniProt_ID']
        query_seq = row[seq_col]
        binding_sites = row[sites_col] if sites_col else None
        
        mapping_result = pair_results.get((uniprot_id, query_seq))
        
        if mapping_result and mapping_result.get('best_overall_alignment'):
            best_alignment = mapping_result['best_overall_alignment']
            
            # Add alignment metadata
            df.at[idx, 'describePROT_ID'] = best_alignment['subject_id']
            df.at[idx, 'describePROT_sequence'] = subject_sequences.get(best_alignment['subject_id'], '')
            df.at[idx, 'alignment_identity'] = best_alignment['identity']
            df.at[idx, 'alignment_evalue'] = best_alignment['evalue']
            df.at[idx, 'alignment_score'] = best_alignment['score']
            
            query_len = len(query_seq)
            aligned_len = best_alignment['query_end'] - best_alignment['query_start'] + 1
            coverage = (aligned_len / query_len) * 100 if query_len > 0 else 0
            df.at[idx, 'alignment_coverage'] = coverage
            
            # Map binding sites
            if binding_sites and not pd.isna(binding_sites):
                site_mapping = map_binding_sites_with_alignment(
                    binding_sites, 
                    best_alignment,
                    subject_sequence=df.at[idx, 'describePROT_sequence'],
                    debug=False,
                    uniprot_id=uniprot_id,
                    query_seq=query_seq
                )
                
                total_original_sites += site_mapping.get('total_sites', 0)
                total_mapped_sites += site_mapping.get('mapped_count', 0)
                total_removed_mismatches += site_mapping.get('removed_mismatches', 0)
                total_removed_out_of_range += site_mapping.get('removed_out_of_range', 0)
                
                df.at[idx, 'renumbered_binding_sites'] = site_mapping.get('mapped_sites_str')
                df.at[idx, 'binding_sites_mapped_count'] = site_mapping.get('mapped_count', 0)
                df.at[idx, 'binding_sites_total_count'] = site_mapping.get('total_sites', 0)
                df.at[idx, 'binding_sites_all_mapped'] = site_mapping.get('all_mapped', False)
                
                if site_mapping.get('all_mapped', False):
                    df.at[idx, 'alignment_status'] = 'all_sites_mapped'
                    success_count += 1
                elif site_mapping.get('mapped_count', 0) > 0:
                    df.at[idx, 'alignment_status'] = 'partial_mapping'
                    success_count += 1
                else:
                    df.at[idx, 'alignment_status'] = 'no_sites_mapped'
                    no_hits_count += 1
            else:
                df.at[idx, 'alignment_status'] = 'aligned_no_sites'
                success_count += 1
        else:
            df.at[idx, 'alignment_status'] = 'no_hit'
            no_hits_count += 1
    
    # =========================================================================
    # STEP 11: Post-alignment statistics
    # =========================================================================
    print(f"\n  {'-'*40}")
    print(f"  AFTER ALIGNMENT (mapped binding sites against describePROT)")
    print(f"  {'-'*40}")
    
    if 'renumbered_binding_sites' in df.columns:
        count_binding_site_mismatches(df, 'renumbered_binding_sites', 'describePROT_sequence', "After mapping (before filtering)")
    
    print(f"\n  Site removal statistics during mapping:")
    print(f"    Total original sites: {total_original_sites:,}")
    if total_original_sites > 0:
        print(f"    Sites kept after mapping: {total_mapped_sites:,} ({total_mapped_sites/total_original_sites*100:.2f}%)")
        print(f"    Removed due to amino acid mismatch: {total_removed_mismatches:,}")
        print(f"    Removed due to out of range: {total_removed_out_of_range:,}")
    else:
        print(f"    No sites to map")
    
    # =========================================================================
    # STEP 12: Filter out failed alignments
    # =========================================================================
    df_filtered = filter_failed_alignments(df, dataset_name)
    
    # =========================================================================
    # STEP 13: Restore preserved columns (affinity, metadata)
    # =========================================================================
    if preserve_cols and '_preserved_data' in df_filtered.columns:
        print(f"\n  Restoring {len(preserve_cols)} preserved columns...")
        for col in preserve_cols:
            df_filtered[col] = df_filtered['_preserved_data'].apply(lambda x: x.get(col) if isinstance(x, dict) else None)
        df_filtered = df_filtered.drop(columns=['_preserved_data'])
    
    # =========================================================================
    # STEP 14: Final validation
    # =========================================================================
    print(f"\n  {'-'*40}")
    print(f"  FINAL VALIDATION (after filtering rows with no hits)")
    print(f"  {'-'*40}")
    if 'Binding_sites' in df_filtered.columns:
        count_binding_site_mismatches(df_filtered, 'Binding_sites', 'describePROT_sequence', "Final Binding_sites column")
    elif 'renumbered_binding_sites' in df_filtered.columns:
        count_binding_site_mismatches(df_filtered, 'renumbered_binding_sites', 'describePROT_sequence', "Final renumbered_binding_sites")
    
    # =========================================================================
    # STEP 15: Print final statistics
    # =========================================================================
    print(f"\n  {'='*50}")
    print(f"  {dataset_name} Alignment Results")
    print(f"  {'='*50}")
    print(f"    Total rows processed: {len(df):,}")
    print(f"    Successfully aligned: {len(df_filtered):,}")
    
    if len(df_filtered) > 0:
        all_mapped = (df_filtered['alignment_status'] == 'all_sites_mapped').sum()
        partial_mapped = (df_filtered['alignment_status'] == 'partial_mapping').sum()
        no_sites = (df_filtered['alignment_status'] == 'aligned_no_sites').sum()
        
        print(f"      - All binding sites mapped: {all_mapped:,}")
        print(f"      - Partial binding sites mapped: {partial_mapped:,}")
        print(f"      - No binding sites to map: {no_sites:,}")
    
    # =========================================================================
    # STEP 16: Drop temporary alignment columns
    # =========================================================================
    columns_to_drop = [
        'describePROT_ID',
        'alignment_identity',
        'alignment_coverage',
        'alignment_score',
        'alignment_evalue',
        'alignment_status',
        'binding_sites_mapped_count',
        'binding_sites_total_count',
        'binding_sites_all_mapped',
    ]
    
    cols_to_drop = [col for col in columns_to_drop if col in df_filtered.columns]
    if cols_to_drop:
        df_filtered = df_filtered.drop(columns=cols_to_drop)
        print(f"    Dropped columns: {cols_to_drop}")
    
    return df_filtered

def validate_and_update_binding_sites(df: pd.DataFrame, dataset_name: str) -> pd.DataFrame:
    """Rename renumbered binding sites to standardized column name."""
    df = df.copy()
    
    renumbered_col = None
    if 'renumbered_binding_sites' in df.columns:
        renumbered_col = 'renumbered_binding_sites'
    elif 'QBioLiP_Renumbered' in df.columns:
        renumbered_col = 'QBioLiP_Renumbered'
    elif 'Binding_site_renumbered' in df.columns:
        renumbered_col = 'Binding_site_renumbered'
    
    if renumbered_col and renumbered_col != 'Binding_sites':
        # Rename instead of duplicate
        df = df.rename(columns={renumbered_col: 'Binding_sites'})
        print(f"\n  {dataset_name}: Renamed '{renumbered_col}' to 'Binding_sites'")
    
    return df

def create_alignment_summary(df: pd.DataFrame, output_path: Path) -> None:
    """
    Create a summary CSV of alignment statistics.
    
    Args:
        df: DataFrame with 'alignment_status' column
        output_path: Path to save CSV file
    """
    summary_data = []
    
    if 'alignment_status' in df.columns:
        status_counts = df['alignment_status'].value_counts()
        for status, count in status_counts.items():
            summary_data.append({
                'category': 'alignment_status',
                'value': status,
                'count': count,
                'percentage': (count / len(df)) * 100 if len(df) > 0 else 0
            })
    
    if summary_data:
        summary_df = pd.DataFrame(summary_data)
        summary_df.to_csv(output_path, index=False)
        print(f"\n  Alignment summary saved to: {output_path}")


if __name__ == "__main__":
    print("align_and_renumber.py loaded successfully")