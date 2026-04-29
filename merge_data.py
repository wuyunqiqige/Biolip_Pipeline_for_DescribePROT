#!/usr/bin/env python3
"""
Part 6: Merge sequences with Q-BioLiP data

This module takes the flattened Q-BioLiP binding site data and merges it with
the amino acid sequences extracted from PDB files. It handles both single PDB
files and bundles, combining results from both sources.
"""

import pandas as pd
import re
from utils import timer


def parse_binding_sites(binding_sites_str):
    """Parse space-separated binding sites into list of (amino_acid, position)"""
    if pd.isna(binding_sites_str) or binding_sites_str == '':
        return []
    
    sites = []
    for part in str(binding_sites_str).strip().split():
        part = part.strip()
        if len(part) >= 2:
            aa = part[0]
            try:
                pos = int(part[1:])
                sites.append((aa, pos))
            except ValueError:
                continue
    return sites


def format_binding_sites(sites):
    """Format list of (amino_acid, position) back to space-separated string"""
    if not sites:
        return None
    return ' '.join([f"{aa}{pos}" for aa, pos in sites])


def validate_single_binding_site(aa, position, sequence):
    """
    Validate a single binding site against the sequence.
    
    Returns:
        tuple: (True, actual_position) if the amino acid matches at the given position,
               (False, None) otherwise
    """
    if not sequence or position < 1 or position > len(sequence):
        return False, None
    
    actual_aa = sequence[position - 1]
    if actual_aa == aa:
        return True, position
    return False, None


def calculate_offsets(binding_sites, sequence):
    """
    Calculate possible offsets for each binding site based on the sequence.
    
    For each binding site (aa, position), the offset is (actual_position - expected_position)
    where actual_position is where the amino acid appears in the sequence.
    
    Returns a list of offsets that work for a continuation of sites.
    """
    if not binding_sites or not sequence:
        return []
    
    # For the first binding site, find all possible positions in the sequence
    first_aa, first_pos = binding_sites[0]
    possible_positions = []
    
    # Find all occurrences of the amino acid in the sequence
    for i, seq_aa in enumerate(sequence):
        if seq_aa == first_aa:
            # Calculate offset: actual position (i+1) - expected position
            offset = (i + 1) - first_pos
            possible_positions.append((i + 1, offset))
    
    return possible_positions


def find_longest_continuation(binding_sites, sequence):
    """
    Find the longest continuation of binding sites that share the same offset.
    
    Returns:
        tuple: (mapped_sites, offset, end_index)
        mapped_sites: list of (aa, actual_position) for mapped sites
        offset: the offset that worked
        end_index: index in binding_sites where continuation ends
    """
    if not binding_sites:
        return [], None, 0
    
    best_mapped_sites = []
    best_offset = None
    best_end_index = 0
    
    # Get possible offsets for the first binding site
    possible_offsets = calculate_offsets(binding_sites[:1], sequence)
    
    # Try each possible offset
    for actual_pos, offset in possible_offsets:
        mapped_sites = []
        valid = True
        end_idx = 0
        
        for i, (aa, expected_pos) in enumerate(binding_sites):
            # Calculate where this site should be with this offset
            target_pos = expected_pos + offset
            
            # Check if target position is within sequence bounds
            if 1 <= target_pos <= len(sequence):
                actual_aa = sequence[target_pos - 1]
                if actual_aa == aa:
                    mapped_sites.append((aa, target_pos))
                    end_idx = i + 1
                else:
                    # Amino acid mismatch - continuation breaks
                    valid = False
                    break
            else:
                # Position out of range - continuation breaks
                valid = False
                break
        
        # Track the best continuation (longest)
        if len(mapped_sites) > len(best_mapped_sites):
            best_mapped_sites = mapped_sites
            best_offset = offset
            best_end_index = end_idx
    
    return best_mapped_sites, best_offset, best_end_index


def align_binding_sites_to_sequence(binding_sites_str, sequence):
    """
    
    Align binding sites to the sequence using offset calculation.
       Process for single site:
        1. Check if the amino acid matches at the given position
        2. Keep if match, discard if not

    Process for multiple sites:
        1. Find longest continuation of binding sites with same offset
        2. Keep the actual sequence positions (not renumbered from 1)
        3. If continuation length > 1, keep them
        4. If continuation length <= 1, discard those sites
        5. Repeat for remaining binding sites after the continuation
    
    Returns:
        tuple: (formatted_sites, kept_count, dropped_count)
        formatted_sites: str of binding sites with actual sequence positions, or None
        kept_count: number of individual sites kept
        dropped_count: number of individual sites dropped
    """
    if pd.isna(binding_sites_str) or binding_sites_str == '':
        return None, 0, 0
    
    if pd.isna(sequence) or sequence == '':
        return None, 0, 0
    
    # Parse binding sites
    original_sites = parse_binding_sites(binding_sites_str)
    total_sites = len(original_sites)
    
    # Handle single binding site case
    if len(original_sites) == 1:
        aa, pos = original_sites[0]
        valid, actual_pos = validate_single_binding_site(aa, pos, sequence)
        if valid:
            return format_binding_sites([(aa, actual_pos)]), 1, 0
        else:
            return None, 0, 1
    
    # Handle multiple binding sites
    all_mapped_sites = []
    remaining_sites = original_sites.copy()
    
    while len(remaining_sites) > 0:
        # If only one site remains, validate it directly
        if len(remaining_sites) == 1:
            aa, pos = remaining_sites[0]
            valid, actual_pos = validate_single_binding_site(aa, pos, sequence)
            if valid:
                all_mapped_sites.append((aa, actual_pos))
            break
        
        # Find longest continuation from the start of remaining sites
        mapped_sites, offset, end_index = find_longest_continuation(remaining_sites, sequence)
        
        # Only keep continuations longer than 1
        if len(mapped_sites) > 1:
            # Keep actual sequence positions (no renumbering)
            all_mapped_sites.extend(mapped_sites)
            
            # Remove processed sites and continue
            remaining_sites = remaining_sites[end_index:]
        else:
            # No valid continuation found, discard the first site and try again
            remaining_sites = remaining_sites[1:]
    
    kept_count = len(all_mapped_sites)
    dropped_count = total_sites - kept_count
    
    if not all_mapped_sites:
        return None, 0, dropped_count
    
    return format_binding_sites(all_mapped_sites), kept_count, dropped_count


@timer
def merge_sequences_with_qbiolip(qbio_flat, easy_seq, diff_seq):
    """
    Merge Q-BioLiP data with extracted sequences.
    
    The Q-BioLiP data has been flattened so each row represents one binding site.
    Each row contains:
        - Assembly ID (PDB identifier)
        - PDB Chain (chain ID where the binding site is located)
    
    We need to attach the corresponding protein sequence for that chain.
    Sequences come from two sources:
        - easy_seq: Sequences from single PDB files
        - diff_seq: Sequences from PDB bundles (after chain mapping)
    
    Args:
        qbio_flat: DataFrame with flattened Q-BioLiP binding site data
                  Contains 'Assembly ID' and 'PDB Chain' columns
        easy_seq: DataFrame from extract_sequences.py
                  Columns: ['Assembly ID', 'Chain', 'Sequence']
        diff_seq: DataFrame from extract_sequences.py
                  Columns: ['Assembly ID', 'Original Chain ID', 'Sequence']
    
    Returns:
        DataFrame with all original columns plus a new 'Sequence' column
    """
    print("Merging sequences with Q-BioLiP data...")
    
    # =========================================================================
    # STEP 1: Merge with easy_seq (single PDB files)
    # =========================================================================
    qbio_merged = pd.merge(
        qbio_flat,
        easy_seq,
        left_on=["Assembly ID", "PDB Chain"],
        right_on=["Assembly ID", "Chain"],
        how="left"
    )
    
    # =========================================================================
    # STEP 2: Merge with diff_seq (bundled PDB files)
    # =========================================================================
    qbio_seq = pd.merge(
        qbio_merged,
        diff_seq,
        left_on=["Assembly ID", "PDB Chain"],
        right_on=["Assembly ID", "Original Chain ID"],
        how="left"
    )
    
    # =========================================================================
    # STEP 3: Combine sequences from both sources
    # =========================================================================
    qbio_seq["Sequence"] = qbio_seq["Sequence_x"].combine_first(qbio_seq["Sequence_y"])
    
    # =========================================================================
    # STEP 4: Clean up intermediate columns
    # =========================================================================
    cols_to_drop = ["Sequence_x", "Sequence_y", "Chain_y", "Original Chain ID"]
    existing_cols = [col for col in cols_to_drop if col in qbio_seq.columns]
    qbio_seq = qbio_seq.drop(columns=existing_cols)
    
    if "Chain_x" in qbio_seq.columns:
        qbio_seq = qbio_seq.rename(columns={"Chain_x": "QbioLiP Chain ID"})
    
    # =========================================================================
    # STEP 5: Remove rows without sequences
    # =========================================================================
    original_count = len(qbio_seq)
    qbio_seq = qbio_seq[~(qbio_seq["Sequence"].isna() | (qbio_seq["Sequence"] == ""))].copy()
    print(f"  Removed {original_count - len(qbio_seq)} rows without sequences")
    
    # =========================================================================
    # STEP 6: Print column headers for debugging
    # =========================================================================
    print(f"\n  Available columns for binding site alignment:")
    print(f"    {qbio_seq.columns.tolist()}")
    
    # Identify binding site column
    binding_site_col = None
    possible_binding_cols = ['Binding Site', 'Binding_site_original', 'Binding_site_pdb', 'Binding sites']
    
    for col in possible_binding_cols:
        if col in qbio_seq.columns:
            binding_site_col = col
            print(f"  Using binding site column: '{binding_site_col}'")
            break
    
    if binding_site_col is None:
        print(f"  WARNING: No binding site column found. Available columns: {qbio_seq.columns.tolist()}")
        return qbio_seq
    
    # Show sample before alignment
    sample_sites = qbio_seq[binding_site_col].dropna().iloc[0] if len(qbio_seq) > 0 else None
    sample_seq = qbio_seq['Sequence'].dropna().iloc[0] if len(qbio_seq) > 0 else None
    print(f"  Sample binding site before alignment: {sample_sites}")
    
    # =========================================================================
    # STEP 7: Align binding sites to sequences
    # =========================================================================
    print("\n  Aligning binding sites to sequences...")
    
    aligned_sites = []
    total_kept_sites = 0
    total_dropped_sites = 0
    rows_with_success = 0
    
    for idx, row in qbio_seq.iterrows():
        original_sites = row.get(binding_site_col)
        sequence = row.get('Sequence')
        
        if pd.notna(original_sites) and original_sites != '' and pd.notna(sequence) and sequence != '':
            aligned, kept, dropped = align_binding_sites_to_sequence(str(original_sites), str(sequence))
            
            if aligned:
                aligned_sites.append(aligned)
                rows_with_success += 1
                total_kept_sites += kept
                total_dropped_sites += dropped
            else:
                aligned_sites.append(None)
                total_dropped_sites += dropped
        else:
            aligned_sites.append(None)
    
    # =========================================================================
    # STEP 8: Create backup and replace binding site column
    # =========================================================================
    print(f"\n  Alignment Statistics (by individual amino acid binding sites):")
    print(f"    Total binding sites kept: {total_kept_sites:,}")
    print(f"    Total binding sites dropped: {total_dropped_sites:,}")
    total_sites = total_kept_sites + total_dropped_sites
    if total_sites > 0:
        print(f"    Percentage kept: {total_kept_sites/total_sites*100:.1f}%")
    print(f"    Rows with at least one aligned site: {rows_with_success:,}")
    print(f"    Rows with no aligned sites: {len(qbio_seq) - rows_with_success:,}")
    
    # Create backup of original column
    qbio_seq[f'{binding_site_col}_original'] = qbio_seq[binding_site_col]
    
    # Replace with aligned versions (or None if alignment failed)
    qbio_seq[binding_site_col] = aligned_sites
    
    # =========================================================================
    # STEP 9: Drop rows with no binding sites after alignment
    # =========================================================================
    rows_before = len(qbio_seq)
    
    # Keep only rows that have non-null binding sites after alignment
    qbio_seq = qbio_seq[qbio_seq[binding_site_col].notna()].copy()
    
    rows_dropped = rows_before - len(qbio_seq)
    print(f"\n  Rows dropped due to no binding sites after alignment: {rows_dropped:,}")
    print(f"  Rows kept with binding sites: {len(qbio_seq):,}")
    
    # Show sample after alignment
    if len(qbio_seq) > 0:
        sample_aligned = qbio_seq[binding_site_col].dropna().iloc[0]
        print(f"\n  Sample binding site after alignment: {sample_aligned}")
    
    return qbio_seq