#!/usr/bin/env python3
"""
Part 6: Merge sequences with Q-BioLiP data

This module takes the flattened Q-BioLiP binding site data and merges it with
the amino acid sequences extracted from PDB files. It handles both single PDB
files and bundles, combining results from both sources.
"""

import pandas as pd
from utils import timer


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
    
    Logic:
        1. First attempt to match with easy_seq (single PDB files)
        2. For unmatched rows, attempt to match with diff_seq (bundles)
        3. Combine the two sequence sources
        4. Clean up intermediate columns
        5. Remove rows that couldn't be matched to any sequence
    """
    print("Merging sequences with Q-BioLiP data...")
    
    # =========================================================================
    # STEP 1: Merge with easy_seq (single PDB files)
    # =========================================================================
    # Match Q-BioLiP rows to sequences from single PDB files
    # Left join ensures we keep all Q-BioLiP rows even if no match found
    #
    # Key mapping:
    #   qbio_flat['Assembly ID']  -> easy_seq['Assembly ID']
    #   qbio_flat['PDB Chain']    -> easy_seq['Chain']
    #
    # Result adds 'Sequence_x' column (from easy_seq) to qbio_merged
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
    # For rows that didn't match a single PDB file, try matching with bundles
    #
    # Key mapping:
    #   qbio_merged['Assembly ID'] -> diff_seq['Assembly ID']
    #   qbio_merged['PDB Chain']   -> diff_seq['Original Chain ID']
    #
    # Result adds 'Sequence_y' column (from diff_seq)
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
    # combine_first() takes Sequence_x if available, otherwise uses Sequence_y
    # This gives priority to single PDB files, falling back to bundles
    #
    # Example:
    #   Row with Sequence_x = 'AGK...', Sequence_y = None  -> uses 'AGK...'
    #   Row with Sequence_x = None, Sequence_y = 'MVT...'  -> uses 'MVT...'
    #   Row with both -> uses Sequence_x (single PDB takes priority)
    qbio_seq["Sequence"] = qbio_seq["Sequence_x"].combine_first(qbio_seq["Sequence_y"])
    
    # =========================================================================
    # STEP 4: Clean up intermediate columns
    # =========================================================================
    # These columns were created by the merges and are no longer needed
    cols_to_drop = ["Sequence_x", "Sequence_y", "Chain_y", "Original Chain ID"]
    existing_cols = [col for col in cols_to_drop if col in qbio_seq.columns]
    qbio_seq = qbio_seq.drop(columns=existing_cols)
    
    # Rename Chain_x to something more descriptive
    # Chain_x came from the easy_seq merge and represents the chain ID
    if "Chain_x" in qbio_seq.columns:
        qbio_seq = qbio_seq.rename(columns={"Chain_x": "QbioLiP Chain ID"})
    
    # =========================================================================
    # STEP 5: Remove rows without sequences
    # =========================================================================
    # If a row couldn't be matched to either source, it has no sequence
    # These rows cannot be aligned and are removed
    original_count = len(qbio_seq)
    qbio_seq = qbio_seq[~(qbio_seq["Sequence"].isna() | (qbio_seq["Sequence"] == ""))].copy()
    print(f"  Removed {original_count - len(qbio_seq)} rows without sequences")
    
    return qbio_seq