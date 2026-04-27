#!/usr/bin/env python3
"""
Parts 2-4: Process Q-BioLiP data, filter by ligands, and flatten binding sites

This module handles the complete Q-BioLiP processing pipeline:
- Loading and filtering Q-BioLiP data
- Filtering by approved ligands
- Flattening binding sites (exploding multi-site entries into individual rows)
- Standardizing column names and values
"""

import pandas as pd
from tqdm import tqdm
from utils import timer


# =============================================================================
# STRING CLEANING UTILITIES
# =============================================================================

def trim_string_columns(df):
    """
    Trim leading/trailing whitespace from all string columns.
    Also converts empty strings and 'nan' strings to proper None/NA values.
    
    Args:
        df: Input DataFrame
    
    Returns:
        DataFrame with cleaned string columns
    
    Logic:
        1. Iterates through all columns
        2. For object (string) columns, converts to string and strips whitespace
        3. Replaces string 'nan' with pandas NA
        4. Replaces empty or whitespace-only strings with None
    """
    df = df.copy()
    
    for col in df.columns:
        if df[col].dtype == 'object':
            # Step 1: Convert to string and strip whitespace
            df[col] = df[col].astype(str).str.strip()
            
            # Step 2: Replace the string 'nan' with actual pandas NA
            df[col] = df[col].replace('nan', pd.NA)
            
            # Step 3: Replace empty or whitespace-only strings with None
            # This ensures empty values become null in final JSON
            df[col] = df[col].apply(lambda x: None if pd.notna(x) and (x == '' or x.isspace()) else x)
    
    return df


def standardize_ligand_names(df):
    """
    Standardize ligand names across Q-BioLiP dataset.
    
    Args:
        df: DataFrame with 'Ligand ID' column
    
    Returns:
        DataFrame with standardized ligand names
    
    Logic:
        1. Capitalizes all ligand IDs (e.g., 'amp' -> 'AMP')
        2. Renames 'DNA+RNA' to 'DRNA' for consistency
    """
    df = df.copy()
    
    # Capitalize all ligand IDs for case-insensitive matching
    df['Ligand ID'] = df['Ligand ID'].astype(str).str.upper().str.strip()
    
    # Special case: rename DNA+RNA to DRNA
    df['Ligand ID'] = df['Ligand ID'].replace('DNA+RNA', 'DRNA')
    
    return df


def standardize_qbiolip_columns(df):
    """
    Standardize Q-BioLiP column names for validation and alignment.
    
    This creates consistent column names across the pipeline:
    - Original Q-BioLiP names → Standardized names
    - Affinity columns are renamed to a common format
    
    Args:
        df: Raw Q-BioLiP DataFrame
    
    Returns:
        DataFrame with renamed columns
    
    Critical: Affinity columns are renamed but NOT dropped!
        'BindingMOAD' → 'Binding_affinity_MOAD'
        'PDBbind-CN' → 'Binding_affinity_PDBbind'
        'BindingDB' → 'Binding_affinity_BindingDB'
    """
    df = df.copy()
    
    rename_map = {
        # Core identifiers
        'Uniprot ID': 'UniProt_ID',
        'Ligand ID': 'Ligand_ID',
        'PDB ID': 'PDB_ID',
        
        # Binding site columns (two versions: original and PDB-formatted)
        'Binding Site': 'Binding_site_original',
        'Binding Site PDB': 'Binding_site_pdb',
        'Site': 'Binding_site_number',
        
        # Sequence column
        'Sequence': 'Receptor_sequence',
        
        # Affinity columns - PRESERVED and renamed
        'BindingMOAD': 'Binding_affinity_MOAD',
        'PDBbind-CN': 'Binding_affinity_PDBbind',
        'BindingDB': 'Binding_affinity_BindingDB'
    }
    
    for old, new in rename_map.items():
        if old in df.columns:
            df[new] = df[old]
    
    return df


# =============================================================================
# MAIN PROCESSING FUNCTIONS
# =============================================================================

@timer
def load_and_filter_qbiolip(qbio_file, describe_df):
    """
    Load Q-BioLiP data and apply initial filters.
    
    Filtering steps:
        1. Load CSV and clean string columns
        2. Drop unnecessary columns (keep affinity columns)
        3. Filter to only UniProt IDs present in describePROT
        4. Filter to only entries with Interaction == 'yes'
    
    Args:
        qbio_file: Path to Q-BioLiP CSV file
        describe_df: DescribePROT DataFrame with 'ACC' column for filtering
    
    Returns:
        Filtered Q-BioLiP DataFrame
    
    Note: Affinity columns (BindingMOAD, PDBbind-CN, BindingDB) are NOT dropped
    """
    print(f"Loading Q-BioLiP data: {qbio_file}")
    qbio = pd.read_csv(qbio_file, low_memory=False)
    
    # Step 1: Clean all string columns (trim whitespace, convert empty to None)
    qbio = trim_string_columns(qbio)
    
    # Step 2: Drop columns we don't need
    # IMPORTANT: Affinity columns are NOT in this list!
    cols_to_drop = [
        "Stoichiometry", "Resolution (Å)", "Assembly Detail", 
        "Pubmed ID", "Relevant"
    ]
    cols_to_drop_existing = [col for col in cols_to_drop if col in qbio.columns]
    qbio_filtered = qbio.drop(columns=cols_to_drop_existing)
    
    # Step 3: Standardize ligand names (capitalize, DNA+RNA → DRNA)
    qbio_filtered = standardize_ligand_names(qbio_filtered)
    
    # Step 4: Filter by UniProt ID (only keep entries in describePROT)
    original_count = len(qbio_filtered)
    qbio_filtered = qbio_filtered[qbio_filtered["Uniprot ID"].isin(describe_df["ACC"])]
    print(f"  After Uniprot filter: {len(qbio_filtered)} rows (removed {original_count - len(qbio_filtered)})")
    
    # Step 5: Filter by Interaction status (only keep 'yes')
    original_count = len(qbio_filtered)
    qbio_filtered = qbio_filtered[qbio_filtered["Interaction"] == "yes"]
    qbio_filtered = qbio_filtered.drop(columns=["Interaction"])
    
    return qbio_filtered


@timer
def filter_by_ligands(qbio_df, approved_ligands_file, ligand_id_col='Ligand ID', drugbank_col='DrugBank'):
    """
    Filter Q-BioLiP data using approved ligands file and add DrugBank column.
    
    Args:
        qbio_df: Q-BioLiP DataFrame
        approved_ligands_file: Excel file with approved ligands
        ligand_id_col: Column name in Excel for ligand IDs (default: 'Ligand ID')
        drugbank_col: Column name in Excel for DrugBank IDs (default: 'DrugBank')
    
    Returns:
        Filtered DataFrame with added 'DrugBank' column
    
    Logic:
        1. Load approved ligands Excel file
        2. Clean and standardize ligand names
        3. Filter qbio_df to only approved ligands
        4. Map DrugBank IDs from approved ligands to qbio_df
    """
    print(f"Filtering by ligands using: {approved_ligands_file}")
    
    # Step 1: Load approved ligands
    approved = pd.read_excel(approved_ligands_file)
    
    # Step 2: Clean string columns in approved ligands
    for col in approved.columns:
        if approved[col].dtype == 'object':
            approved[col] = approved[col].astype(str).str.strip()
    
    # Step 3: Standardize ligand IDs in both DataFrames for matching
    qbio_df["Ligand ID"] = qbio_df["Ligand ID"].astype(str).str.upper().str.strip()
    approved['Ligand_upper'] = approved[ligand_id_col].astype(str).str.upper().str.strip()
    
    # Step 4: Filter to only approved ligands
    original_count = len(qbio_df)
    qbio_filtered = qbio_df[qbio_df['Ligand ID'].isin(approved['Ligand_upper'])].copy()
    print(f"  After ligand filter: {len(qbio_filtered)} rows (removed {original_count - len(qbio_filtered)})")
    
    # Step 5: Add DrugBank column by mapping from approved ligands
    if drugbank_col in approved.columns:
        drugbank_map = approved.set_index('Ligand_upper')[drugbank_col].to_dict()
        qbio_filtered['DrugBank'] = qbio_filtered['Ligand ID'].map(drugbank_map)
    else:
        qbio_filtered['DrugBank'] = None
        print(f"  No DrugBank column found in approved ligands file")
    
    return qbio_filtered


@timer
def flatten_binding_sites(qbio_df):
    """
    Flatten the binding site data by exploding on binding sites.
    
    Original Q-BioLiP rows may have multiple binding sites in one row.
    This function splits them into individual rows (one per binding site).
    
    Example:
        Original: Binding Site = "A:L84 V87 Y88"
        After flattening: 3 rows, one for each site
        
    Args:
        qbio_df: DataFrame with binding sites as space-separated strings
    
    Returns:
        Flattened DataFrame with one row per binding site
    
    Logic:
        1. Split space-separated binding sites into lists
        2. Use pandas explode() to create one row per list item
        3. Extract chain ID and residue from each binding site string
        4. Also process the PDB-formatted binding sites similarly
    """
    print("Flattening binding site data...")
    
    # Step 1: Split binding site strings into lists
    # Regex pattern: split on space only if followed by chain:residue pattern
    # Example: "A:L84 V87 Y88" → ["A:L84", "V87", "Y88"]
    split_patt = r' (?=[A-Za-z0-9]+:)'
    qbio_df["Binding Site"] = qbio_df["Binding Site"].str.split(split_patt)
    qbio_df["Binding Site PDB"] = qbio_df["Binding Site PDB"].str.split(split_patt)
    
    # Step 2: Explode - create one row per binding site
    original_count = len(qbio_df)
    qbio_flat = qbio_df.explode(["Binding Site", "Binding Site PDB"], ignore_index=True)
    print(f"  After exploding: {len(qbio_flat)} rows (from {original_count})")
    
    # Step 3: Extract chain and residue from binding site strings
    # Format: "CHAIN:RESIDUE" (e.g., "A:L84")
    qbio_flat[["Chain", "Binding Site"]] = (
        qbio_flat["Binding Site"]
        .str.extract(r'^([A-Za-z0-9]+):(.*)$')
    )
    
    # Trim extracted binding sites
    qbio_flat["Binding Site"] = qbio_flat["Binding Site"].str.strip()
    
    # Step 4: Extract PDB chain from PDB-formatted binding sites
    # Format: "CHAIN:RESIDUE" (same pattern)
    qbio_flat[["PDB Chain", "Binding Site PDB"]] = (
        qbio_flat["Binding Site PDB"]
        .str.extract(r'^([A-Za-z0-9]+):(.*)$')
    )
    
    # Trim extracted PDB binding sites
    qbio_flat["Binding Site PDB"] = qbio_flat["Binding Site PDB"].str.strip()
    
    print(f"  Final flattened rows: {len(qbio_flat)}")
    
    return qbio_flat