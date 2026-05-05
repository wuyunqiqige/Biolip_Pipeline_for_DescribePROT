#!/usr/bin/env python3
"""
Process BioLiP2.txt file - Filter by Uniprot IDs and approved ligands

BioLiP2 is a tab-separated file with 21 predefined columns.
This module loads, filters, and standardizes BioLiP2 data for merging with Q-BioLiP.
"""

import pandas as pd
from tqdm import tqdm
from utils import timer

# =============================================================================
# COLUMN DEFINITIONS
# =============================================================================

# Complete list of 21 columns in BioLiP2.txt (tab-separated)
COLUMN_NAMES = [
    'PDB_ID',                           # 01 - PDB identifier
    'Receptor_chain',                   # 02 - Chain ID of the receptor
    'Resolution',                       # 03 - Structure resolution in Angstroms
    'Binding_site_number',              # 04 - Sequential number of binding site
    'Ligand_ID_CCD',                    # 05 - Ligand identifier (CCD format)
    'Ligand_chain',                     # 06 - Chain ID of the ligand
    'Ligand_serial_number',             # 07 - Serial number of ligand
    'Binding_site_PDB_numbering',       # 08 - Binding site with PDB numbering
    'Binding_site_renumbered',          # 09 - Binding site with renumbered coordinates
    'Catalytic_site_PDB_numbering',     # 10 - Catalytic site with PDB numbering
    'Catalytic_site_renumbered',        # 11 - Catalytic site with renumbered coordinates
    'EC_number',                        # 12 - Enzyme Commission number
    'GO_terms',                         # 13 - Gene Ontology terms
    'Binding_affinity_manual',          # 14 - Manually curated affinity
    'Binding_affinity_MOAD',            # 15 - Affinity from MOAD database
    'Binding_affinity_PDBbind',         # 16 - Affinity from PDBbind database
    'Binding_affinity_BindingDB',       # 17 - Affinity from BindingDB
    'UniProt_ID',                       # 18 - UniProt accession number
    'PubMed_ID',                        # 19 - PubMed reference ID
    'Residue_sequence_number',          # 20 - Sequential residue numbering
    'Receptor_sequence'                 # 21 - Full receptor amino acid sequence
]

# Columns to keep in final output after filtering
KEEP_COLUMNS = [
    'PDB_ID',
    'Binding_site_number',
    'Ligand_ID_CCD',
    'Binding_site_renumbered',          # Primary binding site column
    'Binding_affinity_MOAD',
    'Binding_affinity_PDBbind',
    'Binding_affinity_BindingDB',
    'UniProt_ID',
    'Receptor_sequence',
    'DrugBank'                          # Added from approved ligands file
]


# =============================================================================
# STRING CLEANING UTILITIES
# =============================================================================

def trim_string_columns(df):
    """
    Trim leading/trailing whitespace from all string columns.
    Converts empty strings and 'nan' strings to proper None/NA values.
    
    Args:
        df: Input DataFrame
    
    Returns:
        DataFrame with cleaned string columns
    """
    df = df.copy()
    
    for col in df.columns:
        if df[col].dtype == 'object':
            # Step 1: Convert to string and strip whitespace
            df[col] = df[col].astype(str).str.strip()
            
            # Step 2: Replace string 'nan' with pandas NA
            df[col] = df[col].replace('nan', pd.NA)
            
            # Step 3: Replace empty or whitespace-only strings with None
            df[col] = df[col].apply(lambda x: None if pd.notna(x) and (x == '' or x.isspace()) else x)
    
    return df


def standardize_ligand_names(df):
    """
    Standardize ligand names across BioLiP2 dataset.
    
    Args:
        df: DataFrame with 'Ligand_ID_CCD' column
    
    Returns:
        DataFrame with standardized ligand names
    
    Logic:
        1. Capitalizes all ligand IDs for case-insensitive matching
        2. Renames 'PEPTIDE' to 'III' (consistent with Q-BioLiP naming)
    """
    df = df.copy()
    
    # Capitalize all ligand IDs
    df['Ligand_ID_CCD'] = df['Ligand_ID_CCD'].astype(str).str.upper().str.strip()
    
    # Rename III to III_
    df['Ligand_ID_CCD'] = df['Ligand_ID_CCD'].replace('III', 'III_')
    # Rename PEPTIDE to III for consistency with Q-BioLiP
    df['Ligand_ID_CCD'] = df['Ligand_ID_CCD'].replace('PEPTIDE', 'III')
    
    return df


def standardize_biolip2_columns(df):
    """
    Standardize BioLiP2 column names for validation and alignment.
    
    Creates a 'Binding_site_original' alias for 'Binding_site_renumbered'
    to maintain consistency with the Q-BioLiP naming convention.
    
    Args:
        df: BioLiP2 DataFrame
    
    Returns:
        DataFrame with standardized column names
    """
    df = df.copy()
    
    # Create alias for binding site column
    if 'Binding_site_renumbered' in df.columns:
        df['Binding_site_original'] = df['Binding_site_renumbered']
    
    # Ensure sequence column exists (some versions may use 'seq')
    if 'Receptor_sequence' not in df.columns and 'seq' in df.columns:
        df['Receptor_sequence'] = df['seq']
    
    return df


def validate_biolip2_binding_sites(df):
    """
    Validate BioLiP2 binding sites against their sequences BEFORE alignment.
    
    This checks if the amino acid at each binding site position matches
    the expected amino acid in the receptor sequence.
    
    Args:
        df: BioLiP2 DataFrame with 'Binding_site_original' and 'Receptor_sequence'
    
    Returns:
        Original DataFrame (validation is for reporting only)
    
    Note: This is diagnostic only - no filtering is applied.
    """
    print("\n  Validating BioLiP2 binding sites (pre-alignment)...")
    
    df_validated = df.copy()
    
    # Find rows that have both binding sites and sequences
    has_binding_sites = df_validated['Binding_site_original'].notna() & (df_validated['Binding_site_original'] != '')
    has_sequence = df_validated['Receptor_sequence'].notna() & (df_validated['Receptor_sequence'] != '')
    
    total_sites = 0
    valid_sites = 0
    
    for idx, row in df_validated[has_binding_sites & has_sequence].iterrows():
        binding_site = str(row['Binding_site_original'])
        sequence = str(row['Receptor_sequence'])
        
        # Parse space-separated binding sites
        sites = binding_site.strip().split()
        total_sites += len(sites)
        
        for site in sites:
            if len(site) >= 2:
                expected_aa = site[0]          # First character is the amino acid
                try:
                    position = int(site[1:])   # Rest is the position number
                    idx_pos = position - 1      # Convert to 0-based index
                    if 0 <= idx_pos < len(sequence) and sequence[idx_pos].upper() == expected_aa.upper():
                        valid_sites += 1
                except ValueError:
                    pass  # Skip malformed site strings
    
    print(f"    Total binding sites: {total_sites:,}")
    print(f"    Valid sites (match sequence): {valid_sites:,}")
    if total_sites > 0:
        print(f"    Validation rate: {valid_sites/total_sites*100:.1f}%")
    
    return df_validated


# =============================================================================
# LOADING AND FILTERING FUNCTIONS
# =============================================================================

@timer
def load_biolip2_file(input_file):
    """
    Load BioLiP2.txt file with proper column names.
    
    BioLiP2 is tab-separated with no header row.
    
    Args:
        input_file: Path to BioLiP2.txt file
    
    Returns:
        DataFrame with named columns
    """
    print(f"Loading BioLiP2 file: {input_file}")
    
    # Read tab-separated file with no header, using predefined column names
    df = pd.read_csv(input_file, sep='\t', header=None, names=COLUMN_NAMES, low_memory=False)
    
    # Clean all string columns
    df = trim_string_columns(df)
    
    print(f"  Loaded {len(df):,} total records from BioLiP2")
    print(f"  Columns: {', '.join(df.columns.tolist()[:5])}...")
    
    return df


@timer
def filter_biolip2_by_uniprot(biolip2_df, describe_df):
    """
    Filter BioLiP2 by Uniprot IDs that exist in describePROT.
    
    Args:
        biolip2_df: BioLiP2 DataFrame
        describe_df: DescribePROT DataFrame with 'ACC' column
    
    Returns:
        Filtered DataFrame (only entries with matching UniProt IDs)
    """
    print("Filtering BioLiP2 by Uniprot IDs...")
    
    original_count = len(biolip2_df)
    
    # Create uppercase version for case-insensitive matching
    biolip2_df['UniProt_ID_upper'] = biolip2_df['UniProt_ID'].astype(str).str.upper().str.strip()
    describe_acc_upper = describe_df['ACC'].astype(str).str.upper().str.strip()
    
    # Filter to only UniProt IDs present in describePROT
    biolip2_filtered = biolip2_df[biolip2_df['UniProt_ID_upper'].isin(describe_acc_upper)].copy()
    biolip2_filtered = biolip2_filtered.drop(columns=['UniProt_ID_upper'])
    
    print(f"  After Uniprot filter: {len(biolip2_filtered):,} rows (removed {original_count - len(biolip2_filtered):,})")
    
    return biolip2_filtered


@timer
def filter_biolip2_by_ligands(biolip2_df, approved_ligands_file, ligand_id_col='Ligand ID', drugbank_col='DrugBank'):
    """
    Filter BioLiP2 by approved ligands and add DrugBank column.
    
    Args:
        biolip2_df: BioLiP2 DataFrame
        approved_ligands_file: Excel file with approved ligands
        ligand_id_col: Column name in Excel for ligand IDs
        drugbank_col: Column name in Excel for DrugBank IDs
    
    Returns:
        Filtered DataFrame with added 'DrugBank' column
    """
    print(f"Filtering BioLiP2 by approved ligands...")
    
    # Step 1: Load approved ligands
    approved = pd.read_excel(approved_ligands_file)
    
    # Step 2: Clean string columns
    for col in approved.columns:
        if approved[col].dtype == 'object':
            approved[col] = approved[col].astype(str).str.strip()
    
    # Step 3: Standardize ligand IDs for matching
    biolip2_df['Ligand_ID_CCD_upper'] = biolip2_df['Ligand_ID_CCD'].astype(str).str.upper().str.strip()
    approved['Ligand_upper'] = approved[ligand_id_col].astype(str).str.upper().str.strip()
    
    # Special case: rename PEPTIDE to III in approved ligands for matching
    approved['Ligand_upper'] = approved['Ligand_upper'].replace('PEPTIDE', 'III')
    
    # Step 4: Filter to only approved ligands
    original_count = len(biolip2_df)
    biolip2_filtered = biolip2_df[biolip2_df['Ligand_ID_CCD_upper'].isin(approved['Ligand_upper'])].copy()
    print(f"  After ligand filter: {len(biolip2_filtered):,} rows (removed {original_count - len(biolip2_filtered):,})")
    
    # Step 5: Add DrugBank column by mapping
    if drugbank_col in approved.columns:
        drugbank_map = approved.set_index('Ligand_upper')[drugbank_col].to_dict()
        biolip2_filtered['DrugBank'] = biolip2_filtered['Ligand_ID_CCD_upper'].map(drugbank_map)
        
        has_drugbank = biolip2_filtered['DrugBank'].notna().sum()
        print(f"  Rows with DrugBank info: {has_drugbank:,} ({has_drugbank/len(biolip2_filtered)*100:.1f}%)")
    else:
        biolip2_filtered['DrugBank'] = None
        print(f"  No DrugBank column found in approved ligands file")
    
    # Remove temporary column
    biolip2_filtered = biolip2_filtered.drop(columns=['Ligand_ID_CCD_upper'])
    
    return biolip2_filtered


# =============================================================================
# MAIN PROCESSING FUNCTION
# =============================================================================

@timer
def process_biolip2(biolip2_file, describe_df, approved_ligands_file, output_file=None):
    """
    Main function to process BioLiP2 file through the entire pipeline.
    
    Processing steps:
        1. Load BioLiP2.txt with proper column names
        2. Standardize ligand names (uppercase, PEPTIDE→III)
        3. Filter by UniProt IDs present in describePROT
        4. Filter by approved ligands
        5. Keep only essential columns
        6. Optionally save to CSV
    
    Args:
        biolip2_file: Path to BioLiP2.txt file
        describe_df: DescribePROT DataFrame for UniProt filtering
        approved_ligands_file: Excel file with approved ligands
        output_file: Optional path to save filtered results
    
    Returns:
        Filtered BioLiP2 DataFrame ready for alignment
    """
    
    # Step 1: Load raw data
    biolip2_df = load_biolip2_file(biolip2_file)
    
    # Step 2: Standardize ligand names
    biolip2_df = standardize_ligand_names(biolip2_df)
    
    # Step 3: Filter by UniProt IDs
    biolip2_df = filter_biolip2_by_uniprot(biolip2_df, describe_df)
    
    if len(biolip2_df) == 0:
        print("  No records remaining after Uniprot filtering. Exiting.")
        return pd.DataFrame()
    
    # Step 4: Filter by approved ligands
    biolip2_df = filter_biolip2_by_ligands(biolip2_df, approved_ligands_file)
    
    if len(biolip2_df) == 0:
        print("  No records remaining after ligand filtering. Exiting.")
        return pd.DataFrame()
    
    # Step 5: Keep only essential columns
    columns_to_keep = [col for col in KEEP_COLUMNS if col in biolip2_df.columns]
    final_df = biolip2_df[columns_to_keep].copy()
    
    print(f"  Final columns: {', '.join(final_df.columns.tolist())}")
    print(f"  Final records: {len(final_df):,}")
    
    # Step 6: Save if output file specified
    if output_file:
        print(f"  Saving to {output_file}...")
        final_df.to_csv(output_file, index=False)
        print(f"  Saved {len(final_df):,} records to {output_file}")
    
    return final_df