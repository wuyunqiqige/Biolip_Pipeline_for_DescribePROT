#!/usr/bin/env python3
"""
Merge Q-BioLiP and BioLiP2 datasets after alignment to describePROT

This module performs the final merging of both datasets and creates:
    1. A merged DataFrame with separate binding site columns for each source
    2. A JSON file with combined PDB IDs for overlapping binding sites
    3. An Excel file with per-ligand binding site counts

Key operations:
    - Merge on (UniProt_ID, Ligand_ID, PDB_ID)
    - Keep binding sites from both sources as separate columns
    - Combine affinity values (prefer BioLiP2)
    - Explode binding sites into individual records
    - Group by (UniProt, Ligand, Binding_site) to combine PDB IDs
"""

import pandas as pd
import json
import re
from pathlib import Path
from utils import timer


# =============================================================================
# BINDING SITE PARSING UTILITIES
# =============================================================================

def parse_binding_sites(binding_sites_str):
    """
    Parse space-separated binding sites into a list.
    
    Args:
        binding_sites_str: String like "L84 V87 Y88" or None
    
    Returns:
        List of site strings: ['L84', 'V87', 'Y88'] or empty list
    """
    if pd.isna(binding_sites_str) or binding_sites_str == '':
        return []
    return binding_sites_str.strip().split()


def parse_binding_sites_to_set(binding_sites_str):
    """
    Parse binding sites string into a set of individual site strings.
    Used when combining sites from multiple rows.
    
    Args:
        binding_sites_str: String like "L84 V87 Y88"
    
    Returns:
        Set of site strings: {'L84', 'V87', 'Y88'}
    """
    if pd.isna(binding_sites_str) or binding_sites_str == '':
        return set()
    
    sites = set()
    for part in str(binding_sites_str).strip().split():
        part = part.strip()
        if part:
            sites.add(part)
    return sites


def get_position(site):
    """
    Extract numeric position from a binding site string.
    
    Args:
        site: String like "L176" or "M78"
    
    Returns:
        Integer position (176) or 0 if not found
    """
    match = re.search(r'\d+', site)
    return int(match.group()) if match else 0

def print_detailed_overlap_statistics(original_df, dataset_name="BINDING SITE OVERLAP"):
    """
    Print comprehensive overlap statistics between Q-BioLiP and BioLiP2 datasets.
    
    Overlap is defined as exact match of:
        - UniProt_ID
        - Ligand_ID  
        - Individual binding site (e.g., 'K57', 'D116')
    
    Args:
        original_df: DataFrame with 'Binding_sites_QBioLiP', 'Binding_sites_BioLiP2',
                    'UniProt_ID', and 'Ligand_ID' columns
        dataset_name: Name for printing header
    """
    if original_df is None or original_df.empty:
        print("  No data for overlap statistics")
        return None
    
    print("\n" + "=" * 80)
    print(f"  {dataset_name}")
    print("=" * 80)
    
    # =========================================================================
    # STEP 1: Build sets of (UniProt_ID, Ligand_ID, Binding_site)
    # =========================================================================
    qbio_sites = set()
    biolip2_sites = set()
    
    # Track counts BEFORE deduplication (raw site occurrences)
    qbio_raw_count = 0
    biolip2_raw_count = 0
    
    for _, row in original_df.iterrows():
        uniprot = row['UniProt_ID']
        ligand = row['Ligand_ID']
        
        # Process Q-BioLiP binding sites
        qbio_sites_str = row.get('Binding_sites_QBioLiP', '')
        if pd.notna(qbio_sites_str) and qbio_sites_str:
            sites = parse_binding_sites(qbio_sites_str)
            qbio_raw_count += len(sites)
            for site in sites:
                qbio_sites.add((uniprot, ligand, site))
        
        # Process BioLiP2 binding sites
        biolip2_sites_str = row.get('Binding_sites_BioLiP2', '')
        if pd.notna(biolip2_sites_str) and biolip2_sites_str:
            sites = parse_binding_sites(biolip2_sites_str)
            biolip2_raw_count += len(sites)
            for site in sites:
                biolip2_sites.add((uniprot, ligand, site))
    
    # =========================================================================
    # STEP 2: Calculate overlap sets
    # =========================================================================
    common = qbio_sites & biolip2_sites
    only_qbio = qbio_sites - biolip2_sites
    only_biolip2 = biolip2_sites - qbio_sites
    total_unique = qbio_sites | biolip2_sites
    
    # =========================================================================
    # STEP 3: Print OVERALL STATISTICS
    # =========================================================================
    print(f"\n  OVERALL BINDING SITE STATISTICS")
    print(f"  {'-' * 60}")
    print(f"  Q-BioLiP raw sites (before dedup):         {qbio_raw_count:,}")
    print(f"  BioLiP2 raw sites (before dedup):          {biolip2_raw_count:,}")
    print(f"  {'-' * 60}")
    print(f"  Q-BioLiP UNIQUE (UniProt+Ligand+Site):     {len(qbio_sites):,}")
    print(f"  BioLiP2 UNIQUE (UniProt+Ligand+Site):      {len(biolip2_sites):,}")
    print(f"  {'-' * 60}")
    print(f"  Common to both datasets (UNIQUE):          {len(common):,}")
    print(f"  Only in Q-BioLiP (UNIQUE):                 {len(only_qbio):,}")
    print(f"  Only in BioLiP2 (UNIQUE):                  {len(only_biolip2):,}")
    print(f"  Total unique binding sites (union):        {len(total_unique):,}")
    
    # =========================================================================
    # STEP 4: Print OVERLAP PERCENTAGES
    # =========================================================================
    if len(total_unique) > 0:
        print(f"\n  OVERLAP ANALYSIS (Based on Unique Sites)")
        print(f"  {'-' * 60}")
        
        overlap_pct_of_total = (len(common) / len(total_unique)) * 100
        print(f"  Shared binding sites (% of total unique):         {overlap_pct_of_total:.1f}%")
        
        if len(qbio_sites) > 0:
            qbio_overlap_pct = (len(common) / len(qbio_sites)) * 100
            print(f"  Q-BioLiP unique sites also in BioLiP2:          {qbio_overlap_pct:.1f}% ({len(common):,}/{len(qbio_sites):,})")
        
        if len(biolip2_sites) > 0:
            biolip2_overlap_pct = (len(common) / len(biolip2_sites)) * 100
            print(f"  BioLiP2 unique sites also in Q-BioLiP:          {biolip2_overlap_pct:.1f}% ({len(common):,}/{len(biolip2_sites):,})")
        
        print(f"  {'-' * 60}")
        
        qbio_dedup_rate = (1 - len(qbio_sites) / qbio_raw_count) * 100 if qbio_raw_count > 0 else 0
        biolip2_dedup_rate = (1 - len(biolip2_sites) / biolip2_raw_count) * 100 if biolip2_raw_count > 0 else 0
        
        print(f"  Q-BioLiP redundancy removed: {qbio_dedup_rate:.1f}% ({qbio_raw_count - len(qbio_sites):,} duplicates)")
        print(f"  BioLiP2 redundancy removed: {biolip2_dedup_rate:.1f}% ({biolip2_raw_count - len(biolip2_sites):,} duplicates)")
    
    # =========================================================================
    # STEP 5: Return statistics for potential further use
    # =========================================================================
    return {
        'qbio_sites': qbio_sites,
        'biolip2_sites': biolip2_sites,
        'common': common,
        'only_qbio': only_qbio,
        'only_biolip2': only_biolip2,
        'total_unique': total_unique,
        'qbio_raw_count': qbio_raw_count,
        'biolip2_raw_count': biolip2_raw_count
    }


# =============================================================================
# GROUPING AND MERGING FUNCTIONS
# =============================================================================

def combine_with_separate_columns(df, group_cols, label):
    """
    Group by specified columns but KEEP Q-BioLiP and BioLiP2 binding sites as separate columns.
    Also adds a merged column combining unique sites from both sources.
    
    This function reduces the data from many rows (one per binding site per source)
    to one row per (UniProt, Ligand, PDB) combination.
    
    Args:
        df: DataFrame with Binding_sites_QBioLiP and Binding_sites_BioLiP2 columns
        group_cols: Columns to group by (e.g., ['UniProt_ID', 'Ligand_ID', 'PDB_ID'])
        label: Label for printing statistics
    
    Returns:
        DataFrame with grouped binding sites as space-separated strings
    
    Example transformation:
        Input (multiple rows):
            Row1: QBioLiP sites: "L84", BioLiP2 sites: None
            Row2: QBioLiP sites: "V87", BioLiP2 sites: "V87"
        
        Output (one row):
            Binding_sites_QBioLiP: "L84 V87"
            Binding_sites_BioLiP2: "V87"
            Binding_sites_Merged: "L84 V87"
    """
    print(f"\n  Grouping by: {group_cols}")
    
    # Identify affinity columns to preserve
    affinity_cols = [col for col in df.columns if 'affinity' in col.lower() or 
                     'BindingMOAD' in col or 'PDBbind' in col or 'BindingDB' in col]
    
    if affinity_cols:
        print(f"  Preserving affinity columns: {affinity_cols}")
    
    grouped = df.groupby(group_cols)
    results = []
    
    total_qbio_count = 0
    total_biolip2_count = 0
    total_merged_count = 0
    
    for keys, group in grouped:
        # Handle different number of keys
        if len(group_cols) == 2:
            uniprot, ligand = keys
            pdb = None
        else:
            uniprot, ligand, pdb = keys
        
        # Collect all binding sites from each source separately
        all_qbio_sites = set()
        all_biolip2_sites = set()
        
        # Collect affinity values (take first non-null, non-empty value)
        # All rows in the group should have the same affinity value
        affinity_values = {}
        for aff_col in affinity_cols:
            # Filter out NaN and empty/whitespace strings
            valid_values = group[aff_col].dropna()
            if len(valid_values) > 0:
                # Convert to string and strip whitespace
                valid_values = valid_values.astype(str).str.strip()
                # Filter out empty strings, "<NA>", and "None"
                valid_values = valid_values[~valid_values.isin(['', '<NA>', 'None'])]
                if len(valid_values) > 0:
                    affinity_values[aff_col] = valid_values.iloc[0]
                else:
                    affinity_values[aff_col] = None
            else:
                affinity_values[aff_col] = None
        
        # Collect binding sites from each row in the group
        for _, row in group.iterrows():
            qbio_sites = parse_binding_sites_to_set(row.get('Binding_sites_QBioLiP', ''))
            biolip2_sites = parse_binding_sites_to_set(row.get('Binding_sites_BioLiP2', ''))
            
            all_qbio_sites.update(qbio_sites)
            all_biolip2_sites.update(biolip2_sites)
        
        # Sort sites by position for consistent ordering
        sorted_qbio = sorted(all_qbio_sites, key=get_position)
        sorted_biolip2 = sorted(all_biolip2_sites, key=get_position)
        
        # Create merged column (union of both)
        all_sites = all_qbio_sites.union(all_biolip2_sites)
        sorted_merged = sorted(all_sites, key=get_position)
        
        # Convert to space-separated strings (or None if empty)
        qbio_str = ' '.join(sorted_qbio) if sorted_qbio else None
        biolip2_str = ' '.join(sorted_biolip2) if sorted_biolip2 else None
        merged_str = ' '.join(sorted_merged) if sorted_merged else None
        
        # Update counts
        total_qbio_count += len(sorted_qbio)
        total_biolip2_count += len(sorted_biolip2)
        total_merged_count += len(sorted_merged)
        
        # Get first row's other columns (they should be identical within group)
        first_row = group.iloc[0]
        
        # Build result dictionary
        if len(group_cols) == 2:
            result = {
                'UniProt_ID': uniprot,
                'Ligand_ID': ligand,
                'Binding_sites_QBioLiP': qbio_str,
                'Binding_sites_BioLiP2': biolip2_str,
                'Binding_sites_Merged': merged_str,
                'describePROT_sequence': first_row.get('describePROT_sequence'),
                'DrugBank': first_row.get('DrugBank')
            }
        else:
            result = {
                'UniProt_ID': uniprot,
                'Ligand_ID': ligand,
                'PDB_ID': pdb,
                'Binding_sites_QBioLiP': qbio_str,
                'Binding_sites_BioLiP2': biolip2_str,
                'Binding_sites_Merged': merged_str,
                'describePROT_sequence': first_row.get('describePROT_sequence'),
                'DrugBank': first_row.get('DrugBank')
            }
        
        # Add affinity values to result
        for aff_col, value in affinity_values.items():
            result[aff_col] = value
        
        results.append(result)
    
    final_df = pd.DataFrame(results)
    
    # Print statistics
    print(f"\n  Original rows: {len(df):,}")
    print(f"  Unique combinations: {len(final_df):,}")
    
    print(f"\n  {'='*50}")
    print(f"  BINDING SITE STATISTICS ({label})")
    print(f"  {'='*50}")
    print(f"  Q-BioLiP total sites: {total_qbio_count:,}")
    print(f"  BioLiP2 total sites: {total_biolip2_count:,}")
    print(f"  Merged total sites: {total_merged_count:,}")
    print(f"  {'='*50}")
    
    return final_df


# =============================================================================
# JSON CREATION FUNCTIONS
# =============================================================================

def create_final_json_with_combined_pdbs(exploded_df, original_df, output_file, ligand_counts_file, save_ligand_counts=True):
    """
    Create final JSON with combined PDB IDs for overlapping binding sites.
    
    This function:
        1. Groups by (UniProt_ID, Ligand_ID, Binding_site)
        2. Combines PDB IDs with underscore (e.g., "7bgn_7bgm")
        3. Takes first non-null affinity value per group
        4. Saves to JSON and Excel
    
    Args:
        exploded_df: DataFrame with exploded binding sites (one row per binding site)
        original_df: Original DataFrame with Binding_sites_QBioLiP and Binding_sites_BioLiP2
        output_file: Path for JSON output
        ligand_counts_file: Path for Excel output
        save_ligand_counts: Whether to save the Excel file
    
    Returns:
        List of JSON records
    """
    
    print("\n" + "=" * 70)
    print("CREATING FINAL JSON WITH COMBINED PDB IDs")
    print("=" * 70)
    
    # Use exploded_df for JSON and counts
    df = exploded_df.copy()
    
    #print(f"\n  DIAGNOSTIC: Input df shape: {df.shape}")
   # print(f"  DIAGNOSTIC: Input df columns: {df.columns.tolist()}")
    #print(f"  DIAGNOSTIC: Input df empty: {df.empty}")
    
    if df.empty:
        print("\n  WARNING: Input DataFrame is empty! No binding sites to process.")
        print("  Creating empty output files...")
        
        if save_ligand_counts:
            empty_counts = pd.DataFrame(columns=['Ligand_ID', 'QBioLiP_Sites', 'BioLiP2_Sites', 'Combined_Unique_Sites'])
            empty_counts.to_excel(ligand_counts_file, index=False)
        
        with open(output_file, 'w') as f:
            json.dump([], f)
        
        return []
    
    # Ensure Ligand_ID column exists
    if 'Ligand_ID' not in df.columns:
        print(f"\n  ERROR: 'Ligand_ID' column not found in DataFrame!")
        for col in df.columns:
            if 'ligand' in col.lower():
                print(f"  Found alternative: '{col}', renaming to 'Ligand_ID'")
                df = df.rename(columns={col: 'Ligand_ID'})
                break
        else:
            raise KeyError("'Ligand_ID' column not found")
    
    # Standardize ligand names
    df['Ligand_ID'] = df['Ligand_ID'].astype(str).str.upper().str.strip()
    df['Ligand_ID'] = df['Ligand_ID'].replace('DNA+RNA', 'DRNA')
    
    # =========================================================================
    # Calculate ligand binding site counts
    # =========================================================================
    print("\n  Calculating ligand binding site counts...")
    
    # Combined counts (from exploded_df)
    ligand_counts = df.groupby('Ligand_ID')['Binding_site'].nunique().reset_index()
    ligand_counts.columns = ['Ligand_ID', 'Combined_Unique_Sites']
    ligand_counts = ligand_counts.sort_values('Combined_Unique_Sites', ascending=False)
    
    # Separate counts from original_df (if available)
    if original_df is not None and not original_df.empty:
        # Q-BioLiP counts
        qbio_sites = []
        for _, row in original_df.iterrows():
            sites = parse_binding_sites(row.get('Binding_sites_QBioLiP', ''))
            for site in sites:
                qbio_sites.append({
                    'Ligand_ID': row['Ligand_ID'],
                    'Binding_site': site
                })
        
        df_qbio = pd.DataFrame(qbio_sites)
        if not df_qbio.empty:
            qbio_counts = df_qbio.groupby('Ligand_ID')['Binding_site'].nunique().reset_index()
            qbio_counts.columns = ['Ligand_ID', 'QBioLiP_Sites']
        else:
            qbio_counts = pd.DataFrame(columns=['Ligand_ID', 'QBioLiP_Sites'])
        
        # BioLiP2 counts
        biolip2_sites = []
        for _, row in original_df.iterrows():
            sites = parse_binding_sites(row.get('Binding_sites_BioLiP2', ''))
            for site in sites:
                biolip2_sites.append({
                    'Ligand_ID': row['Ligand_ID'],
                    'Binding_site': site
                })
        
        df_biolip2 = pd.DataFrame(biolip2_sites)
        if not df_biolip2.empty:
            biolip2_counts = df_biolip2.groupby('Ligand_ID')['Binding_site'].nunique().reset_index()
            biolip2_counts.columns = ['Ligand_ID', 'BioLiP2_Sites']
        else:
            biolip2_counts = pd.DataFrame(columns=['Ligand_ID', 'BioLiP2_Sites'])
        
        # Merge all counts
        all_counts = qbio_counts.merge(biolip2_counts, on='Ligand_ID', how='outer')
        all_counts = all_counts.merge(ligand_counts, on='Ligand_ID', how='outer')
        all_counts = all_counts.fillna(0).astype({'QBioLiP_Sites': int, 'BioLiP2_Sites': int, 'Combined_Unique_Sites': int})
        all_counts = all_counts.sort_values('Combined_Unique_Sites', ascending=False)
    else:
        all_counts = ligand_counts.copy()
        all_counts['QBioLiP_Sites'] = 0
        all_counts['BioLiP2_Sites'] = 0
        all_counts = all_counts[['Ligand_ID', 'QBioLiP_Sites', 'BioLiP2_Sites', 'Combined_Unique_Sites']]


     # =========================================================================
    # ADD OVERLAP STATISTICS HERE - RIGHT AFTER LIGAND COUNTS
    # =========================================================================
    if original_df is not None and not original_df.empty:
        print_detailed_overlap_statistics(original_df, "QBioLiP vs BioLiP2 BINDING SITE OVERLAP")

    # =========================================================================
    # Identify affinity columns for JSON
    # =========================================================================
    affinity_cols = []
    for col in df.columns:
        if 'affinity' in col.lower() or 'BindingMOAD' in col or 'PDBbind' in col or 'BindingDB' in col:
            affinity_cols.append(col)
    
    if affinity_cols:
        print(f"\n  Including affinity columns in JSON: {affinity_cols}")
    
    # =========================================================================
    # Create JSON records with combined PDB IDs
    # =========================================================================
    print("\n  Creating final JSON with combined PDB IDs...")
    
    # Group by (UniProt_ID, Ligand_ID, Binding_site) to combine PDB IDs
    try:
        binding_groups = df.groupby(['UniProt_ID', 'Ligand_ID', 'Binding_site'])
    except KeyError as e:
        print(f"  ERROR: Cannot group by - missing column: {e}")
        return []
    
    json_records = []
    
    for (uniprot, ligand, site), group in binding_groups:
        # Collect all PDB IDs for this binding site
        if 'PDB_ID' not in group.columns:
            continue
            
        pdb_ids = group['PDB_ID'].unique().tolist()
        
        # Combine PDB IDs with underscore if multiple
        if len(pdb_ids) > 1:
            combined_pdb = '_'.join(sorted(pdb_ids))
        else:
            combined_pdb = pdb_ids[0]
        
        # Get DrugBank (should be same across group)
        drugbank = group.iloc[0]['DrugBank'] if 'DrugBank' in group.columns else None
        drugbank = drugbank if pd.notna(drugbank) and drugbank != '' else None
        
        # Build the JSON record
        json_record = {
            'UniProt_ID': uniprot,
            'Ligand_ID': ligand,
            'Binding_site': site,
            'PDB_ID': combined_pdb,
            'DrugBank': drugbank
        }
        
        # Add affinity data (first non-null, non-empty value)
        for aff_col in affinity_cols:
            if aff_col in group.columns:
                non_null_values = group[aff_col].dropna()
                if len(non_null_values) > 0:
                    aff_value = non_null_values.iloc[0]
                    json_key = aff_col.replace('Binding_affinity_', 'Binding_')
                    
                    if pd.notna(aff_value):
                        str_value = str(aff_value).strip()
                        # Filter out empty strings, "<NA>", and "None"
                        if str_value != '' and not str_value.isspace() and str_value != '<NA>' and str_value != 'None':
                            json_record[json_key] = str_value
                        else:
                            json_record[json_key] = None
                    else:
                        json_record[json_key] = None
                else:
                    json_record[aff_col.replace('Binding_affinity_', 'Binding_')] = None
            else:
                json_record[aff_col.replace('Binding_affinity_', 'Binding_')] = None
        
        json_records.append(json_record)
    
    # Sort records for consistent output
    json_records = sorted(json_records, key=lambda x: (x['UniProt_ID'], x['Ligand_ID'], get_position(x['Binding_site'])))
    
    # Save to JSON
    with open(output_file, 'w') as f:
        json.dump(json_records, f, indent=2)
    
    print(f"  Saved JSON to: {output_file} ({len(json_records):,} records)")
    
    # =========================================================================
    # Print overlap statistics
    # =========================================================================
    if original_df is not None and not original_df.empty:
        print("\n" + "=" * 70)
        print("OVERLAP STATISTICS")
        print("=" * 70)
        
        qbio_set = set()
        for _, row in original_df.iterrows():
            sites = parse_binding_sites(row.get('Binding_sites_QBioLiP', ''))
            for site in sites:
                qbio_set.add((row['UniProt_ID'], row['Ligand_ID'], site))
        
        biolip2_set = set()
        for _, row in original_df.iterrows():
            sites = parse_binding_sites(row.get('Binding_sites_BioLiP2', ''))
            for site in sites:
                biolip2_set.add((row['UniProt_ID'], row['Ligand_ID'], site))
        
        common = qbio_set & biolip2_set
        only_qbio = qbio_set - biolip2_set
        only_biolip2 = biolip2_set - qbio_set
        total_unique = qbio_set | biolip2_set
        
        """
        print(f"Q-BioLiP tuples: {len(qbio_set):,}")
        print(f"BioLiP2 tuples: {len(biolip2_set):,}")
        print(f"Common to both: {len(common):,}")
        print(f"Only in Q-BioLiP: {len(only_qbio):,}")
        print(f"Only in BioLiP2: {len(only_biolip2):,}")
        print(f"Total unique tuples: {len(total_unique):,}")
        
        if len(total_unique) > 0:
            print(f"Overlapping Binding AA: {len(common)/len(total_unique)*100:.1f}%")
        """ 

    # Save ligand counts to Excel (conditional)
    if save_ligand_counts:
        all_counts.to_excel(ligand_counts_file, index=False)
        print(f"  Saved ligand counts to: {ligand_counts_file}")
    
    return json_records


# =============================================================================
# MAIN MERGE FUNCTION
# =============================================================================

@timer
def merge_biolip_datasets(qbiolip_df, biolip2_df, remove_duplicates=True, save_merged_csv=True, save_ligand_counts=True):
    """
    Merge Q-BioLiP and BioLiP2 datasets by UniProt_ID, Ligand_ID, and PDB_ID.
    
    This is the main orchestration function that:
        1. Removes duplicate columns from both datasets
        2. Standardizes column names
        3. Merges on (UniProt, Ligand, PDB)
        4. Combines affinity values (preferring BioLiP2)
        5. Groups binding sites into space-separated strings
        6. Explodes into individual records
        7. Creates final JSON with combined PDB IDs
    
    Args:
        qbiolip_df: Processed Q-BioLiP DataFrame (after alignment)
        biolip2_df: Processed BioLiP2 DataFrame (after alignment)
        remove_duplicates: Whether to remove duplicate columns (always True)
        save_merged_csv: Whether to save the merged CSV file
        save_ligand_counts: Whether to save the ligand counts Excel file
    
    Returns:
        Merged DataFrame with separate binding site columns
    """
    
    print(f"\n{'='*60}")
    print(f"MERGING DATASETS")
    print(f"{'='*60}")
    
    print(f"\nQ-BioLiP records: {len(qbiolip_df):,}")
    print(f"BioLiP2 records: {len(biolip2_df):,}")
    
    # Make copies
    qbio = qbiolip_df.copy()
    biolip2 = biolip2_df.copy()
    
    # =========================================================================
    # Remove duplicate columns
    # =========================================================================
    print(f"\n  Q-BioLiP original columns: {len(qbio.columns)}")
    qbio = qbio.loc[:, ~qbio.columns.duplicated()]
    print(f"  Q-BioLiP after removing duplicates: {len(qbio.columns)}")
    
    print(f"  BioLiP2 original columns: {len(biolip2.columns)}")
    biolip2 = biolip2.loc[:, ~biolip2.columns.duplicated()]
    print(f"  BioLiP2 after removing duplicates: {len(biolip2.columns)}")
    
    # =========================================================================
    # Identify affinity columns
    # =========================================================================
    qbio_affinity_cols = [col for col in qbio.columns if 'affinity' in col.lower() or 
                          'BindingMOAD' in col or 'PDBbind' in col or 'BindingDB' in col]
    biolip2_affinity_cols = [col for col in biolip2.columns if 'affinity' in col.lower() or 
                             'BindingMOAD' in col or 'PDBbind' in col or 'BindingDB' in col]
    
    if qbio_affinity_cols:
        print(f"\n  Q-BioLiP affinity columns: {qbio_affinity_cols}")
    if biolip2_affinity_cols:
        print(f"  BioLiP2 affinity columns: {biolip2_affinity_cols}")
    
    # =========================================================================
    # Ensure Ligand_ID column exists
    # =========================================================================
    if 'Ligand_ID' not in qbio.columns:
        if 'Ligand_ID_CCD' in qbio.columns:
            qbio['Ligand_ID'] = qbio['Ligand_ID_CCD']
        elif 'Ligand ID' in qbio.columns:
            qbio['Ligand_ID'] = qbio['Ligand ID']
        elif 'Ligand' in qbio.columns:
            qbio['Ligand_ID'] = qbio['Ligand']
    
    if 'Ligand_ID' not in biolip2.columns:
        if 'Ligand_ID_CCD' in biolip2.columns:
            biolip2['Ligand_ID'] = biolip2['Ligand_ID_CCD']
        elif 'Ligand' in biolip2.columns:
            biolip2['Ligand_ID'] = biolip2['Ligand']
    
    # =========================================================================
    # Rename binding site columns to keep them separate
    # =========================================================================
    if 'Binding_sites' in qbio.columns:
        qbio = qbio.rename(columns={'Binding_sites': 'Binding_sites_QBioLiP'})
    elif 'renumbered_binding_sites' in qbio.columns:
        qbio = qbio.rename(columns={'renumbered_binding_sites': 'Binding_sites_QBioLiP'})
    
    if 'Binding_sites' in biolip2.columns:
        biolip2 = biolip2.rename(columns={'Binding_sites': 'Binding_sites_BioLiP2'})
    elif 'renumbered_binding_sites' in biolip2.columns:
        biolip2 = biolip2.rename(columns={'renumbered_binding_sites': 'Binding_sites_BioLiP2'})
    
    # =========================================================================
    # Standardize affinity column names
    # =========================================================================
    affinity_rename_map = {
        'BindingMOAD': 'Binding_affinity_MOAD',
        'PDBbind-CN': 'Binding_affinity_PDBbind',
        'BindingDB': 'Binding_affinity_BindingDB'
    }
    
    for old, new in affinity_rename_map.items():
        if old in qbio.columns and old != new:
            qbio = qbio.rename(columns={old: new})
        if old in biolip2.columns and old != new:
            biolip2 = biolip2.rename(columns={old: new})
    
    # =========================================================================
    # Build column lists for merge
    # =========================================================================
    base_cols = ['UniProt_ID', 'Ligand_ID', 'PDB_ID', 'Binding_sites_QBioLiP', 
                 'Binding_sites_BioLiP2', 'describePROT_sequence', 'DrugBank']
    
    for aff_col in ['Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB']:
        if aff_col in qbio.columns and aff_col not in base_cols:
            base_cols.append(aff_col)
        if aff_col in biolip2.columns and aff_col not in base_cols:
            base_cols.append(aff_col)
    
    qbio_cols = list(dict.fromkeys(base_cols))
    biolip2_cols = list(dict.fromkeys(base_cols))
    
    # Only keep columns that exist
    qbio = qbio[[col for col in qbio_cols if col in qbio.columns]]
    biolip2 = biolip2[[col for col in biolip2_cols if col in biolip2.columns]]
    
    # =========================================================================
    # Merge datasets
    # =========================================================================
    #print(f"\n  DIAGNOSTIC: Q-BioLiP columns before merge: {qbio.columns.tolist()}")
    #print(f"  DIAGNOSTIC: BioLiP2 columns before merge: {biolip2.columns.tolist()}")
    
    print("\n  Merging datasets on UniProt_ID, Ligand_ID, and PDB_ID...")
    merged_df = pd.merge(qbio, biolip2, on=['UniProt_ID', 'Ligand_ID', 'PDB_ID'], 
                         how='outer', suffixes=('_qbio', '_biolip2'))
    
    #print(f"\n  DIAGNOSTIC: Merged df shape: {merged_df.shape}")
    
    # Remove duplicate columns from merged_df
    if merged_df.columns.duplicated().any():
        dup_count = merged_df.columns.duplicated().sum()
        print(f"  Removing {dup_count} duplicate columns from merged DataFrame")
        merged_df = merged_df.loc[:, ~merged_df.columns.duplicated()]
    
    # =========================================================================
    # Combine columns (describePROT_sequence, DrugBank, affinity)
    # =========================================================================
    # Combine describePROT_sequence
    if 'describePROT_sequence_qbio' in merged_df.columns and 'describePROT_sequence_biolip2' in merged_df.columns:
        merged_df['describePROT_sequence'] = merged_df['describePROT_sequence_qbio'].combine_first(
            merged_df['describePROT_sequence_biolip2']
        )
    elif 'describePROT_sequence_qbio' in merged_df.columns:
        merged_df['describePROT_sequence'] = merged_df['describePROT_sequence_qbio']
    elif 'describePROT_sequence_biolip2' in merged_df.columns:
        merged_df['describePROT_sequence'] = merged_df['describePROT_sequence_biolip2']
    
    # Combine DrugBank
    if 'DrugBank_qbio' in merged_df.columns and 'DrugBank_biolip2' in merged_df.columns:
        merged_df['DrugBank'] = merged_df['DrugBank_qbio'].combine_first(
            merged_df['DrugBank_biolip2']
        )
    elif 'DrugBank_qbio' in merged_df.columns:
        merged_df['DrugBank'] = merged_df['DrugBank_qbio']
    elif 'DrugBank_biolip2' in merged_df.columns:
        merged_df['DrugBank'] = merged_df['DrugBank_biolip2']
    
    # Combine affinity columns (prefer BioLiP2)
    for aff_col in ['Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB']:
        col_qbio = f"{aff_col}_qbio"
        col_biolip2 = f"{aff_col}_biolip2"
        
        # Handle DataFrame to Series conversion if needed
        if col_qbio in merged_df.columns:
            if isinstance(merged_df[col_qbio], pd.DataFrame):
                merged_df[col_qbio] = merged_df[col_qbio].iloc[:, 0] if merged_df[col_qbio].shape[1] > 0 else pd.Series(dtype=object)
        
        if col_biolip2 in merged_df.columns:
            if isinstance(merged_df[col_biolip2], pd.DataFrame):
                merged_df[col_biolip2] = merged_df[col_biolip2].iloc[:, 0] if merged_df[col_biolip2].shape[1] > 0 else pd.Series(dtype=object)
        
        # Combine with preference for BioLiP2
        if col_qbio in merged_df.columns and col_biolip2 in merged_df.columns:
            merged_df[aff_col] = merged_df[col_biolip2].combine_first(merged_df[col_qbio])
            merged_df = merged_df.drop(columns=[col_qbio, col_biolip2])
        elif col_qbio in merged_df.columns:
            merged_df[aff_col] = merged_df[col_qbio]
            merged_df = merged_df.drop(columns=[col_qbio])
        elif col_biolip2 in merged_df.columns:
            merged_df[aff_col] = merged_df[col_biolip2]
            merged_df = merged_df.drop(columns=[col_biolip2])
    
    # Drop intermediate columns
    cols_to_drop = []
    for col in ['describePROT_sequence_qbio', 'describePROT_sequence_biolip2', 
                'DrugBank_qbio', 'DrugBank_biolip2']:
        if col in merged_df.columns:
            cols_to_drop.append(col)
    
    if cols_to_drop:
        merged_df = merged_df.drop(columns=cols_to_drop)
    
    # =========================================================================
    # Filter out rows with no binding sites
    # =========================================================================
    has_qbio = merged_df['Binding_sites_QBioLiP'].notna() if 'Binding_sites_QBioLiP' in merged_df.columns else False
    has_biolip2 = merged_df['Binding_sites_BioLiP2'].notna() if 'Binding_sites_BioLiP2' in merged_df.columns else False
    has_any = has_qbio | has_biolip2
    original_count = len(merged_df)
    merged_df = merged_df[has_any].copy()
    print(f"  Removed {original_count - len(merged_df):,} rows with no binding sites")
    
    print(f"\n  Merge complete: {len(merged_df):,} rows")
    
    # =========================================================================
    # Create final dataset grouped by (UniProt, Ligand, PDB)
    # =========================================================================
    print("\n" + "=" * 70)
    print("FINAL DATASET: GROUPED BY (UNIPROT_ID, LIGAND_ID, PDB_ID) - SEPARATE COLUMNS")
    print("=" * 70)
    final_df = combine_with_separate_columns(
        merged_df,
        ['UniProt_ID', 'Ligand_ID', 'PDB_ID'],
        "BY PDB (SEPARATE)"
    )
    
    # =========================================================================
    # Diagnostic checks
    # =========================================================================
    #print(f"\n  DIAGNOSTIC: final_df shape: {final_df.shape}")
    #print(f"  DIAGNOSTIC: final_df columns: {final_df.columns.tolist()}")
    #print(f"  DIAGNOSTIC: Does final_df have 'Ligand_ID'? {'Ligand_ID' in final_df.columns}")
    
    final_affinity_cols = [col for col in final_df.columns if 'affinity' in col.lower() or 
                           'BindingMOAD' in col or 'PDBbind' in col or 'BindingDB' in col]
    if final_affinity_cols:
        print(f"  Affinity columns in final output: {final_affinity_cols}")
    
    # =========================================================================
    # Create exploded DataFrame for JSON
    # =========================================================================
    json_output = Path("final.json")
    ligand_counts_output = Path("ligand_binding_site_counts.xlsx")
    
    exploded_rows = []
    
    # Identify affinity columns in final_df
    affinity_cols_in_final = [col for col in final_df.columns if 'affinity' in col.lower() or 
                              'BindingMOAD' in col or 'PDBbind' in col or 'BindingDB' in col]
    
    if affinity_cols_in_final:
        print(f"\n  Including affinity columns in exploded DataFrame: {affinity_cols_in_final}")
    
    for idx, row in final_df.iterrows():
        # Get Ligand_ID safely
        ligand_id = row.get('Ligand_ID', None)
        if ligand_id is None:
            for col in row.index:
                if 'ligand' in col.lower():
                    ligand_id = row[col]
                    break
        
        if ligand_id is None:
            continue
        
        # Get affinity values for this row
        affinity_values = {}
        for aff_col in affinity_cols_in_final:
            value = row.get(aff_col)
            if pd.notna(value):
                str_value = str(value).strip()
                if str_value != '' and not str_value.isspace() and str_value != '<NA>' and str_value != 'None':
                    affinity_values[aff_col] = str_value
                else:
                    affinity_values[aff_col] = None
            else:
                affinity_values[aff_col] = None
        
        # Process Q-BioLiP binding sites
        qbio_sites = parse_binding_sites(row.get('Binding_sites_QBioLiP', ''))
        for site in qbio_sites:
            record = {
                'UniProt_ID': row['UniProt_ID'],
                'Ligand_ID': ligand_id,
                'Binding_site': site,
                'PDB_ID': row.get('PDB_ID'),
                'DrugBank': row.get('DrugBank')
            }
            record.update(affinity_values)
            exploded_rows.append(record)
        
        # Process BioLiP2 binding sites
        biolip2_sites = parse_binding_sites(row.get('Binding_sites_BioLiP2', ''))
        for site in biolip2_sites:
            record = {
                'UniProt_ID': row['UniProt_ID'],
                'Ligand_ID': ligand_id,
                'Binding_site': site,
                'PDB_ID': row.get('PDB_ID'),
                'DrugBank': row.get('DrugBank')
            }
            record.update(affinity_values)
            exploded_rows.append(record)
    
    exploded_df = pd.DataFrame(exploded_rows)
    
    #print(f"\n  DIAGNOSTIC: exploded_df shape: {exploded_df.shape}")
    #print(f"  DIAGNOSTIC: exploded_df columns: {exploded_df.columns.tolist()}")
    
    #if not exploded_df.empty:
        #print(f"  DIAGNOSTIC: First row of exploded_df:\n{exploded_df.iloc[0].to_dict()}")
    
    # Create JSON with combined PDB IDs (pass save_ligand_counts flag)
    create_final_json_with_combined_pdbs(exploded_df, final_df, json_output, ligand_counts_output, save_ligand_counts=save_ligand_counts)
    
    # Save merged CSV (conditional)
    if save_merged_csv:
        csv_output = Path("merged_biolip_datasets_by_pdb_separate.csv")
        final_df.to_csv(csv_output, index=False)
        print(f"\n  Saved CSV to: {csv_output}")
    
    return final_df