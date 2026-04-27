# merge_biolip_datasets.py README

## Overview

This module performs the final merging of Q-BioLiP and BioLiP2 datasets after alignment to DescribePROT. It creates a unified dataset with separate binding site columns for each source, combines PDB IDs for overlapping binding sites, and generates JSON output with ligand statistics.

## Purpose

After both datasets have been aligned and renumbered to DescribePROT coordinates, this module:

1. Merges the datasets on (UniProt_ID, Ligand_ID, PDB_ID)
2. Keeps binding sites from both sources as separate columns
3. Combines affinity values (preferring BioLiP2 when available)
4. Groups binding sites to combine PDB IDs (e.g., "7bgn_7bgm")
5. Creates final JSON output and ligand count statistics

## Key Features

| Feature | Description |
|---------|-------------|
| **Separate Source Columns** | Maintains `Binding_sites_QBioLiP` and `Binding_sites_BioLiP2` columns |
| **PDB ID Merging** | Same binding site from multiple PDBs becomes `7bgn_7bgm` |
| **Affinity Priority** | BioLiP2 values preferred over Q-BioLiP when both exist |
| **Ligand Statistics** | Calculates unique binding sites per ligand |
| **Flexible Output** | Optional saving of intermediate CSV and ligand counts |

## Algorithm Flow

1. Load aligned DataFrames for Q-BioLiP and BioLiP2
2. Remove duplicate columns from both datasets
3. Standardize affinity column names
4. Merge on (UniProt_ID, Ligand_ID, PDB_ID) with outer join
5. Combine affinity values (BioLiP2 preferred)
6. Group by (UniProt_ID, Ligand_ID, PDB_ID) to combine binding sites
7. Explode binding sites into individual records (one per site)
8. Group by (UniProt_ID, Ligand_ID, Binding_site) to combine PDB IDs
9. Generate JSON output and ligand statistics

## Function Documentation

### `parse_binding_sites(binding_sites_str)`

Parses space-separated binding sites into a list.

| Parameter | Type | Description |
|-----------|------|-------------|
| `binding_sites_str` | str | String like "L84 V87 Y88" or None |

**Returns:** List of site strings like `['L84', 'V87', 'Y88']` or empty list

### `parse_binding_sites_to_set(binding_sites_str)`

Parses binding sites into a set, used when combining sites from multiple rows.

**Returns:** Set of site strings like `{'L84', 'V87', 'Y88'}`

### `get_position(site)`

Extracts numeric position from a binding site string.

| Parameter | Type | Description |
|-----------|------|-------------|
| `site` | str | String like "L176" or "M78" |

**Returns:** Integer position (176) or 0 if not found

### `combine_with_separate_columns(df, group_cols, label)`

Groups by specified columns while keeping Q-BioLiP and BioLiP2 binding sites as separate columns. Adds a merged column combining unique sites from both sources.

| Parameter | Type | Description |
|-----------|------|-------------|
| `df` | pd.DataFrame | DataFrame with `Binding_sites_QBioLiP` and `Binding_sites_BioLiP2` |
| `group_cols` | list | Columns to group by (e.g., `['UniProt_ID', 'Ligand_ID', 'PDB_ID']`) |
| `label` | str | Label for printing statistics |

**Returns:** DataFrame with grouped binding sites as space-separated strings

**Example transformation:**

Input (multiple rows):

Row1: QBioLiP sites: "L84", BioLiP2 sites: None

Row2: QBioLiP sites: "V87", BioLiP2 sites: "V87"

Output (one row):
Binding_sites_QBioLiP: "L84 V87", 
Binding_sites_BioLiP2: "V87", 
Binding_sites_Merged: "L84 V87"


### `create_final_json_with_combined_pdbs(exploded_df, original_df, output_file, ligand_counts_file, save_ligand_counts=True)`

Creates final JSON with combined PDB IDs for overlapping binding sites.

| Parameter | Type | Description |
|-----------|------|-------------|
| `exploded_df` | pd.DataFrame | DataFrame with exploded binding sites (one row per binding site) |
| `original_df` | pd.DataFrame | Original DataFrame with separate source columns |
| `output_file` | Path | Path for JSON output |
| `ligand_counts_file` | Path | Path for Excel output |
| `save_ligand_counts` | bool | Whether to save the Excel file |

**Returns:** List of JSON records

**Operations:**
1. Groups by (UniProt_ID, Ligand_ID, Binding_site)
2. Combines PDB IDs with underscore (e.g., "7bgn_7bgm")
3. Takes first non-null affinity value per group
4. Saves to JSON and optionally Excel

### `merge_biolip_datasets(qbiolip_df, biolip2_df, remove_duplicates=True, save_merged_csv=True, save_ligand_counts=True)`

Main orchestration function that merges both datasets.

| Parameter | Type | Description |
|-----------|------|-------------|
| `qbiolip_df` | pd.DataFrame | Processed Q-BioLiP DataFrame (after alignment) |
| `biolip2_df` | pd.DataFrame | Processed BioLiP2 DataFrame (after alignment) |
| `remove_duplicates` | bool | Whether to remove duplicate columns |
| `save_merged_csv` | bool | Whether to save the merged CSV file |
| `save_ligand_counts` | bool | Whether to save the ligand counts Excel file |

**Returns:** Merged DataFrame with separate binding site columns

## Input/Output Specifications

### Input DataFrame Requirements

**Q-BioLiP DataFrame requires:**

- `UniProt_ID`
- `Ligand_ID`
- `PDB_ID`
- `Binding_sites_QBioLiP`
- `describePROT_sequence`
- Affinity columns: 
    - `Binding_affinity_MOAD`
    - `Binding_affinity_PDBbind`
    - `Binding_affinity_BindingDB`

**BioLiP2 DataFrame requires:**

- `UniProt_ID`
- `Ligand_ID`
- `PDB_ID`
- `Binding_sites_BioLiP2`
- `describePROT_sequence`
- Affinity columns: 
    - `Binding_affinity_MOAD`
    - `Binding_affinity_PDBbind`
    - `Binding_affinity_BindingDB`

### Output Files

| File | Description |
|------|-------------|
| `final.json` | JSON records with combined PDB IDs |
| `ligand_binding_site_counts.xlsx` | Per-ligand unique binding site counts |
| `merged_biolip_datasets_by_pdb_separate.csv` | Intermediate merged CSV (if `save_merged_csv=True`) |

### JSON Output Format Example

```json
{
  "UniProt_ID": "A0A072U2X9",
  "Ligand_ID": "AMP",
  "Binding_site": "M78",
  "PDB_ID": "7bgn_7bgm",
  "DrugBank": "DB00131",
  "Binding_MOAD": null,
  "Binding_PDBbind": null,
  "Binding_BindingDB": null
}
```
## Column Preservation Logic

Affinity and metadata columns are preserved by:

1. Identifying columns to preserve - Scans for patterns like BindingMOAD, PDBbind, BindingDB, DrugBank, etc.

2. Taking first non-null value - When grouping, takes the first valid affinity value (all rows in group should have same value)

3. Preferring BioLiP2 - When merging, BioLiP2 affinity values take precedence over Q-BioLiP

4. Filtering empty values - Converts empty strings, <NA>, and "None" to actual null in JSON

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| `'Ligand_ID' column not found` | DataFrame missing `Ligand_ID` | Ensure column exists or rename alternative |
| Missing `PDB_ID` column | PDB identifier not present | Verify input DataFrames have `PDB_ID` |
| Duplicate columns | Multiple columns with same name | Automatically removed |