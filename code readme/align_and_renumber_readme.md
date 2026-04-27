# align_and_renumber.py README

## Overview

This module performs the task of mapping binding site positions from original PDB numbering to DescribePROT reference sequence coordinates using BLAST alignment.

## Purpose

Q-BioLiP and BioLiP2 databases provide binding site positions based on PDB file numbering, but DescribePROT uses different reference sequences. This module:

1. Aligns PDB sequences to DescribePROT references using BLAST

2. Creates position maps between the two numbering systems

3. Renumbers binding sites to DescribePROT coordinates

4. Validates amino acid matches (only keeps correct mappings)

5. Preserves affinity data and metadata through the process

## Key Features

| Feature | Description |
|---------|-------------|
| **Pairwise BLAST** | Each sequence aligned directly to its matching DescribePROT reference using Uniprot ID |
| **Position Mapping** | Creates dictionaries mapping PDB positions → DescribePROT positions |
| **Amino Acid Validation** | Only keeps binding sites where the amino acid matches the reference sequence |
| **Result Caching** | Stores alignment results in `.pkl` files to avoid redundant BLAST runs |
| **Column Preservation** | Maintains affinity and metadata through the alignment process using dictionary storage |
| **Batch Processing** | Groups identical (UniProt_ID, Sequence) pairs to minimize BLAST calls |

## Algorithm Flow
1. Input DataFrame + DescribePROT reference
2. Group by (UniProt_ID, Sequence)
3. Check Cache for existing alignments
4. For each unique pair:
   - If cached: load position map
   - If sequences identical: create 1:1 map
   - Else:
     - Run BLAST
     - Parse alignment
     - Create position map
5. Apply mapping to all rows in group
6. Validate each binding site (amino acid match)
7. Restore preserved columns (affinity, metadata)
8. Filter out failed alignments
9. Output: DataFrame with renumbered binding sites


## Function Documentation

### Core Alignment Functions

#### `align_and_renumber_dataset(df, describeprot_df, dataset_name, use_cache=True)`\

Main orchestration function that processes entire datasets.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `df` | pd.DataFrame | Input DataFrame with sequences and binding sites |
| `describeprot_df` | pd.DataFrame | DescribePROT reference with 'ACC' and 'seq' columns |
| `dataset_name` | str | "QBioLiP" or "BioLiP2" (affects column detection) |
| `use_cache` | bool | Whether to use/save alignment cache (.pkl files) |

**Returns:** DataFrame with renumbered binding sites and describePROT sequences

### BLAST Functions

#### `run_pairwise_blast(query_seq, query_id, target_seq, target_id)`
Runs BLAST between query and target sequences.

**Parameters:**
- `query_seq`: Amino acid sequence from PDB
- `query_id`: UniProt ID for query
- `target_seq`: DescribePROT reference sequence
- `target_id`: UniProt ID for target

**Returns:** List of hit dictionaries with alignment information

#### `parse_blast_alignment(hit, query_seq, subject_seq)`
Parses BLAST alignment to create position mapping dictionary.

**Example:**

```Text
Query:   M A V T G - - K R

Subject: M A V T G G G K R

Mapping: 1→1, 2→2, 3→3, 4→4, 5→5, 6→8 (shifted due to gaps)
```


#### `create_perfect_match_result(uniprot_id, target_seq, query_seq)`
Optimization for identical sequences - creates 1:1 position map without BLAST.

### Binding Site Mapping Functions

#### `map_binding_sites_with_alignment(binding_sites_str, alignment, subject_sequence=None, ...)`
Core mapping function that renumbers binding sites.

**Process:**
1. Parse original binding sites (e.g., "L84 V87" → [('L',84), ('V',87)])
2. Map each position using alignment position_map
3. Validate amino acid matches reference
4. Only keep sites that pass validation

**Returns:**
```python
{
    'mapped_sites': [('L',84), ('V',87)],
    'mapped_sites_str': "L84 V87",
    'all_mapped': True/False,
    'mapped_count': 2,
    'total_sites': 2,
    'removed_mismatches': 0,
    'removed_out_of_range': 0
}
```

## Validation Functions

### `count_binding_site_mismatches(df, binding_col, seq_col, label)`

Diagnostic function that validates binding sites against reference sequences. Prints statistics but does not filter.

### `filter_failed_alignments(df, dataset_name)`

Removes records where BLAST found no alignment.

**alignment_status values kept:**

- `all_sites_mapped` - All binding sites mapped successfully
- `partial_mapping` - Some sites mapped
- `aligned_no_sites` - Alignment succeeded but no binding sites to map
- `no_sites_mapped` - Alignment succeeded but no sites mapped successfully

**Removed:**

- `no_hit` - BLAST found no alignment
- `no_sequence` - No sequence available for alignment

## Preservation Functions

### `validate_and_update_binding_sites(df, dataset_name)`

Creates unified 'Binding_sites' column from renumbered sites for downstream compatibility.

### Column Preservation Logic

Affinity and metadata columns are preserved by:

1. **Identify columns to preserve** - Scans for patterns like `affinity`, `BindingMOAD`, `PDBbind`, `BindingDB`, `DrugBank`, `PDB_ID`, etc.

2. **Store as dictionary** - Converts all preserved columns for each row into a dictionary and saves in temporary `_preserved_data` column. Each row keeps its own dictionary with its specific affinity values.

3. **Survive grouping** - When rows are grouped by `(UniProt_ID, Sequence)` for BLAST, each dictionary stays attached to its original row.

4. **Restore after alignment** - After BLAST completes, extracts values from each row's dictionary back into individual columns.

5. **Clean up** - Drops the temporary `_preserved_data` column.

## Cache System

### Cache Files

| File | Description |
|------|-------------|
| `QBioLiP_alignment_cache.pkl` | Cached alignments for Q-BioLiP dataset |
| `BioLiP2_alignment_cache.pkl` | Cached alignments for BioLiP2 dataset |

```python
{
    ('UniProt_ID', 'Sequence'): {
        'best_overall_alignment': {
            'subject_id': 'P00720',
            'identity': 98.5,
            'evalue': 1.2e-45,
            'score': 1234,
            'position_map': {1: 1, 2: 2, 84: 84, 85: 86, ...}
        }
    }
}
```

## Input/Output Specifications

### Input DataFrame Requirements

**Required columns:**

- `UniProt_ID` - UniProt accession number
- Sequence column (one of): `Sequence`, `Receptor_sequence`, or `seq`
- Binding sites column (dataset-specific):
  - **QBioLiP**: `Binding_site_pdb` (preferred) or `Binding_site_original`
  - **BioLiP2**: `Binding_site_renumbered`

**Preserved columns:**

- Affinity: `Binding_affinity_MOAD`, `Binding_affinity_PDBbind`, `Binding_affinity_BindingDB`
- Metadata: `DrugBank`, `PDB_ID`, `Ligand_ID`, etc.

### Output DataFrame Columns

| Column | Description |
|--------|-------------|
| `describePROT_sequence` | Reference sequence from DescribePROT |
| `renumbered_binding_sites` | Binding sites in DescribePROT coordinates |
| `alignment_identity` | BLAST percent identity |
| `alignment_coverage` | Fraction of query aligned |
| `alignment_status` | Status code for filtering |
| Preserved columns | Original affinity and metadata |


## Issues with Q-BioLiP

```Text 
============================================================
  ALIGNMENT FOR QBioLiP
============================================================

  ----------------------------------------
  BEFORE ALIGNMENT
  ----------------------------------------
    Binding sites column 'Binding_site_pdb' has 328,765 non-null values
    Sample binding site: K1056 E1058 E1117
    Original validation: 445,537/2,056,753 matches (21.66%)

  Processing stats:
    From cache: 55,202
    Direct matches (no BLAST): 701
    Pairwise BLAST runs: 54,501

  ----------------------------------------
  AFTER ALIGNMENT
  ----------------------------------------
    After mapping: 431,682/431,682 matches (100.00%)

  Site removal statistics:
    Total original sites: 2,026,448
    Sites kept after mapping: 431,682 (21.30%)
    Removed due to amino acid mismatch: 1,234,785

  QBioLiP Alignment Results
  ==================================================
    Total rows processed: 328,765
    Successfully aligned: 325,203
      - All binding sites mapped: 55,905
      - Partial binding sites mapped: 61,684
```
Logic will need to be updated to correct Q-BioLiP parsing losing so much Binding AA information. 

## Error Handling

| Error | Cause | Solution |
|-------|-------|----------|
| No sequence column found | DataFrame missing `Sequence`/`Receptor_sequence`/`seq` | Check column names |
| No binding sites column found | Dataset name doesn't match columns | Verify `dataset_name` parameter |
| BLAST fails | BLAST+ not installed or not in PATH | Install BLAST+ and verify PATH |