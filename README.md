# BioLiP Binding Site Database Pipeline For DescribePROT

A comprehensive pipeline for processing, aligning, and merging binding site data from Q-BioLiP and BioLiP2 databases into a unified, searchable SQLite database.

## Overview

This pipeline processes millions of binding site records from two major protein-ligand interaction databases, aligns sequences to DescribePROT reference sequences, renumbers binding sites to standard coordinates, and creates a fully normalized database.

### Key Features

- **BLAST alignment** to DescribePROT reference sequences
- **PDB ID merging** - combines same binding sites from multiple PDBs (e.g., `7bgn_7bgm`)
- **Affinity data preservation** - keeps binding affinity values from multiple sources (MOAD, PDBbind, BindingDB)
- **Normalized SQLite database** - fully relational with foreign keys

## Pipeline Steps

1. Load and filter DescribePROT reference data
2. Process Q-BioLiP data (filter, flatten, extract sequences)
3. Process BioLiP2 data (filter, standardize)
4. Align sequences to DescribePROT and renumber binding sites via BLAST
5. Merge both datasets and generate final JSON output
6. Create normalized SQLite database

## File Descriptions

| File | Description |
|------|-------------|
| `main.py` | Main orchestration script - runs the entire pipeline |
| `config.py` | Configuration settings (file paths, batch sizes, column names) |
| `utils.py` | Utility functions (logging, file checking, timing decorator) |
| `filter_describeprot.py` | Streams and filters DescribePROT JSON file using ijson |
| `process_qbiolip.py` | Loads Q-BioLiP data, filters by UniProt/ligands, flattens binding sites by chain |
| `process_biolip2.py` | Loads BioLiP2 data, filters by UniProt/ligands, standardizes columns |
| `extract_sequences.py` | Extracts amino acid sequences from single PDB files and gzipped bundles for Q-BioLiP |
| `merge_data.py` | Merges extracted sequences with Q-BioLiP |
| `align_and_renumber.py` | Performs pairwise BLAST alignment to DescribePROT, maps binding site positions, validates amino acid matches |
| `merge_biolip_datasets.py` | Merges Q-BioLiP and BioLiP2 datasets, combines PDB IDs, creates final JSON |
| `create_database.py` | Creates normalized SQLite database from final.json with proper relationships |

## Input Files Required

Place these files in the data directory:

| File | Description | Source | Link |
|------|-------------|--------|------|
| `entire_database_AF.json` | DescribePROT sequences with UniProt IDs | DescribePROT database |
| `Q-BioLiP_all.csv` | Q-BioLiP binding site data | Q-BioLiP database |
| `BioLiP2.txt` | BioLiP2 binding site data (tab-separated) | BioLiP2 database |
| `approved_ligands.xlsx` | Approved ligands with DrugBank IDs | User created |
| `rec_pdb/` | Folder containing PDB files and bundles | PDB database |
| `filtered_describePROT.csv` | DescribePROT cache filtered by ["ACC", "seq"] | created by user |

Please note that while filtered_describePROT is not technically a required file to run this pipeline, and is instead a cache, the filter describePROT python file takes a very long time to run, more than overnight - for speed, I'd recommend using the cache until a more efficient way of removing descriptor data from the original describePROT database json is developed. 

### PDB Folder Structure

- **rec_pdb**
  - `{pdb_id}.pdb`
  - `{pdb_id}-pdb-bundle.tar.gz`

## Output Files

| File | Description |
|------|-------------|
| `final.json` | Final binding site records with combined PDB IDs |
| `ligand_binding_site_counts.xlsx` | Per-ligand unique binding site counts |
| `qbiolip_aligned.csv` | Q-BioLiP data with renumbered binding sites |
| `biolip2_aligned.csv` | BioLiP2 data with renumbered binding sites |
| `merged_biolip_datasets_by_pdb_separate.csv` | Full merged dataset |
| `biolip.db` | SQLite database with normalized tables |

### Database Tables

- **proteins**
  - `protein_id` (INTEGER, PRIMARY KEY)
  - `uniprot_id` (TEXT, UNIQUE, NOT NULL)

- **ligands**
  - `ligand_id` (INTEGER, PRIMARY KEY)
  - `ligand_code` (TEXT, UNIQUE, NOT NULL)

- **drugs**
  - `drug_id` (INTEGER, PRIMARY KEY)
  - `drugbank_id` (TEXT, UNIQUE, NOT NULL)
  - `ligand_id` (INTEGER, FOREIGN KEY → ligands.ligand_id)

- **pdb_structures**
  - `pdb_id` (INTEGER, PRIMARY KEY)
  - `pdb_code` (TEXT, UNIQUE, NOT NULL)

- **binding_sites**
  - `binding_site_id` (INTEGER, PRIMARY KEY)
  - `protein_id` (INTEGER, FOREIGN KEY → proteins.protein_id)
  - `ligand_id` (INTEGER, FOREIGN KEY → ligands.ligand_id)
  - `binding_site` (TEXT, NOT NULL)
  - `pdb_id` (INTEGER, FOREIGN KEY → pdb_structures.pdb_id)
  - `drugbank_id` (TEXT)
  - `binding_moad` (TEXT)
  - `binding_pdbbind` (TEXT)
  - `binding_bindingdb` (TEXT)

## Installation

### Prerequisites

- Python 3.8+
- BLAST+ installed and in PATH

### Install Dependencies

```bash
pip install pandas numpy biopython tqdm ijson openpyxl tabulate
```
### BLAST Installation

Download BLAST+ from: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/


## Usage

### Command Line Options

| Option | Description |
|--------|-------------|
| `--data-dir DIR` | Directory containing input files (default: current directory) |
| `--only-final` | Only save final.json and database (skip all intermediate files) |
| `--skip-db` | Skip creating database entirely |
| `--no-qbio-aligned` | Do not save qbiolip_aligned.csv |
| `--no-biolip2-aligned` | Do not save biolip2_aligned.csv |
| `--no-merged-csv` | Do not save merged_biolip_datasets_by_pdb_separate.csv |
| `--no-ligand-counts` | Do not save ligand_binding_site_counts.xlsx |
| `--skip-alignment` | Skip BLAST alignment step (use pre-aligned data) |
| `--no-cache` | Disable caching of alignment results |
| `--ligand-json` | Path to ligand.json file with ligand names (optional) |
| `--verbose` | Print verbose output |
| `--show-requirements` | Show required input files and exit |


### Examples

```Bash
# Full pipeline with default outputs
python main.py

# Only save final.json and database
python main.py --only-final

# Only final.json (no database)
python main.py --only-final --skip-db

# Skip intermediate files but keep CSV
python main.py --no-qbio-aligned --no-biolip2-aligned

# Skip alignment (use cached results)
python main.py --skip-alignment --no-cache

# Skip everything except final.json
python main.py --only-final --skip-db

# Custom data directory
python main.py --data-dir ./my_data

# Verbose output for debugging
python main.py --verbose

# Show required files without running
python main.py --show-requirements
```

## Cache Files

| Cache File | Purpose |
|------------|---------|
| `filtered_describePROT.csv` | Filtered DescribePROT sequences |
| `QBioLiP_alignment_cache.pkl` | Q-BioLiP BLAST results with alignments stored |
| `BioLiP2_alignment_cache.pkl` | BioLiP2 BLAST results with alignments stored |
| `pdb_cache/single_pdb_sequences.csv` | Extracted single PDB sequences |
| `pdb_cache/bundle_pdb_sequences.csv` | Extracted bundle PDB sequences |

# Future Improvements 

- filter_describeprot:

    - The code to filter the entire describeprot database takes a long time to clean the data to only include ("ACC", "seq"). It is encouraged that filtered_describePROT.csv  or filtered_describePROT.json is used until a more efficient method of removing predictor information is created. 

- align_and_renumber:

    - Blast results only use top hit (one with highest identity) - needs to also include other hits for binding AA's not found in top hit. A lot of Binding information from QBioLiP is lost in the process of aligning to describeprot sequences, while biolip2's alignment to describeprot is near perfect. This is because Biolip2 gived binding amino acid positions renumbered from 1, while qbiolips sequences need to be aligned to pdb sequences, then realigned after BLASTP. 

    - Blast alignment also slows down significantly when connected to the internet, at least for me. Disconnecting sped the process up. 

- process_qbiolip & process_biolip2

    - The column names are hardwired, and not input by user - any changes to feature headers in the future will mean features won't be found by the pipeline. Biolip2 had to have their feature names added manually to the dataset from annotation information. 

- create_database:

    - This function creates tables based solely on the output final.json. Therefore, it doesn't contain sequences as DescribePROT likely contains its own normalized table for sequences. The foreign key for sequences will need to be added in fitting these tables for input into DescribePROT.  

