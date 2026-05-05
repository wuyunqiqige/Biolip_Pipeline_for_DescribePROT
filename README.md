# BioLiP Binding Site Database Pipeline For DescribePROT

A comprehensive pipeline for processing, aligning, and merging binding site data from Q-BioLiP and BioLiP2 databases into a unified, searchable SQLite database.

## Overview

This pipeline processes millions of binding site records from two major protein-ligand interaction databases, aligns sequences to DescribePROT reference sequences, renumbers binding sites to standard coordinates, and creates a fully normalized database.

### Key Features

- **BLAST alignment** to DescribePROT reference sequences
- **PDB ID merging** - combines same binding sites from multiple PDBs (e.g., `7bgn_7bgm`)
- **Affinity data preservation** - keeps binding affinity values from multiple sources (MOAD, PDBbind, BindingDB)
- **Normalized SQLite database** - fully relational with foreign keys
- **Overlap statistics** - automatic comparison between Q-BioLiP and BioLiP2 datasets

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
| `merge_biolip_datasets.py` | Merges Q-BioLiP and BioLiP2 datasets, combines PDB IDs, creates final JSON, prints overlap statistics |
| `create_database.py` | Creates normalized SQLite database from final.json with proper relationships |

## Input Files Required

Place these files in the data directory:

| File | Description | Source | Link |
|------|-------------|--------|------|
| `entire_database_AF.json`* | DescribePROT sequences with UniProt IDs | DescribePROT database | http://biomine.cs.vcu.edu/servers/DESCRIBEPROT/ |
| `Q-BioLiP_all.csv` | Q-BioLiP binding site data | Q-BioLiP database | https://yanglab.qd.sdu.edu.cn/Q-BioLiP/DATA/Q-BioLiP_all.csv |
| `BioLiP2.txt` | BioLiP2 binding site data (tab-separated and extracted) | BioLiP2 database | https://aideepmed.com/BioLiP/download/BioLiP.txt.gz |
| `approved_ligands.xlsx` | Approved ligands excel with DrugBank IDs | User created | In github |
| `rec_pdb/` | Folder containing PDB files and bundles | Q-BioLiP database | https://yanglab.qd.sdu.edu.cn/Q-BioLiP/Download/download_auto.html |
| `filtered_describePROT.csv`** | DescribePROT cache filtered by ["ACC", "seq"] | created by user | https://github.com/sagegint/Biolip_Pipeline_for_DescribePROT/blob/main/filtered_describePROT.csv |

*If filtered_describePROT already exists within the directory, you do not need to have entire_database_AF.json in the directory. 

**Please note that while filtered_describePROT is not technically a required file to run this pipeline, and is instead a cache, the filter describePROT python file takes a very long time to run, more than overnight - for speed, I'd recommend using the cache until a more efficient way of removing descriptor data from the original describePROT database json is developed. 

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

## Overlap Statistics

The pipeline automatically prints detailed overlap statistics between Q-BioLiP and BioLiP2 datasets, showing:

- Raw sites before deduplication
- Unique sites after deduplication
- Common binding sites between datasets
- Dataset-specific binding sites
- Overlap percentages
- Redundancy/duplication rates

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
| `overlap_statistics.pkl` | Cached overlap statistics between datasets |

# Future Improvements 

- filter_describeprot:

    - The code to filter the entire describeprot database takes a long time to clean the data to only include ("ACC", "seq"). It is encouraged that filtered_describePROT.csv  or filtered_describePROT.json is used until a more efficient method of removing predictor information is created. 

- align_and_renumber:

    - Blast results only use top hit (one with highest identity) - needs to also include other hits for binding AA's not found in top hit. A lot of Binding information from QBioLiP is lost in the process of aligning to describeprot sequences, while biolip2's alignment to describeprot is near perfect. This is because Biolip2 gived binding amino acid positions renumbered from 1, while qbiolips sequences need to be aligned to pdb sequences, then realigned after BLASTP. 

    - Blast alignment also slows down significantly when connected to the internet, at least for me. Disconnecting sped the process up. 

- process_qbiolip & process_biolip2

    - The column names are taken from the source databases, and not determined by user at start - any changes to feature headers in the future will mean features won't be found by the pipeline. Biolip2 had to have their feature names added manually to the dataset from annotation information. 

- create_database:

    - This function creates tables based solely on the output final.json. Therefore, it doesn't contain sequences as DescribePROT likely contains its own normalized table for sequences. The foreign key for sequences will need to be added in fitting these tables for input into DescribePROT.  


## Example Output 

```Text
PS C:\Users\sageg\Desktop\VCU\SR_design_project\QBioLiP\Create_Dataset\Biolip_Pipeline_for_DescribePROT> python3 .\main.py
C:\Users\sageg\.pyenv\pyenv-win\versions\3.13.2\Lib\site-packages\Bio\Application\__init__.py:39: BiopythonDeprecationWarning: The Bio.Application modules and modules relying on it have been deprecated.

Due to the on going maintenance burden of keeping command line application
wrappers up to date, we have decided to deprecate and eventually remove these
modules.

We instead now recommend building your command line and invoking it directly
with the subprocess module.
  warnings.warn(
============================================================
Creating Final Dataset for DescribePROT from BioLiP databases
============================================================
Data directory: .
Skip alignment: False
Use cache: True
Create database: True
----------------------------------------
Output options:
  - qbiolip_aligned.csv: YES
  - biolip2_aligned.csv: YES
  - merged_biolip_datasets_by_pdb_separate.csv: YES
  - ligand_binding_site_counts.xlsx: YES
  - final.json: YES (always)
  - biolip.db: YES
------------------------------------------------------------

Checking input files...
All required files found!

============================================================
STEP 1: Loading DescribePROT data
============================================================
  Found cached filtered file: filtered_describePROT.json
  Loading directly from cache (skipping JSON processing)...
  Loaded 2,276,602 records from cache

  [DescribePROT] After loading - Columns (3):
    1. ACC
    2. ACC_entry
    3. seq

  [DescribePROT] After loading - Sample entry (first row):
    ACC: A0A0A7EPL0
    ACC_entry: PIAL1_ARATH
    seq: MVIPATSRFGFRAEFNTKEFQASCISLANEIDAAIGRNEVPGNIQELALILNNVCRRKCDDYQTRAVVMALMISVKS...
----------------------------------------

============================================================
STEP 2: Loading and processing Q-BioLiP data
============================================================
Loading Q-BioLiP data: Q-BioLiP_all.csv
  After Uniprot filter: 999493 rows (removed 2223369)
  load_and_filter_qbiolip took 43.90 seconds

  [Q-BioLiP] After load_and_filter_qbiolip - Columns (12):
    1. Q-BioLiP ID
    2. PDB ID
    3. Assembly ID
    4. Uniprot ID
    5. Ligand ID
    6. Site
    7. Ligand Detail
    8. Binding Site
    9. Binding Site PDB
    10. BindingMOAD
    11. PDBbind-CN
    12. BindingDB

  [Q-BioLiP] After load_and_filter_qbiolip - Sample entry (first row):
    Q-BioLiP ID: BL4
    PDB ID: 102l
    Assembly ID: 102l_1
    Uniprot ID: P00720
    Ligand ID: CL
    Site: BS01
    Ligand Detail: 102l_1_CL_B
    Binding Site: A:T143 N145 R146
    Binding Site PDB: A:T142 N144 R145
    BindingMOAD: None
    PDBbind-CN: None
    BindingDB: None
----------------------------------------

  [Q-BioLiP] After trim and ligand name standardization - Columns (12):
    1. Q-BioLiP ID
    2. PDB ID
    3. Assembly ID
    4. Uniprot ID
    5. Ligand ID
    6. Site
    7. Ligand Detail
    8. Binding Site
    9. Binding Site PDB
    10. BindingMOAD
    11. PDBbind-CN
    12. BindingDB

  [Q-BioLiP] After trim and ligand name standardization - Sample entry (first row):
    Q-BioLiP ID: BL4
    PDB ID: 102l
    Assembly ID: 102l_1
    Uniprot ID: P00720
    Ligand ID: CL
    Site: BS01
    Ligand Detail: 102l_1_CL_B
    Binding Site: A:T143 N145 R146
    Binding Site PDB: A:T142 N144 R145
    BindingMOAD: None
    PDBbind-CN: None
    BindingDB: None
----------------------------------------

  [Q-BioLiP] After standardize_qbiolip_columns - Columns (11):
    1. Q-BioLiP ID
    2. Assembly ID
    3. Ligand Detail
    4. UniProt_ID
    5. Ligand_ID
    6. PDB_ID
    7. Binding_site_pdb
    8. Binding_site_number
    9. Binding_affinity_MOAD
    10. Binding_affinity_PDBbind
    11. Binding_affinity_BindingDB

  [Q-BioLiP] After standardize_qbiolip_columns - Sample entry (first row):
    Q-BioLiP ID: BL4
    Assembly ID: 102l_1
    Ligand Detail: 102l_1_CL_B
    UniProt_ID: P00720
    Ligand_ID: CL
    PDB_ID: 102l
    Binding_site_pdb: A:T142 N144 R145
    Binding_site_number: BS01
    Binding_affinity_MOAD: None
    Binding_affinity_PDBbind: None
    Binding_affinity_BindingDB: None
----------------------------------------

----------------------------------------
Filtering by approved ligands...
Filtering by ligands using: approved_ligands.xlsx
  After ligand filter: 262581 rows (removed 669479)
  filter_by_ligands took 1.28 seconds

  [Q-BioLiP] After filter_by_ligands - Columns (12):
    1. Q-BioLiP ID
    2. Assembly ID
    3. Ligand Detail
    4. UniProt_ID
    5. Ligand_ID
    6. PDB_ID
    7. Binding_site_pdb
    8. Binding_site_number
    9. Binding_affinity_MOAD
    10. Binding_affinity_PDBbind
    11. Binding_affinity_BindingDB
    12. DrugBank

  [Q-BioLiP] After filter_by_ligands - Sample entry (first row):
    Q-BioLiP ID: BL83
    Assembly ID: 117e_1
    Ligand Detail: 117e_1_MN_K
    UniProt_ID: P00817
    Ligand_ID: MN
    PDB_ID: 117e
    Binding_site_pdb: B:K1056 E1058 E1117
    Binding_site_number: BS01
    Binding_affinity_MOAD: None
    Binding_affinity_PDBbind: None
    Binding_affinity_BindingDB: None
    DrugBank: nan
----------------------------------------

----------------------------------------
Flattening binding sites...
Processing binding site data (keeping grouped)...
  Processed 262581 grouped rows
  flatten_binding_sites took 0.60 seconds

  [Q-BioLiP] After flatten_binding_sites - Columns (14):
    1. Q-BioLiP ID
    2. Assembly ID
    3. Ligand Detail
    4. UniProt_ID
    5. Ligand_ID
    6. PDB_ID
    7. Binding_site_pdb
    8. Binding_site_number
    9. Binding_affinity_MOAD
    10. Binding_affinity_PDBbind
    11. Binding_affinity_BindingDB
    12. DrugBank
    13. Chain
    14. PDB_Chain

  [Q-BioLiP] After flatten_binding_sites - Sample entry (first row):
    Q-BioLiP ID: BL83
    Assembly ID: 117e_1
    Ligand Detail: 117e_1_MN_K
    UniProt_ID: P00817
    Ligand_ID: MN
    PDB_ID: 117e
    Binding_site_pdb: K1056 E1058 E1117
    Binding_site_number: BS01
    Binding_affinity_MOAD: None
    Binding_affinity_PDBbind: None
    Binding_affinity_BindingDB: None
    DrugBank: nan
    Chain: B
    PDB_Chain: B
----------------------------------------

----------------------------------------
Extracting sequences from PDB files...
Extracting sequences from PDB files in: rec_pdb
  Found cached sequence data!
  Loading single PDB sequences from: pdb_cache\single_pdb_sequences.csv
  Loading bundle sequences from: pdb_cache\bundle_pdb_sequences.csv
  Loaded 145883 single and 4307 bundle sequences from cache
  extract_all_sequences took 0.37 seconds

  [PDB] Extracted sequences (single PDBs) - Columns (3):
    1. Assembly ID
    2. Chain
    3. Sequence

  [PDB] Extracted sequences (single PDBs) - Sample entry (first row):
    Assembly ID: 6ml4_1
    Chain: A
    Sequence: GSKSFTCDQCGKYFSQKRQLKSHYRVHTGHSLPECSHCHRKFMDVSQLKKHLRTHTGEKPFTCEICGKSFTAKSSLQ...
----------------------------------------

  [PDB] Extracted sequences (bundles) - Columns (3):
    1. Assembly ID
    2. Original Chain ID
    3. Sequence

  [PDB] Extracted sequences (bundles) - Sample entry (first row):
    Assembly ID: 3j9l_1
    Original Chain ID: A
    Sequence: YQYKDILSVFEDAFVDNFDCKDVQDMPKSILSKEEIDHIIMSKDAVSGTLRLFWTLLSKQEEMVQKFVEEVLRINYK...
----------------------------------------

----------------------------------------
Merging sequences with Q-BioLiP data...
Merging sequences with Q-BioLiP data...
  Removed 2 rows without sequences

  Available columns for binding site alignment:
    ['Q-BioLiP ID', 'Assembly ID', 'Ligand Detail', 'UniProt_ID', 'Ligand_ID', 'PDB_ID', 'Binding_site_pdb', 'Binding_site_number', 'Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB', 'DrugBank', 'PDB_Chain', 'Sequence']
  Using binding site column: 'Binding_site_pdb'
  Sample binding site before alignment: K1056 E1058 E1117

  Aligning binding sites to sequences...

  Alignment Statistics (by individual amino acid binding sites):
    Total binding sites kept: 2,013,816
    Total binding sites dropped: 42,959
    Percentage kept: 97.9%
    Rows with at least one aligned site: 245,299
    Rows with no aligned sites: 17,280

  Rows dropped due to no binding sites after alignment: 17,280
  Rows kept with binding sites: 245,299

  Sample binding site after alignment: K56 E58 E117
  merge_sequences_with_qbiolip took 16.32 seconds

  [Q-BioLiP] Final Q-BioLiP (after merge) - Columns (15):
    1. Q-BioLiP ID
    2. Assembly ID
    3. Ligand Detail
    4. UniProt_ID
    5. Ligand_ID
    6. PDB_ID
    7. Binding_site_pdb
    8. Binding_site_number
    9. Binding_affinity_MOAD
    10. Binding_affinity_PDBbind
    11. Binding_affinity_BindingDB
    12. DrugBank
    13. PDB_Chain
    14. Binding_site_pdb_original
    15. Receptor_sequence

  [Q-BioLiP] Final Q-BioLiP (after merge) - Sample entry (first row):
    Q-BioLiP ID: BL83
    Assembly ID: 117e_1
    Ligand Detail: 117e_1_MN_K
    UniProt_ID: P00817
    Ligand_ID: MN
    PDB_ID: 117e
    Binding_site_pdb: K56 E58 E117
    Binding_site_number: BS01
    Binding_affinity_MOAD: None
    Binding_affinity_PDBbind: None
    Binding_affinity_BindingDB: None
    DrugBank: nan
    PDB_Chain: B
    Binding_site_pdb_original: K1056 E1058 E1117
    Receptor_sequence: TYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGKL...
----------------------------------------

  Q-BioLiP final records: 245,299

============================================================
STEP 3: Loading and processing BioLiP2 data
============================================================
Loading BioLiP2 file: BioLiP2.txt
  Loaded 989,058 total records from BioLiP2
  Columns: PDB_ID, Receptor_chain, Resolution, Binding_site_number, Ligand_ID_CCD...
  load_biolip2_file took 13.61 seconds
Filtering BioLiP2 by Uniprot IDs...
  After Uniprot filter: 605,518 rows (removed 383,540)
  filter_biolip2_by_uniprot took 1.92 seconds
Filtering BioLiP2 by approved ligands...
  After ligand filter: 453,676 rows (removed 151,842)
  Rows with DrugBank info: 453,676 (100.0%)
  filter_biolip2_by_ligands took 0.85 seconds
  Final columns: PDB_ID, Binding_site_number, Ligand_ID_CCD, Binding_site_renumbered, Binding_affinity_MOAD, Binding_affinity_PDBbind, Binding_affinity_BindingDB, UniProt_ID, Receptor_sequence, DrugBank
  Final records: 453,676
  Saving to biolip2_filtered.csv...
  Saved 453,676 records to biolip2_filtered.csv
  process_biolip2 took 20.16 seconds

  [BioLiP2] After process_biolip2 - Columns (10):
    1. PDB_ID
    2. Binding_site_number
    3. Ligand_ID_CCD
    4. Binding_site_renumbered
    5. Binding_affinity_MOAD
    6. Binding_affinity_PDBbind
    7. Binding_affinity_BindingDB
    8. UniProt_ID
    9. Receptor_sequence
    10. DrugBank

  [BioLiP2] After process_biolip2 - Sample entry (first row):
    PDB_ID: 117e
    Binding_site_number: BS01
    Ligand_ID_CCD: MN
    Binding_site_renumbered: D115 E117 D120 D152
    Binding_affinity_MOAD: <NA>
    Binding_affinity_PDBbind: nan
    Binding_affinity_BindingDB: <NA>
    UniProt_ID: P00817
    Receptor_sequence: TYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGKL...
    DrugBank: nan
----------------------------------------

  [BioLiP2] After standardize_biolip2_columns - Columns (11):
    1. PDB_ID
    2. Binding_site_number
    3. Ligand_ID_CCD
    4. Binding_site_renumbered
    5. Binding_affinity_MOAD
    6. Binding_affinity_PDBbind
    7. Binding_affinity_BindingDB
    8. UniProt_ID
    9. Receptor_sequence
    10. DrugBank
    11. Binding_site_original

  [BioLiP2] After standardize_biolip2_columns - Sample entry (first row):
    PDB_ID: 117e
    Binding_site_number: BS01
    Ligand_ID_CCD: MN
    Binding_site_renumbered: D115 E117 D120 D152
    Binding_affinity_MOAD: <NA>
    Binding_affinity_PDBbind: nan
    Binding_affinity_BindingDB: <NA>
    UniProt_ID: P00817
    Receptor_sequence: TYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGKL...
    DrugBank: nan
    Binding_site_original: D115 E117 D120 D152
----------------------------------------

  [BioLiP2] After trim and ligand standardization - Columns (11):
    1. PDB_ID
    2. Binding_site_number
    3. Ligand_ID_CCD
    4. Binding_site_renumbered
    5. Binding_affinity_MOAD
    6. Binding_affinity_PDBbind
    7. Binding_affinity_BindingDB
    8. UniProt_ID
    9. Receptor_sequence
    10. DrugBank
    11. Binding_site_original

  [BioLiP2] After trim and ligand standardization - Sample entry (first row):
    PDB_ID: 117e
    Binding_site_number: BS01
    Ligand_ID_CCD: MN
    Binding_site_renumbered: D115 E117 D120 D152
    Binding_affinity_MOAD: <NA>
    Binding_affinity_PDBbind: nan
    Binding_affinity_BindingDB: <NA>
    UniProt_ID: P00817
    Receptor_sequence: TYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGKL...
    DrugBank: <NA>
    Binding_site_original: D115 E117 D120 D152
----------------------------------------

  BioLiP2 final records: 453,676

============================================================
STEP 4: Aligning sequences to describePROT and renumbering binding sites
============================================================
This will:
  1. Align each sequence to describePROT reference
  2. Map binding sites to describePROT coordinates
  3. Use best hit for each binding site independently
  4. Drop only sequences with NO BLAST hits
------------------------------------------------------------

  Loaded 2,276,602 describePROT sequences

----------------------------------------
Aligning Q-BioLiP dataset...

============================================================
  ALIGNMENT FOR QBioLiP
============================================================

  Initializing alignment for QBioLiP...
  Loaded 2,276,602 describePROT sequences
  Loaded 51257 cached alignments
  Using sequence column: 'Receptor_sequence'
  Using Q-BioLiP binding sites column: 'Binding_site_pdb'
  Preserving 8 additional columns: ['Assembly ID', 'Ligand Detail', 'Ligand_ID', 'PDB_ID', 'Binding_affinity_MOAD']...

  ----------------------------------------
  BEFORE ALIGNMENT (original binding sites against original sequences)
  ----------------------------------------
    Binding sites column 'Binding_site_pdb' has 245,299 non-null values
    Sample binding site: K56 E58 E117

  Finding unique (UniProt_ID, Sequence) pairs...
  Found 51,257 unique pairs out of 245,299 total rows
  Aligning unique pairs:   0%|                                 | 1/51257 [00:00<8:27:23,  1.68it/s]
  Checkpoint saved: 1000/51257 pairs
  Aligning unique pairs:   2%|▌                              | 1000/51257 [00:01<00:52, 949.81it/s]
  Checkpoint saved: 2000/51257 pairs
  Aligning unique pairs:   4%|█▏                            | 2000/51257 [00:01<00:41, 1185.12it/s]
  Checkpoint saved: 3000/51257 pairs
  Aligning unique pairs:   6%|█▊                            | 3000/51257 [00:02<00:39, 1235.66it/s]
  Checkpoint saved: 4000/51257 pairs
  Aligning unique pairs:   8%|██▎                           | 4000/51257 [00:03<00:37, 1265.21it/s]
  Checkpoint saved: 5000/51257 pairs
  Aligning unique pairs:  10%|██▉                           | 5000/51257 [00:04<00:35, 1286.72it/s]
  Checkpoint saved: 6000/51257 pairs
  Aligning unique pairs:  12%|███▌                          | 6000/51257 [00:04<00:33, 1334.57it/s]
  Checkpoint saved: 7000/51257 pairs
  Aligning unique pairs:  14%|████                          | 7000/51257 [00:05<00:33, 1316.29it/s]
  Checkpoint saved: 8000/51257 pairs
  Aligning unique pairs:  16%|████▋                         | 8000/51257 [00:06<00:32, 1320.05it/s]
  Checkpoint saved: 9000/51257 pairs
  Aligning unique pairs:  18%|█████▎                        | 9000/51257 [00:07<00:31, 1338.50it/s]
  Checkpoint saved: 10000/51257 pairs
  Aligning unique pairs:  20%|█████▋                       | 10000/51257 [00:07<00:31, 1316.77it/s]
  Checkpoint saved: 11000/51257 pairs
  Aligning unique pairs:  21%|██████▏                      | 11000/51257 [00:08<00:30, 1327.92it/s]
  Checkpoint saved: 12000/51257 pairs
  Aligning unique pairs:  23%|██████▊                      | 12000/51257 [00:09<00:29, 1313.52it/s]
  Checkpoint saved: 13000/51257 pairs
  Aligning unique pairs:  25%|███████▎                     | 13000/51257 [00:10<00:29, 1311.52it/s]
  Checkpoint saved: 14000/51257 pairs
  Aligning unique pairs:  27%|███████▉                     | 14000/51257 [00:10<00:27, 1342.81it/s]
  Checkpoint saved: 15000/51257 pairs
  Aligning unique pairs:  29%|████████▍                    | 15000/51257 [00:11<00:26, 1349.32it/s]
  Checkpoint saved: 16000/51257 pairs
  Aligning unique pairs:  31%|█████████                    | 16000/51257 [00:12<00:25, 1359.39it/s]
  Checkpoint saved: 17000/51257 pairs
  Aligning unique pairs:  33%|█████████▌                   | 17000/51257 [00:13<00:25, 1340.91it/s]
  Checkpoint saved: 18000/51257 pairs
  Aligning unique pairs:  35%|██████████▏                  | 18000/51257 [00:13<00:25, 1326.19it/s]
  Checkpoint saved: 19000/51257 pairs
  Aligning unique pairs:  37%|██████████▋                  | 19000/51257 [00:14<00:23, 1345.02it/s]
  Checkpoint saved: 20000/51257 pairs
  Aligning unique pairs:  39%|███████████▎                 | 20000/51257 [00:15<00:23, 1353.21it/s]
  Checkpoint saved: 21000/51257 pairs
  Aligning unique pairs:  41%|███████████▉                 | 21000/51257 [00:16<00:22, 1337.35it/s]
  Checkpoint saved: 22000/51257 pairs
  Aligning unique pairs:  43%|████████████▍                | 22000/51257 [00:16<00:22, 1326.66it/s]
  Checkpoint saved: 23000/51257 pairs
  Aligning unique pairs:  45%|█████████████                | 23000/51257 [00:17<00:21, 1324.14it/s]
  Checkpoint saved: 24000/51257 pairs
  Aligning unique pairs:  47%|█████████████▌               | 24000/51257 [00:18<00:20, 1327.09it/s]
  Checkpoint saved: 25000/51257 pairs
  Aligning unique pairs:  49%|██████████████▏              | 25000/51257 [00:19<00:20, 1294.01it/s]
  Checkpoint saved: 26000/51257 pairs
  Aligning unique pairs:  51%|██████████████▋              | 26000/51257 [00:20<00:19, 1277.81it/s]
  Checkpoint saved: 27000/51257 pairs
  Aligning unique pairs:  53%|███████████████▎             | 27000/51257 [00:20<00:18, 1300.87it/s]
  Checkpoint saved: 28000/51257 pairs
  Aligning unique pairs:  55%|███████████████▊             | 28000/51257 [00:21<00:17, 1311.14it/s]
  Checkpoint saved: 29000/51257 pairs
  Aligning unique pairs:  57%|████████████████▍            | 29000/51257 [00:22<00:16, 1315.00it/s]
  Checkpoint saved: 30000/51257 pairs
  Aligning unique pairs:  59%|████████████████▉            | 30000/51257 [00:23<00:16, 1310.26it/s]
  Checkpoint saved: 31000/51257 pairs
  Aligning unique pairs:  60%|█████████████████▌           | 31000/51257 [00:23<00:15, 1326.26it/s]
  Checkpoint saved: 32000/51257 pairs
  Aligning unique pairs:  62%|██████████████████           | 32000/51257 [00:24<00:14, 1324.48it/s]
  Checkpoint saved: 33000/51257 pairs
  Aligning unique pairs:  64%|██████████████████▋          | 33000/51257 [00:25<00:13, 1325.75it/s]
  Checkpoint saved: 34000/51257 pairs
  Aligning unique pairs:  66%|███████████████████▏         | 34000/51257 [00:26<00:12, 1329.73it/s]
  Checkpoint saved: 35000/51257 pairs
  Aligning unique pairs:  68%|███████████████████▊         | 35000/51257 [00:26<00:12, 1348.93it/s]
  Checkpoint saved: 36000/51257 pairs
  Aligning unique pairs:  70%|████████████████████▎        | 36000/51257 [00:27<00:11, 1298.81it/s]
  Checkpoint saved: 37000/51257 pairs
  Aligning unique pairs:  72%|████████████████████▉        | 37000/51257 [00:28<00:10, 1325.88it/s]
  Checkpoint saved: 38000/51257 pairs
  Aligning unique pairs:  74%|█████████████████████▍       | 38000/51257 [00:29<00:09, 1355.23it/s]
  Checkpoint saved: 39000/51257 pairs
  Aligning unique pairs:  76%|██████████████████████       | 39000/51257 [00:29<00:09, 1339.76it/s]
  Checkpoint saved: 40000/51257 pairs
  Aligning unique pairs:  78%|██████████████████████▋      | 40000/51257 [00:30<00:08, 1331.48it/s]
  Checkpoint saved: 41000/51257 pairs
  Aligning unique pairs:  80%|███████████████████████▏     | 41000/51257 [00:31<00:07, 1352.14it/s]
  Checkpoint saved: 42000/51257 pairs
  Aligning unique pairs:  82%|███████████████████████▊     | 42000/51257 [00:32<00:06, 1326.12it/s]
  Checkpoint saved: 43000/51257 pairs
  Aligning unique pairs:  84%|████████████████████████▎    | 43000/51257 [00:32<00:06, 1338.86it/s]
  Checkpoint saved: 44000/51257 pairs
  Aligning unique pairs:  86%|████████████████████████▉    | 44000/51257 [00:33<00:05, 1339.45it/s]
  Checkpoint saved: 45000/51257 pairs
  Aligning unique pairs:  88%|█████████████████████████▍   | 45000/51257 [00:34<00:04, 1351.37it/s]
  Checkpoint saved: 46000/51257 pairs
  Aligning unique pairs:  90%|██████████████████████████   | 46000/51257 [00:35<00:03, 1335.58it/s]
  Checkpoint saved: 47000/51257 pairs
  Aligning unique pairs:  92%|██████████████████████████▌  | 47000/51257 [00:35<00:03, 1380.48it/s]
  Checkpoint saved: 48000/51257 pairs
  Aligning unique pairs:  94%|███████████████████████████▏ | 48000/51257 [00:36<00:02, 1364.49it/s]
  Checkpoint saved: 49000/51257 pairs
  Aligning unique pairs:  96%|███████████████████████████▋ | 49000/51257 [00:37<00:01, 1377.14it/s]
  Checkpoint saved: 50000/51257 pairs
  Aligning unique pairs:  98%|████████████████████████████▎| 50000/51257 [00:37<00:00, 1343.33it/s]
  Checkpoint saved: 51000/51257 pairs
  Aligning unique pairs: 100%|█████████████████████████████| 51257/51257 [00:38<00:00, 1325.45it/s]

  Final cache saved: 51257 alignments to QBioLiP_alignment_cache.pkl

  Processing stats:
    From cache: 51,257
    Direct matches (no BLAST): 0
    Pairwise BLAST runs: 0
    No target sequence found: 0

  Applying results to 245,299 rows...
  Applying results: 100%|████████████████████████████████| 245299/245299 [00:38<00:00, 6340.82it/s]

  ----------------------------------------
  AFTER ALIGNMENT (mapped binding sites against describePROT)
  ----------------------------------------
    After mapping (before filtering): 1,958,293/1,958,293 matches (100.00%)
      Mismatches: 0 (0.00%)
      Out of range: 0 (0.00%)

  Site removal statistics during mapping:
    Total original sites: 1,991,818
    Sites kept after mapping: 1,958,293 (98.32%)
    Removed due to amino acid mismatch: 12,802
    Removed due to out of range: 0

  QBioLiP: Removed 2,079 records with NO BLAST hit, kept 243,220
    Kept breakdown: {'all_sites_mapped': 229303, 'partial_mapping': 11590, 'no_sites_mapped': 2327}

  Restoring 8 preserved columns...

  ----------------------------------------
  FINAL VALIDATION (after filtering rows with no hits)
  ----------------------------------------
    Final renumbered_binding_sites: 1,958,293/1,958,293 matches (100.00%)
      Mismatches: 0 (0.00%)
      Out of range: 0 (0.00%)

  ==================================================
  QBioLiP Alignment Results
  ==================================================
    Total rows processed: 245,299
    Successfully aligned: 243,220
      - All binding sites mapped: 229,303
      - Partial binding sites mapped: 11,590
      - No binding sites to map: 0
    Dropped columns: ['describePROT_ID', 'alignment_identity', 'alignment_coverage', 'alignment_score', 'alignment_evalue', 'alignment_status', 'binding_sites_mapped_count', 'binding_sites_total_count', 'binding_sites_all_mapped']

  [Q-BioLiP] After alignment - Columns (17):
    1. Q-BioLiP ID
    2. Assembly ID
    3. Ligand Detail
    4. UniProt_ID
    5. Ligand_ID
    6. PDB_ID
    7. Binding_site_pdb
    8. Binding_site_number
    9. Binding_affinity_MOAD
    10. Binding_affinity_PDBbind
    11. Binding_affinity_BindingDB
    12. DrugBank
    13. PDB_Chain
    14. Binding_site_pdb_original
    15. Receptor_sequence
    16. describePROT_sequence
    17. renumbered_binding_sites

  [Q-BioLiP] After alignment - Sample entry (first row):
    Q-BioLiP ID: BL83
    Assembly ID: 117e_1
    Ligand Detail: 117e_1_MN_K
    UniProt_ID: P00817
    Ligand_ID: MN
    PDB_ID: 117e
    Binding_site_pdb: K56 E58 E117
    Binding_site_number: BS01
    Binding_affinity_MOAD: None
    Binding_affinity_PDBbind: None
    Binding_affinity_BindingDB: None
    DrugBank: nan
    PDB_Chain: B
    Binding_site_pdb_original: K1056 E1058 E1117
    Receptor_sequence: TYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGKL...
    describePROT_sequence: MTYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGK...
    renumbered_binding_sites: K57 E59
----------------------------------------

  Q-BioLiP: Renamed 'renumbered_binding_sites' to 'Binding_sites'

  [Q-BioLiP] After validate_and_update_binding_sites - Columns (17):
    1. Q-BioLiP ID
    2. Assembly ID
    3. Ligand Detail
    4. UniProt_ID
    5. Ligand_ID
    6. PDB_ID
    7. Binding_site_pdb
    8. Binding_site_number
    9. Binding_affinity_MOAD
    10. Binding_affinity_PDBbind
    11. Binding_affinity_BindingDB
    12. DrugBank
    13. PDB_Chain
    14. Binding_site_pdb_original
    15. Receptor_sequence
    16. describePROT_sequence
    17. Binding_sites

  [Q-BioLiP] After validate_and_update_binding_sites - Sample entry (first row):
    Q-BioLiP ID: BL83
    Assembly ID: 117e_1
    Ligand Detail: 117e_1_MN_K
    UniProt_ID: P00817
    Ligand_ID: MN
    PDB_ID: 117e
    Binding_site_pdb: K56 E58 E117
    Binding_site_number: BS01
    Binding_affinity_MOAD: None
    Binding_affinity_PDBbind: None
    Binding_affinity_BindingDB: None
    DrugBank: nan
    PDB_Chain: B
    Binding_site_pdb_original: K1056 E1058 E1117
    Receptor_sequence: TYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGKL...
    describePROT_sequence: MTYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGK...
    Binding_sites: K57 E59
----------------------------------------

----------------------------------------
Aligning BioLiP2 dataset...

============================================================
  ALIGNMENT FOR BioLiP2
============================================================

  Initializing alignment for BioLiP2...
  Loaded 2,276,602 describePROT sequences
  Loaded 84945 cached alignments
  Using sequence column: 'Receptor_sequence'
  Using BioLiP2 binding sites column: 'Binding_site_renumbered'
  Preserving 6 additional columns: ['PDB_ID', 'Ligand_ID_CCD', 'Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB']...

  ----------------------------------------
  BEFORE ALIGNMENT (original binding sites against original sequences)
  ----------------------------------------
    Binding sites column 'Binding_site_renumbered' has 453,676 non-null values
    Sample binding site: D115 E117 D120 D152

  Finding unique (UniProt_ID, Sequence) pairs...
  Found 84,942 unique pairs out of 453,676 total rows
  Aligning unique pairs:   0%|                                | 1/84942 [00:00<17:38:46,  1.34it/s]
  Checkpoint saved: 1000/84942 pairs
  Aligning unique pairs:   1%|▎                              | 1000/84942 [00:01<02:04, 674.78it/s]
  Checkpoint saved: 2000/84942 pairs
  Aligning unique pairs:   3%|▊                              | 2367/84942 [00:03<01:33, 885.34it/s]
  Checkpoint saved: 3000/84942 pairs
  Aligning unique pairs:   4%|█                              | 3000/84942 [00:04<01:45, 773.77it/s]
  Checkpoint saved: 4000/84942 pairs
  Aligning unique pairs:   5%|█▍                             | 4000/84942 [00:05<01:36, 840.11it/s]
  Checkpoint saved: 5000/84942 pairs
  Aligning unique pairs:   6%|█▊                             | 5000/84942 [00:06<01:31, 878.41it/s]
  Checkpoint saved: 6000/84942 pairs
  Aligning unique pairs:   7%|██▏                            | 6000/84942 [00:07<01:30, 876.53it/s]
  Checkpoint saved: 7000/84942 pairs
  Aligning unique pairs:   8%|██▌                            | 7000/84942 [00:08<01:29, 870.66it/s]
  Checkpoint saved: 8000/84942 pairs
  Aligning unique pairs:   9%|██▉                            | 8000/84942 [00:09<01:28, 870.54it/s]
  Checkpoint saved: 9000/84942 pairs
  Aligning unique pairs:  11%|███▎                           | 9000/84942 [00:10<01:26, 874.95it/s]
  Checkpoint saved: 10000/84942 pairs
  Aligning unique pairs:  12%|███▌                          | 10000/84942 [00:11<01:25, 876.86it/s]
  Checkpoint saved: 11000/84942 pairs
  Aligning unique pairs:  13%|███▉                          | 11000/84942 [00:13<01:24, 872.29it/s]
  Checkpoint saved: 12000/84942 pairs
  Aligning unique pairs:  14%|████▏                         | 12000/84942 [00:14<01:23, 878.30it/s]
  Checkpoint saved: 13000/84942 pairs
  Aligning unique pairs:  15%|████▌                         | 13000/84942 [00:15<01:20, 889.54it/s]
  Checkpoint saved: 14000/84942 pairs
  Aligning unique pairs:  16%|████▉                         | 14000/84942 [00:16<01:19, 889.01it/s]
  Checkpoint saved: 15000/84942 pairs
  Aligning unique pairs:  18%|█████▎                        | 15000/84942 [00:17<01:19, 883.11it/s]
  Checkpoint saved: 16000/84942 pairs
  Aligning unique pairs:  19%|█████▋                        | 16000/84942 [00:18<01:17, 890.31it/s]
  Checkpoint saved: 17000/84942 pairs
  Aligning unique pairs:  20%|██████                        | 17000/84942 [00:19<01:16, 883.68it/s]
  Checkpoint saved: 18000/84942 pairs
  Aligning unique pairs:  21%|██████▎                       | 18000/84942 [00:20<01:15, 884.27it/s]
  Checkpoint saved: 19000/84942 pairs
  Aligning unique pairs:  22%|██████▋                       | 19000/84942 [00:22<01:15, 872.49it/s]
  Checkpoint saved: 20000/84942 pairs
  Aligning unique pairs:  24%|███████                       | 20000/84942 [00:23<01:14, 866.42it/s]
  Checkpoint saved: 21000/84942 pairs
  Aligning unique pairs:  25%|███████▍                      | 21000/84942 [00:24<01:13, 868.71it/s]
  Checkpoint saved: 22000/84942 pairs
  Aligning unique pairs:  26%|███████▊                      | 22000/84942 [00:25<01:12, 866.20it/s]
  Checkpoint saved: 23000/84942 pairs
  Aligning unique pairs:  27%|████████                      | 23000/84942 [00:26<01:10, 876.93it/s]
  Checkpoint saved: 24000/84942 pairs
  Aligning unique pairs:  28%|████████▍                     | 24000/84942 [00:27<01:10, 866.88it/s]
  Checkpoint saved: 25000/84942 pairs
  Aligning unique pairs:  29%|████████▊                     | 25000/84942 [00:29<01:09, 864.36it/s]
  Checkpoint saved: 26000/84942 pairs
  Aligning unique pairs:  31%|█████████▏                    | 26000/84942 [00:30<01:08, 862.34it/s]
  Checkpoint saved: 27000/84942 pairs
  Aligning unique pairs:  32%|█████████▌                    | 27000/84942 [00:31<01:06, 866.62it/s]
  Checkpoint saved: 28000/84942 pairs
  Aligning unique pairs:  33%|█████████▉                    | 28000/84942 [00:32<01:05, 866.52it/s]
  Checkpoint saved: 29000/84942 pairs
  Aligning unique pairs:  34%|██████████▏                   | 29000/84942 [00:33<01:03, 878.84it/s]
  Checkpoint saved: 30000/84942 pairs
  Aligning unique pairs:  35%|██████████▌                   | 30000/84942 [00:34<01:02, 877.82it/s]
  Checkpoint saved: 31000/84942 pairs
  Aligning unique pairs:  36%|██████████▉                   | 31000/84942 [00:35<00:59, 899.46it/s]
  Checkpoint saved: 32000/84942 pairs
  Aligning unique pairs:  38%|███████████▎                  | 32000/84942 [00:36<00:58, 904.09it/s]
  Checkpoint saved: 33000/84942 pairs
  Aligning unique pairs:  39%|███████████▋                  | 33000/84942 [00:38<00:57, 902.30it/s]
  Checkpoint saved: 34000/84942 pairs
  Aligning unique pairs:  40%|████████████                  | 34000/84942 [00:39<00:57, 888.18it/s]
  Checkpoint saved: 35000/84942 pairs
  Aligning unique pairs:  41%|████████████▎                 | 35000/84942 [00:40<00:56, 887.76it/s]
  Checkpoint saved: 36000/84942 pairs
  Aligning unique pairs:  42%|████████████▋                 | 36000/84942 [00:41<00:55, 875.50it/s]
  Checkpoint saved: 37000/84942 pairs
  Aligning unique pairs:  44%|█████████████                 | 37000/84942 [00:42<00:54, 877.69it/s]
  Checkpoint saved: 38000/84942 pairs
  Aligning unique pairs:  45%|█████████████▍                | 38000/84942 [00:43<00:53, 883.19it/s]
  Checkpoint saved: 39000/84942 pairs
  Aligning unique pairs:  46%|█████████████▊                | 39000/84942 [00:44<00:50, 901.67it/s]
  Checkpoint saved: 40000/84942 pairs
  Aligning unique pairs:  47%|██████████████▏               | 40000/84942 [00:45<00:49, 911.34it/s]
  Checkpoint saved: 41000/84942 pairs
  Aligning unique pairs:  49%|██████████████▌               | 41335/84942 [00:47<00:44, 972.25it/s]
  Checkpoint saved: 42000/84942 pairs
  Aligning unique pairs:  49%|██████████████▊               | 42000/84942 [00:48<00:48, 883.21it/s]
  Checkpoint saved: 43000/84942 pairs
  Aligning unique pairs:  51%|███████████████▏              | 43000/84942 [00:49<00:46, 897.07it/s]
  Checkpoint saved: 44000/84942 pairs
  Aligning unique pairs:  52%|███████████████▌              | 44000/84942 [00:50<00:46, 889.06it/s]
  Checkpoint saved: 45000/84942 pairs
  Aligning unique pairs:  53%|███████████████▉              | 45000/84942 [00:51<00:45, 882.92it/s]
  Checkpoint saved: 46000/84942 pairs
  Aligning unique pairs:  54%|████████████████▏             | 46000/84942 [00:52<00:44, 879.49it/s]
  Checkpoint saved: 47000/84942 pairs
  Aligning unique pairs:  55%|████████████████▌             | 47000/84942 [00:53<00:42, 882.95it/s]
  Checkpoint saved: 48000/84942 pairs
  Aligning unique pairs:  57%|████████████████▉             | 48000/84942 [00:54<00:41, 896.81it/s]
  Checkpoint saved: 49000/84942 pairs
  Aligning unique pairs:  58%|█████████████████▎            | 49000/84942 [00:55<00:39, 911.38it/s]
  Checkpoint saved: 50000/84942 pairs
  Aligning unique pairs:  59%|█████████████████▋            | 50000/84942 [00:57<00:38, 906.32it/s]
  Checkpoint saved: 51000/84942 pairs
  Aligning unique pairs:  60%|██████████████████            | 51000/84942 [00:58<00:38, 884.22it/s]
  Checkpoint saved: 52000/84942 pairs
  Aligning unique pairs:  61%|██████████████████▎           | 52000/84942 [00:59<00:37, 886.92it/s]
  Checkpoint saved: 53000/84942 pairs
  Aligning unique pairs:  62%|██████████████████▋           | 53000/84942 [01:00<00:36, 875.71it/s]
  Checkpoint saved: 54000/84942 pairs
  Aligning unique pairs:  64%|███████████████████           | 54000/84942 [01:01<00:35, 868.31it/s]
  Checkpoint saved: 55000/84942 pairs
  Aligning unique pairs:  65%|███████████████████▍          | 55000/84942 [01:02<00:34, 873.69it/s]
  Checkpoint saved: 56000/84942 pairs
  Aligning unique pairs:  66%|███████████████████▊          | 56000/84942 [01:03<00:33, 874.83it/s]
  Checkpoint saved: 57000/84942 pairs
  Aligning unique pairs:  67%|████████████████████▏         | 57000/84942 [01:05<00:31, 876.09it/s]
  Checkpoint saved: 58000/84942 pairs
  Aligning unique pairs:  68%|████████████████████▍         | 58000/84942 [01:06<00:31, 852.92it/s]
  Checkpoint saved: 59000/84942 pairs
  Aligning unique pairs:  69%|████████████████████▊         | 59000/84942 [01:07<00:30, 859.14it/s]
  Checkpoint saved: 60000/84942 pairs
  Aligning unique pairs:  71%|█████████████████████▏        | 60000/84942 [01:08<00:29, 841.28it/s]
  Checkpoint saved: 61000/84942 pairs
  Aligning unique pairs:  72%|█████████████████████▌        | 61000/84942 [01:09<00:28, 847.44it/s]
  Checkpoint saved: 62000/84942 pairs
  Aligning unique pairs:  73%|█████████████████████▉        | 62000/84942 [01:11<00:26, 850.58it/s]
  Checkpoint saved: 63000/84942 pairs
  Aligning unique pairs:  74%|██████████████████████▎       | 63000/84942 [01:12<00:25, 844.05it/s]
  Checkpoint saved: 64000/84942 pairs
  Aligning unique pairs:  75%|██████████████████████▌       | 64000/84942 [01:13<00:25, 832.25it/s]
  Checkpoint saved: 65000/84942 pairs
  Aligning unique pairs:  77%|██████████████████████▉       | 65000/84942 [01:14<00:23, 853.57it/s]
  Checkpoint saved: 66000/84942 pairs
  Aligning unique pairs:  78%|███████████████████████▎      | 66000/84942 [01:15<00:21, 865.00it/s]
  Checkpoint saved: 67000/84942 pairs
  Aligning unique pairs:  79%|███████████████████████▋      | 67000/84942 [01:16<00:20, 863.70it/s]
  Checkpoint saved: 68000/84942 pairs
  Aligning unique pairs:  80%|████████████████████████      | 68000/84942 [01:18<00:19, 869.70it/s]
  Checkpoint saved: 69000/84942 pairs
  Aligning unique pairs:  81%|████████████████████████▎     | 69000/84942 [01:19<00:18, 846.16it/s]
  Checkpoint saved: 70000/84942 pairs
  Aligning unique pairs:  82%|████████████████████████▋     | 70000/84942 [01:20<00:17, 848.85it/s]
  Checkpoint saved: 71000/84942 pairs
  Aligning unique pairs:  84%|█████████████████████████▏    | 71341/84942 [01:21<00:14, 909.26it/s]
  Checkpoint saved: 72000/84942 pairs
  Aligning unique pairs:  85%|█████████████████████████▍    | 72000/84942 [01:22<00:15, 830.73it/s]
  Checkpoint saved: 73000/84942 pairs
  Aligning unique pairs:  86%|█████████████████████████▊    | 73000/84942 [01:23<00:14, 841.74it/s]
  Checkpoint saved: 74000/84942 pairs
  Aligning unique pairs:  87%|██████████████████████████▏   | 74000/84942 [01:25<00:12, 861.97it/s]
  Checkpoint saved: 75000/84942 pairs
  Aligning unique pairs:  88%|██████████████████████████▍   | 75000/84942 [01:26<00:11, 867.63it/s]
  Checkpoint saved: 76000/84942 pairs
  Aligning unique pairs:  89%|██████████████████████████▊   | 76000/84942 [01:27<00:10, 845.94it/s]
  Checkpoint saved: 77000/84942 pairs
  Aligning unique pairs:  91%|███████████████████████████▏  | 77000/84942 [01:28<00:09, 860.22it/s]
  Checkpoint saved: 78000/84942 pairs
  Aligning unique pairs:  92%|███████████████████████████▌  | 78000/84942 [01:29<00:08, 851.67it/s]
  Checkpoint saved: 79000/84942 pairs
  Aligning unique pairs:  93%|███████████████████████████▉  | 79000/84942 [01:30<00:07, 845.89it/s]
  Checkpoint saved: 80000/84942 pairs
  Aligning unique pairs:  94%|████████████████████████████▎ | 80000/84942 [01:32<00:05, 845.76it/s]
  Checkpoint saved: 81000/84942 pairs
  Aligning unique pairs:  95%|████████████████████████████▌ | 81000/84942 [01:33<00:04, 855.36it/s]
  Checkpoint saved: 82000/84942 pairs
  Aligning unique pairs:  97%|████████████████████████████▉ | 82000/84942 [01:34<00:03, 857.12it/s]
  Checkpoint saved: 83000/84942 pairs
  Aligning unique pairs:  98%|█████████████████████████████▎| 83000/84942 [01:35<00:02, 850.31it/s]
  Checkpoint saved: 84000/84942 pairs
  Aligning unique pairs: 100%|██████████████████████████████| 84942/84942 [01:36<00:00, 877.53it/s]

  Final cache saved: 84945 alignments to BioLiP2_alignment_cache.pkl

  Processing stats:
    From cache: 84,942
    Direct matches (no BLAST): 0
    Pairwise BLAST runs: 0
    No target sequence found: 0

  Applying results to 453,676 rows...
  Applying results: 100%|████████████████████████████████| 453676/453676 [01:10<00:00, 6469.98it/s]

  ----------------------------------------
  AFTER ALIGNMENT (mapped binding sites against describePROT)
  ----------------------------------------
    After mapping (before filtering): 5,354,760/5,354,760 matches (100.00%)
      Mismatches: 0 (0.00%)
      Out of range: 0 (0.00%)

  Site removal statistics during mapping:
    Total original sites: 5,404,506
    Sites kept after mapping: 5,354,760 (99.08%)
    Removed due to amino acid mismatch: 12,091
    Removed due to out of range: 0

  BioLiP2: Removed 602 records with NO BLAST hit, kept 453,074
    Kept breakdown: {'all_sites_mapped': 436209, 'partial_mapping': 13102, 'no_sites_mapped': 3763}

  Restoring 6 preserved columns...

  ----------------------------------------
  FINAL VALIDATION (after filtering rows with no hits)
  ----------------------------------------
    Final renumbered_binding_sites: 5,354,760/5,354,760 matches (100.00%)
      Mismatches: 0 (0.00%)
      Out of range: 0 (0.00%)

  ==================================================
  BioLiP2 Alignment Results
  ==================================================
    Total rows processed: 453,676
    Successfully aligned: 453,074
      - All binding sites mapped: 436,209
      - Partial binding sites mapped: 13,102
      - No binding sites to map: 0
    Dropped columns: ['describePROT_ID', 'alignment_identity', 'alignment_coverage', 'alignment_score', 'alignment_evalue', 'alignment_status', 'binding_sites_mapped_count', 'binding_sites_total_count', 'binding_sites_all_mapped']

  [BioLiP2] After alignment - Columns (13):
    1. PDB_ID
    2. Binding_site_number
    3. Ligand_ID_CCD
    4. Binding_site_renumbered
    5. Binding_affinity_MOAD
    6. Binding_affinity_PDBbind
    7. Binding_affinity_BindingDB
    8. UniProt_ID
    9. Receptor_sequence
    10. DrugBank
    11. Binding_site_original
    12. describePROT_sequence
    13. renumbered_binding_sites

  [BioLiP2] After alignment - Sample entry (first row):
    PDB_ID: 117e
    Binding_site_number: BS01
    Ligand_ID_CCD: MN
    Binding_site_renumbered: D115 E117 D120 D152
    Binding_affinity_MOAD: <NA>
    Binding_affinity_PDBbind: nan
    Binding_affinity_BindingDB: <NA>
    UniProt_ID: P00817
    Receptor_sequence: TYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGKL...
    DrugBank: None
    Binding_site_original: D115 E117 D120 D152
    describePROT_sequence: MTYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGK...
    renumbered_binding_sites: D116 D121 D153
----------------------------------------

  BioLiP2: Renamed 'renumbered_binding_sites' to 'Binding_sites'

  [BioLiP2] After validate_and_update_binding_sites - Columns (13):
    1. PDB_ID
    2. Binding_site_number
    3. Ligand_ID_CCD
    4. Binding_site_renumbered
    5. Binding_affinity_MOAD
    6. Binding_affinity_PDBbind
    7. Binding_affinity_BindingDB
    8. UniProt_ID
    9. Receptor_sequence
    10. DrugBank
    11. Binding_site_original
    12. describePROT_sequence
    13. Binding_sites

  [BioLiP2] After validate_and_update_binding_sites - Sample entry (first row):
    PDB_ID: 117e
    Binding_site_number: BS01
    Ligand_ID_CCD: MN
    Binding_site_renumbered: D115 E117 D120 D152
    Binding_affinity_MOAD: <NA>
    Binding_affinity_PDBbind: nan
    Binding_affinity_BindingDB: <NA>
    UniProt_ID: P00817
    Receptor_sequence: TYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGKL...
    DrugBank: None
    Binding_site_original: D115 E117 D120 D152
    describePROT_sequence: MTYTTRQIGAKNTLEYKVYIEKDGKPVSAFHDIPLYADKENNIFNMVVEIPRWTNAKLEITKEETLNPIIQDTKKGK...
    Binding_sites: D116 D121 D153
----------------------------------------

  Q-BioLiP aligned data saved to: qbiolip_aligned.csv
  BioLiP2 aligned data saved to: biolip2_aligned.csv

============================================================
STEP 5: Merging Q-BioLiP and BioLiP2 datasets
============================================================

============================================================
MERGING DATASETS
============================================================

Q-BioLiP records: 243,220
BioLiP2 records: 453,074

  Q-BioLiP original columns: 17
  Q-BioLiP after removing duplicates: 17
  BioLiP2 original columns: 13
  BioLiP2 after removing duplicates: 13

  Q-BioLiP affinity columns: ['Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB']
  BioLiP2 affinity columns: ['Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB']

  Merging datasets on UniProt_ID, Ligand_ID, and PDB_ID...
  Removed 13,258 rows with no binding sites

  Merge complete: 1,860,072 rows

======================================================================
FINAL DATASET: GROUPED BY (UNIPROT_ID, LIGAND_ID, PDB_ID) - SEPARATE COLUMNS
======================================================================

  Grouping by: ['UniProt_ID', 'Ligand_ID', 'PDB_ID']
  Preserving affinity columns: ['Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB']

  Original rows: 1,860,072
  Unique combinations: 225,888

  ==================================================
  BINDING SITE STATISTICS (BY PDB (SEPARATE))
  ==================================================
  Q-BioLiP total sites: 796,278
  BioLiP2 total sites: 3,971,026
  Merged total sites: 4,218,094
  ==================================================
  Affinity columns in final output: ['Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB']

  Including affinity columns in exploded DataFrame: ['Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB']

======================================================================
CREATING FINAL JSON WITH COMBINED PDB IDs
======================================================================

  Calculating ligand binding site counts...

================================================================================
  QBioLiP vs BioLiP2 BINDING SITE OVERLAP
================================================================================

  OVERALL BINDING SITE STATISTICS
  ------------------------------------------------------------
  Q-BioLiP raw sites (before dedup):         796,278
  BioLiP2 raw sites (before dedup):          3,971,026
  ------------------------------------------------------------
  Q-BioLiP UNIQUE (UniProt+Ligand+Site):     290,177
  BioLiP2 UNIQUE (UniProt+Ligand+Site):      492,140
  ------------------------------------------------------------
  Common to both datasets (UNIQUE):          200,299
  Only in Q-BioLiP (UNIQUE):                 89,878
  Only in BioLiP2 (UNIQUE):                  291,841
  Total unique binding sites (union):        582,018

  OVERLAP ANALYSIS (Based on Unique Sites)
  ------------------------------------------------------------
  Shared binding sites (% of total unique):         34.4%
  Q-BioLiP unique sites also in BioLiP2:          69.0% (200,299/290,177)
  BioLiP2 unique sites also in Q-BioLiP:          40.7% (200,299/492,140)
  ------------------------------------------------------------
  Q-BioLiP redundancy removed: 63.6% (506,101 duplicates)
  BioLiP2 redundancy removed: 87.6% (3,478,886 duplicates)

  Including affinity columns in JSON: ['Binding_affinity_MOAD', 'Binding_affinity_PDBbind', 'Binding_affinity_BindingDB']

  Creating final JSON with combined PDB IDs...
  Saved JSON to: final.json (582,018 records)

======================================================================
OVERLAP STATISTICS
======================================================================
  Saved ligand counts to: ligand_binding_site_counts.xlsx

  Saved CSV to: merged_biolip_datasets_by_pdb_separate.csv
  merge_biolip_datasets took 465.06 seconds

  Merged dataset size: 225,888 records

============================================================
STEP 6: Creating SQLite database
============================================================
  Creating database: biolip.db
Loading final.json...
  Total binding site records: 582,018

Collecting unique entities...
  Scanning records: 100%|█████████████████████████████| 582018/582018 [00:00<00:00, 1046197.86it/s]

Inserting 20,324 proteins...
  Proteins: 100%|████████████████████████████████████████| 20324/20324 [00:00<00:00, 294801.29it/s]

Inserting 770 ligands...
  Ligands: 100%|█████████████████████████████████████████████| 770/770 [00:00<00:00, 210631.58it/s]

Inserting 125,408 PDB structures...
  PDB structures: 100%|█████████████████████████████████| 125408/125408 [00:03<00:00, 35460.37it/s]

Inserting 700 drugs...
  Drugs: 100%|███████████████████████████████████████████████| 700/700 [00:00<00:00, 200862.89it/s]

Inserting 582,018 binding sites...
  Binding sites: 100%|█████████████████████████████████| 582018/582018 [00:04<00:00, 135591.83it/s]

============================================================
DATABASE STATISTICS
============================================================
  Proteins: 20,324
  Ligands: 770
  Drugs: 700
  PDB structures: 125,408
  Binding sites: 582,018

  Database saved to: biolip.db

============================================================
PROCESSING COMPLETE
============================================================

Output files generated:
  - qbiolip_aligned.csv (Q-BioLiP with renumbered sites)
  - biolip2_aligned.csv (BioLiP2 with renumbered sites)
  - final.json (Final JSON with combined PDB IDs)
  - ligand_binding_site_counts.xlsx (Ligand counts)
  - biolip.db (SQLite database) 
```
