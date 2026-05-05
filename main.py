#!/usr/bin/env python3
"""
Main orchestration script for the BioLiP processing pipeline.

This script coordinates all modules to:
    1. Load and filter DescribePROT reference data
    2. Process Q-BioLiP data (filter, flatten, extract sequences)
    3. Process BioLiP2 data (filter, standardize)
    4. Align sequences to DescribePROT and renumber binding sites
    5. Merge both datasets and generate final JSON output
    6. Create normalized tables

The pipeline is designed to handle large datasets efficiently with:
    - Streaming JSON parsing for DescribePROT
    - Caching of intermediate results
    - Batch processing for memory efficiency
    - Checkpoint saves for long-running BLAST operations

Usage:
    python main.py --data-dir /path/to/data --only-final --skip-db --no-qbio-aligned 
    --no-biolip2-aligned --no-merged-csv --no-ligand-counts --skip-alignment --no-cache 
    --ligand-json custom_ligand.json
"""

import sys
import argparse
import pandas as pd
from pathlib import Path

# Import all modules
from config import Config
from utils import setup_logging, check_input_files, timer
from filter_describeprot import filter_describeprot_data
from process_qbiolip import (
    load_and_filter_qbiolip, 
    filter_by_ligands, 
    flatten_binding_sites,
    standardize_qbiolip_columns,
    trim_string_columns,
    standardize_ligand_names
)
from process_biolip2 import (
    process_biolip2,
    standardize_biolip2_columns,
    trim_string_columns as trim_biolip2_strings,
    standardize_ligand_names as standardize_biolip2_ligands
)
from extract_sequences import extract_all_sequences
from merge_data import merge_sequences_with_qbiolip
from merge_biolip_datasets import merge_biolip_datasets
from align_and_renumber import align_and_renumber_dataset, validate_and_update_binding_sites
from create_database import create_database, load_data_to_database


def print_columns(df, step_name, dataset_name):
    """Print column names for a DataFrame at a specific step."""
    print(f"\n  [{dataset_name}] {step_name} - Columns ({len(df.columns)}):")
    for i, col in enumerate(df.columns):
        print(f"    {i+1}. {col}")
    
    # Print first row sample
    if len(df) > 0:
        print(f"\n  [{dataset_name}] {step_name} - Sample entry (first row):")
        first_row = df.iloc[0]
        for col in df.columns:
            # Truncate long values for readability
            value = str(first_row[col])
            if len(value) > 80:
                value = value[:77] + "..."
            print(f"    {col}: {value}")
        print("-" * 40)


def main():
    """
    Main function - orchestrates the entire pipeline.
    """
    
    # =========================================================================
    # COMMAND LINE ARGUMENT PARSING
    # =========================================================================
    parser = argparse.ArgumentParser(
        description='Process Q-BioLiP and BioLiP2 data: validate, align to describePROT, and merge',
        epilog='Place all required files in the data directory (see config.py for details)'
    )
    parser.add_argument('--data-dir', default='.', 
                       help='Directory containing input files (default: current directory)')
    parser.add_argument('--verbose', '-v', action='store_true', 
                       help='Print verbose output')
    parser.add_argument('--show-requirements', action='store_true',
                       help='Show required input files and exit')
    parser.add_argument('--skip-alignment', action='store_true',
                       help='Skip sequence alignment to describePROT (not recommended)')
    parser.add_argument('--no-cache', action='store_true',
                       help='Disable caching of alignment results')
    
    # =========================================================================
    # OUTPUT FILE OPTIONS (all optional, default True)
    # =========================================================================
    parser.add_argument('--no-qbio-aligned', action='store_true',
                       help='Do not save qbiolip_aligned.csv (intermediate)')
    parser.add_argument('--no-biolip2-aligned', action='store_true',
                       help='Do not save biolip2_aligned.csv (intermediate)')
    parser.add_argument('--no-merged-csv', action='store_true',
                       help='Do not save merged_biolip_datasets_by_pdb_separate.csv (intermediate)')
    parser.add_argument('--no-ligand-counts', action='store_true',
                       help='Do not save ligand_binding_site_counts.xlsx')
    parser.add_argument('--only-final', action='store_true',
                       help='Only save final.json and SQLite database (skip all intermediate files)')
    
    # =========================================================================
    # DATABASE OPTIONS
    # =========================================================================
    parser.add_argument('--skip-db', action='store_true',
                       help='Skip creating SQLite database')
    parser.add_argument('--ligand-json', default='ligand.json',
                       help='Path to ligand.json file with ligand names (optional)')
    
    args = parser.parse_args()
    
    # Handle --only-final flag (sets all other --no-* flags to True)
    if args.only_final:
        args.no_qbio_aligned = True
        args.no_biolip2_aligned = True
        args.no_merged_csv = True
        args.no_ligand_counts = True
    
    # Handle --show-requirements flag
    if args.show_requirements:
        from utils import print_file_requirements
        print_file_requirements(Config)
        return
    
    # Setup logging
    logger = setup_logging(args.verbose)
    
    # Print pipeline header
    print("=" * 60)
    print("Creating Final Dataset for DescribePROT from BioLiP databases")
    print("=" * 60)
    print(f"Data directory: {args.data_dir}")
    print(f"Skip alignment: {args.skip_alignment}")
    print(f"Use cache: {not args.no_cache}")
    print(f"Create database: {not args.skip_db}")
    print("-" * 40)
    print("Output options:")
    print(f"  - qbiolip_aligned.csv: {'NO' if args.no_qbio_aligned else 'YES'}")
    print(f"  - biolip2_aligned.csv: {'NO' if args.no_biolip2_aligned else 'YES'}")
    print(f"  - merged_biolip_datasets_by_pdb_separate.csv: {'NO' if args.no_merged_csv else 'YES'}")
    print(f"  - ligand_binding_site_counts.xlsx: {'NO' if args.no_ligand_counts else 'YES'}")
    print(f"  - final.json: YES (always)")
    print(f"  - biolip.db: {'YES' if not args.skip_db else 'NO'}")
    print("-" * 60)
    
   # =========================================================================
    # INPUT FILE VALIDATION (Skip original DescribePROT if cache exists)
    # =========================================================================
    print("\nChecking input files...")

    # Check for cached DescribePROT file first
    data_dir = Path(args.data_dir)
    describeprot_cache_file = data_dir / Config.INPUT_FILES['describeprot'].get('cache', 'filtered_describePROT.json')

    # If cache exists, we don't need the original file
    if describeprot_cache_file.exists():
        print(f"  Found cached DescribePROT file: {describeprot_cache_file}")
        print("  Will use cache instead of original entire_database_AF.json")
        print("  Original entire_database_AF.json is NOT required")
    else:
        print("  No cached DescribePROT file found. Will need original entire_database_AF.json")

    # Check for missing files
    missing_files = []

    # Check each required file
    for file_key, file_info in Config.INPUT_FILES.items():
        # Skip describeprot entirely if cache exists
        if file_key == 'describeprot' and describeprot_cache_file.exists():
            continue
        
        file_path = Config.get_input_path(file_key, args.data_dir)
        if not file_path.exists():
            missing_files.append({
                'path': str(file_path),
                'description': file_info.get('description', f'Required {file_key} file')
            })

    # Also check for PDB folder
    pdb_path = Config.get_pdb_path(args.data_dir)
    if not pdb_path.exists():
        missing_files.append({
            'path': str(pdb_path),
            'description': Config.get_pdb_description()
        })

    if missing_files:
        print("ERROR: Missing required files/folders:")
        print("-" * 60)
        for file_info in missing_files:
            print(f"\n  {file_info['path']}")
            print(f"     {file_info['description']}")
        print("\n" + "=" * 60)
        print(f"Please place all required files in: {Path(args.data_dir).absolute()}")
        sys.exit(1)

    print("All required files found!")
    
    try:
        # =====================================================================
        # STEP 1: LOAD DESCRIBEPROT DATA
        # =====================================================================
        print("\n" + "=" * 60)
        print("STEP 1: Loading DescribePROT data")
        print("=" * 60)

        # Check for cached file first
        describeprot_cache_file = Path(args.data_dir) / Config.INPUT_FILES['describeprot'].get('cache', 'filtered_describePROT.json')

        if describeprot_cache_file.exists():
            print(f"  Found cached filtered file: {describeprot_cache_file}")
            import json
            with open(describeprot_cache_file, 'r') as f:
                cached_data = json.load(f)
            describe_df = pd.DataFrame(cached_data)
            print(f"  Loaded {len(describe_df):,} records from cache")
        else:
            # No cache - need original file
            describeprot_file = Config.get_input_path('describeprot', args.data_dir)
            if not describeprot_file.exists():
                print(f"  ERROR: No DescribePROT file found!")
                print(f"  Please provide either:")
                print(f"    1. Cached file: {describeprot_cache_file}")
                print(f"    2. Original file: {describeprot_file}")
                sys.exit(1)
            
            describeprot_cache = Config.get_cache_path('describeprot', args.data_dir)
            describe_df = filter_describeprot_data(
                describeprot_file,
                keep_columns=Config.DESCRIBEPROT_KEEP_COLUMNS,
                batch_size=Config.BATCH_SIZE,
                cache_file=describeprot_cache
            )

        print_columns(describe_df, "After loading", "DescribePROT")
        
        # =====================================================================
        # STEP 2: LOAD AND PROCESS Q-BIOLIP DATA
        # =====================================================================
        print("\n" + "=" * 60)
        print("STEP 2: Loading and processing Q-BioLiP data")
        print("=" * 60)
        qbiolip_file = Config.get_input_path('qbiolip', args.data_dir)

        # 2a: Load and filter
        qbio_filtered = load_and_filter_qbiolip(qbiolip_file, describe_df)
        print_columns(qbio_filtered, "After load_and_filter_qbiolip", "Q-BioLiP")

        # 2b: Trim strings and standardize ligand names (BEFORE standardizing columns)
        qbio_filtered = trim_string_columns(qbio_filtered)
        qbio_filtered = standardize_ligand_names(qbio_filtered)
        print_columns(qbio_filtered, "After trim and ligand name standardization", "Q-BioLiP")

        # 2c: Standardize columns (renames and drops originals)
        qbio_filtered = standardize_qbiolip_columns(qbio_filtered)
        print_columns(qbio_filtered, "After standardize_qbiolip_columns", "Q-BioLiP")

        # 2d: Filter by approved ligands
        print("\n" + "-" * 40)
        print("Filtering by approved ligands...")
        approved_ligands_file = Config.get_input_path('approved_ligands', args.data_dir)

        qbio_with_ligands = filter_by_ligands(
            qbio_filtered, 
            approved_ligands_file,
            ligand_id_col='Ligand ID',
            drugbank_col='DrugBank'
        )
        print_columns(qbio_with_ligands, "After filter_by_ligands", "Q-BioLiP")

        # 2e: Flatten binding sites
        print("\n" + "-" * 40)
        print("Flattening binding sites...")
        qbio_flat = flatten_binding_sites(qbio_with_ligands)
        print_columns(qbio_flat, "After flatten_binding_sites", "Q-BioLiP")

        # 2f: Extract sequences from PDB files
        print("\n" + "-" * 40)
        print("Extracting sequences from PDB files...")
        pdb_folder = Config.get_pdb_path(args.data_dir)
        pdb_cache_dir = Config.get_pdb_cache_dir(args.data_dir)
        pdb_cache_dir.mkdir(exist_ok=True)

        easy_seq, diff_seq = extract_all_sequences(
            pdb_folder, 
            qbio_flat,
            cache_dir=pdb_cache_dir
        )

        print_columns(easy_seq, "Extracted sequences (single PDBs)", "PDB")
        print_columns(diff_seq, "Extracted sequences (bundles)", "PDB")

        # 2g: Merge sequences with Q-BioLiP binding data
        print("\n" + "-" * 40)
        print("Merging sequences with Q-BioLiP data...")
        qbio_final = merge_sequences_with_qbiolip(qbio_flat, easy_seq, diff_seq)
        qbio_final = standardize_qbiolip_columns(qbio_final)

        print_columns(qbio_final, "Final Q-BioLiP (after merge)", "Q-BioLiP")
        print(f"\n  Q-BioLiP final records: {len(qbio_final):,}")
        
        # =====================================================================
        # STEP 3: LOAD AND PROCESS BIOLIP2 DATA
        # =====================================================================
        print("\n" + "=" * 60)
        print("STEP 3: Loading and processing BioLiP2 data")
        print("=" * 60)

        biolip2_file = Config.get_input_path('biolip2', args.data_dir)
        biolip2_output = Config.get_output_path('biolip2', args.data_dir)

        biolip2_df = None
        if biolip2_file.exists():
            biolip2_df = process_biolip2(
                biolip2_file,
                describe_df,
                approved_ligands_file,
                output_file=biolip2_output
            )
            print_columns(biolip2_df, "After process_biolip2", "BioLiP2")
            
            if biolip2_df is not None and len(biolip2_df) > 0:
                biolip2_df = standardize_biolip2_columns(biolip2_df)
                print_columns(biolip2_df, "After standardize_biolip2_columns", "BioLiP2")
                
                biolip2_df = trim_biolip2_strings(biolip2_df)
                biolip2_df = standardize_biolip2_ligands(biolip2_df)
                print_columns(biolip2_df, "After trim and ligand standardization", "BioLiP2")
                
                print(f"\n  BioLiP2 final records: {len(biolip2_df):,}")
            else:
                print("  BioLiP2 processing returned no records")
                biolip2_df = None
        else:
            print(f"  BioLiP2 file not found: {biolip2_file}")
            print("  Skipping BioLiP2 processing...")
        
        # =====================================================================
        # STEP 4: ALIGN TO DESCRIBEPROT AND RENUMBER
        # =====================================================================
        if not args.skip_alignment:
            print("\n" + "=" * 60)
            print("STEP 4: Aligning sequences to describePROT and renumbering binding sites")
            print("=" * 60)
            print("This will:")
            print("  1. Align each sequence to describePROT reference")
            print("  2. Map binding sites to describePROT coordinates")
            print("  3. Use best hit for each binding site independently")
            print("  4. Drop only sequences with NO BLAST hits")
            print("-" * 60)
            
            full_describe_df = describe_df[['ACC', 'seq']].copy()
            print(f"\n  Loaded {len(full_describe_df):,} describePROT sequences")
            
            use_cache = not args.no_cache
            
            # Align Q-BioLiP dataset
            print("\n" + "-" * 40)
            print("Aligning Q-BioLiP dataset...")
            qbio_aligned = align_and_renumber_dataset(
                qbio_final,
                full_describe_df,
                dataset_name="QBioLiP",
                use_cache=use_cache
            )
            print_columns(qbio_aligned, "After alignment", "Q-BioLiP")
            
            qbio_aligned = validate_and_update_binding_sites(qbio_aligned, "Q-BioLiP")
            print_columns(qbio_aligned, "After validate_and_update_binding_sites", "Q-BioLiP")
            
            # Align BioLiP2 dataset (if available)
            biolip2_aligned = None
            if biolip2_df is not None and len(biolip2_df) > 0:
                print("\n" + "-" * 40)
                print("Aligning BioLiP2 dataset...")
                biolip2_aligned = align_and_renumber_dataset(
                    biolip2_df,
                    full_describe_df,
                    dataset_name="BioLiP2",
                    use_cache=use_cache
                )
                print_columns(biolip2_aligned, "After alignment", "BioLiP2")
                
                biolip2_aligned = validate_and_update_binding_sites(biolip2_aligned, "BioLiP2")
                print_columns(biolip2_aligned, "After validate_and_update_binding_sites", "BioLiP2")
            
            # Save intermediate aligned datasets (optional)
            if not args.no_qbio_aligned:
                qbio_aligned_output = Path(args.data_dir) / "qbiolip_aligned.csv"
                qbio_aligned.to_csv(qbio_aligned_output, index=False)
                print(f"\n  Q-BioLiP aligned data saved to: {qbio_aligned_output}")
            
            if biolip2_aligned is not None and not args.no_biolip2_aligned:
                biolip2_aligned_output = Path(args.data_dir) / "biolip2_aligned.csv"
                biolip2_aligned.to_csv(biolip2_aligned_output, index=False)
                print(f"  BioLiP2 aligned data saved to: {biolip2_aligned_output}")
            
            final_qbio = qbio_aligned
            final_biolip2 = biolip2_aligned
        else:
            print("\n" + "=" * 60)
            print("SKIPPING SEQUENCE ALIGNMENT (NOT RECOMMENDED)")
            print("=" * 60)
            final_qbio = qbio_final
            final_biolip2 = biolip2_df
        
        # =====================================================================
        # STEP 5: MERGE DATASETS AND GENERATE OUTPUTS
        # =====================================================================
        if final_biolip2 is not None and len(final_biolip2) > 0:
            print("\n" + "=" * 60)
            print("STEP 5: Merging Q-BioLiP and BioLiP2 datasets")
            print("=" * 60)
            
            merged_df = merge_biolip_datasets(
                final_qbio,
                final_biolip2,
                remove_duplicates=True,
                save_merged_csv=not args.no_merged_csv,
                save_ligand_counts=not args.no_ligand_counts
            )
            
            print(f"\n  Merged dataset size: {len(merged_df):,} records"
                  if len(merged_df) > 0 else "\n  Merged dataset is empty")
            
        else:
            print("\n" + "=" * 60)
            print("NO BIOLIP2 DATA TO MERGE")
            print("=" * 60)
            
            if not args.skip_alignment:
                final_output = Path(args.data_dir) / "qbiolip_aligned_final.csv"
            else:
                final_output = Path(args.data_dir) / "qbiolip_filtered.csv"
                
            final_qbio.to_csv(final_output, index=False)
            print(f"\n  Q-BioLiP results saved to: {final_output}")
            print(f"  Total records: {len(final_qbio):,}")
        
        # =====================================================================
        # STEP 6: CREATE SQLITE DATABASE
        # =====================================================================
        if not args.skip_db:
            print("\n" + "=" * 60)
            print("STEP 6: Creating SQLite database")
            print("=" * 60)
            
            json_file = Path("final.json")
            if not json_file.exists():
                print(f"  Warning: {json_file} not found. Skipping database creation.")
            else:
                print(f"  Creating database: biolip.db")
                conn = create_database("biolip.db")
                load_data_to_database(conn, str(json_file))
                conn.close()
                print(f"\n  Database saved to: biolip.db")
                
        # =====================================================================
        # FINAL SUMMARY
        # =====================================================================
        print("\n" + "=" * 60)
        print("PROCESSING COMPLETE")
        print("=" * 60)
        print("\nOutput files generated:")
        
        if not args.skip_alignment and not args.no_qbio_aligned:
            print(f"  - qbiolip_aligned.csv (Q-BioLiP with renumbered sites)")
        if not args.skip_alignment and final_biolip2 is not None and len(final_biolip2) > 0 and not args.no_biolip2_aligned:
            print(f"  - biolip2_aligned.csv (BioLiP2 with renumbered sites)")
        print(f"  - final.json (Final JSON with combined PDB IDs)")
        if not args.no_ligand_counts:
            print(f"  - ligand_binding_site_counts.xlsx (Ligand counts)")
        if not args.skip_db:
            print(f"  - biolip.db (SQLite database)")
        
        print("\n" + "=" * 60)
        
    except KeyboardInterrupt:
        print("\n\nUser interrupted. Exiting.")
        sys.exit(1)
    except Exception as e:
        logger.error(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()