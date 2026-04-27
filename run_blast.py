#!/usr/bin/env python3
"""
Run BLAST for all sequences in merged dataset against describePROT sequences
Caches results so it only needs to run once
"""

import pandas as pd
import os
import subprocess
import tempfile
import time
import pickle
import re
from pathlib import Path
from tqdm import tqdm

def load_datasets(data_dir='.'):
    """Load the merged dataset and describePROT data"""
    
    # Load merged dataset
    merged_file = Path(data_dir) / "merged_biolip_datasets.csv"
    if not merged_file.exists():
        raise FileNotFoundError(f"Merged dataset not found: {merged_file}")
    
    print(f"Loading merged dataset: {merged_file}")
    merged_df = pd.read_csv(merged_file, low_memory=False)
    
    # Load describePROT data (use filtered version if available)
    describeprot_file = Path(data_dir) / "filtered_describePROT.json"
    if not describeprot_file.exists():
        describeprot_file = Path(data_dir) / "entire_database_AF.json"
    
    print(f"Loading describePROT data: {describeprot_file}")
    describeprot_df = pd.read_json(describeprot_file)
    
    return merged_df, describeprot_df

def combine_sequences(merged_df):
    """
    Create a single sequence column by combining sequences from both sources
    Prioritizes BioLiP2 sequence if available, otherwise uses Q-BioLiP sequence
    """
    df = merged_df.copy()
    
    # Create combined sequence column - BioLiP2 first
    if 'Receptor_sequence_BioLiP2' in df.columns and 'Receptor_sequence_QBioLiP' in df.columns:
        df['Receptor_sequence'] = df['Receptor_sequence_BioLiP2'].combine_first(
            df['Receptor_sequence_QBioLiP']
        )
    elif 'Receptor_sequence_BioLiP2' in df.columns:
        df['Receptor_sequence'] = df['Receptor_sequence_BioLiP2']
    elif 'Receptor_sequence_QBioLiP' in df.columns:
        df['Receptor_sequence'] = df['Receptor_sequence_QBioLiP']
    else:
        print("ERROR: No sequence column found!")
        return df
    
    return df

def prepare_sequences(df, id_col, seq_col):
    """Prepare sequences for BLAST by removing rows with missing data"""
    # Check if columns exist
    if id_col not in df.columns:
        print(f"  Warning: Column '{id_col}' not found.")
        return pd.DataFrame()
    if seq_col not in df.columns:
        print(f"  Warning: Column '{seq_col}' not found.")
        return pd.DataFrame()
    
    df_clean = df.dropna(subset=[id_col, seq_col]).copy()
    df_clean[id_col] = df_clean[id_col].astype(str).str.strip()
    df_clean[seq_col] = df_clean[seq_col].astype(str).str.strip()
    
    # Filter out empty sequences
    df_clean = df_clean[df_clean[seq_col] != '']
    df_clean = df_clean[df_clean[seq_col] != 'nan']
    
    return df_clean

def run_blast_for_pair(query_id, query_seq, subject_id, subject_seq, output_dir):
    """Run BLAST for a single pair and return the output file path"""
    
    # Create temporary FASTA files
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as query_f:
        query_f.write(f">{query_id}\n{query_seq}\n")
        query_fasta = query_f.name
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as subject_f:
        subject_f.write(f">{subject_id}\n{subject_seq}\n")
        subject_fasta = subject_f.name
    
    # Output file path
    safe_id = str(query_id).replace('/', '_').replace('|', '_').replace(' ', '_')
    blast_output = os.path.join(output_dir, f"{safe_id}_blast.txt")
    
    # Run BLAST
    cmd = [
        'blastp',
        '-query', query_fasta,
        '-subject', subject_fasta,
        '-out', blast_output,
        '-outfmt', '0'
    ]
    
    success = False
    try:
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        
        if result.returncode == 0 and os.path.exists(blast_output) and os.path.getsize(blast_output) > 0:
            success = True
        else:
            print(f"  BLAST failed for {query_id}: {result.stderr[:100] if result.stderr else 'Unknown error'}")
            
    except subprocess.TimeoutExpired:
        print(f"  BLAST timeout for {query_id}")
    except Exception as e:
        print(f"  Error processing {query_id}: {str(e)[:50]}")
    
    # Clean up temp files
    try:
        os.unlink(query_fasta)
        os.unlink(subject_fasta)
    except:
        pass
    
    return success, blast_output

def run_blast_with_cache(merged_df, describeprot_df, output_dir="blast_results", cache_file="blast_cache.pkl"):
    """
    Run BLAST for all unique UniProt IDs in the merged dataset against describePROT
    Uses caching to avoid re-running BLAST for already processed IDs
    """
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Combine sequences from both sources
    print("\nCombining sequences from Q-BioLiP and BioLiP2...")
    merged_df = combine_sequences(merged_df)
    
    # Prepare sequences
    print("\nPreparing sequences...")
    
    # For merged dataset, use UniProt_ID and combined Receptor_sequence
    merged_clean = prepare_sequences(merged_df, 'UniProt_ID', 'Receptor_sequence')
    
    if len(merged_clean) == 0:
        print("No valid sequences found in merged dataset!")
        return
    
    print(f"  Merged dataset: {len(merged_clean):,} rows with valid sequences")
    print(f"  Unique UniProt IDs in merged dataset: {merged_clean['UniProt_ID'].nunique():,}")
    
    # For describePROT, use ACC column
    describeprot_clean = prepare_sequences(describeprot_df, 'ACC', 'seq')
    print(f"  DescribePROT: {len(describeprot_clean):,} rows with valid sequences")
    
    # Get unique IDs from merged dataset
    unique_ids = merged_clean['UniProt_ID'].unique().tolist()
    print(f"\nUnique UniProt IDs to process: {len(unique_ids):,}")
    
    # Load cache if exists
    cache_path = Path(output_dir) / cache_file
    processed_ids = set()
    failed_ids = set()
    
    if cache_path.exists():
        print(f"\nLoading cache from {cache_path}")
        try:
            with open(cache_path, 'rb') as f:
                cache_data = pickle.load(f)
                processed_ids = set(cache_data.get('processed_ids', []))
                failed_ids = set(cache_data.get('failed_ids', []))
            print(f"  Previously processed: {len(processed_ids):,} IDs")
            print(f"  Previously failed: {len(failed_ids):,} IDs")
        except Exception as e:
            print(f"  Error loading cache: {e}")
    
    # Filter IDs that need processing
    ids_to_process = [uid for uid in unique_ids if uid not in processed_ids and uid not in failed_ids]
    print(f"\nIDs to process: {len(ids_to_process):,}")
    
    if len(ids_to_process) == 0:
        print("All IDs already processed!")
        return
    
    # Create mapping from UniProt ID to sequence
    id_to_seq = dict(zip(merged_clean['UniProt_ID'], merged_clean['Receptor_sequence']))
    
    # Create mapping for describePROT
    dp_id_to_seq = dict(zip(describeprot_clean['ACC'], describeprot_clean['seq']))
    
    # Run BLAST for each ID
    print("\nRunning BLAST...")
    start_time = time.time()
    success_count = 0
    new_processed = []
    new_failed = []
    
    for idx, acc_id in enumerate(tqdm(ids_to_process, desc="Processing IDs"), 1):
        # Get sequences
        final_seq = id_to_seq.get(acc_id, '')
        dp_seq = dp_id_to_seq.get(acc_id, '')
        
        if not final_seq or not dp_seq:
            tqdm.write(f"  Missing sequence for {acc_id}")
            new_failed.append(acc_id)
            continue
        
        # Run BLAST
        success, blast_file = run_blast_for_pair(
            acc_id, final_seq, 
            acc_id, dp_seq, 
            output_dir
        )
        
        if success:
            success_count += 1
            new_processed.append(acc_id)
        else:
            new_failed.append(acc_id)
        
        # Save cache periodically
        if idx % 100 == 0:
            # Update cache
            all_processed = processed_ids.union(set(new_processed))
            all_failed = failed_ids.union(set(new_failed))
            cache_data = {
                'processed_ids': list(all_processed),
                'failed_ids': list(all_failed),
                'last_update': time.time()
            }
            with open(cache_path, 'wb') as f:
                pickle.dump(cache_data, f)
            
            elapsed = time.time() - start_time
            tqdm.write(f"\n  Progress: {idx}/{len(ids_to_process)}")
            tqdm.write(f"  Success: {success_count}, Failed: {len(new_failed)}")
            tqdm.write(f"  Elapsed: {elapsed/60:.1f} minutes")
    
    # Final cache update
    all_processed = processed_ids.union(set(new_processed))
    all_failed = failed_ids.union(set(new_failed))
    cache_data = {
        'processed_ids': list(all_processed),
        'failed_ids': list(all_failed),
        'last_update': time.time()
    }
    with open(cache_path, 'wb') as f:
        pickle.dump(cache_data, f)
    
    # Final summary
    total_time = time.time() - start_time
    print("\n" + "="*60)
    print("BLAST RUN COMPLETE")
    print("="*60)
    print(f"Total unique IDs in merged dataset: {len(unique_ids):,}")
    print(f"Successfully processed (this run): {success_count}")
    print(f"Failed (this run): {len(new_failed)}")
    print(f"Total processed (all runs): {len(all_processed):,}")
    print(f"Total failed (all runs): {len(all_failed):,}")
    print(f"Total time this run: {total_time/60:.1f} minutes")
    print(f"BLAST output files saved to: {output_dir}")
    print(f"Cache saved to: {cache_path}")
    print("="*60)

if __name__ == "__main__":
    print("run_blast.py loaded successfully")