#!/usr/bin/env python3
"""
Part 1: Filter describePROT JSON file - True streaming with caching support

This module processes the large describePROT JSON file without loading the entire
file into memory. It uses streaming to read JSON objects one at a time, filters
to keep only necessary columns, and caches the result for faster subsequent runs.
"""

import ijson
import pandas as pd
import io
import json
from pathlib import Path


def filter_describeprot_data(input_file, keep_columns=None, batch_size=10000, cache_file=None):
    """
    Filter describePROT JSON file using true streaming.
    Never loads entire file into memory - processes record by record.
    Uses cached version if available for faster subsequent runs.
    
    Args:
        input_file: Path to the describePROT JSON file
        keep_columns: List of column names to keep (default: ["ACC", "ACC_entry", "seq"])
        batch_size: Number of records to accumulate before writing to CSV buffer
        cache_file: Optional custom cache file path
    
    Returns:
        pandas DataFrame with filtered describePROT data
    
    Logic Flow:
        1. Check for existing cache file (multiple possible locations)
        2. If cache exists, load it directly (fast path)
        3. If no cache, stream through JSON file and filter records
        4. Save filtered result to cache for future runs
    """
    
    if keep_columns is None:
        keep_columns = ["ACC", "ACC_entry", "seq"]
    
    input_path = Path(input_file)
    
    # =========================================================================
    # STEP 1: Define possible cache file locations to check
    # =========================================================================
    # Check multiple possible cache locations in priority order
    possible_caches = [
        cache_file,  # User-specified cache
        input_path.parent / "filtered_describePROT.json",  
        input_path.parent / "filtered_describePROT.csv",
        input_path.parent / f"{input_path.stem}_filtered.json",
        input_path.parent / f"{input_path.stem}_filtered.csv",
    ]
    
    # Remove None values from the list
    possible_caches = [c for c in possible_caches if c is not None]
    
    # =========================================================================
    # STEP 2: Check for existing cache file
    # =========================================================================
    existing_cache = None
    for cache_path in possible_caches:
        if cache_path.exists():
            existing_cache = cache_path
            break
    
    # =========================================================================
    # STEP 3: If cache found, load it directly (FAST PATH)
    # =========================================================================
    if existing_cache:
        print(f"  Found cached filtered file: {existing_cache}")
        print(f"  Loading directly from cache (skipping JSON processing)...")
        
        # Handle different file formats (CSV vs JSON)
        if existing_cache.suffix == '.csv':
            # CSV format - read in chunks to avoid memory issues
            chunks = []
            for chunk in pd.read_csv(existing_cache, chunksize=batch_size):
                chunks.append(chunk)
            describe_df = pd.concat(chunks, ignore_index=True)
            
        elif existing_cache.suffix == '.json':
            # JSON format - use ijson to stream read
            csv_buffer = io.StringIO()
            header_written = False
            batch = []
            processed = 0
            
            with open(existing_cache, 'rb') as f:
                # ijson.items streams JSON objects one at a time
                for record in ijson.items(f, 'item'):
                    batch.append(record)
                    processed += 1
                    
                    if processed % 1000 == 0:
                        print(f"  Loading cache: {processed:,} records...", end='\r')
                    
                    # Write batch when it reaches batch_size
                    if len(batch) >= batch_size:
                        df = pd.DataFrame(batch)
                        df.to_csv(csv_buffer, index=False, header=not header_written, mode="a")
                        header_written = True
                        batch = []
                
                # Write remaining records
                if batch:
                    df = pd.DataFrame(batch)
                    df.to_csv(csv_buffer, index=False, header=not header_written, mode="a")
            
            csv_buffer.seek(0)  # Reset buffer position to beginning
            describe_df = pd.read_csv(csv_buffer)
        else:
            # Unknown extension - try reading as CSV
            describe_df = pd.read_csv(existing_cache)
        
        # Ensure ACC column is uppercase for case-insensitive matching
        if 'ACC' in describe_df.columns:
            describe_df['ACC'] = describe_df['ACC'].astype(str).str.upper().str.strip()
        
        print(f"  Loaded {len(describe_df):,} records from cache")
        return describe_df
    
    # =========================================================================
    # STEP 4: NO CACHE FOUND - Process JSON file from scratch (SLOW PATH)
    # =========================================================================
    print(f"  No cache found. Processing JSON file: {input_file}")
    print(f"  Cache will be saved to: {possible_caches[0] if possible_caches else 'filtered_describePROT.json'}")
    print(f"  This may take a while...")
    
    # StringIO buffer to accumulate CSV data before creating DataFrame
    csv_buffer = io.StringIO()
    header_written = False
    batch = []          # Batch of records to write
    processed = 0       # Total records processed
    
    # =========================================================================
    # STEP 4a: Stream through JSON file manually (without ijson for large files)
    # =========================================================================
    # This manual parsing handles the specific JSON format of describePROT
    # It reads 1MB chunks and extracts individual JSON objects by counting braces
    
    with open(input_file, 'rb') as f:
        buffer = b''           # Buffer for current JSON object
        in_object = False      # Whether we're currently inside a JSON object
        brace_count = 0        # Count of nested braces (handles nested objects)
        
        while True:
            # Read 1MB chunks for memory efficiency
            chunk = f.read(1024 * 1024)
            if not chunk:
                break
            
            # Process each character in the chunk
            for byte in chunk:
                char = chr(byte)
                buffer += bytes([byte])
                
                if char == '{':
                    if not in_object:
                        # Start of a new JSON object
                        in_object = True
                        buffer = b'{'
                        brace_count = 1
                    else:
                        # Nested object opening brace
                        brace_count += 1
                        
                elif char == '}':
                    if in_object:
                        brace_count -= 1
                        
                        if brace_count == 0:
                            # End of JSON object - we have a complete record
                            try:
                                obj_str = buffer.decode('utf-8')
                                record = json.loads(obj_str)
                                
                                # Filter to keep only specified columns
                                filtered = {k: record.get(k, '') for k in keep_columns}
                                batch.append(filtered)
                                processed += 1
                                
                                if processed % 1000 == 0:
                                    print(f"  Processed {processed:,} records...", end='\r')
                                
                                # Write batch when it reaches batch_size
                                if len(batch) >= batch_size:
                                    df = pd.DataFrame(batch)
                                    # Append to CSV buffer (write header only once)
                                    df.to_csv(csv_buffer, index=False, header=not header_written, mode="a")
                                    header_written = True
                                    batch = []
                                    
                            except json.JSONDecodeError:
                                # Skip malformed JSON objects
                                pass
                            
                            # Reset for next object
                            in_object = False
                            buffer = b''
                            
                elif in_object:
                    # Continue building the current JSON object
                    pass
    
    # =========================================================================
    # STEP 4b: Write any remaining records in the final batch
    # =========================================================================
    if batch:
        df = pd.DataFrame(batch)
        df.to_csv(csv_buffer, index=False, header=not header_written, mode="a")
    
    print(f"\n  Processed {processed:,} total records")
    
    # Check if any records were found
    if processed == 0:
        print("  No records found. Check file format.")
        return pd.DataFrame()
    
    # =========================================================================
    # STEP 4c: Load filtered data into DataFrame
    # =========================================================================
    csv_buffer.seek(0)  # Reset buffer position to beginning
    describe_df = pd.read_csv(csv_buffer)
    
    # Convert ACC column to uppercase for case-insensitive matching with UniProt IDs
    if 'ACC' in describe_df.columns:
        describe_df['ACC'] = describe_df['ACC'].astype(str).str.upper().str.strip()
        print(f"  Converted ACC column to uppercase")
    
    # =========================================================================
    # STEP 4d: Save to cache for future runs
    # =========================================================================
    cache_path = possible_caches[0] if possible_caches else input_path.parent / "filtered_describePROT.csv"
    print(f"  Saving filtered data to cache: {cache_path}")
    describe_df.to_csv(cache_path, index=False)
    
    print(f"  Loaded {len(describe_df):,} records from describePROT")
    
    return describe_df