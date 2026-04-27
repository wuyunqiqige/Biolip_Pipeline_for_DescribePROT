#!/usr/bin/env python3
"""
Shared utility functions

This module provides common utility functions used throughout the pipeline:
- Logging setup
- File validation
- Progress display
- Function timing decorator
"""

import sys
from pathlib import Path
from datetime import datetime
import logging


def setup_logging(verbose=False):
    """
    Setup logging configuration for the pipeline.
    
    Args:
        verbose: If True, set logging level to DEBUG for detailed output.
                 If False, set to INFO for only important messages.
    
    Returns:
        logging.Logger: Configured logger instance
    
    Usage:
        logger = setup_logging(args.verbose)
        logger.info("Processing started")
        logger.debug("Detailed debug info")  # Only appears if verbose=True
    """
    level = logging.DEBUG if verbose else logging.INFO
    
    # Configure basic logging to output to console (stderr via stdout)
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[logging.StreamHandler(sys.stdout)]
    )
    return logging.getLogger(__name__)


def check_input_files(data_dir='.', config=None):
    """
    Check if all required input files and folders exist.
    
    Args:
        data_dir: Base directory containing input files (default: current directory)
        config: Config class instance (loads Config if None)
    
    Returns:
        list: List of dictionaries with missing file information.
              Each dict contains: 'path', 'description', 'key'
              Returns empty list if all files found.
    
    Usage:
        missing = check_input_files('./data', Config)
        if missing:
            print("Missing files:")
            for f in missing:
                print(f"  {f['path']} - {f['description']}")
    
    Logic:
        1. Iterates through all keys in config.INPUT_FILES
        2. For each, constructs full path and checks existence
        3. Also checks for PDB folder existence
        4. Collects missing items in a list for reporting
    """
    if config is None:
        from config import Config
        config = Config
    
    data_path = Path(data_dir)
    missing_files = []
    
    # Step 1: Check each input file defined in config
    for file_key in config.INPUT_FILES:
        file_path = config.get_input_path(file_key, data_dir)
        if not file_path.exists():
            missing_files.append({
                'path': str(file_path),
                'description': config.get_input_description(file_key),
                'key': file_key
            })
    
    # Step 2: Check PDB folder existence
    pdb_path = config.get_pdb_path(data_dir)
    if not pdb_path.exists():
        missing_files.append({
            'path': str(pdb_path),
            'description': config.get_pdb_description(),
            'key': 'pdb_folder'
        })
    
    return missing_files


def print_file_requirements(config=None):
    """
    Print all required input files with their descriptions.
    Used by --show-requirements command line flag.
    
    Args:
        config: Config class instance (loads Config if None)
    
    Usage:
        python main.py --show-requirements
    
    Output:
        ======================================================================
        REQUIRED INPUT FILES
        ======================================================================
        
          entire_database_AF.json
             Main describePROT database file
        
          Q-BioLiP_all.csv
             Q-BioLiP database with protein-ligand binding site information
        ...
    """
    if config is None:
        from config import Config
        config = Config
    
    print("\n" + "=" * 70)
    print("REQUIRED INPUT FILES")
    print("=" * 70)
    
    # Print each input file with its description
    for key, info in config.INPUT_FILES.items():
        print(f"\n  {info['filename']}")
        print(f"     {info['description']}")
    
    # Print PDB folder requirement
    print(f"\n  {config.PDB_FOLDER_NAME}/")
    print(f"     {config.PDB_FOLDER['description']}")
    print("\n" + "=" * 70)


def print_progress(current, total, prefix='Progress:', suffix='Complete', bar_length=50):
    """
    Print a text-based progress bar to the console.
    
    Args:
        current: Current progress value (0 to total)
        total: Total value to reach (100% complete)
        prefix: Text to display before the progress bar
        suffix: Text to display after the progress bar
        bar_length: Length of the progress bar in characters
    
    Usage:
        for i in range(100):
            print_progress(i+1, 100, prefix='Processing:', suffix='Done')
            # do some work
    
    Output Example:
        Progress: [==========          ] 50% Complete
    
    Logic:
        1. Calculates completion percentage (current / total)
        2. Creates arrow string of '=' characters based on percentage
        3. Uses carriage return '\r' to overwrite same line
        4. Prints newline when complete
    """
    if total == 0:
        return
    
    percent = current / total
    # Calculate how many '=' characters to display
    arrow = '=' * int(round(percent * bar_length))
    # Fill remaining space with spaces
    spaces = ' ' * (bar_length - len(arrow))
    
    # \r returns cursor to start of line for overwriting
    sys.stdout.write(f'\r{prefix} [{arrow}{spaces}] {percent:.1%} {suffix}')
    sys.stdout.flush()  # Force immediate output
    
    # Print newline when complete
    if current >= total:
        print()


def timer(func):
    """
    Decorator that times function execution and prints elapsed time.
    
    Args:
        func: The function to be timed
    
    Returns:
        wrapper: Wrapped function that prints execution time
    
    Usage:
        @timer
        def my_function():
            # do something
            return result
    
    Output:
        my_function took 2.34 seconds
    
    Logic:
        1. Records start time before function call
        2. Executes the original function
        3. Records end time after function completes
        4. Calculates and prints elapsed time
        5. Returns the function's result unchanged
    
    Note:
        This decorator preserves the original function's return value.
        Time is printed to stdout but not returned.
    """
    def wrapper(*args, **kwargs):
        start = datetime.now()          # Start timer
        result = func(*args, **kwargs)  # Execute original function
        end = datetime.now()            # End timer
        elapsed = (end - start).total_seconds()  # Calculate elapsed
        print(f"  {func.__name__} took {elapsed:.2f} seconds")
        return result                    # Return original function result
    return wrapper