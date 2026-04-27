#!/usr/bin/env python3
"""
Configuration settings for Q-BioLiP processing

This module centralizes all configuration parameters including:
- File paths for input/output files
- Processing parameters (batch sizes, caching)
- Column names to keep/drop
- Directory structures
"""

from pathlib import Path
import pandas as pd

class Config:
    """
    Configuration class for all processing parameters
    
    Uses class methods so configuration can be accessed without instantiation.
    All paths are relative to a base directory (default: current directory).
    """
    
    # =========================================================================
    # INPUT FILE CONFIGURATION
    # =========================================================================
    # Each entry contains:
    #   - filename: actual file name in the data directory
    #   - description: human-readable description of the file
    #   - cache: (optional) cache file name for faster subsequent runs
    
    INPUT_FILES = {
        'describeprot': {
            'filename': 'entire_database_AF.json',     # Changed from describePROT.csv
            'description': 'Main describePROT database file',
            'cache': 'filtered_describePROT.json'
        },
        'qbiolip': {
            'filename': 'Q-BioLiP_all.csv',
            'description': 'Q-BioLiP database with protein-ligand binding site information'
            # Note: No cache specified - Q-BioLiP is processed fresh each time
        },
        'biolip2': {
            'filename': 'BioLiP2.txt',
            'description': 'BioLiP2 database with pre-computed sequences',
            'cache': 'BioLiP2_filtered.csv'          # Filtered BioLiP2 data cache
        },
        'approved_ligands': {
            'filename': 'approved_ligands.xlsx',
            'description': 'Single Excel file with all approved ligands'
            # Note: Excel format (.xlsx) requires openpyxl
        }
    }
    
    # =========================================================================
    # PDB DIRECTORY CONFIGURATION
    # =========================================================================
    
    PDB_FOLDER = {
        'name': 'rec_pdb',
        'description': 'Directory containing PDB structure files'
    }
    
    # =========================================================================
    # COLUMN FILTERING
    # =========================================================================
    
    # Columns to keep from describePROT database
    # ACC: UniProt accession number (primary key for matching)
    # ACC_entry: Alternative accession
    # seq: Protein sequence (for alignment)
    DESCRIBEPROT_KEEP_COLUMNS = ["ACC", "seq"] 
    
    # Columns to drop from Q-BioLiP (affinity columns are NOT in this list)
    # Note: BindingMOAD, PDBbind-CN, BindingDB are preserved automatically
    QBIOLIP_DROP_COLUMNS = [
        "Stoichiometry",      
        "Resolution (Å)",     
        "Assembly Detail",   
        "Pubmed ID",          
        "Relevant"            
    ]
    
    # =========================================================================
    # LIGAND FILE CONFIGURATION
    # =========================================================================
    
    LIGAND_FILE_COLUMNS = {
        'ligand_id': 'Ligand ID',    # Column name for ligand identifiers
        'drugbank': 'DrugBank'       # Column name for DrugBank IDs
    }
    
    # =========================================================================
    # PROCESSING PARAMETERS
    # =========================================================================
    
    BATCH_SIZE = 10000               # Rows per batch for large file processing
    PDB_FOLDER_NAME = "rec_pdb"      # Directory name for PDB files
    PDB_CACHE_DIR = "pdb_cache"      # Directory for cached PDB sequences
    
    # Cache settings - speeds up subsequent runs significantly
    USE_CACHE = True                 # Enable/disable all caching
    
    # =========================================================================
    # OUTPUT FILE CONFIGURATION
    # =========================================================================
    # All output files are written to the base directory
    
    OUTPUT_FILES = {
        # Q-BioLiP outputs
        'qbiolip': 'qbiowithsequences.csv',           # Q-BioLiP with extracted sequences
        'qbiolip_aligned': 'qbiolip_aligned.csv',     # Q-BioLiP after alignment to describePROT
        
        # BioLiP2 outputs
        'biolip2': 'biolip2_filtered.csv',            # Filtered BioLiP2 data
        'biolip2_aligned': 'biolip2_aligned.csv',     # BioLiP2 after alignment
        
        # Final outputs (created by merge_biolip_datasets.py)
        'final_json': 'final.json',                   # Main output - combined binding sites
        'ligand_counts': 'ligand_binding_site_counts.xlsx',  # Per ligand binding AA counts 
    }
    
    @classmethod
    def get_input_path(cls, file_key, base_dir='.'):
        """
        Get full path for an input file.
        
        Args:
            file_key: Key in INPUT_FILES dictionary (e.g., 'qbiolip', 'biolip2')
            base_dir: Base directory containing the data files
        
        Returns:
            Path object pointing to the input file
        """
        return Path(base_dir) / cls.INPUT_FILES[file_key]['filename']
    
    @classmethod
    def get_cache_path(cls, file_key, base_dir='.'):
        """
        Get cache path for an input file.
        
        Args:
            file_key: Key in INPUT_FILES dictionary
            base_dir: Base directory containing cache files
        
        Returns:
            Path object for cache file, or None if no cache defined
        """
        if 'cache' in cls.INPUT_FILES[file_key]:
            return Path(base_dir) / cls.INPUT_FILES[file_key]['cache']
        return None
    
    @classmethod
    def get_output_path(cls, file_key, base_dir='.'):
        """
        Get output path for a result file.
        
        Args:
            file_key: Key in OUTPUT_FILES dictionary
            base_dir: Base directory for output files
        
        Returns:
            Path object for output file
        """
        return Path(base_dir) / cls.OUTPUT_FILES.get(file_key, f"{file_key}_output.csv")
    
    @classmethod
    def get_pdb_path(cls, base_dir='.'):
        """
        Get path to PDB folder.
        
        Args:
            base_dir: Base directory containing the rec_pdb folder
        
        Returns:
            Path object pointing to PDB directory
        """
        return Path(base_dir) / cls.PDB_FOLDER_NAME
    
    @classmethod
    def get_pdb_cache_dir(cls, base_dir='.'):
        """
        Get path to PDB cache directory.
        
        Args:
            base_dir: Base directory for cache
        
        Returns:
            Path object for PDB cache directory
        """
        return Path(base_dir) / cls.PDB_CACHE_DIR
    
    @classmethod
    def get_input_description(cls, file_key):
        """
        Get human-readable description for an input file.
        
        Args:
            file_key: Key in INPUT_FILES dictionary
        
        Returns:
            Description string
        """
        return cls.INPUT_FILES[file_key]['description']
    
    @classmethod
    def get_pdb_description(cls):
        """
        Get description for PDB folder.
        
        Returns:
            Description string for PDB folder
        """
        return cls.PDB_FOLDER['description']