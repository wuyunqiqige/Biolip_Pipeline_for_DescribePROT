#!/usr/bin/env python3
"""
Create normalized SQLite database from final.json only
"""

import sqlite3
import json
from pathlib import Path
import argparse
from tqdm import tqdm


def create_database(db_path="biolip.db"):
    """Create normalized SQLite database schema"""
    
    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()
    
    # Drop existing tables
    cursor.execute("DROP TABLE IF EXISTS binding_sites")
    cursor.execute("DROP TABLE IF EXISTS proteins")
    cursor.execute("DROP TABLE IF EXISTS ligands")
    cursor.execute("DROP TABLE IF EXISTS drugs")
    cursor.execute("DROP TABLE IF EXISTS pdb_structures")
    
    # Create proteins table
    cursor.execute("""
        CREATE TABLE proteins (
            protein_id INTEGER PRIMARY KEY AUTOINCREMENT,
            uniprot_id TEXT UNIQUE NOT NULL
        )
    """)
    
    # Create ligands table
    cursor.execute("""
        CREATE TABLE ligands (
            ligand_id INTEGER PRIMARY KEY AUTOINCREMENT,
            ligand_code TEXT UNIQUE NOT NULL
        )
    """)
    
    # Create drugs table
    cursor.execute("""
        CREATE TABLE drugs (
            drug_id INTEGER PRIMARY KEY AUTOINCREMENT,
            drugbank_id TEXT UNIQUE,
            ligand_id INTEGER NOT NULL,
            FOREIGN KEY (ligand_id) REFERENCES ligands(ligand_id)
        )
    """)
    
    # Create PDB structures table
    cursor.execute("""
        CREATE TABLE pdb_structures (
            pdb_id INTEGER PRIMARY KEY AUTOINCREMENT,
            pdb_code TEXT UNIQUE NOT NULL
        )
    """)
    
    # Create binding sites table
    cursor.execute("""
        CREATE TABLE binding_sites (
            binding_site_id INTEGER PRIMARY KEY AUTOINCREMENT,
            protein_id INTEGER NOT NULL,
            ligand_id INTEGER NOT NULL,
            binding_site TEXT NOT NULL,
            pdb_id INTEGER NOT NULL,
            drugbank_id TEXT,
            binding_moad TEXT,
            binding_pdbbind TEXT,
            binding_bindingdb TEXT,
            FOREIGN KEY (protein_id) REFERENCES proteins(protein_id),
            FOREIGN KEY (ligand_id) REFERENCES ligands(ligand_id),
            FOREIGN KEY (pdb_id) REFERENCES pdb_structures(pdb_id)
        )
    """)
    
    # Create indexes
    cursor.execute("CREATE INDEX idx_proteins_uniprot ON proteins(uniprot_id)")
    cursor.execute("CREATE INDEX idx_ligands_code ON ligands(ligand_code)")
    cursor.execute("CREATE INDEX idx_drugs_drugbank ON drugs(drugbank_id)")
    cursor.execute("CREATE INDEX idx_pdb_code ON pdb_structures(pdb_code)")
    cursor.execute("CREATE INDEX idx_binding_protein ON binding_sites(protein_id)")
    cursor.execute("CREATE INDEX idx_binding_ligand ON binding_sites(ligand_id)")
    cursor.execute("CREATE INDEX idx_binding_pdb ON binding_sites(pdb_id)")
    
    conn.commit()
    return conn


def load_data_to_database(conn, final_json):
    """Load data from final.json into database"""
    
    cursor = conn.cursor()
    
    print(f"Loading {final_json}...")
    with open(final_json, 'r') as f:
        data = json.load(f)
    
    print(f"  Total binding site records: {len(data):,}")
    
    # Collect unique entities from JSON
    print("\nCollecting unique entities...")
    proteins = set()
    ligands = set()
    drugs = {}  # drugbank_id -> ligand_code
    pdb_codes = set()
    binding_records = []
    
    for record in tqdm(data, desc="  Scanning records"):
        uniprot = record['UniProt_ID']
        ligand = record['Ligand_ID']
        drugbank = record.get('DrugBank')
        pdb_code = record['PDB_ID']
        binding_site = record['Binding_site']
        binding_moad = record.get('Binding_MOAD')
        binding_pdbbind = record.get('Binding_PDBbind')
        binding_bindingdb = record.get('Binding_BindingDB')
        
        proteins.add(uniprot)
        ligands.add(ligand)
        
        if drugbank and drugbank not in drugs:
            drugs[drugbank] = ligand
        
        pdb_codes.add(pdb_code)
        
        binding_records.append({
            'uniprot': uniprot,
            'ligand': ligand,
            'binding_site': binding_site,
            'pdb_code': pdb_code,
            'drugbank': drugbank,
            'binding_moad': binding_moad,
            'binding_pdbbind': binding_pdbbind,
            'binding_bindingdb': binding_bindingdb
        })
    
    # Insert proteins
    print(f"\nInserting {len(proteins):,} proteins...")
    for uniprot in tqdm(proteins, desc="  Proteins"):
        cursor.execute(
            "INSERT OR IGNORE INTO proteins (uniprot_id) VALUES (?)",
            (uniprot,)
        )
    conn.commit()
    
    # Get protein IDs
    protein_ids = {}
    cursor.execute("SELECT protein_id, uniprot_id FROM proteins")
    for row in cursor.fetchall():
        protein_ids[row[1]] = row[0]
    
    # Insert ligands
    print(f"\nInserting {len(ligands):,} ligands...")
    for ligand_code in tqdm(ligands, desc="  Ligands"):
        cursor.execute(
            "INSERT OR IGNORE INTO ligands (ligand_code) VALUES (?)",
            (ligand_code,)
        )
    conn.commit()
    
    # Get ligand IDs
    ligand_ids = {}
    cursor.execute("SELECT ligand_id, ligand_code FROM ligands")
    for row in cursor.fetchall():
        ligand_ids[row[1]] = row[0]
    
    # Insert PDB structures
    print(f"\nInserting {len(pdb_codes):,} PDB structures...")
    for pdb_code in tqdm(pdb_codes, desc="  PDB structures"):
        cursor.execute(
            "INSERT OR IGNORE INTO pdb_structures (pdb_code) VALUES (?)",
            (pdb_code,)
        )
    conn.commit()
    
    # Get PDB IDs
    pdb_ids = {}
    cursor.execute("SELECT pdb_id, pdb_code FROM pdb_structures")
    for row in cursor.fetchall():
        pdb_ids[row[1]] = row[0]
    
    # Insert drugs
    print(f"\nInserting {len(drugs):,} drugs...")
    for drugbank, ligand in tqdm(drugs.items(), desc="  Drugs"):
        if ligand in ligand_ids:
            cursor.execute(
                "INSERT OR IGNORE INTO drugs (drugbank_id, ligand_id) VALUES (?, ?)",
                (drugbank, ligand_ids[ligand])
            )
    conn.commit()
    
    # Insert binding sites
    print(f"\nInserting {len(binding_records):,} binding sites...")
    for record in tqdm(binding_records, desc="  Binding sites"):
        cursor.execute("""
            INSERT INTO binding_sites 
            (protein_id, ligand_id, binding_site, pdb_id, drugbank_id, 
             binding_moad, binding_pdbbind, binding_bindingdb)
            VALUES (?, ?, ?, ?, ?, ?, ?, ?)
        """, (
            protein_ids[record['uniprot']],
            ligand_ids[record['ligand']],
            record['binding_site'],
            pdb_ids[record['pdb_code']],
            record['drugbank'],
            record['binding_moad'],
            record['binding_pdbbind'],
            record['binding_bindingdb']
        ))
    conn.commit()
    
    # Print statistics
    print("\n" + "=" * 60)
    print("DATABASE STATISTICS")
    print("=" * 60)
    cursor.execute("SELECT COUNT(*) FROM proteins")
    print(f"  Proteins: {cursor.fetchone()[0]:,}")
    cursor.execute("SELECT COUNT(*) FROM ligands")
    print(f"  Ligands: {cursor.fetchone()[0]:,}")
    cursor.execute("SELECT COUNT(*) FROM drugs")
    print(f"  Drugs: {cursor.fetchone()[0]:,}")
    cursor.execute("SELECT COUNT(*) FROM pdb_structures")
    print(f"  PDB structures: {cursor.fetchone()[0]:,}")
    cursor.execute("SELECT COUNT(*) FROM binding_sites")
    print(f"  Binding sites: {cursor.fetchone()[0]:,}")


def main():
    parser = argparse.ArgumentParser(
        description='Create SQLite database from final.json'
    )
    parser.add_argument('--json', default='final.json',
                       help='Path to final.json file')
    parser.add_argument('--db', default='biolip.db',
                       help='Output SQLite database file')
    
    args = parser.parse_args()
    
    print("=" * 60)
    print("CREATING DATABASE FROM final.json")
    print("=" * 60)
    
    if not Path(args.json).exists():
        print(f"ERROR: {args.json} not found!")
        return
    
    conn = create_database(args.db)
    load_data_to_database(conn, args.json)
    conn.close()
    
    print(f"\nDatabase saved to: {args.db}")


if __name__ == "__main__":
    main()