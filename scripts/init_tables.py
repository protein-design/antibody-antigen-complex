#!/usr/bin/python
'''
build table pdb and insert data into table pdb
'''
import os
import sys
project_dir = os.path.dirname(os.getcwd())
print('project direcotry is ', project_dir)
if project_dir not in sys.path:
    sys.path.append(project_dir)

from src.build_db import BuildDb
from src.utils import Utils

create_pdb_query = """
    CREATE TABLE pdb(
        pdb_id VARCHAR(10),
        pdb_url VARCHAR(50),
        local_path VARCHAR(100),
        release_date DATE,
        resolution float,
        name VARCHAR(500),
        CONSTRAINT pk_pdb PRIMARY KEY (pdb_id)
    );
"""

create_compound_query = """
    CREATE TABLE compound(
        pdb_id VARCHAR(10),
        compound_no INT,
        chain VARCHAR(200),
        compound_text Text,
        CONSTRAINT pk_compound PRIMARY KEY (pdb_id, compound_no),
        CONSTRAINT fk_compound_pdb FOREIGN KEY (pdb_id)
            REFERENCES pdb(pdb_id)  ON DELETE CASCADE
    );
"""

create_complex_query = """
    CREATE TABLE complex(
        pdb_id VARCHAR(10),
        complex_type VARCHAR(50),
        CONSTRAINT pk_complex PRIMARY KEY (pdb_id, complex_type),
        CONSTRAINT fk_complex_pdb FOREIGN KEY (pdb_id)
            REFERENCES pdb(pdb_id) ON DELETE CASCADE
    );
"""

create_abag_query = """
    CREATE TABLE abag (
        pdb_id VARCHAR(10),
        compound_no INT,
        chain_type VARCHAR(30),
        fragment VARCHAR(30),
        CONSTRAINT pk_abag PRIMARY KEY (pdb_id, compound_no),
        CONSTRAINT fk_abag_compound FOREIGN KEY (pdb_id, compound_no)
            REFERENCES compound(pdb_id, compound_no) ON DELETE CASCADE
    );
"""

def main():
    b = BuildDb()

    meta = [
        ('pdb', create_pdb_query),
        ('complex', create_complex_query),
        ('compound', create_compound_query),
        ('abag', create_abag_query),
    ]
    # drop tables
    for table_name, create_query in meta[::-1]:
        res = b.drop_table(table_name)
        print(f'Try to drop table {table_name}:', res)
    
    # create table
    for table_name, create_query in meta:
        res = b.create_table(create_query)
        print(f'Try to create table {table_name}:', res)

if __name__ == "__main__":
    main()
