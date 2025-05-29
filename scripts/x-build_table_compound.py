'''
the script must be executed after build_db_pdb.py
insert data into table compound and abag
'''
import os
import sys
project_dir = os.path.dirname(os.getcwd())
print('project direcotry is ', project_dir)
if project_dir not in sys.path:
    sys.path.append(project_dir)
    
from src.build_db import BuildDb
from src.utils import Utils

abag_create_query = f"""
    CREATE TABLE abag (
        pdb_id VARCHAR(10),
        compound_no INT,
        chain_type VARCHAR(30),
        fragment VARCHAR(30),
        CONSTRAINT pk_abag PRIMARY KEY (pdb_id, compound_no),
        CONSTRAINT fk_abag_compound FOREIGN KEY (pdb_id, compound_no)
            REFERENCES compound(pdb_id, compound_no)
            ON DELETE CASCADE
    );
"""


def main(pdb_dir):
    pdb_path_dict = Utils.get_raw_pdb(pdb_dir)
    b = BuildDb(pdb_path_dict)

    # drop
    table_name = 'compound'
    res = b.drop_table(table_name)
    print(f'try to drop table {table_name}: {res}')

    # create table
    query = f"""
        CREATE TABLE {table_name}(
            pdb_id VARCHAR(10),
            compound_no INT,
            chain VARCHAR(200),
            compound_text Text,
            CONSTRAINT pk_compound PRIMARY KEY (pdb_id, compound_no),
            CONSTRAINT fk_compound_pdb FOREIGN KEY (pdb_id)
                REFERENCES pdb(pdb_id)
                ON DELETE CASCADE
        );
    """
    res = b.create_table(query)
    print(f'try to create table {table_name}: {res}')


    # get records
    records_iter = b.record_compound()
    insert_query = f"""
        INSERT INTO {table_name} (pdb_id, compound_no, chain, compound_text)
        VALUES (%s, %s, %s, %s)
    """
    succeed, failed = 0, []
    for pdb_id, record in records_iter:
        # insert data into table compound
        m = b.insert_data(insert_query, record)
        if m:
            succeed += 1
        else:
            failed.append(pdb_id)
    print(f"succeed={succeed}, failed={len(failed)}")
    print(failed)

if __name__ == "__main__":
    pdb_dir = sys.argv[1]
    main(pdb_dir)
