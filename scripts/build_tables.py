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


def build_pdb(pdb_dir):
    '''
    retrieve data from pdb files and insert recrods into table pdb
    '''
    print("\nTry to retrieve data and insert recrods into table pdb...")
    pdb_path_dict = Utils.get_raw_pdb(pdb_dir)
    b = BuildDb(pdb_path_dict)
    # get records
    records_iter = b.record_pdb()

    # insert data into table pdb
    insert_query = """
        INSERT INTO pdb (pdb_id, pdb_url, local_path, release_date, resolution, name)
        VALUES (%s, %s, %s, %s, %s, %s)
    """

    succeed = 0
    failed = []
    for pdb_id, record in records_iter:
        m = b.insert_data(insert_query, record)
        if m:
            succeed += 1
        else:
            failed.append(pdb_id)
    print(f"failed:{len(failed)}, successful={succeed}")
    print(failed)

def build_compound(pdb_dir):
    '''
    retrieve data from pdb and insert recrods into table compound
    '''
    print("\nTry to retrieve data and insert recrods into table compound...")
    pdb_path_dict = Utils.get_raw_pdb(pdb_dir)
    b = BuildDb(pdb_path_dict)
    # get records
    records_iter = b.record_compound()

    # insert data into table compound
    insert_query = """
        INSERT INTO compound (pdb_id, compound_no, chain, compound_text)
        VALUES (%s, %s, %s, %s)
    """

    succeed, failed = 0, []
    for pdb_id, record in records_iter:
        m = b.insert_data(insert_query, record)
        if m:
            succeed += 1
        else:
            failed.append(pdb_id)
    print(f"succeed={succeed}, failed={len(failed)}")
    print(failed)

if __name__ == "__main__":
    pdb_dir = sys.argv[1]

    # table pdb
    build_pdb(pdb_dir)

    # table compound
    build_compound(pdb_dir)

    print("Done. Great\n\n")
