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
from src.query_db import QueryDb
from src.utils import Utils

pdb_query = """
    INSERT INTO pdb (pdb_id, pdb_url, local_path, release_date, resolution, name)
    VALUES (%s, %s, %s, %s, %s, %s)
"""

compound_query = """
    INSERT INTO compound (pdb_id, compound_no, chain, compound_text)
    VALUES (%s, %s, %s, %s)
"""

def synchronize_pdb(pdb_dir):
    url = "rsync.rcsb.org"
    endpoint = "ftp_data/structures/divided/pdb/"
    cmd = f"rsync -rlpt -v -z --delete --port=33444 {url}::{endpoint} {pdb_dir}"
    res = os.system(cmd)
    return res

def update_pdb(pdb_path_dict):
    # get new records
    print('Try to update records of table pdb: ', len(pdb_path_dict))
    b = BuildDb(pdb_path_dict, verbose=True)
    records_iter = b.record_pdb()

    # insert data into table pdb
    succeed, failed = 0, []
    for pdb_id, record in records_iter:
        m = b.insert_data(pdb_query, record)
        if m:
            succeed += 1
        else:
            failed.append(pdb_id)
    print(f"table pdb: succeed ={succeed}, failed ={len(failed)}")
    print(failed)

def update_compound(pdb_path_dict):
    print('Try to update records of table compound: ', len(pdb_path_dict))
    b = BuildDb(pdb_path_dict, verbose=True)
    # get records
    records_iter = b.record_compound()

    # insert data into table compound and abag
    succeed, failed = 0, []
    for pdb_id, compound_record in records_iter:
        # insert data into table compound
        m = b.insert_data(compound_query, compound_record)
        if m:
            succeed += 1
        else:
            failed.append(pdb_id)
    print(f"table compound: succeed ={succeed}, failed ={len(failed)}")
    print(failed)


if __name__ == "__main__":
    pdb_dir = sys.argv[1]
    
    # synchronize remote ~ local
    synchronize_pdb(pdb_dir)

    # update table pdb
    db_pdb_ids = QueryDb().get_pdb_ids('pdb')
    print(f"Number of pbd in pdb is {len(db_pdb_ids)}.")
    # screen updated pdb_ids
    pdb_path_dict = Utils.filter_out_raw_pdb(pdb_dir, db_pdb_ids)
    # insert new records
    if pdb_path_dict:
        update_pdb(pdb_path_dict)
    else:
        print('No updates on table pdb.')

    # update table compound
    db_pdb_ids = QueryDb().get_pdb_ids('compound')
    print(f"Number of pbd in compound is {len(db_pdb_ids)}.")
    pdb_path_dict = Utils.filter_out_raw_pdb(pdb_dir, db_pdb_ids)
    if pdb_path_dict:
        update_compound(pdb_path_dict)
    else:
        print('No updates on table compound.')
