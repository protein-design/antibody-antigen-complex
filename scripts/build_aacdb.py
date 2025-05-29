'''
build table pdb and insert data into table pdb
'''
import pandas as pd
import os
import sys
project_dir = os.path.dirname(os.getcwd())
print('project direcotry is ', project_dir)
if project_dir not in sys.path:
    sys.path.append(project_dir)

from src.build_db import BuildDb
from src.utils import Utils

table_name = 'aacdb'
create_query = f"""
    CREATE TABLE {table_name}(
        pdb_id VARCHAR(10),
        antibody_chains VARCHAR(20),
        antigen_chains VARCHAR(20),
        antibody VARCHAR(100),
        protein	 VARCHAR(200),
        resolution float
    );
"""
insert_query = f"""
    INSERT INTO {table_name}(pdb_id, antibody_chains, 
        antigen_chains, antibody, protein, resolution)
    VALUES (%s, %s, %s, %s, %s, %s)
"""

def main(data_dir):
    # download
    url = 'https://i.uestc.edu.cn/AACDB/data_zip'
    os.system(f"wget -c {url}/protein_table.txt -P {data_dir}")

    # get data
    infile = os.path.join(data_dir, 'protein_table.txt')
    aacdb_summary = pd.read_csv(infile, sep='\t')
    print(aacdb_summary.shape)

    # get records
    b = BuildDb(aacdb_summary)
    records_iter = b.record_aacdb()

    # create table
    print(f'try to drop table {table_name}')
    res = b.drop_table(table_name)
    print(f'Try to create table {table_name}')
    res = b.create_table(create_query)

    # insert data into table pdb
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

if __name__ == "__main__":
    data_dir = sys.argv[1]
    main(data_dir)
