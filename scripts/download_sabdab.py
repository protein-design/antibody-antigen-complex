'''
Download data from database SAbDab
https://opig.stats.ox.ac.uk
'''
import numpy as np
import pandas as pd
import os
import sys
project_dir = os.path.dirname(os.getcwd())
print('project direcotry is ', project_dir)
if project_dir not in sys.path:
    sys.path.append(project_dir)

import requests
import subprocess

from src.dir import Dir
from src.collect import Collect
from src.utils import Utils
from src.build_db import BuildDb


table_name = 'sabdab'
create_query = f"""
    CREATE TABLE {table_name}(
        pdb_id VARCHAR(10),
        model INT,
        heavy_chain VARCHAR(20),
        light_chain VARCHAR(20),
        antigen_chain VARCHAR(30),
        antigen_type VARCHAR(50),
        resolution VARCHAR(50)
    );
"""
insert_query = f"""
    INSERT INTO {table_name}(pdb_id, model, heavy_chain, light_chain,
        antigen_chain, antigen_type, resolution)
    VALUES (%s, %s, %s, %s, %s, %s, %s)
"""

class download:
    url = 'https://opig.stats.ox.ac.uk/webapps/abdb/entries'

    def __init__(self, outdir):
        self.outdir = outdir
        self.summary_dir = None
    
    def __call__(self, pdb_dir:str=None):
        # retrieve list of pdb_id from sabdab_pdbs_list.txt
        # Warning: the list pdb in this file doesn't cover all pdb in sabdab.
        self.summary_dir = os.path.join(self.outdir, 'summary')
        Dir(self.summary_dir).init_dir()

        # summary/*.tsv
        sabdab_ids = self.retrieve_sabdab_pdb_ids()
        self.download_summary(sabdab_ids)

        if pdb_dir:
            # retrieve pdb_ids from PDB
            pdb_ids = Utils.filter_out_raw_pdb(pdb_dir, sabdab_ids)
            # download summary/*.tsv
            self.download_summary(pdb_ids)

        # database
        b = BuildDb()
        # create table
        print(f'try to drop table {table_name}')
        res = b.drop_table(table_name)
        print(f'Try to create table {table_name}')
        res = b.create_table(create_query)

        # insert data into table pdb
        succeed = 0
        failed = []
        records_iter = self.record_sabdab()
        for pdb_id, record in records_iter:
            m = b.insert_data(insert_query, record)
            if m:
                succeed += 1
            else:
                failed.append(pdb_id)
        print(f"failed:{len(failed)}, successful={succeed}")
        print(failed)


    def retrieve_sabdab_pdb_ids(self):
        '''
        retrieve list of pdb_id
        update self.pdb_ids
        '''
        try:
            # https://opig.stats.ox.ac.uk/webapps/abdb/entries/sabdab_pdbs_list.txt
            res = requests.get(f'{self.url}/sabdab_pdbs_list.txt')
            pdb_ids = res.text.split('\n')
            return pdb_ids
        except Exception as e:
            print(f"Error: can't retrieve pdb ids from {self.url}. error={e}")
        return []

    def download_summary(self, pdb_ids:list):
        if not pdb_ids:
            print(f"WARNING: summary can not be downloaded.")
            return None

        n, m, k = 0, 0, 0
        with open(os.path.join(self.outdir, 'error_summary.txt'), 'w') as f:
            for pdb_id in pdb_ids:
                outfile = os.path.join(self.summary_dir, f'{pdb_id}.tsv')
                if not os.path.isfile(outfile):
                    endpoint = f'{self.url}/{pdb_id}/summary/{pdb_id}.tsv'
                    cmd = ['wget', '-c', endpoint, '-P', self.summary_dir]
                    try:
                        subprocess.run(cmd, check=True)
                        n += 1
                    except Exception as e:
                        err = f"{pdb_id} | {e}\n"
                        f.write(err)
                        m += 1
                else:
                    k += 1
        print(f"download summary: succeed={n}, failed={m}, skipped={k}")

    def record_sabdab(self):
        sabdab = []
        file_iter = Dir(self.summary_dir).recrusive_files()
        for path in file_iter:
            if path.endswith('.tsv'):
                df = pd.read_csv(path, sep='\t')
                record = []
                pdb_id = None
                if len(df) > 0:
                    for ix, row in df.iterrows():
                        row = row.replace({
                            np.nan: None,
                            'NOT': None,
                        })
                        pdb_id = str(row['pdb']).upper()
                        rec = (
                            pdb_id,
                            row['model'],
                            row['Hchain'],
                            row['Lchain'],
                            row['antigen_chain'],
                            row['antigen_type'],
                            row['resolution'],
                        )
                        record.append(rec)
                    yield pdb_id, record
        


if __name__ == '__main__':
    outdir = sys.argv[1]
    pdb_dir = sys.argv[2] if len(sys.argv) > 2 else None
    download(outdir)(pdb_dir)