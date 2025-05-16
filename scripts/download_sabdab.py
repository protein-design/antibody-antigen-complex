'''
Download data from database SAbDab
https://opig.stats.ox.ac.uk
'''
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

class download:
    url = 'https://opig.stats.ox.ac.uk/webapps/abdb/entries'

    def __init__(self, outdir):
        self.outdir = outdir
    
    def __call__(self, pdb_dir:str=None):
        # retrieve list of pdb_id from sabdab_pdbs_list.txt
        # Warning: the list pdb in this file doesn't cover all pdb in sabdab.
        sabdab_ids = self.retrieve_sabdab_pdb_ids()
        # summary/*.tsv
        self.download_summary(sabdab_ids)

        if pdb_dir:
            # retrieve pdb_ids from PDB
            pdb_ids = []
            raw_iter = Collect(pdb_dir).get_raw_pdb()
            for pdb_id, _ in raw_iter:
                pdb_id = pdb_id.lower()
                if pdb_id not in sabdab_ids:
                    pdb_ids.append(pdb_id)
            # download summary/*.tsv
            self.download_summary(pdb_ids)


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
        outdir = os.path.join(self.outdir, 'summary')
        Dir(outdir).init_dir()
        
        n, m, k = 0, 0, 0
        with open(os.path.join(self.outdir, 'error_summary.txt'), 'w') as f:
            for pdb_id in pdb_ids:
                outfile = os.path.join(outdir, f'{pdb_id}.tsv')
                if not os.path.isfile(outfile):
                    endpoint = f'{self.url}/{pdb_id}/summary/{pdb_id}.tsv'
                    cmd = ['wget', '-c', endpoint, '-P', outdir]
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

if __name__ == '__main__':
    outdir = sys.argv[1]
    pdb_dir = sys.argv[2] if len(sys.argv) > 2 else None
    download(outdir)(pdb_dir)