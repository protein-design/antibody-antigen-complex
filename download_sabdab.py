'''
Download data from database SAbDab
https://opig.stats.ox.ac.uk
'''
import os
import sys
import requests
import subprocess

from src.dir import Dir

class download:
    url = 'https://opig.stats.ox.ac.uk/webapps/abdb/entries'

    def __init__(self, outdir):
        self.outdir = outdir
        # default objects
        self.pdb_ids = []
    
    def __call__(self):
        # retrieve list of pdb_id 
        self.retrieve_pdb_ids()
        # summary/*.tsv
        self.download_summary()

    def retrieve_pdb_ids(self):
        '''
        retrieve list of pdb_id
        update self.pdb_ids
        '''
        try:
            res = requests.get(f'{self.url}/sabdab_pdbs_list.txt')
            self.pdb_ids = res.text.split('\n')
        except Exception as e:
            print(f"Error: can't retrieve pdb ids from {self.url}. error={e}")

    def download_summary(self):
        outdir = os.path.join(self.outdir, 'summary')
        Dir(outdir).init_dir()
        
        n, m, k = 0, 0, 0
        with open(os.path.join(self.outdir, 'error_summary.txt'), 'w') as f:
            for pdb_id in self.pdb_ids:
                outfile = os.path.join(outdir, f'{pdb_id}.tsv')
                if not os.path.isfile(outfile):
                    endpoint = f'{self.url}/{pdb_id}/summary/{pdb_id}.tsv'
                    cmd = ['wget', '-c', endpoint, '-P', outdir]
                    try:
                        subprocess.run(cmd, check=True)
                        n += 1
                    except Exception as e:
                        f.write(pdb_id + '\n')
                        m += 1
                else:
                    k += 1
        print(f"download summary: succeed={n}, failed={m}, skipped={k}")

if __name__ == '__main__':
    outdir = sys.argv[1]
    download(outdir)()