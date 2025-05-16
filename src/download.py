'''
download data from APIs
'''
import os
import json
import pickle
import requests
import urllib

from src.dir import Dir

class Download:

    @staticmethod
    def download_pdb(pdb_id:str, pdb_dir:str=None) -> str:
        '''
        download *.pdb from PDB
        the local dir is ../pdb in default
        '''
        if pdb_dir is None:
            pdb_dir = '../pdb'
        pdb_file = os.path.join(pdb_dir, f'{pdb_id}.pdb')
        if not os.path.isfile(pdb_file):
            print(f"download {pdb_id} to {pdb_file}")
            url = 'https://files.rcsb.org/view'
            os.system(f"wget -qnc {url}/{pdb_id}.pdb -P {pdb_dir}")
        print(f"{pdb_file} is ready.")
        return pdb_file

    def getpdb(pdb_entry, out_path):
        """
        Get the PDB file from SAbDab.
        Check that it has successfully downloaded.
        """
        out_file = os.path.join(out_path, f"{pdb_entry}.pdb")
        endpoint = "https://opig.stats.ox.ac.uk/webapps/abdb/entries"
        url = f"{endpoint}/{pdb_entry}/structure/{pdb_entry}.pdb"
        urllib.request.urlretrieve(url, out_file)
        if os.path.isfile(out_file):
            Retrieved = open(out_file).read()
            if not Retrieved.count("ATOM"):
                print("Failed to retrieve PDB file from SAbDab")
                os.remove(out_file)
                return False
            else:
                return True
        else:
            return False

    @staticmethod
    def pdb_to_uniprot(pdb_id) -> dict:
        '''
        Uniprot-SIFTS
        mapping of pdb -> uniprot
        '''
        pdb_id = pdb_id.lower()
        url = f"https://www.ebi.ac.uk/pdbe/api/mappings/{pdb_id}"
        print(url)
        try:
            res = requests.get(url)
            data = json.loads(res.text)
            return data[pdb_id]
        except Exception as e:
            print(e)
        return None

