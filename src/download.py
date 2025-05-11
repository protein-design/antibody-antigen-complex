'''
download data from APIs
'''
import os
import json
import pickle
import requests
from dir import Dir
import urllib

class Download:

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

    def getpdb(pdb_entry, out_path):
        """
        Get the PDB file from sabdab.
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