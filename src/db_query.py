import os
import pandas as pd
import mysql.connector
from copy import deepcopy

class DbQuery:
    config = {
        'host': 'localhost',
        'user': 'admin',
        'password': 'strong_password',
        'database': 'COMPLEX'
    }

    def __init__(self, verbose:bool=False):
        self.verbose = verbose

    def list_data(self, query:str, is_df:bool=False):
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query)
            rows = cursor.fetchall()
            if is_df:
                df = pd.DataFrame(rows)
                return df
            return rows
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None

    def table_names(self):
        query = "SHOW TABLES;"
        res = self.list_data(query, False)
        res = pd.DataFrame(res).iloc[:,0]
        return tuple(res)

    def get_latest_date(self, table_name):
        query = f"""
            SELECT MAX(release_date) AS date FROM {table_name};
        """
        res = self.list_data(query, False)
        return res[0]['date']

    def get_pdb_ids(self, table_name:str):
        query = f"""
            SELECT DISTINCT pdb_id FROM {table_name};
        """
        pdb_ids = self.list_data(query, True)
        if len(pdb_ids) > 0:
            return tuple(pdb_ids['pdb_id'])
        return tuple([])

    def delete_pdb_ids(self, table_name:str):
        query = f"""
            DELETE FROM {table_name}
            WHERE chain_type = "x";
        """
        pdb_ids = self.list_data(query)
        return tuple(pdb_ids['pdb_id'])

    def count_complex(self):
        query = f"""
            SELECT complex_type, COUNT(DISTINCT pdb_id) AS pdb_ids
            FROM complex
            GROUP BY complex_type
            ORDER BY complex_type;
        """
        res = self.list_data(query)
        return pd.DataFrame(res)

    def scan_abag(self):
        query = f"""
            SELECT A.pdb_id, A.compound_no, A.chain_type, B.chain, C.local_path
            FROM abag A, compound B, pdb C
            WHERE A.pdb_id = B.pdb_id
                AND A.compound_no = B.compound_no
                AND A.pdb_id = C.pdb_id
            ORDER BY A.pdb_id, A.compound_no;
        """
        res = self.list_data(query)
        g = pd.DataFrame(res).groupby(['pdb_id', 'local_path'])
        for (pdb_id, local_path), sub in g:
            rec = {
                'pdb_id': pdb_id,
                'local_path': local_path,
                'chains': {},
            }
            for item in sub.to_dict(orient='records'):
                chain_type = item['chain_type']
                if chain_type not in rec['chains']:
                    rec['chains'][chain_type] = []
                rec['chains'][chain_type] += item['chain'].split(',')
            yield rec
