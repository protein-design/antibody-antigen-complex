import os
import numpy as np
import pandas as pd
import mysql.connector

from src.utils import Utils
from src.process_pdb import ProcessPdb

class DbBuild:
    config = {
        'host': 'localhost',
        'user': 'admin',
        'password': 'strong_password',
        'database': 'COMPLEX'
    }

    def __init__(self, data=None, verbose:bool=False):
        self.data = data
        self.verbose = verbose

    def drop_table(self, table_name):
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            query = f"""
                DROP TABLE {table_name};
            """
            cursor.execute(query)
            connection.commit()
            return cursor
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None

    def create_table(self, query):
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query)
            connection.commit()
            return cursor
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None

    def insert_data(self, insert_query, records):
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.executemany(insert_query, records)
            connection.commit()
            return cursor.rowcount
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}. data={records}")
            pass
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None

    def insert_complex(self):
        '''
        table complex
        '''
        cols = ['pdb_id', 'complex_type']
        df = self.data[cols].drop_duplicates()
        data = df.to_records(index=False)
        insert_query = """
            INSERT INTO complex (pdb_id, complex_type)
            VALUES (%s, %s)
        """

        n, m = 0, 0
        for rec in data:
            records = [tuple(rec),]
            r = self.insert_data(insert_query, records)
            if r:
                n+= 1
            else:
                m += 1
        print(f"table complex: total={len(data)}, succeed={n}, unchanged={m}")

    def insert_abag(self):
        '''
        table abag
        '''
        # avoid error of int64
        self.data['compound_no'] = self.data['compound_no'].astype(object)
        cols = ['pdb_id', 'compound_no', 'chain_type', 'fragment']
        self.data = self.data.where(pd.notnull(self.data), None)
        data = self.data[cols].to_records(index=False)
        insert_query = """
            INSERT INTO abag(pdb_id, compound_no, chain_type, fragment)
            VALUES (%s, %s, %s, %s)
        """

        n, m = 0, 0
        for rec in data:
            records = [tuple(rec),]
            r = self.insert_data(insert_query, records)
            if r:
                n+= 1
            else:
                m += 1
        print(f"table abag: total={len(self.data)}, succeed={n}, unchanged={m}")

    def record_pdb(self) -> tuple:
        '''
        self.data is pdb paths
        record for insertion into table pdb
        '''
        for pdb_id, local_path in self.data.items():
            pdb_url = f"https://www.rcsb.org/structure/{pdb_id}"
            p = ProcessPdb(local_path, None, False)
            release_date = p.structure.header['release_date']
            name = p.structure.header['name']
            resolution = p.structure.header['resolution']
            rec = (pdb_id, pdb_url, local_path, release_date, resolution, name)
            yield pdb_id, [rec,]

    def record_compound(self) -> tuple:
        '''
        self.data is pdb paths
        record for insertion into table compound
        '''
        for pdb_id, local_path in self.data.items():
            p = ProcessPdb(local_path, None, False)
            compound = p.structure.header['compound']
            record = []
            for compound_no, cp in compound.items():
                compound_no = int(compound_no)
                chain = cp.get('chain').replace(' ', '')
                _text = []
                for _key, _value in cp.items():
                    if _key not in ('chain',) and _value:
                        _text.append(f"{_key}: {_value}")
                rec = (pdb_id, compound_no, chain, ' ; '.join(_text))
                record.append(rec)
            if record:
                yield pdb_id, record

    def record_aacdb(self):
        '''
        self.data is df of aacdb summary
        '''
        for ix, row in self.data.iterrows():
            pdb_id = row['pdb'].upper()
            chains = row['chains'].split('_')
            rec = (
                pdb_id,
                chains[0],
                chains[1] if len(chains) > 1 else None,
                row['antibody'] ,
                row['protein'],
                None if np.isnan(row['resolution']) else row['resolution'],
            )
            yield pdb_id, [rec,]

