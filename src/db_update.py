import os
import pandas as pd
import mysql.connector
from copy import deepcopy

class DbUpdate:
    config = {
        'host': 'localhost',
        'user': 'admin',
        'password': 'strong_password',
        'database': 'COMPLEX'
    }

    def __init__(self, verbose:bool=False):
        self.verbose = verbose

    def execute_data(self, query, record):
        record = deepcopy(record)
        try:
            connection = mysql.connector.connect(**self.config)
            cursor = connection.cursor(dictionary=True)
            cursor.execute(query, record)
            connection.commit()
            return cursor.rowcount
        except Exception as e:
            if self.verbose:
                print(f"Error: {e}")
        finally:
            if 'connection' in locals() and connection.is_connected():
                cursor.close()
                connection.close()
        return None

    def update_compound(self, df, col_name:str):
        '''
        table is compound
        '''
        update_query = f"""
            UPDATE compound
            SET {col_name} = %({col_name})s
            WHERE pdb_id = %(pdb_id)s AND compound_no = %(compound_no)s
        """
        n,m=0,0
        cols = [col_name, 'pdb_id', 'compound_no']
        data = df[cols].to_dict(orient='records')
        for rec in data:
            r = self.execute_data(update_query, rec)
            if r:
                n+= 1
            else:
                m += 1
        print(f"{col_name}: succeed={n}, unchanged={m}")
    
    def update_compound_table(self, table_name, col_name, df):
        '''
        The compound tables could be compountd, abag
        where primary key is )pdb_id, compound_no)
        '''
        update_query = f"""
            UPDATE {table_name}
            SET {col_name} = %({col_name})s
            WHERE pdb_id = %(pdb_id)s AND compound_no = %(compound_no)s
        """
        n,m=0,0
        cols = ['pdb_id', 'compound_no', col_name]
        data = df[cols].to_dict(orient='records')
        for rec in data:
            # records = [rec,]
            r = self.execute_data(update_query, rec)
            if r:
                n+= 1
            else:
                m += 1
        print(f"{table_name}.{col_name}: succeed={n}, unchanged={m}")

