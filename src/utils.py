'''
'''
import os
import json
import pickle
from datetime import date

from src.dir import Dir

class Utils:

    @staticmethod
    def iter_raw_pdb(indir):
        file_iter = Dir(indir).recrusive_files()
        for path in file_iter:
            file_name = os.path.basename(path)
            pdb_id = file_name.replace('pdb', '').replace('.ent.gz', '').upper()
            if pdb_id:
                yield pdb_id, path

    @staticmethod
    def get_raw_pdb(indir) -> dict:
        res = {}
        file_iter = Utils.iter_raw_pdb(indir)
        for pdb_id, path in file_iter:
            res[pdb_id] = path
        return res

    @staticmethod
    def filter_out_raw_pdb(indir, pdb_ids):
        res = {}
        file_iter = Utils.iter_raw_pdb(indir)
        for pdb_id, path in file_iter:
            if pdb_id not in pdb_ids: 
                res[pdb_id] = path
        return res

    @staticmethod
    def filter_in_raw_pdb(indir, pdb_ids):
        file_iter = Utils.iter_raw_pdb(indir)
        for pdb_id, path in file_iter:
            if pdb_id in pdb_ids:
                yield pdb_id, path

    @staticmethod
    def get_new_raw_pdb(indir, date_stamp):
        res = {}
        file_iter = Utils.iter_raw_pdb(indir)
        for pdb_id, path in file_iter:
            mod_time = os.path.getmtime(path)
            local_time = date.fromtimestamp(mod_time)
            if local_time >= date_stamp:
                res[pdb_id] = path
        return res

    @staticmethod
    def get_meta(outputs_dir):
        ready = {}
        file_iter = Dir(outputs_dir).recrusive_files()
        for path in file_iter:
            if path.endswith('meta.json'):
                pdb_id = os.path.basename(os.path.dirname(path))
                ready[pdb_id] = path
        return ready

    @staticmethod
    def iter_meta(outputs_dir):
        file_iter = Dir(outputs_dir).recrusive_files()
        for path in file_iter:
            if path.endswith('meta.json'):
                pdb_id = os.path.basename(os.path.dirname(path))
                with open(path, 'r') as f:
                    data = json.load(f)
                    yield pdb_id, data

    @staticmethod
    def iter_files(outputs_dir, file_name):
        file_iter = Dir(outputs_dir).recrusive_files()
        for path in file_iter:
            if path.endswith(file_name):
                pdb_id = os.path.basename(os.path.dirname(path))
                yield pdb_id, path

    @staticmethod
    def match_chains(df1, df2):
        '''
        intersection by column names
        '''
        match = df1.columns.intersection(df2.columns)
        df1_match = df1[match].T
        df2_match = df2[match].T
        print(df1_match.shape, df2_match.shape)
        return df1_match, df2_match

    @staticmethod
    def iter_to_txt(_iter, outprefix:str, chunk_size:int=10_000):
        pool, file_list = [], []
        n, ix = 0, 0
        for pdb_id, path in _iter:
            pool.append(path)
            n += 1
            if len(pool) >= chunk_size:
                outfile = f"{outprefix}_{ix}.txt"
                with open(outfile, 'w') as f:
                    f.write('\n'.join(pool))
                pool = []
                ix += 1
                file_list.append(outfile)
        else:
            if pool:
                outfile = f"{outprefix}_{ix}.txt"
                with open(outfile, 'w') as f:
                    f.write('\n'.join(pool))
                file_list.append(outfile)
        print('total', n)
        return file_list

    @staticmethod
    def iter_to_json(_iter, outprefix:str, chunk_size:int=10_000):
        pool, file_list = [], []
        n, ix = 0, 0
        for record in _iter:
            pool.append(record)
            n += 1
            if len(pool) >= chunk_size:
                outfile = f"{outprefix}_{ix}.json"
                with open(outfile, 'w') as f:
                    json.dump(pool, f, indent=4)
                pool = []
                ix += 1
                file_list.append(outfile)
        else:
            if pool:
                outfile = f"{outprefix}_{ix}.json"
                with open(outfile, 'w') as f:
                    json.dump(pool, f, indent=4)
                file_list.append(outfile)
        print('total', n)
        return file_list