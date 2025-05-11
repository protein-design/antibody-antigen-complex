'''
'''
import numpy as np
import pandas as pd
import os
import json
import pickle
import requests

from dir import Dir
from utils import Utils
from analyze_dist import AnalyzeDist
from analyze_sasa import AnalyzeSasa

class Collect:
    def __init__(self, outputs_dir:str, type_name:str=None, chain_len:int=None):
        self.outputs_dir = outputs_dir
        self.type_name = type_name
        self.chain_len = chain_len

    def meta(self):
        file_iter = Dir(self.outputs_dir).recrusive_files()
        for path in file_iter:
            if path.endswith('meta.json'):
                with open(path, 'r') as f:
                    data = json.load(f)
                    # print(data)
                    meta = []
                    for rec in data['chains']:
                        df = pd.DataFrame(rec['chains'])
                        df = df.assign(
                            model = rec['model_id'],
                            pdb_id = os.path.basename(os.path.dirname(path)),
                        )
                        meta.append(df)
                    yield meta


    def pair_sasa(self):
        names, antibody_sasa, antigen_sasa = [], [], []
        _iter = Utils.iter_files(self.outputs_dir, 'freesasa.p')
        for pdb_id, pfile in _iter:
            p = AnalyzeSasa(pfile)
            sasa = p.get_chain_data(self.type_name, self.chain_len)
            for chain_id in sasa:
                names.append((pdb_id, chain_id[0], chain_id[1]))
                antibody, antigen = sasa[chain_id]
                antibody_sasa.append(antibody['value'])
                antigen_sasa.append(antigen['value'])
        
        antibody_sasa = pd.concat(antibody_sasa, axis=1)
        antibody_sasa.columns = names

        antigen_sasa = pd.concat(antigen_sasa, axis=1)
        antigen_sasa.columns = names
        
        print(antibody_sasa.shape, antigen_sasa.shape)
        return antibody_sasa, antigen_sasa

    def pair_dist(self):
        names, antibody_dist, antigen_dist = [], [], []
        _iter = Utils.iter_files(self.outputs_dir, 'dist.p')
        for pdb_id, pfile in _iter:
            p = AnalyzeDist(pfile)
            dist = p.get_chain_data(self.type_name, self.chain_len)
            for chain_id in dist:
                names.append((pdb_id, chain_id[0], chain_id[1]))
                antibody, antigen = dist[chain_id]
                antibody_dist.append(antibody['value'])
                antigen_dist.append(antigen['value'])
        
        antibody_dist = pd.concat(antibody_dist, axis=1)
        antibody_dist.columns = names
        antibody_dist = -antibody_dist
        
        antigen_dist = pd.concat(antigen_dist, axis=1)
        antigen_dist.columns = names
        antigen_dist = -antigen_dist
        
        print(antibody_dist.shape, antigen_dist.shape)
        return antibody_dist, antigen_dist

    def receptor_sasa(self):
        names, antigen = [], []
        _iter = Utils.iter_files(self.outputs_dir, 'freesasa.p')
        for pdb_id, pfile in _iter:
            p = AnalyzeSasa(pfile)
            sasa = p.get_receptor_data()
            for chain_id in sasa:
                names.append((pdb_id, chain_id[0], chain_id[1]))
                _antigen = sasa[chain_id]
                antigen.append(_antigen['value'].iloc[:self.chain_len])
        antigen = pd.concat(antigen, axis=1)
        antigen = antigen.fillna(0)
        antigen.columns = names
        print(antigen.shape)
        return antigen

    def receptor_dist(self):
        names, antigen = [], []
        _iter = Utils.iter_files(self.outputs_dir, 'dist.p')
        for pdb_id, pfile in _iter:
            p = AnalyzeDist(pfile)
            dist = p.get_receptor_data()
            for chain_id in dist:
                names.append((pdb_id, chain_id[0], chain_id[1]))
                _antigen = dist[chain_id]
                antigen.append(_antigen['value'].iloc[:self.chain_len])
        antigen = pd.concat(antigen, axis=1)
        antigen = antigen.fillna(0)
        antigen.columns = names
        antigen = -antigen
        print(antigen.shape)
        return antigen