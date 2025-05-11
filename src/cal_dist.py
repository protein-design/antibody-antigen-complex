'''
distance between two atoms
'''
import os
import numpy as np
import pandas as pd
import subprocess

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

from process_pdb import ProcessPdb
from amino_acids import EXTRA_LONGER_NAMES

class CalDist:

    def __init__(self, pdb_file:str, dist_threshold:float=None):
        self.pdb_file = pdb_file
        self.dist_threshold = 5 if dist_threshold is None else dist_threshold
        # load
        self.load_structure()
    
    def load_structure(self):
        '''
        '''
        self.pdb_dir = os.path.dirname(self.pdb_file)
        self.pdb_name = os.path.splitext(os.path.basename(self.pdb_file))[0]
        self.u = mda.Universe(self.pdb_file)
    
    def cal_dist(self, chain1_id:str, chain2_id:str):
        '''
        '''
        # Select chains
        chain1 = self.u.select_atoms(f"chainID {chain1_id}")
        chain2 = self.u.select_atoms(f"chainID {chain2_id}")

        # Get the coordinates
        coords1 = chain1.positions
        coords2 = chain2.positions
        # print(len(coords1), len(coords2))
        # Calculate the distances between each pair of atoms
        dists = distance_array(coords1, coords2)

        # value is minimum distance
        df1 = pd.DataFrame({
            'res': chain1.resnames,
            'res_no': chain1.resids,
            'atom': [a.name for a in chain1.atoms],
            'value': np.apply_along_axis(np.min, 1, dists),
        })
        df1['aa'] = [EXTRA_LONGER_NAMES.get(i, '') for i in df1['res']]
        g1 = df1.groupby(['res_no']).agg('min').reset_index()
        # 
        df2 = pd.DataFrame({
            'res': chain2.resnames,
            'res_no': chain2.resids,
            'atom': [a.name for a in chain2.atoms],
            'value': np.apply_along_axis(np.min, 0, dists),
        })
        df2['aa'] = [EXTRA_LONGER_NAMES.get(i, '') for i in df2['res']]
        g2 = df2.groupby(['res_no']).agg('min').reset_index()
        return g1, g2

    def print_dist(self, chain1_id:str, chain2_id:str):
        '''
        value is minimum distance
        '''
        g1, g2 = self.cal_dist(chain1_id, chain2_id)
        g1 = g1[g1['value']<=self.dist_threshold]
        print(g1)
        g2 = g2[g2['value']<=self.dist_threshold]
        print(g2)