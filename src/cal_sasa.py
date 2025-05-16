'''
SASA solvent accessible surface area
“A” (Atom), “R” (Residue), “C” (Chain), “M” (Model), or “S” (Structure).
'''
import os
import numpy as np
import pandas as pd
import re
import subprocess

import freesasa
from Bio.PDB.SASA import ShrakeRupley

from src.process_pdb import ProcessPdb
from src.amino_acids import EXTRA_LONGER_NAMES

class CalSasa:

    def __init__(self, pdb_file:str, chain_ids:str=None, prob_radious:float=None):
        self.pdb_file = pdb_file
        self.pdb_dir = os.path.dirname(pdb_file)
        self.pdb_name = os.path.splitext(os.path.basename(pdb_file))[0]
        self.chain_ids = chain_ids
        self.prob_radious = 1.4 if prob_radious is None else prob_radious
        # updated by self.cal_freesasa()
        self.df = None
    
    def cal_freesasa(self, method=None):
        try:
            hetatm = True if re.findall(r'\-|\+', self.chain_ids) else False
            print(self.chain_ids, hetatm)
            structure = freesasa.Structure(
                self.pdb_file,
                options = {'hetatm': hetatm,}
            )
            method = freesasa.LeeRichards if method == 'L' \
                else freesasa.ShrakeRupley
            calc_args = freesasa.Parameters({
                'algorithm': method,
                'probe-radius': self.prob_radious,
                'n-slices': 100,
            })
            result = freesasa.calc(structure)
        except Exception as e:
            print(e)
            return None

        sasa = []
        res_areas = result.residueAreas()
        for chain_id, chain_sasa in res_areas.items():
            for res_area in chain_sasa.values():
                _sasa = {
                    'chain_id': chain_id,
                    'res_no': res_area.residueNumber,
                    'res': res_area.residueType,
                    'total': res_area.total,
                    'relative_total': res_area.relativeTotal,
                    'apolar': res_area.apolar,
                    'relative_apolar': res_area.relativeApolar,
                    'polar': res_area.polar,
                    'relative_polar': res_area.relativePolar,
                    'side_chain': res_area.sideChain,
                    'relative_side_chain': res_area.relativeSideChain,
                    'main_chain': res_area.mainChain,
                    'relative_main_chain': res_area.relativeMainChain,
                }
                _sasa['aa'] = EXTRA_LONGER_NAMES.get(_sasa['res'], '')
                sasa.append(_sasa)
        self.df = pd.DataFrame(sasa)
        self.df = self.df.set_index(['chain_id', 'res_no'])
        return self.df

    @staticmethod
    def cal_delta_sasa(bounded_df, unbounded_df):
        '''
        calculate delta SASA
        '''
        matched = bounded_df.index.intersection(unbounded_df.index)
        delta = bounded_df.loc[matched, ['res', 'aa']]
        delta['value'] = unbounded_df.loc[matched, 'total'] - bounded_df.loc[matched, 'total']
        delta = delta.reset_index()

        # check if tight connection between two chains
        #  at least tightly connected 3 residues
        sig = delta[delta['value'] > 0]
        if len(sig) >= 3:
            return delta
        if len(sig) > 0:
            print(sig)
        return None

    def cal_total_sasa(self) -> dict:
        """
        Calculate total SASA using biopython
        """
        p = ProcessPdb(self.pdb_file)
        p.load_structure()

        sr = ShrakeRupley(probe_radius=self.prob_radious, n_points=100)
        sr.compute(p.structure, level="S")
        return p.structure.sasa

    def cal_total_freesasa(self):
        '''
        calculate SASA using freesasa python module
        '''
        structure = freesasa.Structure(self.pdb_file)
        result = freesasa.calc(structure)
        res = freesasa.classifyResults(result, structure)
        res['total'] = result.totalArea()
        return res

    def cal_sasa_dict(self) -> dict:
        """
        Calculate SASA for all residues
        """
        p = ProcessPdb(self.pdb_file)
        p.load_structure()

        sr = ShrakeRupley(probe_radius=self.prob_radious, n_points=100)
        sr.compute(p.structure, level="R")

        sasa_dict = {}
        for model in p.structure:
            for chain in model:
                for residue in chain:
                    res_id = (chain.id, residue.id[1])
                    sasa_dict[res_id] = 0
                    for atom in residue:
                        sasa = result.atomArea(atom.serial_number)
                        sasa_dict[res_id] += sasa
        return sasa_dict

    def cal_sasa_df(self, model_id, chain_id):
        '''
        return data frame
        '''
        p = ProcessPdb(self.pdb_file)
        p.load_structure()

        sr = ShrakeRupley(probe_radius=self.prob_radious, n_points=100)
        sr.compute(p.structure, level="C")
    
        for model in p.structure:
            for chain in model:
                if model.id == model_id and chain.id == chain_id:
                    sasa = []
                    for res in chain:
                        sasa.append((res.id, res.sasa))
                    df = pd.DataFrame(sasa, columns=['res_id', 'sasa'])
                    return df
    
    def run_freesasa_cmd(self, outdir=None, method=None):
        '''
        freesasa should be installed
        '''
        if outdir is None:
            outdir = self.pdb_dir
        outfile = os.path.join(outdir, self.pdb_name + '.txt')
        method = '-S' if method == 'S&R' else '-L'
        # run command
        cmd = ['freesasa', self.pdb_file, f"--output={outfile}",\
            '--format=seq', f'-p {self.prob_radious}', method]
        subprocess.call(cmd)
