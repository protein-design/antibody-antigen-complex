'''
retrieve antibody-antigen complex for fine-tune LLM
'''
from copy import deepcopy
import itertools
import pickle
import os
import re
from Bio.PDB import PDBParser, PDBIO, Structure, Model

from process_pdb import ProcessPdb
from src.cal_sasa import CalSasa
from src.cal_dist import CalDist

class ParseAbAg(ProcessPdb):

    def __init__(self, pdb_file, outdir=None):
        super().__init__(pdb_file, outdir)
        self.bounded = []
    
    def filter_antibody_antigen(self, outdir:str=None):
        '''
        update self.info
        '''
        # resolution
        resolution = self.structure.header.get('resolution')
        if resolution is None or int(resolution) > 5:
            return None
        
        #Initialize outdir only if complex are detected
        # 3 chains complex: heavy + light ~ receptor
        self.detect_heavy_light()
        if self.info:
            self.init_outdir(outdir)
            return self.info
        
        # 2 chains complex: other cases
        self.detect_antibody()
        if self.info:
            self.init_outdir(outdir)
            return self.info
        
    def detect_antibody(self):
        '''
        - parse single domain antibody (nanobody) - receptor
        '''
        antibody, receptor = [], []
        for v in self.structure.header['compound'].values():
            chain_id = v['chain'].upper().replace(' ', '')
            ids = chain_id.split(',')
            mol_str = v['molecule'].lower()
            #other case 1DZB: single-chain variable fragment 
            # 4HEP: VHH17 domain
            if re.findall(r'antibody|scfv|immunoglobulin|vhh17', mol_str):
                antibody.extend(ids)
            else:
                receptor.extend(ids)
        
        # two chains are existing
        if antibody == [] or receptor == []:
            return None

        # update bounded pairs using castescian product
        for _a, _r in itertools.product(antibody, receptor):
            _info = {
                _a: 'antibody',
                _r: 'receptor',
            }
            if _info not in self.info:
                self.info.append(_info)
            _bounded = {
                'type': 'antibody-receptor',
                'chains': [[_a,], [_r,], [_a, _r]],
                'pdb_files': [],
            }
            if _bounded not in self.bounded:
                self.bounded.append(_bounded)

    def detect_heavy_light(self):
        '''
        parse heavy, light, receptor (antigen) chain
        '''
        heavy, light, receptor = [], [], []
        for v in self.structure.header['compound'].values():
            chain_id = v['chain'].upper().replace(' ', '')
            ids = chain_id.split(',')
            if 'heavy chain' in v['molecule'].lower():
                heavy.extend(ids)
            elif 'light chain' in v['molecule'].lower():
                light.extend(ids)
            else:
                receptor.extend(ids)
        # print(heavy, light, receptor)

        # 3 chains are existing
        if heavy == [] and light == [] and receptor == []:
            return None

        if self.structure_id == '3BQU':
            heavy, light, receptor = ['D',], ['C',], ['B',]

        # heavy/light chain
        for _name, donor in [('heavy', heavy), ('light', light)]:
            for _d, _r in itertools.product(donor, receptor):
                _info = {
                    _d: _name,
                    _r: 'receptor',
                }
                if _info not in self.info:
                    self.info.append(_info)
                _bounded = {
                    'type': _name + '-receptor',
                    'chains': [[_d,], [_r,], [_d, _r]],
                    'pdb_files': [],
                }
                if _bounded not in self.bounded:
                    self.bounded.append(_bounded)

    def chains_to_pdb(self):
        """
        split chains into *pdb
        """
        pool = {}
        for rec in self.bounded:
            for ids in rec['chains']:
                ids_str = ''.join(ids)
                if ids_str not in pool:
                    outfile = self.split_by_chain(ids)
                    pool[ids_str] = outfile
                if pool[ids_str] not in rec['pdb_files']:
                    rec['pdb_files'].append(pool[ids_str])

    def save_df(self, df, file_name):
        if self.outdir:
            file_name = file_name + '.csv'
            outfile = os.path.join(self.outdir, file_name)
            df.to_csv(outfile)
            return outfile
        return None

    def build_freesasa(self):
        """
        update self.bounded
        Note: althernative approach
            example: CalSasa('./pdb/1A14_H_L_N.pdb').run_freesasa_cmd()
        """
        data = []
        for rec in self.bounded:
            rec_data = {}
            c1, c2, c12 = rec['chains']
            c1s, c2s, c12s = ''.join(c1), ''.join(c2), ''.join(c12)
            f1, f2, f12 = rec['pdb_files']
            # calculate SASA
            bounded = CalSasa(f12)
            df12 = bounded.cal_freesasa()
            # delta SASA: light/heavy ~ bounded
            df1 = CalSasa(f1).cal_freesasa()
            delta1 = bounded.cal_delta_sasa(df1)
            # delta SASA: receptor ~ bounded
            df2 = CalSasa(f2).cal_freesasa()
            delta2 = bounded.cal_delta_sasa(df2)

            if delta1 is not None and delta2 is not None:
                k = rec['type']
                s1 = self.save_df(df1, f'{k}_{c1s}_sasa')
                s12 = self.save_df(df12, f'{k}_{c12s}_sasa')
                s2 = self.save_df(df2, f'{k}_{c2s}_sasa')
    
                name1 = f'{k}_{c12s}_{c1s}_delta_total_sasa'
                out_delta1 = self.save_df(delta1, name1)
                name2 = f'{k}_{c12s}_{c2s}_delta_total_sasa'
                out_delta2 = self.save_df(delta1, name2)

                # update self.bounded
                rec_data = deepcopy(rec)
                rec_data.update({
                    'sasa_files': [s1, s2, s12],
                    'delta_total_sasa': {
                        'delta1': delta1,
                        'delta1_file': out_delta1,
                        'delta2': delta2,
                        'delta2_file': out_delta2,
                    },
                })
            if rec_data:
                data.append(rec_data)
        # to pickle
        outfile = os.path.join(self.outdir, 'freesasa.p')
        with open(outfile, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        return outfile

    def build_dist(self):
        '''
        distance of residues between two chains
        Note: must be two chains
        '''
        data = []
        for rec in self.bounded:
            rec_data = {}
            c1, c2, c12 = rec['chains']
            c1s, c2s, c12s = ''.join(c1), ''.join(c2), ''.join(c12)
            # calculate distance of any two residues
            f12 = rec['pdb_files'][-1]
            dist1, dist2 = CalDist(f12).cal_dist(c1s, c2s)

            if dist1 is not None and dist2 is not None:
                k = rec['type']
                name1 = f'{k}_{c1s}_{c2s}_dist'
                out_dist1 = self.save_df(dist1, name1)
                name2 = f'{k}_{c2s}_{c1s}_dist'
                out_dist2 = self.save_df(dist1, name2)

                # update self.bounded
                rec_data = deepcopy(rec)
                rec_data.update({
                    'dist': {
                        'dist1': dist1,
                        'dist1_file': out_dist1,
                        'dist2': dist2,
                        'dist2_file': out_dist2,
                    },
                })
            if rec_data:
                data.append(rec_data)
        # to pickle
        outfile = os.path.join(self.outdir, 'dist.p')
        with open(outfile, 'wb') as f:
            pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
        return outfile