'''
retrieve antibody-antigen complex for fine-tune LLM
'''
from copy import deepcopy
import itertools
import pickle
import os
import re
from Bio.PDB import PDBParser, PDBIO, Structure, Model

from src.process_pdb import ProcessPdb
from src.cal_sasa import CalSasa
from src.cal_dist import CalDist
from src.parse_header import ParseHeader


class ParseAbAg(ProcessPdb):

    def __init__(self, pdb_file, outdir=None, verbose:bool=True):
        super().__init__(pdb_file, outdir, verbose)
        self.bounded = []
        self.complex = {}
    
    def filter_antibody_antigen(self, outdir:str=None):
        '''
        update self.info
        '''
        self.complex = {
            'pdb_id': self.structure_id,
            'resolution': self.structure.header['resolution'],
            'release_date': self.structure.header['release_date'],
        }
        #Initialize outdir only if complex are detected
        # 1 chains complex: heavy + light ~ receptor
        self.detect_heavy_light()
        if self.info:
            self.init_outdir(outdir)
            return self.info
        
        # 2 chains complex: single-domain antibody
        self.detect_antibody()
        if self.info:
            self.init_outdir(outdir)
            return self.info
        
        # 3 chain complex: antibody as antigen
        self.detect_antibody_antibody()
        if self.info:
            self.init_outdir(outdir)
            return self.info

    def detect_antibody_antibody(self):
        '''
        - parse antibody - antibody (antigen)
        update self.complex, self.info, self.bounded
        '''
        p = ParseHeader(self.structure)
        donor, receptor = p.antibody_x()
        if donor and receptor == []:
            # antibody antigen
            if len(donor) > 1:
                self.complex['complex_type'] = 'antibody_antibody'
                self.complex['compound'] = p.get_compound()
                for _d, _r in itertools.combinations(donor, 2):
                    self.update_info_bounded(_d, _r)
            # antigen is ligand
            else:
                self.complex['complex_type'] = 'antibody_ligand'
                self.complex['compound'] = p.get_compound()
                for _d in donor:
                    ligand = ('ligand', _d[1] + '-')
                    self.update_info_bounded(_d, ligand)

    def detect_antibody(self):
        '''
        - parse single domain antibody (nanobody) - receptor
        update self.complex, self.info, self.bounded
        '''
        # manual type: match single or antigen
        MANUAL_DEFINE = {
            '6YE3': ['single', 'single', 'antigen'],
            '7R23': ['single', 'antigen'],
        }
        chain_types = MANUAL_DEFINE.get(self.structure_id)
        p = ParseHeader(self.structure)
        antibody, receptor = p.singledomain_antigen(chain_types)
        if antibody == [] or receptor == []:
            return None

        self.complex['complex_type'] = 'single-domain-antibody_antigen'
        self.complex['compound'] = p.get_compound()
        # update bounded pairs using castescian product
        for _a, _r in itertools.product(antibody, receptor):
            self.update_info_bounded(_a, _r)

    def detect_heavy_light(self):
        '''
        parse heavy, light, receptor (antigen) chain
        update self.complex, self.info, self.bounded
        '''
        p = ParseHeader(self.structure)
        heavy_light, receptor = p.heavylight_antigen()

        # type chains should be existing
        if heavy_light == [] or receptor == []:
            return None

        self.complex['complex_type'] = 'heavy-light-antibody_antigen'
        self.complex['compound'] = p.get_compound()
        # update heavy/light chains
        for _d, _r in itertools.product(heavy_light, receptor):
            self.update_info_bounded(_d, _r)

    def update_info_bounded(self, antibody:tuple, receptor:tuple):
        '''
        '''
        antibody_name, antibody_chain = antibody
        receptor_name, receptor_chain = receptor
        # self.info
        _info = {
            antibody_chain: antibody_name,
            receptor_chain: receptor_name,
        }
        if _info not in self.info:
            self.info.append(_info)

        # self.bounded
        bound_chain = [antibody_chain, receptor_chain]
        if receptor_chain[-1] == '-' and antibody_chain == receptor_chain[:-1]:
            bound_chain = [antibody_chain + '+',]
        _bounded = {
            'type': antibody_name + '-' + receptor_name,
            'chains': [[antibody_chain,], [receptor_chain,], bound_chain],
            'pdb_files': [],
        }
        if _bounded not in self.bounded:
            self.bounded.append(_bounded)

    def update_chains(self, new_key:str, info:dict=None):
        '''
        integrate pairwise chains from self.info  to self.chains
        '''
        if info is None:
            info = self.info
        if info:
            print(f"Pairwise chains of {self.structure_id} are updated.")
        for rec in info:
            for chain_id, val in rec.items():
                for _model in self.chains:
                    for _chain in _model['chains']:
                        if _chain['chain_id'] == chain_id:
                            _chain[new_key] = val

    def chains_to_pdb(self):
        """
        split chains into *pdb
        """
        pool = {}
        for rec in self.bounded:
            for ids in rec['chains']:
                ids_str = ''.join(ids)
                if ids_str not in pool:
                    outfile = self.split_by_chain(ids, ids_str)
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
            c1, c2, c12 = rec['chains']
            c1s, c2s, c12s = ''.join(c1), ''.join(c2), ''.join(c12)
            f1, f2, f12 = rec['pdb_files']
            k = rec['type']

            s1, s12, s2 = None, None, None
            delta1, delta2 = None, None
            out_delta1, out_delta2 = None, None
            # calculate SASA of bounded
            df12 = CalSasa(f12, c12s).cal_freesasa()
            if df12 is None:
                print(f"warning: Bounded-SASA for {self.structure_id} failed.")
                return None

            s12 = self.save_df(df12, f'{k}_{c12s}_sasa')
            # delta SASA: light/heavy ~ bounded
            df1 = CalSasa(f1, c1s).cal_freesasa()
            if df1 is not None:
                delta1 = CalSasa.cal_delta_sasa(df12, df1)
                if delta1 is not None:
                    s1 = self.save_df(df1, f'{k}_{c1s}_sasa')
                    name1 = f'{k}_{c12s}_{c1s}_delta_total_sasa'
                    out_delta1 = self.save_df(delta1, name1)
            # delta SASA: receptor ~ bounded
            df2 = CalSasa(f2, c2s).cal_freesasa()
            if df2 is not None:
                delta2 = CalSasa.cal_delta_sasa(df12, df2)
                if delta2 is not None:
                    s2 = self.save_df(df2, f'{k}_{c2s}_sasa')
                    name2 = f'{k}_{c12s}_{c2s}_delta_total_sasa'
                    out_delta2 = self.save_df(delta1, name2)
            print(df12.shape, df1.shape)

            if (delta1 is not None) or (delta2 is not None):
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
                data.append(rec_data)
        if data:
            # to pickle
            outfile = os.path.join(self.outdir, 'freesasa.p')
            with open(outfile, 'wb') as f:
                pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
            return outfile
        print(f"WARNING: Total SASA for {self.structure_id} failed.")
        return None

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
        if data:
            outfile = os.path.join(self.outdir, 'dist.p')
            with open(outfile, 'wb') as f:
                pickle.dump(data, f, protocol=pickle.HIGHEST_PROTOCOL)
            return outfile
        print(f"WARNING: Distance of residues for {self.structure_id} failed.")
        return None