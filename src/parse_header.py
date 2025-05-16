'''
Analyze annotations in PDB
detect chains of antibody or antigen
'''
import re

class ParseHeader:
    def __init__(self, structure):
        self.structure = structure
        self.init_compound()

    def init_compound(self):
        '''
        self.compound
        '''
        self.compound = []
        _compound = self.structure.header['compound']
        for i in _compound:
            rec = dict(_compound[i])
            rec['compound_id'] = i
            for k in list(rec):
                rec[k] = rec[k].lower()
            self.compound.append(rec)
    
    def set_ids(self, rec:dict, chain_type:str=None) -> list:
        '''
        standardize format if 'ids'
        '''
        if 'chain' in rec:
            ids = rec['chain'].upper().replace(' ', '').split(',')
            if chain_type is None:
                chain_type = 'other'
            return [(chain_type, i) for i in ids]
        return []

    def get_compound(self):
        for rec in self.compound:
            for k in ('misc', 'compound_id', 'engineered'):
                if k in rec:
                    del rec[k]
        return self.compound

    def heavylight_antigen(self):
        '''
        antibody: heavy/light chain
        antigen: protein/peptide chain
        '''
        antibody, antigen = [], []
        for rec in self.compound:
            ids = self.heavy_chain(rec)
            if ids:
                rec['chain_type'] = 'heavy'
                antibody.extend(ids)
                continue

            ids = self.light_chain(rec)
            if ids:
                rec['chain_type'] = 'light'
                antibody.extend(ids)
                continue

            rec['chain_type'] = 'antigen'
            ids = self.set_ids(rec, 'antigen')
            antigen.extend(ids)
        return antibody, antigen

    def singledomain_antigen(self, chain_types:list=None):
        '''
        antibody: single-domain chain
        antigen: protein/peptide chain
        '''
        single, antigen = [], []
        if chain_types:
            for chain_type, rec in zip(chain_types, self.compound):
                rec['chain_type'] = chain_type
                ids = self.set_ids(rec, chain_type)
                if chain_type == 'single':
                    single.extend(ids)
                else:
                    antigen.extend(ids)
        else:
            for rec in self.compound:
                ids = self.single_domain_chain(rec)
                if ids:
                    rec['chain_type'] = 'single'
                    single.extend(ids)
                else:
                    rec['chain_type'] = 'antigen'
                    ids = self.set_ids(rec, 'antigen')
                    antigen.extend(ids)
        return single, antigen

    def antibody_x(self):
        '''
        antibody: heavy/light/single-domain chain
        x could be antibody chain as antigen or ligand/small molecule
        '''
        antibody, antigen = [], []
        for rec in self.compound:
            ids = self.heavy_chain(rec)
            if ids:
                rec['chain_type'] = 'heavy'
                antibody.extend(ids)
                continue

            ids = self.light_chain(rec)
            if ids:
                rec['chain_type'] = 'light'
                antibody.extend(ids)
                continue
            
            ids = self.single_domain_chain(rec)
            if ids:
                rec['chain_type'] = 'single'
                antibody.extend(ids)
                continue

            rec['chain_type'] = 'antigen'
            ids = self.set_ids(rec, 'antigen')
            antigen.extend(ids)
        return antibody, antigen

    def _source_str(self, rec):
        _source = ['molecule', 'fragment', 'synonym', 'other_details',]
        _str = ' '.join([rec.get(i, '') for i in _source])
        return _str

    def heavy_chain(self, rec:dict):
        '''
        particular keywords:
            7D00: _heavy_, _light_
            5GRJ: h chain, l chain
            7BEM: vh domain, vl domain
            7NWL: ts2/16 vh(..., ts2/16 vl-...
            4KFZ: 'anti-lmo2 vh'
            3U0T: lc fab, hc fab
            7E8F: 'n9 h', '368-2 h', '604 h'
            5NJ6: 'fab3939 h', 'fab3939 l'
            7Y0W: 'bd55-5514h', 'bd55-5840h'
        '''
        pattern_mol = [
            'heavy ', 'heavy.?chain', 'vh domain',
            'h[-| |_]chain', '_heavy_', 'chain[ |-|_]h',
            'vh[\-|\(| ]', 'vhch', 'hc fab',
            'fabh', 'fab.*[-|_]hc', 'fab\d+ h',
            'bd55-\d+h', 'bd-\d+h', '^\d+h', 'sw[a|c]\d+h',
            'h\d+', '\d+ h', 'ab\d*-h', '\db-h',
        ]
        pattern_mol = '|'.join(pattern_mol)
        exact_mol = ('anti-lmo2 vh', 'n9 h', '368-2 h',)
        if re.findall(pattern_mol, self._source_str(rec)) or \
            rec.get('molecule', '') in exact_mol:
            ids = self.set_ids(rec, 'heavy')
            return ids
        return []

    def light_chain(self, rec:dict):
        '''
        7E8F: 'n9 l', '368-2 l', '604 l'
        7Y0W: 'bd55-5514l', 'bd55-5840l'
        '''
        pattern_mol = [
            'light ', 'light.?chain', 'vl domain',
            'l[-| |_]chain', '_light_', 'chain[ |-|_]l',
            'vl[\-|\(| ]', 'vlch', 'lc fab',
            'fabl', 'fab-lc', 'fab\d+ l',
            'bd55-\d+l', 'bd-\d+l', '^\d+l', 'sw[a|c]\d+l',
            'l\d+', '\d+ l', 'ab\d*-l', '\db-l',
        ]
        pattern_mol = '|'.join(pattern_mol)
        exact_mol = ('n9 l', '368-2 l', )
        if re.findall(pattern_mol, self._source_str(rec)) or \
            rec.get('molecule', '') in exact_mol:
            ids = self.set_ids(rec, 'light')
            return ids
        return []

    def single_domain_chain(self, rec:dict):
        '''
        particular keywords
            4HEP: vhh17
            4PIR: vhh15
            5FUC: vhh6
            1DZB: scfv fragment
            6Y6C: single chain variable
            6TOU: 'single-chain fv'
            1NMB: fab
            6UI1: 'cia-h7', 'cia-d12', 'cia-b5'
            4UV7: 'gc1118a'
            5GRU: diabody
            7ZMN: gluebody
            5F97: nanbody  (typo error)
            6QGY: 'nanob12'
            6QGW: 'nanoe6'
            6QGX: 'nanof7'
            4BEL: 'xa4813',
            6UHT: 'jli-g10' is VHH domain
            7K84: 'jle-e5'
            7D2Z: sr31 is sybody
            7AEJ: VHH
            3MAC: 'fab8062'
            7K7Y: 'jle-e9' is vhh domain
            6UKT: 'nb8109', 'nb8117', 'nb8119', 'nb9156'
            4W6W: 'nbfedf6'
            6RU5: 'hc3nb1'
            6EQI: 'nb696'
            6EY6: 'nb130'
            5O0W: 'nb474'
            7EOW: 'caplacizumab'
            7RNN: nanobodies
            5VNW: 'nb.b201'
            5A2I: 'ig lambda-1 chain v region s43',
            7LX5: 'wnb 2', 'wnb 10'
            5GS2: 'anti-mbp', 'repebody'
            6I53: 'megabody38'
            6UC6: 'jli-h11'
            6OQ8: 7f is vhh domain
            7D5P: 'icab'
            3G9A: 'minimizer'
            7D8B: 'vh-s4'
            7WD1: 'r14'
            7FBJ: 'new antigen receptor variable domain'
            7NA9: 'jsg-c1'
            7FAU: 'nb_1b11'
            7S83: vnar
            7X2M: '1-2c7'
            7WVF: 'mab12'
            7R1Z: nbarc-h11, nbarc-c11
            7SAK: 'lam4'
            7WHI: 'bn03_nano1', 'bn03_nano2'
            7UI1: '230al-37'
            9FVC: 'vhh_vcp#2'
            8XQJ: 'jn241'
            7XQV: 'rh57'
            7VAB: nonobody
            8VA0: 'jgfn4'
            8COE: 'lcp0195'
            7V61: '3e8'
            7UA2: 'lmiv230-01', '230al-18'
            8HBG: 'm678f nab'
        '''
        # 'nb8109', 'nb8117', 'nb8119', 'nb9156', 'nb696', 'nb130', 'nb474'
        #  'vhh1', 'vhh15',
        # 'nanoe6', 'nanof7',
        # 'jle-e9', 'jli-g10', 'jli-h11', 'jle-e5',
        pattern = [
            'antibody', 'diabody', 'sybody', 'gluebody',
            'introbody', 'intrabody', 'megabody', 'repebody',
            'nanobody', 'nanbody', 'nanobodies', 'nanoboy', 'nonobody',
            'antigen receptor', 'immunoglobulin', 'promacrobody',
            'igg ', 'igh ', 'ig lambda', 'ig kappa', 'ig gamma', 
            'fab ', 'fab[\-|\_][a-z]', 'fv ', 'fab\d+',
            'scfv', 'sr31', 'sdab\d*', 'sumo-',
            'vhh\d+', 'vhh[ |\-|\_]', '\(vh\)',
            'single.chain (?:variable|fv|fragment)',
            'nb\.[b|x]\d+', '^nb[-|_]', '^nbarc', '^nbroco',
            'n\d+', 'nb\d+', 'nba\d+', 'm\d+f nab',
            'nano.*[0-9]', '_nano\d+', 'nanosota',
            'wnb \d+', 'vnar', '^anti-', 'lag\d+',
            '^230', '^lam\d+', '^mb\d*', 's.*-[h|k]c',
            'jl[a-z]-[a-z]\d+', 'bd55-', 'cia-h\d+',
        ]
        pattern = '|'.join(pattern)
        exact_mol = (
            'gc1118a', 'hc3nb1', 'xa4813', '230al-37',
            'nbfedf6', 'jn241', 'jgfn4', 'm3f', 'svf16',
            'caplacizumab', 'megabody38', 'anti-mbp', 
            'syb_nl5', 'ca14', 'lmiv230-01', '230al-18',
            '7f', '5d', 'e3', 'icab', 'minimizer', 'vh-s4',
            'r14', 'jsg-c1', '1-2c7', 'mab12', 'lam4', 'rh57',
            'lcp0195', '3e8', 'nbe201', 'sia',
            'sy45', 'tnf30', 'alb8', 'unidab f11a',
        )
        if re.findall(pattern, self._source_str(rec)) or \
            rec['molecule'] in  exact_mol:
            ids = self.set_ids(rec, 'single')
            return ids
        return []
    
    def antibody_chain(self, rec:dict):
        pattern_mol = r'antibody|scfv|immunoglobulin|nanobody|vhh17|fv4155|fab'
        if re.findall(pattern_mol, rec['molecule']) or \
            'fv fragment' in rec.get('fragment', ''):
            return True
        return False
