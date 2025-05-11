'''
'''
import os
import re
import pickle

class AnalyzeComplex:
    def __init__(self, pickle_file, max_span:int=None, verbose:bool=True):
        self.max_span = max_span if max_span else 30
        self.verbose = True if verbose else False
        # load data
        self.load_data(pickle_file)

    def load_data(self, infile):
        '''
        self.data
        self.outdir
        '''
        self.outdir = os.path.dirname(infile)
        with open(infile, 'rb') as f:
            self.data = pickle.load(f)
        if self.verbose:
            print('number of complex', len(self.data))

    def print_meta(self):
        '''
        input pickle: freesasa.p, or dis.p
        print pairwise chain type and chain id
        '''
        for _data in self.data:
            print('key1', list(_data))
            k1, k2 = _data['type'].split('-')
            c1, c2 = _data['chains'][:2]
            print(f'{k1} chain {c1[0]} ~ {k2} chain {c2[0]}')
            print('---'*10)
    
    def cal_span(self, a, b):
        '''
        input pickle: freesasa.p, or dis.p
        '''
        if isinstance(a, str):
            a = int(re.findall(r'\d+', a)[0])
        if isinstance(b, str):
            b = int(re.findall(r'\d+', b)[0])
        return b - a

    def build_dataset(self, min_fragment:int=1):
        '''
        input pickle: sasa_seq.p, or dis_seq.p
        retrieve dataset from machine learning
        '''
        ds_data = []
        for _data in self.data:
            rec = {}
            for key, seq_rec in _data.items():
                seq_list = [i['seq'] for i in seq_rec if \
                    len(i['seq']) >= min_fragment]
                # take antigen or receptor as input
                if key == 'receptor':
                    rec['input'] = ' | '.join(seq_list)
                    rec['input_label'] = key
                else:
                    # take heavy or light chain or single domain 
                    #  chain of antibody as output
                    rec['output'] = ' | '.join(seq_list)
                    rec['output_label'] = key
            # both input and output are existing
            if rec['input'] and rec['output']:
                ds_data.append(rec)
        return ds_data
