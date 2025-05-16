'''
'''
import os
import re
import pickle

from src.analyze_complex import AnalyzeComplex

class AnalyzeSasa(AnalyzeComplex):
    def __init__(self, pickle_file, delta_threshold:float=None, \
        max_span:int=None, verbose:bool=None):
        super().__init__(pickle_file, max_span, verbose)
        self.delta_threshold = delta_threshold if delta_threshold else 0

    def print_delta_sasa(self, ix:int):
        '''
        print delta SASA
        '''
        k1, k2 = self.data[ix]['type'].split('-')
        c1, c2 = self.data[ix]['chains'][:2]
        print(f'complex: {k1} chain {c1[0]} ~ {k2} chain {c2[0]}')


        print('---'*10)
        print(f'{k1} chain:')
        delta1 = self.data[ix]['delta_total_sasa']['delta1']
        sig1 = delta1[delta1['value'] > self.delta_threshold]
        print(sig1)

        print('---'*10)
        print(f'{k2} chain:')
        delta2 = self.data[ix]['delta_total_sasa']['delta2']
        sig2 = delta2[delta2['value'] > self.delta_threshold]
        print(sig2)
    
    def retrieve_seq(self):
        pairs = []
        for _data in self.data:
            k1, k2 = _data['type'].split('-')
            # delta df
            delta1 = _data['delta_total_sasa']['delta1']
            delta2 = _data['delta_total_sasa']['delta2']
            # seq
            motif1 = self.parse_seq(delta1)
            motif2 = self.parse_seq(delta2)
            pairs.append((k1, motif1, k2, motif2))
        # export
        cond = f"delta-{self.delta_threshold}_span-{self.max_span}"
        outfile = os.path.join(self.outdir, f'seq_sasa_{cond}.p')
        if self.verbose:
            print(f"Sequence-SASA data: {outfile}.")
        with open(outfile, 'wb') as f:
            pickle.dump(pairs, f, protocol=pickle.HIGHEST_PROTOCOL)
        return pairs, outfile

    def parse_seq(self, delta):
        '''
        '''
        sig = delta[delta['value'] > self.delta_threshold]

        # group significant point
        group, curr = [], []
        for ix, res_no in sig['res_no'].items():
            if curr: 
                span = self.cal_span(curr[-1][1], res_no)
                if span <= self.max_span:
                    curr.append((ix, res_no))
                else:
                    group.append(curr)
                    curr = [(ix, res_no),]
            else:
                curr.append((ix, res_no))
        else:
            if curr:
                group.append(curr)

        # convert to seq
        motif = []
        for g in group:
            start, end = g[0][0], g[-1][0]
            sub = delta.loc[start:end,:]
            seq = ''.join(sub['aa'])
            motif.append({
                'seq': seq,
                'sig_res': len(g),
                'start': g[0],
                'end': g[-1],
                'delta_sasa_cutoff': self.delta_threshold,
                'max_span': self.max_span,
            })
        return motif

    def get_spans(self):
        res1, res2 = [], []
        for _data in self.data:
            # delta df
            delta1 = _data['delta_total_sasa']['delta1']
            delta2 = _data['delta_total_sasa']['delta2']
            # 
            sig1 = delta1[delta1['value']>0]['res_no']
            sig2 = delta2[delta2['value']>0]['res_no']
            #spans
            span1 = [self.cal_span(a, b) for a, b in zip(sig1[:-1], sig1[1:])]
            span2 = [self.cal_span(a, b) for a, b in zip(sig2[:-1], sig2[1:])]
            # filter
            span1 = [i for i in span1 if i > 1]
            span2 = [i for i in span2 if i > 1]
            res1.append(span1)
            res2.append(span2)
        return res1, res2
    
    def get_fragments(self):
        res1, res2 = [], []
        for _data in self.data:
            # delta df
            delta1 = _data['delta_total_sasa']['delta1']
            delta2 = _data['delta_total_sasa']['delta2']
            # 
            sig1 = list(delta1[delta1['value']>0]['res_no'])
            sig2 = list(delta2[delta2['value']>0]['res_no'])
            # 
            frag1, frag1_len = [], 1
            for ix in range(0, len(sig1)-1):
                span = self.cal_span(sig1[ix],sig1[ix+1])
                if span > 1:
                    frag1.append(frag1_len)
                    frag1_len = 1
                else:
                    frag1_len += 1
            else:
                frag1.append(frag1_len)
                frag1_len = 1
            # 
            frag2, frag2_len = [], 1
            for ix in range(0, len(sig2)-1):
                span = self.cal_span(sig2[ix],sig2[ix+1])
                if span > 1:
                    frag2.append(frag2_len)
                    frag2_len = 1
                else:
                    frag2_len += 1
            else:
                frag2.append(frag2_len)
                frag2_len = 1
            # 
            frag1 = [i for i in frag1 if i>1]
            frag2 = [i for i in frag2 if i>1]
            res1.append(frag1)
            res2.append(frag2)
        return res1, res2

    def get_chain_data(self, type_name:str, chain_len:int=None):
        '''
        '''
        res = {}
        for _data in self.data:
            df1, df2 = None, None
            chain_id = (_data['chains'][0][0], _data['chains'][1][0])
            if _data['type'] == type_name:
                df1 = _data['delta_total_sasa']['delta1']
                df2 = _data['delta_total_sasa']['delta2']
            if df1 is not None and df2 is not None:
                if chain_len is None or (chain_len and chain_len == len(df1)):
                    res[chain_id] = (df1, df2)
        return res

    def get_receptor_data(self, chain_len:int=None):
        '''
        '''
        res = {}
        for _data in self.data:
            chain_id = (_data['chains'][0][0], _data['chains'][1][0])
            df = _data['delta_total_sasa']['delta2']
            if chain_len is None or (chain_len and chain_len == len(df)):
                res[chain_id] = df
        return res