'''
'''
import os
import re
import pickle

from src.analyze_complex import AnalyzeComplex

class AnalyzeDist(AnalyzeComplex):
    def __init__(self, pickle_file, max_span:int=None, \
        dist_cutoff:int=None, verbose:bool=None):
        super().__init__(pickle_file, max_span, verbose)
        self.cutoff = dist_cutoff if dist_cutoff else 5

    def print_min_dist(self, ix:int):
        '''
        '''
        dist_data = self.data[ix]['dist']
        # print(dist_data['type'])

        g1, g2 = dist_data['dist1'], dist_data['dist2']
        print('---'*10)
        sig1 = g1[g1['value']<=self.cutoff]
        print(g1)
        print('---'*10)
        sig2 = g2[g2['value']<=self.cutoff]
        print(sig2)

    def retrieve_seq(self):
        pairs = []
        for _data in self.data:
            k1, k2 = _data['type'].split('-')
            # dist df
            dist1 = _data['dist']['dist1']
            dist2 = _data['dist']['dist2']
            # seq
            motif1 = self.parse_seq(dist1)
            motif2 = self.parse_seq(dist2)
            if motif1 and motif2:
                pairs.append((k1, motif1, k2, motif2))
        # export
        cond = f"dist-{self.cutoff}_span-{self.max_span}"
        outfile = os.path.join(self.outdir, f'seq_dist_{cond}.p')
        if self.verbose:
            print(f"Sequence-Distance data: {outfile}.")
        with open(outfile, 'wb') as f:
            pickle.dump(pairs, f, protocol=pickle.HIGHEST_PROTOCOL)
        return pairs, outfile

    def parse_seq(self, df):
        '''
        '''
        # values less than 5 in default
        sig = df[df['value'] <= self.cutoff]

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
            sub = df.loc[start:end,:]
            seq = ''.join(sub['aa'])
            motif.append({
                'seq': seq,
                'sig_res': len(g),
                'start': g[0],
                'end': g[-1],
                'dist_cutoff': self.cutoff,
                'max_span': self.max_span,
            })
        return motif

    def iter_data(self):
        for _data in self.data:
            if 'dist' in _data:
                # distance df
                dist1 = _data['dist']['dist1']
                dist2 = _data['dist']['dist2']
                yield dist1, dist2

    def get_spans(self):
        res1, res2 = [], []
        for dist1, dist2 in self.iter_data():
            sig1 = dist1[dist1['value']<=self.cutoff]['res_no']
            sig2 = dist2[dist2['value']<=self.cutoff]['res_no']
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
        for dist1, dist2 in self.iter_data():
            sig1 = list(dist1[dist1['value']<=self.cutoff]['res_no'])
            sig2 = list(dist2[dist2['value']<=self.cutoff]['res_no'])

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
        res = {}
        for _data in self.data:
            df1, df2 = None, None
            chain_id = (_data['chains'][0][0], _data['chains'][1][0])
            if _data['type'] == type_name:
                df1 = _data['dist']['dist1']
                df2 = _data['dist']['dist2']
            if (df1 is not None) and (df2 is not None):
                if chain_len is None or (chain_len and chain_len == len(df1)):
                    res[chain_id] = (df1, df2)
        return res
    
    def get_receptor_data(self, chain_len:int=None):
        '''
        '''
        res = {}
        for _data in self.data:
            chain_id = (_data['chains'][0][0], _data['chains'][1][0])
            df = _data['dist']['dist2']
            if chain_len is None or (chain_len and chain_len == len(df)):
                res[chain_id] = df
        return res