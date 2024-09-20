#Class for TCR data

import pandas as pd
import re

def parse_imgt(gene):
    #Convert TRB1-1*00(xx) into TRB1-1 as per IMGT format
    gene = re.sub('\*00.+$','',gene)
    gene = re.sub('TCR','TR',gene)
    gene = re.sub('0(?=[1-9])','',gene)


    return gene

class TcrData:
    def __init__(self):
        self.raw = pd.DataFrame(columns=['tcrseq', 'cdr3', 'count','vgene', 'freq'])
        self.df = pd.DataFrame(columns=['tcrseq', 'cdr3', 'count', 'freq'])
        self.cdr3hash = dict()

    def read_mixcr(self,file, minFreq = 0):
        data = pd.read_csv(file, sep='\t', low_memory=False)

        data = data[data['aaSeqCDR3'].str.contains('^C[A-Z]{3,25}F$', regex=True, na=False)]

        data['vGeneName']=data['allVHitsWithScore'].apply(lambda vgene: parse_imgt(vgene))

        data['jGeneName']=data['allJHitsWithScore'].apply(lambda jgene: parse_imgt(jgene))


        self.raw['tcrseq'] = data['vGeneName'] + "\t" + data['aaSeqCDR3'] + "\t" + data['jGeneName']
        self.raw['vgene'] = data['vGeneName']
        self.raw['jgene'] = data['jGeneName']
        self.raw['cdr3'] = data['aaSeqCDR3']
        self.raw['count'] = data['cloneCount']
        self.raw['freq'] = data['cloneFraction']
        self.raw['freq'] = self.raw['freq']/self.raw['freq'].sum()

        self.raw = self.raw.dropna(axis=0)

        self.raw = self.raw[self.raw['count'] >= minFreq]
