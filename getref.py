#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import pandas as pd

'''
Obtains reference base for a list of positions from a reference genome and
fills in blank spaces with this information - thus providing expected base for
each clade. 
'''

def getref(inputfasta, uniqCSS):

    refseq = SeqIO.read(inputfasta, "fasta")
    print(refseq.id)

    RefPos = []
    RefSNP = []

    CSStable = pd.read_csv(uniqCSS, index_col='POS')
    snppos = CSStable.index.to_list()
    
    for x in snppos:
        pos = int(x)
        pypos = pos-1
        RefPos.append(pos)
        RefSNP.append(refseq[pypos])
    
    output_df = pd.DataFrame({'POS': RefPos, 'Ref': RefSNP}).set_index('POS')
    output_df.to_csv('refcalls.csv')

    #CSStable = pd.read_csv(uniqCSS, index_col='POS')
    CSSwithref = CSStable.apply(lambda x: x.fillna(output_df['Ref']))
    CSSwithref.to_csv('CSSwithref.csv')


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('inputfasta')
    parser.add_argument('uniqCSS')

    args = parser.parse_args()

    getref(**vars(args))
    