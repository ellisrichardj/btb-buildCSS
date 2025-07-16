#!/usr/bin/env python3

import pandas as pd
import glob
import os
import argparse
from functools import reduce


def defineCSS(pathtocladeSNPs):

    clade_files = glob.glob(os.path.join(pathtocladeSNPs, '*cladesigSNPs.csv'))
    list_of_dfs = []

    for file in clade_files:
        df = pd.read_csv(file, index_col='POS')
        list_of_dfs.append(df)

    allCSS_df = reduce(lambda left,right: pd.merge(left,right,on='POS',how='outer'), list_of_dfs)
    allCSS_df.sort_index(inplace=True)
    allCSS_df.sort_index(axis=1, inplace=True)
    monom = allCSS_df.nunique(axis=1) # outputs '1' if the polymorphism is consistent  
    numclades = allCSS_df.count(axis=1) # counts number of clades with SNP
    allCSS_df['Monomorphic'] = monom
    allCSS_df['Numclades'] = numclades
    allCSS_df.to_csv('allCSS.csv')

    # TODO need to add in somthing to capture SNPs in all except B6-16
    refCSS_df = allCSS_df.loc[(allCSS_df['Numclades'].ge(41)) & (allCSS_df['B6-16'].isna())]
    refCSS_df.to_csv('refCSS.csv')


    uniqCSS_df = allCSS_df[(allCSS_df['Numclades'] == 1) | 
                           ((allCSS_df['Numclades'].ge(41)) & (allCSS_df['B6-16'].isna()))]
    uniqCSS_df.to_csv('uniqCSS.csv')
    #pos_df = uniqCSS_df['POS']
    uniqCSS_df.to_csv('snppos.txt', columns=[], header=False)

    clade_dfs = [uniqCSS_df[[col]] for col in uniqCSS_df.columns]
    for x in clade_dfs:
        #x.dropna(inplace=True)
        clade = x.columns[-1]
        x.to_csv('{}.csv'.format(clade))

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pathtocladeSNPs')

    args = parser.parse_args()

    defineCSS(**vars(args))