#!/usr/bin/env python3

import pandas as pd
import argparse
import os
from functools import reduce
import glob


def extractSNPs(pathtoSNPtables):

    # takes multiple SNP tables to generate a list of all SNPs across samples and common SNPs only

    all_files = glob.glob(os.path.join(pathtoSNPtables, '*.tab'))
    clade = os.path.basename(os.path.dirname(pathtoSNPtables).split("_")[-2])
    list_of_dfs = []

    for file in all_files:
        df = pd.read_csv(file, sep='\t', usecols=['POS', 'ALT'], index_col='POS')
        filename = os.path.basename(file)
        colname = filename.split("_")[0]
        df.rename({'ALT': colname}, axis=1, inplace=True)
        list_of_dfs.append(df)

    allSNP_df = reduce(lambda left,right: pd.merge(left,right,on='POS',how='outer'), list_of_dfs)
    allSNP_df.sort_index(inplace=True)
    variable = allSNP_df.apply(lambda row: len(row.unique()) != 1, axis=1)
    allSNP_df['Variable'] = variable
    numsamples = len(allSNP_df.columns)
    varCount = allSNP_df.count(axis=1)
    allSNP_df['Majority'] = varCount/numsamples

    # output SNPs which vary within a clade    
    varSNP_df = allSNP_df[allSNP_df['Variable'] == True]
    varSNP_df.to_csv('{}_varSNPs.csv'.format(clade))

    # output SNPs which are constant within a clade
    commonSNP_df = allSNP_df[allSNP_df['Variable'] == False]
    commonSNP_df.columns.values[0] = clade
    cladeSNP_df = commonSNP_df.iloc[:, 0]
    cladeSNP_df.to_csv('{}_cladeSNPs.csv'.format(clade))

    # output SNPs which are represented by more than 99% of samples in a clade
    # change to 0.95 (95%) for B6-11, B6-53 and Microti
    signifSNP_df = allSNP_df[allSNP_df['Majority'] >= 0.99]
    signifSNP_df.fillna(signifSNP_df.mode(axis=1))
    signifSNP_df.columns.values[0] = clade
    cladesigSNP_df = signifSNP_df.iloc[:, 0]
    cladesigSNP_df.to_csv('{}_cladesigSNPs.csv'.format(clade))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('pathtoSNPtables')

    args = parser.parse_args()

    extractSNPs(**vars(args))
