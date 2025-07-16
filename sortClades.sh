#!/bin/bash

# get list of samples in clades

set -eo pipefail

cladelist=$1
clade=$2

# Separate sample list based on clade and outcome
#awk -F, '{print >> ($9"_"$7".csv")}' $clean

# pull snp tables into folder
mkdir ${clade}_snps

while IFS=, read -r Submission Sample GenomeCov MeanDepth pcMapped group Ncount Path;
do
    aws s3 cp "${Path}snpTables/${Sample}_snps.tab" "${clade}_snps/${Sample}_snps.tab";
done < <(tail -n +2 $cladelist)
