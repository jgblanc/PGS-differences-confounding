#!/bin/bash
pfile_path=$1
outfile=$2
overlap_snps=$3
id=$4

plink2 \
    --pfile $pfile_path \
    --extract $overlap_snps \
    --geno-counts \
    --keep $id \
    --memory 100000 \
    --threads 16 \
    --out $outfile
