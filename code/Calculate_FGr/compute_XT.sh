#!/bin/bash
pfile_path=$1
pheno_path=$2
outfile=$3
overlap_snps=$4
ids=$5

plink2 \
  --pfile $pfile_path \
  --extract $overlap_snps \
  --glm omit-ref \
  --pheno  $pheno_path \
  --pheno-name Tvec \
  --geno-counts \
  --keep $ids \
  --memory 100000 \
  --threads 16 \
  --out $outfile
