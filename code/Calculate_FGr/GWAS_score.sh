#!/bin/bash
pfile_path=$1
beta_path=$2
outfile=$3
overlap_snps=$4
id=$5

plink2 \
  --pfile $pfile_path \
  --extract $overlap_snps \
  --score $beta_path center header-read cols=dosagesum,scoresums \
  --keep $id \
  --memory 50000 \
  --threads 8 \
  --out $outfile
