#!/bin/bash
#SBATCH --job-name=main
#SBATCH --output=logs/main.out
#SBATCH --error=logs/main.err
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=2000
#SBATCH --partition=tier1q

module load gcc/12.1.0
module load python/3.10.5
module load plink/2.0
module load gcta/1.94.1
module load R
module load bcftools
module load samtools

echo "SLURM_JOBID="$SLURM_JOBID
cat snakefile_4PopSplit
snakemake -s snakefile_4PopSplit --unlock
snakemake --profile cluster-setup/ -s snakefile_4PopSplit  --rerun-incomplete
#snakemake -j1 --rerun-incomplete -s snakefile_main 

#snakemake --cores all --rerun-incomplete -s snakefile_UKBB    
#snakemake -j1 --rerun-incomplete 
