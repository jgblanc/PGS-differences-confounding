#!/bin/bash
#SBATCH --job-name=grids
#SBATCH --output=logs/grids.out
#SBATCH --error=logs/grids.err
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
cat snakefile_SimpleGrid
snakemake -s snakefile_SimpleGrid --unlock
#snakemake -np -s snakefile_SimpleGrid 
snakemake --profile cluster-setup/ -s snakefile_SimpleGrid 

