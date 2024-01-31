# This script takes in a specified demographic model and generates genotype data

import argparse

parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

req_grp.add_argument("--sample_size","-s",dest="ss",help="diploid sample size within each deme",type=int,required=True)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
req_grp.add_argument("--chr","-c",dest="chr",help="chromosome number",type=str,required=True)
parser.add_argument("--Ne","-N",dest="npop",help="effective pop size (def:1e4)",type=int,default=10000,nargs="?")
parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:1e7)",type=int,default=10000000,nargs="?")
parser.add_argument("--rho","-r",dest="rho",help="recombination rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--migrate","-m",dest="migrate",help="migration rate (def:0.05)",type=float,default=0.05,nargs="?")
parser.add_argument("--ndemes","-d",dest="ndemes",help="number of demes - must be a squared number (e.g. 25 or 100",type=int,default=36,nargs="?")
parser.add_argument("--tmove","-t",dest="tmove",help="time (g) to panmixia. can be -9 (inf) or any positive integer",type=int,default=100,nargs="?")
parser.add_argument("--seed","-sd",dest="seed",help="seed for simulation",type=int,default=1,nargs="?")
parser.add_argument("--target_n_snps","-nsnp",dest="nsnp",help="how many independent snps do you want",type=int,default=1,nargs="?")
args=parser.parse_args()

print(args)

import msprime
import numpy as np
import statistics
import math
import allel
import matplotlib.pyplot as plt
import statsmodels.api as sm
import pandas as pd
from scipy import (stats,ndimage)
import os
import concurrent.futures

# move time 
tmove = args.tmove

#number of demes
d=args.ndemes

#dimension of square grid
dim=int(np.sqrt(d))
print(dim)

#define 2d grid to with deme identity
pmat=np.arange(0,d).reshape(dim,dim)

#define function to generate adjacency matrix
#arguments:
#m = migration rate in one direction
#nd = number of demes
def step_mig_mat(m,nd):
    #m is the uni-directional symmetric migration rate
    #NOTE: nd MUST be a squared number
    if(math.sqrt(nd).is_integer()==False):
        raise ValueError("nd must be a squared number (e.g. 4, 9, 16 ...) for the 2D model")
    else:
        nd2=int(math.sqrt(nd))
        #create matrix which will be used to determine which cells are adjacent in 2-dimensions
        #diagonals not considered for now but can be incorporated later if needed
        pmat=np.arange(0,nd).reshape(nd2,nd2)

        #create empty migration matrix to be filled in. This will be the output
        mmat=np.zeros(shape=[nd,nd])

        #go through each cell in pmat and find out which cells are adjacent
        #first define functions to take care of corners and sides
        def contain(ix,max_ix):
            if ix<0:
                return(0)
            if ix>(max_ix-1):
                return(max_ix-1)
            else:
                return(ix)

        for ii in range(nd):
            center_ix=np.where(pmat==ii)
            top_ix=pmat[contain(center_ix[0]-1,nd2),contain(center_ix[1],nd2)]
            bottom_ix=pmat[contain(center_ix[0]+1,nd2),contain(center_ix[1],nd2)]
            left_ix=pmat[contain(center_ix[0],nd2),contain(center_ix[1]-1,nd2)]
            right_ix=pmat[contain(center_ix[0],nd2),contain(center_ix[1]+1,nd2)]

            mmat[ii,top_ix]=mmat[ii,bottom_ix]=mmat[ii,left_ix]=mmat[ii,right_ix]=m
            mmat[top_ix,ii]=mmat[bottom_ix,ii]=mmat[left_ix,ii]=mmat[right_ix,ii]=m

            mmat[ii,ii]=0

    return(mmat)

#generate migration matrix with migration rate provided by user
mig_mat=step_mig_mat(m=args.migrate,nd=d)


#diploid sample size within each deme
ss=args.ss


#N is the population size for each deme
#ss_each is the haploid sample size for each deme
#l is the length of the chromosome
#tmove is the number of generations past which all lineages are moved into one deme.
#The is to reduce computational time when the no. of lineages << ndemes
#also to mimic migration of an ancient population after which structure is established
#set to 1000 generations by default

ss_each=2*ss
sample_sizes=[ss_each]*d
population_configurations = [msprime.PopulationConfiguration(sample_size=k) for k in sample_sizes]


def make_tree(population_configurations, tmove, mmat):
    
    if tmove==-9:
         ts=msprime.simulate(Ne=args.npop,
                          population_configurations=population_configurations,
                          migration_matrix=mmat,
                          length=args.length)
         ts=msprime.mutate(ts, rate=args.mu)                  
    else:
        
        # add an extra deme to move lineages to after tmove generations
        population_configurations.append(msprime.PopulationConfiguration(sample_size = 0))

        #no migration to or from this deme until tmove generations
        mmat = np.append(mmat, np.zeros( (1,d) ), axis = 0)
        mmat = np.append(mmat, np.zeros(( (d+1) ,1)), axis = 1)

        #specify demographic event - move all lineages to deme d+1 after tmove generations
        demog=[
            msprime.MassMigration(
                time=tmove,
                source=i,
                destination=d,
                proportion=1.0) for i in range(d)]

        demog.append(#change migration rate among demes to be 0
            msprime.MigrationRateChange(
                time=tmove,
                rate=0))

        ts=msprime.simulate(Ne=args.npop,
                              population_configurations=population_configurations,
                              migration_matrix=mmat,
                              length=args.length,
                           demographic_events=demog)
        ts=msprime.mutate(ts, rate=args.mu)                    

    return(ts)


#simulate!
# Set up output containers 
output = []
header = []

# Function to Simulate SNPs 
def simulate_snp(i):
    
    #if i % 100 == 0:
    print("Simulating ~SNP {}".format(i))
    
    # Make Tree 
    tree_sequence = make_tree(population_configurations, tmove, mig_mat)

    idx=0
    
    # Save to VCF
   # print("writing genotype to vcf file")
    with open(args.outpre + "_" + str(i) + ".vcf", "w") as vcf_file:
        tree_sequence.write_vcf(vcf_file, ploidy=2, contig_id=args.chr)
        
    # Read in VCF
    #print("reading vcf")
    #print(args.outpre + "_" + str(i) + ".vcf")
    with open(args.outpre + "_" + str(i) + ".vcf", "r") as file:
        
        if i == 0:
            # Read the first six lines of the file
            first_six_lines = [file.readline().strip() for _ in range(6)]
            header.extend(first_six_lines)
        
        # Read the ith line and append it to the list
        target_line=None
        for line_number, line in enumerate(file, start=1):
            if line_number == idx + 7:
                target_line = line.strip()
                fields = target_line.strip().split("\t")
                fields[1] = str(i)
                fields[0] = str(args.chr)
                target_line = "\t".join(fields)
                break
        
    # Remove VCF file
    os.remove(args.outpre + "_" + str(i) + ".vcf")
    
    return target_line


# Process snps conccurently 
def process_snp_chunk(chunk):
    with concurrent.futures.ThreadPoolExecutor() as executor:
        return list(executor.map(simulate_snp, chunk))

# Chunk size (100 SNPs at a time)
chunk_size = 1

# Process SNPs in chunks
for chunk_start in range(0, args.nsnp, chunk_size):
    chunk_end = min(chunk_start + chunk_size, args.nsnp)
    snp_chunk = range(chunk_start, chunk_end)

    # Set up concurrent execution for each chunk
    chunk_results = process_snp_chunk(snp_chunk)
    output.extend(chunk_results)


print("Writing VCF")
# Write outout vcf, removing trees with no snps 
output = list(filter(lambda x: x is not None, output))
header.extend(output)
with open(args.outpre+".vcf","w") as vcf_file:
     for item in header:
            vcf_file.write(str(item) + '\n')



#write deme id/subpopulation for each individual
deme_id=[[i]*ss for i in range(0,d)]
#flatten
deme_id=[item for sublist in deme_id for item in sublist]

#write longitude and latitude for each individual
bplace_x=[]
bplace_y=[]
for i in range(0,dim):
    bplace_y.extend( [i] * ss * dim )
    bplace_x.extend([item for item in range(0,dim) for i in range(ss)])

#fid and iid
fid=["tsk_"+str(i) for i in range(0,(ss*d))]
iid=["tsk_"+str(i) for i in range(0,(ss*d))]

popdf=pd.DataFrame({"FID":fid,
                  "IID":iid,
                  "POP":deme_id,
                  "Latitude":bplace_y,
                  "Longitude":bplace_x})

popdf.to_csv(args.outpre+".pop",sep="\t",header=False,index=False)
