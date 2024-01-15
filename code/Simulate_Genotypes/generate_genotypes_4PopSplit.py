# This script takes in a specified demographic model and generates genotype data
## Adapted from: https://elifesciences.org/articles/61548 


# Import modules 
import msprime
import argparse
import pandas as pd
import numpy as np
import os
import concurrent.futures

# Parse Inputs 
parser=argparse.ArgumentParser()
req_grp=parser.add_argument_group(title="Required arguments")

parser.add_argument("--NA","-A",dest="NA",help="diploid effective population size in population A",type=int,default=10000)
parser.add_argument("--NB","-B",dest="NB",help="diploid effective population size in population B",type=int,default=10000)
parser.add_argument("--NC","-C",dest="NC",help="diploid effective population size in population C",type=int,default=10000)
parser.add_argument("--ND","-D",dest="ND",help="diploid effective population size in population A",type=int,default=10000)
parser.add_argument("--Nanc","-anc",dest="Nanc",help="diploid effective ancestral population size before the first split",type=int,default=10000)
parser.add_argument("--sample_size_A","-a",dest="sample_A",help="haploid sample size in population A (must be divisible by ploidy)",type=int,default=50)
parser.add_argument("--sample_size_B","-b",dest="sample_B",help="haploid sample size in population B (must be divisible by ploidy)",type=int,default=50)
parser.add_argument("--sample_size_C","-c",dest="sample_C",help="haploid sample size in population C (must be divisible by ploidy)",type=int,default=50)
parser.add_argument("--sample_size_D","-d",dest="sample_D",help="haploid sample size in population D (must be divisible by ploidy)",type=int,default=50)
req_grp.add_argument("--outpre","-o",dest="outpre",help="output file prefix",type=str,required=True)
parser.add_argument("--length","-L",dest="length",help="length of chromosome (bp) (def:1e7)",type=int,default=10000000,nargs="?")
parser.add_argument("--rho","-r",dest="rho",help="recombination rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--mu","-u",dest="mu",help="mutation rate (def:1e-08)",type=float,default=1e-08,nargs="?")
parser.add_argument("--split1","-s1",dest="split_time1",help="time when 1 pop splits to 2",type=int,default=140e3,nargs="?")
parser.add_argument("--split2","-s2",dest="split_time2",help="time when 2 pops split to 4",type=int,default=70e3,nargs="?")
parser.add_argument("--ploidy","-p",dest="ploidy",help="ploidy of individuals",type=int,default=2,nargs="?")
req_grp.add_argument("--chr","-ch",dest="chrom_num",help="number of chromosome",type=int,required=True)
parser.add_argument("--target_n_snps","-nsnp",dest="nsnp",help="how many independent snps do you want",type=int,default=1,nargs="?")
args=parser.parse_args()

print(args)


# Define Function to generate tree sequency under 4 population split model - single population splits into 2 at `split_time_1` and each of those populations splits into 2 at `split_time_2`
def split(N_A, N_B, N_C, N_D, split_time1, split_time2, sample_A, sample_B, sample_C, sample_D, seg_length, recomb_rate, mut_rate, N_anc, target_n_snps):

    # Times are provided in years, so we convert into generations.
    generation_time = 25
    T_S1 = split_time1 / generation_time
    T_S2 = split_time2 / generation_time


    # Population IDs correspond to their indexes in the population
    # configuration array. Therefore, we have 0=A, 1=B, 2=C, 3=D initially.
    population_configurations = [
        msprime.PopulationConfiguration(
            sample_size=sample_A, initial_size=N_A),
        msprime.PopulationConfiguration(
            sample_size=sample_B, initial_size=N_B), 
        msprime.PopulationConfiguration(
            sample_size=sample_C, initial_size=N_C),
        msprime.PopulationConfiguration(
            sample_size=sample_D, initial_size=N_D)
    ]

    demographic_events = [
        msprime.MassMigration(
            time=T_S2, source=3, destination=2, proportion=1.0),
        msprime.MassMigration(
            time=T_S2, source=1, destination=0, proportion=1.0),
        msprime.MassMigration(
            time=T_S1, source=2, destination = 0, proportion=1.0),
        msprime.PopulationParametersChange(
            time=T_S1, initial_size=N_anc, growth_rate=0, population_id=0)
    ]
    
    # Use the demography debugger to print out the demographic history
    # that we have just described.
    #dd = msprime.DemographyDebugger(
    #    population_configurations=population_configurations,
    #    demographic_events=demographic_events)
    #dd.print_history()
    
    ts = msprime.simulate(population_configurations=population_configurations,
                         demographic_events=demographic_events, length=seg_length, 
                         num_replicates=target_n_snps, mutation_rate=mut_rate)
                         
    return ts

# Generate Tree Sequences and Save to VCF

print("simulating genotypes under demographic model")
ts = split(args.NA, args.NB, args.NC, args.ND, args.split_time1, args.split_time2, args.sample_A, args.sample_B, args.sample_C, args.sample_D, 
args.length, args.rho, args.mu, args.Nanc, args.nsnp)


output = []
header = []

def process_tree_sequence(i, tree_sequence):

    if i % 1000 == 0:
        print("Simulating ~SNP {}".format(i))
    
    # Extract Haplotype
    H = tree_sequence.genotype_matrix()
    p, n = H.shape
    
    
    # Get SNP global frequencies
    snp_frequencies = np.sum(H, axis=1) / (n)
    
    # Find singleton frequencies 
    threshold = 1 / n
    selected_indices = np.where(snp_frequencies > threshold)[0]
     
    # Select random SNP
    if len(selected_indices) > 0:
        idx = np.random.choice(selected_indices, 1)

    # Save to VCF
    #print("writing genotype to vcf file")
    with open(args.outpre + "_" + str(i) + ".vcf", "w") as vcf_file:
        tree_sequence.write_vcf(vcf_file, ploidy=args.ploidy, contig_id=1 + i)
        
    # Read in VCF
    #print("reading vcf")
    with open(args.outpre + "_" + str(i) + ".vcf", "r") as file:
        
        if i == 0:
            # Read the first six lines of the file
            first_six_lines = [file.readline().strip() for _ in range(6)]
            header.extend(first_six_lines)
        
        # Read the ith line and append it to the list
        for line_number, line in enumerate(file, start=1):
            if line_number == idx + 7:
                target_line = line.strip()
                output.append(target_line)
        
    # Remove VCF file
    os.remove(args.outpre + "_" + str(i) + ".vcf")


# Function to process a chunk of the loop
def process_chunk(chunk):
    for i, tree_sequence in chunk:
        process_tree_sequence(i, tree_sequence)


# Chunk the tree sequence for parallel processing
chunk_size = 1000  # You can adjust this based on the optimal size for your task
chunks = [(i, tree_sequence) for i, tree_sequence in enumerate(ts)]
chunked_data = [chunks[i:i + chunk_size] for i in range(0, len(chunks), chunk_size)]

# Use ThreadPoolExecutor for parallel execution
with concurrent.futures.ThreadPoolExecutor() as executor:
    all_results = executor.map(process_chunk, chunked_data)

# Write outout vcf 
header.extend(output)
with open(args.outpre+".vcf","w") as vcf_file:
     for item in header:
            vcf_file.write(str(item) + '\n')
            
# Write population information file (population identity) 

#write population for each individual
deme_id=[]
for i in range(0, int(float(args.sample_A)/float(args.ploidy))):
    deme_id.append("A")
for i in range(0, int(float(args.sample_B)/float(args.ploidy))):
    deme_id.append("B")
for i in range(0, int(float(args.sample_C)/float(args.ploidy))):
    deme_id.append("C")
for i in range(0, int(float(args.sample_D)/float(args.ploidy))):
    deme_id.append("D")
    
#flatten
deme_id=[item for sublist in deme_id for item in sublist]

#fid and iid
fid=["tsk_"+str(i) for i in range(0,(int(float(args.sample_A)/float(args.ploidy)) + int(float(args.sample_B)/float(args.ploidy)) + int(float(args.sample_C)/float(args.ploidy)) + int(float(args.sample_D)/float(args.ploidy))))]
iid=["tsk_"+str(i) for i in range(0,(int(float(args.sample_A)/float(args.ploidy)) + int(float(args.sample_B)/float(args.ploidy)) + int(float(args.sample_C)/float(args.ploidy)) + int(float(args.sample_D)/float(args.ploidy))))]

popdf=pd.DataFrame({"FID":fid,
                  "IID":iid,
                  "POP":deme_id})

popdf.to_csv(args.outpre+".pop",sep="\t",header=True,index=False)
