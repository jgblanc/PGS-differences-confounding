#script takes in genetic values and generates phenotypes with a given amount of stratification

args=commandArgs(TRUE)

if(length(args)<5){stop("Rscript simulate_phenotypes_4PopSplit.R <genetic value file> <pop file> <output_file> <hertaibility> <environmental shift>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

gvalue_file = args[1] #genetic values
popfile = args[2] # pop file from msprime simulation
output_file = args[3] #name of output file
h2 = as.numeric(args[4]) #target heritability
env_s=as.numeric(args[5]) #total environmental shift

print(paste0("The stratification shift is ",as.numeric(env_s)))
print(paste0("The heritability is " , as.numeric(h2)))

# Load true genetic values
prs=fread(gvalue_file)
colnames(prs)<-c("IID","dosage","prs")
sample_size=nrow(prs)

#load file containing the pop ID for each inidividual
pop=fread(popfile,header=F)
colnames(pop)<-c("IID","FID","pop")

#add this info to prs file
prs=merge(prs, pop, by="IID", sort=F)

# Draw random environment
pops <- unique(prs$pop)
prs$env = rnorm(sample_size,0, sqrt(1 - h2))


# Add stratification effect to environment
prs <- prs %>% group_by(pop) %>% mutate(env = ifelse(pop == pops[1], env + env_s, env ))

# Add environmental effect to
prs = prs %>%
  mutate(pheno_strat = prs + env)

# Print heritability
print(paste("h2 :",
            round(cor(prs$prs,
                      prs$pheno_strat)^2,2)))

# Print mean phenotype in each pop
print(prs %>% group_by(pop) %>% summarise(mean(pheno_strat)))

# Write to output file
fwrite(
  prs%>%
    mutate(FID=IID)%>%
    ungroup() %>%
    select(FID,IID,pheno_strat),
  output_file,
  col.names=T,row.names=F,quote=F,sep="\t")
