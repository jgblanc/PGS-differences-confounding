#!/usr/bin/env Rscript

args=commandArgs(TRUE)

if( length(args) != 5){stop("Usage: <causal effects file> <glm.linear file> <number of snps> <output file prefix causal> <output file prefix noncausal>") }

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(R.utils)
}))

geffects_file=args[1] # Causal SNP IDs
gwas_file=args[2] # GWAS file
num_snps = as.numeric(args[3])
output_file_causal=args[4]
output_file_noncausal=args[5]


# Function to get the effect size for the T allele
flip_effect = function(gwas_df,beta_colname){
  gwas_df = gwas_df[A1=="A", beta_colname := -BETA]
  gwas_df = gwas_df[A1=="T", beta_colname := BETA]
  gwas_df$A1="T"
  gwas_df = gwas_df[,.(CHROM,POS,ID,A1,beta_colname,P)]
  colnames(gwas_df)[5] = beta_colname
  return(gwas_df)
}



# Read in list of causal variants and clean data frame
causal=fread(geffects_file)
colnames(causal)=c("rsid","allele","esize")
causal=causal%>%
  separate(rsid,into=c("CHROM","POS","ref","alt"),sep="_",remove=F)
causal$POS=as.numeric(causal$POS)
causal$CHROM=as.numeric(causal$CHROM)

# Read in GWAS results
gwas1=fread(gwas_file,fill=T)
gwas1 <- gwas1 %>% select("Chr", "SNP", "bp", "A1","b","p")
colnames(gwas1) <- c("CHROM", "ID", "POS", "A1", "BETA1", "P")
print(head(gwas1))

# Flip to get correct effect sizes
gwas1 = flip_effect(gwas1,beta_colname = "BETA1")

# Select effect sizes for causal variants
gwas1.1 = gwas1[ID%in%causal$rsid]
gwas.causal=gwas1.1[,c("ID","A1","BETA1")]
colnames(gwas.causal) <- c("ID", "A1", "BETA_Strat")

# Write Beta hat for causal effects
fwrite(gwas.causal, output_file_causal,
       col.names=T,row.names=F,quote=F,sep="\t")


# Select the top associated variants
gwas.red <- gwas1 %>%
  slice_min(order_by = P, n = num_snps)
gwas.red=gwas.red[,c("ID","A1","BETA1")]
colnames(gwas.red) <- c("ID", "A1", "BETA_Strat")


# Write Beta hat from p-value
fwrite(gwas.red,
       output_file_noncausal,
       col.names=T,row.names=F,quote=F,sep="\t")


