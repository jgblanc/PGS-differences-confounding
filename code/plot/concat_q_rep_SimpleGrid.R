# This script takes a setof effect sizes and test vector and computes Q

args=commandArgs(TRUE)

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

nfile = length(args)
print(nfile)

# Read in first file and collect all parameters
out <- fread(args[1])[,2:10]
tmp_name <- strsplit(args[1], "/")[[1]]
out$rep <- tmp_name[4]
out$gwas_size <- as.numeric(strsplit(tmp_name[5], "-")[[1]][2])
out$test_size <- as.numeric(strsplit(tmp_name[6], "-")[[1]][2])
out$h2 <- as.numeric(strsplit(tmp_name[7], "-")[[1]][2])
out$num_causal <- as.numeric(strsplit(tmp_name[8], "-")[[1]][2])
out$env <- as.numeric(strsplit(tmp_name[9], "_")[[1]][2])
out$pheno <- tmp_name[10]
out$test <- tmp_name[11]

# Loops through all other files and add onto output
for (i in 2:(nfile-1)) {

  tmp <- fread(args[i])[, 2:10]
  tmp_name <- strsplit(args[i], "/")[[1]]
  tmp$rep <- tmp_name[4]
  tmp$gwas_size <- as.numeric(strsplit(tmp_name[5], "-")[[1]][2])
  tmp$test_size <- as.numeric(strsplit(tmp_name[6], "-")[[1]][2])
  tmp$h2 <- as.numeric(strsplit(tmp_name[7], "-")[[1]][2])
  tmp$num_causal <- as.numeric(strsplit(tmp_name[8], "-")[[1]][2])
  tmp$env <- as.numeric(strsplit(tmp_name[9], "_")[[1]][2])
  tmp$pheno <- tmp_name[10]
  tmp$test <- tmp_name[11]
  out <- rbind(out, tmp)

}

# Save output
fwrite(out, args[nfile] ,row.names=T,quote=F,sep="\t", col.names = T)


