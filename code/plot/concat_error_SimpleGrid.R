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
out <- fread(args[1])
tmp_name <- strsplit(args[1], "/")[[1]]
out$rep <- tmp_name[4]
out$gwas_size <- as.numeric(strsplit(tmp_name[5], "-")[[1]][2])
out$test_size <- as.numeric(strsplit(tmp_name[6], "-")[[1]][2])
out$test <- tmp_name[7]
tmp1 <- strsplit(tmp_name[8], "-")[[1]][2]
out$L <- strsplit(tmp1, ".txt")[[1]][1]

# Loops through all other files and add onto output
for (i in 2:(nfile-1)) {

  tmp <- fread(args[i])
  tmp_name <- strsplit(args[i], "/")[[1]]
  tmp$rep <- tmp_name[4]
  tmp$gwas_size <- as.numeric(strsplit(tmp_name[5], "-")[[1]][2])
  tmp$test_size <- as.numeric(strsplit(tmp_name[6], "-")[[1]][2])
  tmp$test <- tmp_name[7]
  tmp1 <- strsplit(tmp_name[8], "-")[[1]][2]
  tmp$L <- strsplit(tmp1, ".txt")[[1]][1]
  out <- rbind(out, tmp)

}

# Save output
fwrite(out, args[nfile] ,row.names=F,quote=F,sep="\t", col.names = T)


