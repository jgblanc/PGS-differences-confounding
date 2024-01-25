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
name <- colnames(out)
tmp_name <- strsplit(args[1], "/")[[1]]
out$rep <- tmp_name[4]
out$config <- tmp_name[5]
out$gwas_size <- as.numeric(strsplit(tmp_name[6], "-")[[1]][2])
out$test_size <- as.numeric(strsplit(tmp_name[7], "-")[[1]][2])
tmp <- strsplit(tmp_name[8], "-")[[1]][3]
print(tmp)
out$L <- as.numeric(strsplit(tmp, ".txt")[[1]][1])
out$type <- name
colnames(out)[1] <- "error"


# Loops through all other files and add onto output
for (i in 2:(nfile-1)) {

  tmp <- fread(args[i])
  name <- colnames(tmp)
  tmp_name <- strsplit(args[1], "/")[[1]]
  tmp$rep <- tmp_name[4]
  tmp$config <- tmp_name[5]
  tmp$gwas_size <- as.numeric(strsplit(tmp_name[6], "-")[[1]][2])
  tmp$test_size <- as.numeric(strsplit(tmp_name[7], "-")[[1]][2])
  tmp2 <- strsplit(tmp_name[8], "-")[[1]][3]
  tmp$L <- as.numeric(strsplit(tmp2, ".txt")[[1]][1])
  tmp$type <- name
  colnames(tmp)[1] <- "error"

  out <- rbind(out, tmp)

}

# Save output
fwrite(out, args[nfile] ,row.names=T,quote=F,sep="\t", col.names = T)


