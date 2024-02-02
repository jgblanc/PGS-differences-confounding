# This script computes FGr by doing GX'T

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript calc_FGr.R <plink prefix> <tvec file> <test ids> <gwas ids> <common snps> <out prefix>")}

suppressWarnings(suppressMessages({
  library(pgenlibr)
  library(data.table)
  library(dplyr)
  library(pracma)
}))

plink_prefix = args[1] # Prefix to plink files
tvec_file = args[2] # Path to test vectors
testIDs_file = args[3]
gwasIDs_file = args[4]
snp_file = args[5]
out_prefix = args[6] # Path to output directory
outfile = args[7]

####################
## Functions #######
####################

# Compute G %*% t(X) %*% T
compute_b <- function(path_to_genos, path_to_testvec, outpath) {

  # Compute t(X)T
  outfile_XT <- paste0(outpath, "xt_temp")
  cmd_XT <- paste("sh code/Calculate_FGr/compute_XT.sh", path_to_genos, path_to_testvec, outfile_XT, snp_file, testIDs_file, sep = " ")
  system(cmd_XT)

  # Adjust Betas to account for variance in x

  # Read in betas and genotype counts
  beta_plink <- fread(paste0(outpath, "xt_temp.Tvec.glm.linear"))
  count_plink <- fread(paste0(outpath, "xt_temp.gcount"))

  # Calculate length of mean centered genotypes from counts
  nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
  counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
  mean_gc <- counts / nOBS
  length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)

  # Fix betas
  betas_plink_norm <- beta_plink$BETA * length_mc_genos * (1/(n-1))

  # Compute GWAS genotype counts
  outfile_count <- paste0(outpath, "G_count")
  cmd_count <- paste("sh code/Calculate_FGr/compute_GWAS_count.sh", path_to_genos, outfile_count, snp_file, gwasIDs_file, sep = " ")
  system(cmd_count)

  # Calculate length of mean centered genotypes from counts
  count_plink <- fread(paste0(outpath, "G_count.gcount"))
  nOBS <- (count_plink$HOM_REF_CT + count_plink$HET_REF_ALT_CTS + count_plink$TWO_ALT_GENO_CTS)
  counts <- (count_plink$HOM_REF_CT * 0) + (count_plink$HET_REF_ALT_CTS * 1) + (count_plink$TWO_ALT_GENO_CTS * 2)
  mean_gc <- counts / nOBS
  length_mc_genos <- (count_plink$HOM_REF_CT * (-1 * mean_gc)^2) + (count_plink$HET_REF_ALT_CTS * (1 - mean_gc)^2) +  (count_plink$TWO_ALT_GENO_CTS * (2 - mean_gc)^2)
  length_mc_genos <- length_mc_genos * (1/(m-1))

  #  Re-write .linear file with correct betas
  beta_plink$BETA <- betas_plink_norm * (1/length_mc_genos)
  beta_reformat <- beta_plink %>% dplyr::select(ID, A1, BETA)
  fwrite(beta_reformat, paste0(outpath, "xt_temp.Tvec.glm.linear"), sep = "\t")

  # Compute b
  outfile_b <- paste0(outpath, "b")
  cmd_b <- paste("sh code/Calculate_FGr/GWAS_score.sh", path_to_genos, paste0(outpath, "xt_temp.Tvec.glm.linear"), outfile_b, snp_file, gwasIDs_file, sep = " ")
  system(cmd_b)

  # Read in and return b
  b = fread(paste0(outpath, "b.sscore"))
  b = as.matrix(b$BETA_SUM)

  # Remove temporary files
  fn <- paste0(outpath, "xt_temp*")
  cmd <- paste("rm", fn, sep = " ")
  system(cmd)
  fn <- paste0(outpath, "b*")
  cmd <- paste("rm", fn, sep = " ")
  system(cmd)
  fn <- paste0(outpath, "G_count*")
  cmd <- paste("rm", fn, sep = " ")
  system(cmd)

  return(b)
}

#####################
##     Main       ###
#####################

# Gather parameters
gwasID <- fread(gwasIDs_file)[, 1:2]
m <- nrow(gwasID)
testID <- fread(testIDs_file)
n <- nrow(testID)


# Compute b
b = compute_b(path_to_genos = plink_prefix, path_to_testvec = tvec_file, outpath = out_prefix)
b = as.data.frame(scale(b))
gwasID$POP <- b
colnames(gwasID) <- c("FID", "IID", "FGr")

fwrite(gwasID, outfile, row.names = F, col.names = T, quote = F, sep = "\t")




