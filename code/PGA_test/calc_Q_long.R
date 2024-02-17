# This script takes a setof effect sizes and test vector and computes Q

args=commandArgs(TRUE)

if(length(args)!=7){stop("Rscript calc_Q.R <prefix to genos> <causal betas> <non causal betas> <test vec> <true effects>
                         <number of times to resample in empirical null> <outfile>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

genos_prefix = args[1]
cbetas_file = args[2]
ncbetas_file = args[3]
pop_file = args[4]
true_file = args[5]
num_resample = as.numeric(args[6])
outfile = args[7]


# Function to read in genotype matrix for a set of variants
read_genos <- function(geno_prefix, betas) {

  pvar <- pgenlibr::NewPvar(paste0(geno_prefix, ".pvar"))
  d1 <- pgenlibr::NewPgen(paste0(geno_prefix, ".pgen"))
  var.ids <- betas$ID
  var.indx <- rep(0, length(var.ids))
  for (i in 1:length(var.indx)) {
    var.indx[i] <- pgenlibr::GetVariantsById(pvar,var.ids[i])
  }
  X <- ReadList(d1,var.indx, meanimpute=F)
  colnames(X) <- var.ids

  return(X)
}

# Function to compute PGS
pgs <- function(X, betas) {

  Bhat_strat <- as.matrix(betas)
  Z_strat <- X %*% Bhat_strat
  out <- cbind(Z_strat)
  colnames(out) <- c("STRAT")
  return(out)
}

# Function to calculate Q
calc_q <- function(sscore, tvec) {

  tvec <- as.matrix(tvec)

  qhat <-  (1/N) * t(tvec) %*% sscore

  return(qhat)
}

# Function to flip effect sizes
flip <- function(betas) {
  new_betas <- sample(c(-1,1), length(betas),  replace = T) * betas
  return(new_betas)
}

# Function to flip effect sizes and recompute Qx
en <- function(X, betas, tvec) {

  # Flip effect sizes
  betas <- flip(betas)

  # Compute PGS
  prs <- pgs(X, betas)

  # Calc Qx - Test
  q <- t(calc_q(prs, tvec))[1,1]

  return(q)
}



# Read in betas
cbetas <- fread(cbetas_file)
ncbetas <- fread(ncbetas_file)

# Read in all IDs
fam <- fread(paste0(genos_prefix, ".psam"))

# Make Long test vector
dfPop <- fread(pop_file)
colnames(dfPop) <- c("#FID","IID","Pop", "Lat", "Long")
dfTvec <- inner_join(dfPop, fam)
Tvec <- dfTvec$Long - mean(dfTvec$Long)

# count number of TP individuals
N <- length(Tvec)








# Read in genotype matrix and select the correct rows
cX <- read_genos(genos_prefix, cbetas)
rownames(cX) <- fam$IID
cX <- cX[rownames(cX) %in% dfTvec$IID, ]

ncX <- read_genos(genos_prefix, ncbetas)
rownames(ncX) <- fam$IID
ncX <- ncX[rownames(ncX) %in% dfTvec$IID, ]


# Compute q
cQ <- calc_q(pgs(cX,cbetas$BETA_Strat), Tvec)
ncQ <- calc_q(pgs(ncX,ncbetas$BETA_Strat), Tvec)

# Compute true values of q
true_betas <- fread(true_file)
colnames(true_betas) <- c("ID", "A1", "BETA_Strat")
trueQ <- calc_q(pgs(cX,true_betas$BETA_Strat), Tvec)

## Compute empirical p-values

# causal
redraws <- matrix(0, ncol = 1, nrow = num_resample)
for (i in 1:num_resample){
  redraws[i,] <- en(cX, cbetas$BETA_Strat, Tvec)
}
all_strat <- redraws[,1]
cP <- length(all_strat[abs(all_strat) > abs(cQ[1,1])]) /length(all_strat)
cStan <- cQ / sd(all_strat)

# noncausal
redraws <- matrix(0, ncol = 1, nrow = num_resample)
for (i in 1:num_resample){
  redraws[i,] <- en(ncX, ncbetas$BETA_Strat, Tvec)
}
all_strat <- redraws[,1]
ncP <- length(all_strat[abs(all_strat) > abs(ncQ[1,1])]) /length(all_strat)
ncStan <- ncQ / sd(all_strat)

# Get names
tmp <- strsplit(outfile, "/")[[1]][12]
tmp2 <- strsplit(tmp, "_")
type <- strsplit(tmp2[[1]][1], "-")[[1]][2]
tmp3 <- strsplit(tmp2[[1]][2], ".txt")[[1]][1]
snps <- as.numeric(strsplit(tmp3, "-")[[1]][2])


# Assemble output
out <- as.data.frame(matrix(NA, nrow= 1, ncol = 11))
colnames(out) <- c("c-q", "nc-q", "true-q", "c-p", "nc-p", "c-bias", "nc-bias", "type","L", "stan-c", "stan-nc")
out[1,] <- c(cQ, ncQ, trueQ, cP, ncP, (cQ - trueQ), (ncQ - trueQ), type, snps, cStan, ncStan)

# Save output
print(out)
fwrite(out, outfile ,row.names=T,quote=F,sep="\t", col.names = T)


