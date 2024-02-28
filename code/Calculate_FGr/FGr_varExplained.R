## This script calculates the relative weights of GWAS PCs on T^GWAS

## Requires: gwas_pca.eigenvec, Tm.txt

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_Tm.R <eigenvecs> <Tm> <outfile name>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
}))


vecs_file = args[1] # eigenvectors
Tm_file = args[2] # Tm
out_file = args[3]

# Load gwas eigen vecs
vecs <- fread(vecs_file)
vecs <- vecs[,3:ncol(vecs)]
vecs <- apply(vecs, 2, as.numeric)
vecs <- scale(vecs)


# Load T^GWAS
Tm = fread(Tm_file)
print(head(Tm))
Tm = Tm$Fgr

# Normalize T to have var 1
Tm <- scale(Tm)

# Get number of PCs
nPC <- ncol(vecs)

# Construct data fram to collate results
dfOut <- matrix(NA, nrow = ncol(vecs), ncol = 2)
dfOut[,1] <- seq(1,ncol(vecs))

# Compute multiple R^2 and rho(PC, FGr)
for (i in 1:nrow(dfOut)) {


  # Compute variance explained
  mod <- lm(Tm ~ vecs[,1:i])
  r2 <- cor(Tm, fitted(mod))^2
  print(r2)


  # Collect output
  dfOut[i,2] <- r2

}


# Save weights
fwrite(dfOut, out_file,row.names=F,quote=F,sep="\t", col.names = T)
