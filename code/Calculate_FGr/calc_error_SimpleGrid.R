# This script computes the error in either FGr or PCs for Simple Grid

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_FGr.R <tilde Fgr file> <covars file> <outfile> ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

tFGr_file = args[1]
covar_file = args[2]
out_file = args[3]


# Read in Tilde FGr
dfTilde <- fread(tFGr_file)

# Read in covar file
dfCovar <- fread(covar_file)

# Find (1 - b2)
find_b2_FGr <- function(x,y){

  y <- as.matrix(y)
  x <- as.matrix(x)

  # Compute variance explained
  r <- cov(x,y) / sqrt(var(x), var(y))
  r2 <- r^2

  error = 1-r2
  return(error)
}

# Find (1 - b2)
find_b2_PCs <- function(x,y){

  y <- as.matrix(y)
  x <- as.matrix(x)
  # Compute variance explained
  mod <- lm(x ~ y)
  r2 <- cor(x, fitted(mod))^2
  b2 = 1-r2
  return(b2)
}



# Find angle
#errorFGr <- find_b2(dfTilde$tildeFGr, dfCovar$FGr)
#error10 <- find_b2(dfTilde$tildeFGr, dfCovar[,6:15])
#error35 <- find_b2(dfTilde$tildeFGr, dfCovar[,6:40])
#error100 <- find_b2(dfTilde$tildeFGr, dfCovar[,6:105])

errorFGr <- find_b2(dfTilde$FGr, dfCovar$FGr)
error10 <- find_b2(dfTilde$FGr, dfCovar[,6:15])
error35 <- find_b2(dfTilde$FGr, dfCovar[,6:40])
error100 <- find_b2(dfTilde$FGr, dfCovar[,6:105])



# Save output
out <- matrix(NA, nrow =1, ncol=4)
out[1,] <- c(errorFGr, error10, error35, error100)
out <- as.data.frame(out)
colnames(out) <- c("FGr", "10", "35", "100")
print(out)
fwrite(out, out_file, row.names = F, col.names = T, quote = F, sep = "\t")




