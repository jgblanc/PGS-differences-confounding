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
find_b2_one <- function(x,y){
  dot.prod <- t(x)%*%y
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  b2 = 1-(dot.prod / (norm.x * norm.y))^2
  return(b2)
}

# Find (1 - b2) - mulitple vectors
find_b2_many <- function(x,y){
  tmp <- rep(0, ncol(y))
  y <- as.matrix(y)
  for (i in 1:ncol(y)) {
    k <- y[,i]
    dot.prod <- t(x)%*%k
    norm.x <- norm(x,type="2")
    norm.k <- norm(k,type="2")
    tmp[i] = (dot.prod / (norm.x * norm.k))^2
  }
  b2 = 1-sum(tmp)
  return(as.numeric(b2))
}

# Find angle
errorFGr <- find_b2_one(dfTilde$tildeFGr, dfCovar$FGr)
error10 <- find_b2_many(dfTilde$tildeFGr, dfCovar[,6:15])
error35 <- find_b2_many(dfTilde$tildeFGr, dfCovar[,6:40])
error100 <- find_b2_many(dfTilde$tildeFGr, dfCovar[,6:150])




# Save output
out <- matrix(NA, nrow =1, ncol=4)
out[1,] <- c(errorFGr, error10, error35, error100)
out <- as.data.frame(out)
colnames(out) <- c("FGr", "10", "35", "100")
print(out)
fwrite(out, out_file, row.names = F, col.names = T, quote = F, sep = "\t")




