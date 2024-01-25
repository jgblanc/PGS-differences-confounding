# This script computes the error in either FGr or PC1 for the 4PopSplit example

args=commandArgs(TRUE)

if(length(args)<3){stop("Rscript calc_FGr.R <op file> <vector file> <outfile> ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

pop_file = args[1]
vec_file = args[2]
out_file = args[3]


# Read in popfile
pops <- fread(pop_file)

# Read in vector file
vecs <- fread(vec_file)
vec <- vecs[,3]
name <- colnames(vec)
vec <- as.matrix(vec)



# merge to get GWAS popIDs
ID <- inner_join(pops, vecs)
pops <- unique(ID$POP)
ID <- ID %>% mutate(ID = case_when(POP == pops[1] ~ 0, POP == pops[2] ~ 1))
idvec <- as.matrix(ID$ID - mean(ID$ID))

# Function to find angle between vectors
find_angle <- function(x,y){
  dot.prod <- t(x)%*%y
  norm.x <- norm(x,type="2")
  norm.y <- norm(y,type="2")
  theta <- acos(dot.prod / (norm.x * norm.y))
  return(as.numeric(theta))
}

# Find angle
theta <- find_angle(vec, idvec)

# Find error
error <- theta / pi

# Save output
out <- as.data.frame(error)
colnames(out) <- name
fwrite(out, out_file, row.names = F, col.names = T, quote = F, sep = "\t")




