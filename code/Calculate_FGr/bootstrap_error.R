# This script computes the error in either FGr using bootstrap resampling

args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript bootstrap_error.R <FGr file> <outfile> ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

fgr_file = args[1]
out_file = args[2]


# Read in FGr dataframe
df <- fread(fgr_file)[,1:4]

# Compute number of individuals per deme
n <- nrow(df %>% filter(POP == 0))
print(n)

# Compute bootstrap samples
sigma_j <- rep(0, 36)
for (j in 0:35) {

  # Select only individuals from a single deme
  dfTmp <- df %>% filter(POP == j)
  FGr <- dfTmp$FGr
  FGrBar <- mean(FGr)

  # For each deme compute variance
  sigma_j[j+1] <- sum(FGr - FGrBar) * (1/(nrow(dfTmp) - 1))

}


print(sigma_j)
# Compute error
error <- ((1/36) * sum(sigma_j)) / var(df$FGr)

# Save output
out <- as.data.frame(error)
colnames(out) <- "Error"
fwrite(out, out_file, row.names = F, col.names = T, quote = F, sep = "\t")




