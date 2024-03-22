# This script computes an estimate of \tilde{F}_{Gr}
args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_LOCO_FGr.R <FGr file> <outfile> ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

rep = args[1]
outfile = args[2]

# Get rep number
repNum = as.numeric(strsplit(rep, "")[[1]][2])

# Read in target FGr
covars <- fread(args[2 + repNum])[,1:4]
M <- nrow(covars)

# Compute average FGr per deme over other reps
dfAvg <- matrix(NA,nrow = 36, ncol = 99)
x <- seq(1, 100)
x <- x[-repNum]
for(r in 1:99) {

  print(paste0("r ", i))

  i <- x[r]

  # Read in covar file and take mean per deme
  tmp <- fread(args[2 + i])[,1:4]
  tmpAvg <- tmp %>% group_by(POP) %>% summarise(avgFGr = mean(FGr))

  # add to output list
  dfAvg[,r] <- tmpAvg

}

# Compute average across rows
tildeFGr <- rowMeans(tmpAvg)
print(tildeFGr)

# Add column to covars file
covars$tildeFGr <- 0
for (i in 1:nrow(covars)) {

  popID <- covars[i,3]
  covars[,5] <- tildeFGr[popID + 1]

}


# Save output
fwrite(covars, outfile, row.names = F, col.names = T, quote = F, sep = "\t")




