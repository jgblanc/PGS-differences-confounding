# This script computes an estimate of \tilde{F}_{Gr}
args=commandArgs(TRUE)

if(length(args)<2){stop("Rscript compute_LOCO_FGr.R <FGr file> <outfile> ")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))

pc_num <- 100
outfile = args[1]
print(length(args))

### Read in all covar files ###
rep_files <- list()
for (rep in 2:length(args)) {
  covar <- fread(args[rep])
  rep_files[[rep-1]] <- covar
}

### Compute Error in FGr ###

## Compute variance within demes
avg_sigma_per_rep <- rep(0, 100)
for (r in 1:100) {

  print(r)
  covar <- rep_files[[r]]

  # Calc variance per deme
  sigma_j <- rep(0, 36)
  for (j in 0:35) {

    # Select only individuals from a single deme
    dfTmp <- covar %>% filter(POP == j)
    FGr <- dfTmp$FGr
    n <- nrow(dfTmp)

    # Compute mean
    FGrBar <- mean(FGr)

    # For each deme compute sample variance
    sigma_j[j+1] <- sum((FGr - FGrBar)^2) * (1/(n - 1))

  }
  # Take the average variance within demes
  avg_sigma_per_rep[r] <- mean(sigma_j)
}

## Compute variance of means

# Find mean per deme across reps
deme_means <- matrix(NA, nrow = 36, ncol = 100)
for (r in 1:100) {

  # Find mean FGr per deme per replicate
  tmp <- rep_files[[r]]
  tmp_sum <- tmp %>% group_by(POP) %>% summarise(meanDeme = mean(FGr))
  deme_means[, r] <- tmp_sum$meanDeme

}
meanDeme <- rowMeans(deme_means)

# Compute variance of deme means across replcates
deme_vars <- matrix(NA, nrow = 36, ncol = 100)
for (r in 1:100) {

  tmp <- rep_files[[r]]
  tmp_sum <- tmp %>% group_by(POP) %>% summarise(meanDeme = mean(FGr))
  deme_vars[,r] <-  (tmp_sum$meanDeme - meanDeme)^2

}
variance_deme_means <- rowSums(deme_vars) * (1/35) # check 35 or 36

## Compute emprical variance of target vector
empricalVar <- rep(0,100)
for (r in 1:100) {
  tmp <- rep_files[[r]]
  empricalVar[r] <- var(tmp$FGr)
}


## Compute Error per replicate
Error <- rep(0,100)
for (r in 1:100) {
  Error[r] <- (mean(avg_sigma_per_rep) + mean(variance_deme_means)) / empricalVar[r]
}
FGrError <- mean(Error)

### Compute Lower Bound Error in PCs ###

## Compute variance within demes
avg_sigma_per_rep <- matrix(NA, nrow = pc_num, ncol = 100)
for (pc in 1:pc_num) {

  print(paste0("PC ", pc))

  for (r in 1:100) {

    covar <- rep_files[[r]]

    # Calc variance per deme
    sigma_j <- rep(0, 36)
    for (j in 0:35) {

      # Select only individuals from a single deme
      dfTmp <- covar %>% filter(POP == j)
      colName <- paste0("PC",pc)
      hatVec <- as.matrix(dfTmp[,..colName])
      n <- nrow(dfTmp)

      # Compute mean
      hatVecBar <- mean(hatVec)

      # For each deme compute sample variance
      sigma_j[j+1] <- sum((hatVec - hatVecBar)^2) * (1/(n - 1))

    }
    # Take the average variance within demes
    avg_sigma_per_rep[pc, r] <- mean(sigma_j)
  }
}

## Compute emprical variance of target vector
empricalVar <- matrix(NA, nrow=pc_num, ncol=100)
for (pc in 1:pc_num) {
  for (r in 1:100) {
    tmp <- rep_files[[r]]
    colName <- paste0("PC",pc)
    hatVec <- as.matrix(tmp[,..colName])
    empricalVar[pc, r] <- var(hatVec)[1,1]
  }
}

## Compute Error
Error_PC <- matrix(NA, nrow=pc_num, ncol=100)
for (pc in 1:pc_num) {
  for (r in 1:100) {
    Error_PC[pc, r] <- (rowMeans(avg_sigma_per_rep)[pc]) / empricalVar[r]
  }
}
rowMeans(Error_PC)


# Save output
dfOut <- as.data.frame(c(FGrError, rowMeans(Error_PC)))
colnames(dfOut) <- "Error"
dfOut$Type <- c("FGr", paste0("PC", seq(1,pc_num)))
fwrite(dfOut,outfile, col.names = T, row.names = F, quote = F)




