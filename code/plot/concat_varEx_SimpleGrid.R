
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
out$rep <- as.character(1)

# Loops through all other files and add onto output
for (i in 2:(nfile-1)) {

  tmp <- fread(args[i])
  tmp$rep <- as.character(i)
  out <- rbind(out, tmp)

}

# Compute mean and confidence intervals
out <- out %>% group_by(rep) %>% summarise(num = n(), avg_r2 = mean(V2),
                                           lowerci = avg_r2 - (1.96 * (sd(V2)  /  sqrt(n()))),
                                           upperci = avg_r2 + (1.96 * (sd(V2)  /  sqrt(n()))))

# Save output
fwrite(out, args[nfile] ,row.names=T,quote=F,sep="\t", col.names = T)


