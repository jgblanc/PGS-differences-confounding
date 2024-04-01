# This script averages polygenic scores over replications

args=commandArgs(TRUE)

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


outfile = args[1]
popfile = args[2]
nfile = length(args) - 2
print(nfile)

# Read in popfile
pops <- fread(popfile)
colnames(pops) <- c("IID","FID","Pop", "Lat", "Long")

# Read in all PRS and combine with popfile
df <- fread(args[3])
df <- inner_join(df, pops)

# Add the rest of the files
for (i in 4:length(args)) {

  tmp <- fread(args[i])
  tmp <- inner_join(tmp, pops)

  df <- rbind(df, tmp)
}
df <- df %>% select(FGr, Pop, Lat, Long)

# Compute the average
dfOut <- df %>% group_by(Pop) %>% summarise(avg_FGr = mean(FGr), LAT = mean(Lat), LONG = mean(Long))

# Save output
fwrite(dfOut, outfile ,row.names=F,quote=F,sep="\t", col.names = T)

