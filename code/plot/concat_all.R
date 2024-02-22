args=commandArgs(TRUE)

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
}))

nfile = length(args)
print(nfile)

# Add all reps together
df <- fread(args[1])

if (nfile > 2) {

  for (i in 2:(nfile - 1)) {
     
    tmp <- fread(args[i])
    df <- rbind(df, tmp)
  }


}


# Save output
fwrite(df, args[nfile] ,row.names=T,quote=F,sep="\t", col.names = T)


