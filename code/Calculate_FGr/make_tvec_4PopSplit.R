# Make Test Vector

library(dplyr)
library(data.table)

args=commandArgs(TRUE)

if( length(args) != 2){stop("Usage: <pop file> <fam file> <output file> ") }

# Parse args
id_file = args[1]
output_file = args[2]


# Read in ID file
ids <- fread(id_file)

# Make mean centered Test vec (PopID)
test_pops <- unique(ids$POP)
pop <- ids %>% mutate(Tvec = case_when(POP == test_pops[1] ~ 1, POP == test_pops[2] ~ 0)) %>% mutate(Tvec = Tvec - mean(Tvec))

# Write to file
write.table(pop, output_file,row.names=F,quote=F,sep="\t", col.names = T)
