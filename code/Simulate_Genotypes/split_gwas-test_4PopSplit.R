## This script uses the 4 Pop split metadata file to split IDs into GWAS and Test Panels

args=commandArgs(TRUE)

# Split Population IDs into two sets of test/gwas pairs

library(data.table)
library(dplyr)

# Assign all arguments
popfile = args[1]
outC1GWAS = args[2]
outC1Test = args[3]
outC2GWAS = args[4]
outC2Test = args[5]
gSize = as.numeric(args[6])
tSize = as.numeric(args[7])

# Read in pop info
pop <- fread(popfile, header = T)


# Case 1: A,C gwas B,D test
gwas <- subset(pop, pop$POP == "A" | pop$POP == "C")
test <- subset(pop, pop$POP == "B" | pop$POP == "D")

# Down sample test and GWAS set
test <- test %>% group_by(POP) %>% sample_n(tSize/2)
gwas <- gwas %>% group_by(POP) %>% sample_n(gSize/2)

write.table(gwas, outC1GWAS, quote = F, col.names = T, row.names = F)
write.table(test, outC1Test, quote = F, col.names = T, row.names = F)

# Case 2: A,B gwas, C,D test
gwas <- subset(pop, pop$POP == "A" | pop$POP == "B")
test <- subset(pop, pop$POP == "C" | pop$POP == "D")

# Down sample test and GWAS set
test <- test %>% group_by(POP) %>% sample_n(tSize/2)
gwas <- gwas %>% group_by(POP) %>% sample_n(gSize/2)

write.table(gwas, outC2GWAS, quote = F, col.names = T, row.names = F)
write.table(test, outC2Test, quote = F, col.names = T, row.names = F)
