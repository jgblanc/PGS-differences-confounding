## This script uses the 4 Pop split metadata file to split IDs into GWAS and Test Panels

args=commandArgs(TRUE)

# Split Population IDs into two sets of test/gwas pairs

library(data.table)
library(dplyr)

set.seed(as.numeric(args[7]))


# Read in pop info
pop <- fread(args[1], header = T)

# Get test set size
size = as.numeric(args[2])
print(size)

# Case 1: A,C gwas B,D test
gwas <- subset(pop, pop$POP == "A" | pop$POP == "C")
test <- subset(pop, pop$POP == "B" | pop$POP == "D")

# Down sample test set
test <- test %>% group_by(POP) %>% sample_n(size/2)

write.table(gwas, args[3], quote = F, col.names = T, row.names = F)
write.table(test, args[4], quote = F, col.names = T, row.names = F)

# Case 2: A,B gwas, C,D test
gwas <- subset(pop, pop$POP == "A" | pop$POP == "B")
test <- subset(pop, pop$POP == "C" | pop$POP == "D")

# Down sample test set
test <- test %>% group_by(POP) %>% sample_n(size/2)

write.table(gwas, args[5], quote = F, col.names = T, row.names = F)
write.table(test, args[6], quote = F, col.names = T, row.names = F)
