# Split individuals into GWAS and test panel

args=commandArgs(TRUE)

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


ss_test = as.numeric(args[1]) / 36
ss_gwas = as.numeric(args[2]) / 36
input_file = args[3]
gwas_file = args[4]
test_file = args[5]

# Read in pop info
pop <- fread(input_file, header = F)
colnames(pop) <- c("FID", "IID", "POP", "LAT", "LONG" )

# Sample ss_test individuals per deme
df_test <- pop %>% group_by(POP) %>% sample_n(ss_test)
x <- paste0("A", as.character(df_test$POP))
df_test$pop <- x
#test <- df_test[,1:2]

# Use the rest of the individuals for the GWAS panel
df_gwas <- anti_join(pop, df_test) %>% group_by(POP) %>% sample_n(ss_gwas)
x <- paste0("A", as.character(df_gwas$POP))
df_gwas$pop <- x
#gwas <- df_gwas[,1:2]

# Write panel IDs to file
write.table(df_gwas, gwas_file, quote = F, col.names = T, row.names = F)
write.table(df_test, test_file, quote = F, col.names = T, row.names = F)

