# Format Tm and other covariates into a txt file that will work with Plink 2

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


args=commandArgs(TRUE)

if( length(args) != 4){stop("Usage: <pop file> <FGr> <Eigenvec file> <output file>") }

# Parse args
pop_file = args[1]
FGr_file = args[2]
evec_file = args[3]
output_file = args[4]

# Read in files
pops <- fread(pop_file, header = T)

# Read in FGr
FGr <- fread(FGr_file, header = T)


# Format datafile
df <- dplyr::inner_join(pops, FGr) %>% select("FID", "IID", "POP", "FGr")
pops <- unique(df$POP)
df <- df %>% mutate(PopID = case_when((POP == pops[1]) ~ 1, (POP == pops[2]) ~ 0))
df$PopID <- df$PopID - mean(df$PopID)

# Read in eigenvec files
vecs <- fread(evec_file)
colnames(vecs)[1] <- "FID"
df <- inner_join(df, vecs)

# Write Tm to file
fwrite(df, output_file,row.names=F,quote=F,sep="\t", col.names = T)
