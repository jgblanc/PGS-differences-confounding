# Format Tm and other covariates into a txt file that will work with Plink 2

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
}))


args=commandArgs(TRUE)

if( length(args) != 5){stop("Usage: <pop file> <FGr> <Eigenvec file> <output file>") }

# Parse args
pop_file = args[1]
FGr_file = args[2]
evec_file = args[3]
output_file = args[4]
test_type = args[5]

# Read in files
pops <- fread(pop_file, header = T)
colnames(pops)<-c("IID","FID","POP", "Lat", "Long")

# Read in FGr
FGr <- fread(FGr_file, header = T)

# Format datafile
df <- dplyr::inner_join(pops, FGr) %>% select("FID", "IID", "POP", "FGr", "Lat")

# Add column with test pattern
if (test_type == "LAT") {

  print(test_type)

  # Test vector is latitude
  df$PopID <- df$Lat
  df$PopID <- df$PopID - mean(df$PopID)



} else if (test_type == "PS") {

  print(test_type)

  # Test vector is deme 25 vs everyone else (mean center)
  df <- df %>% mutate(PopID = case_when(POP == 25 ~ 1, POP != 25 ~ 0)) %>% mutate(PopID = PopID - mean(PopID))

} else {
  stop("Please enter acceptable test type: LAT, PS")
}
df <- df %>% select("FID", "IID", "POP", "FGr", "PopID")

# Read in eigenvec files
vecs <- fread(evec_file)
colnames(vecs)[1] <- "FID"
df <- inner_join(df, vecs)

# Write Tm to file
fwrite(df, output_file,row.names=F,quote=F,sep="\t", col.names = T)
