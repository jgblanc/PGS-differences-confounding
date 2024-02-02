# Make Test Vector

library(dplyr)
library(data.table)

args=commandArgs(TRUE)

if( length(args) != 3){stop("Usage: <pop file> <fam file> <type of test> <output file> ") }

# Parse args
id_file = args[1]
test_type = args[2]
output_file = args[3]

# Read in Fam file
pop <- fread(id_file, header = T)[,1:5]

if (test_type == "LAT") {

  print(test_type)

  # Test vector is latitude
  Tvec <- pop$LAT

  # Mean center
  Tvec <- (Tvec-mean(Tvec))

} else if (test_type == "PS") {

  print(test_type)

  # Test vector is deme 25 vs everyone else (mean center)
  pop <- pop %>% mutate(T1 = case_when(POP == 25 ~ 1, POP != 25 ~ 0)) %>% mutate(T1 = T1 - mean(T1))
  Tvec <- pop$T1

} else {
  stop("Please enter acceptable test type: LAT, PS")
}


# Write to file
pop$Tvec <- Tvec
pop <- pop %>% select("FID", "IID", "POP", "Tvec")
write.table(pop, output_file,row.names=F,quote=F,sep="\t", col.names = T)
