args=commandArgs(TRUE)


library(data.table)
library(dplyr)
library(Matrix)

print(args[1])
print(args[2])
print(args[3])
print(args[4])

vecs_file = args[1] # eigenvectors
popfile = args[2]
gvalue_file = args[3]
out_file = args[4]

pheno_pattern = strsplit(out_file, "/")[[1]][10]
print(pheno_pattern)

env_s = as.numeric(strsplit(strsplit(out_file, "/")[[1]][9], "_")[[1]][2])
print(env_s)

# Load gwas eigen vecs
vecs <- fread(vecs_file)
vecs <- vecs[,3:ncol(vecs)]
vecs <- apply(vecs, 2, as.numeric)
vecs <- scale(vecs)

## Get confounder

# Read in genetic values
prs=fread(gvalue_file)
colnames(prs)<-c("IID","dosage","prs")

# Get simulation sample size
sample_size=nrow(prs)


# Load file containing the pop ID for each inidividual
pop=fread(popfile,header=F)
print(head(pop))
colnames(pop)<-c("IID","FID","Pop", "Lat", "Long")

# Add this info to genetic value file
prs=merge(prs, pop, by="IID", sort=F)
prs$env = 0

if (pheno_pattern == "LAT") {

  print(pheno_pattern)

  # Shift the environmental variable by env_s/lat
  num_lat_demes <- length(unique(prs$Lat))
  if (env_s == 0) {
    shifts <- rep(0, num_lat_demes)
  } else {
    shifts <- seq(0, env_s, env_s / (num_lat_demes-1))
  }

  # Add stratification along Latitude
  for (i in 1:num_lat_demes) {
    prs <- prs %>% group_by(Lat) %>% mutate(env = ifelse(Lat == (i-1), env + shifts[i], env))
  }

} else if (pheno_pattern == "DIAG") {

  print(pheno_pattern)

  # Shift the environmental variable by env_s/# of diagonal squares
  num_diag_demes <- length(unique(prs$Lat))
  if (env_s == 0) {
    shifts <- rep(0, num_diag_demes)
  } else {
    shifts <- seq(0, env_s, env_s / (num_diag_demes-1))
  }

  # Add stratification along diagonal Latitude
  id_diag <- c(0,7,14, 21, 28, 35)
  for (i in 1:length(id_diag)) {
    prs <- prs %>% group_by(Pop) %>% mutate(env = ifelse(Pop == id_diag[i], env + shifts[i], env))
  }
} else if (pheno_pattern == "PS") {

  print(pheno_pattern)

  # Add stratification on deme 25
  id_diag <- c(25)
  for (i in 1:length(id_diag)) {
    prs <- prs %>% group_by(Pop) %>% mutate(env = ifelse(Pop == id_diag[i], env + env_s, env))
  }

} else {
  stop("Please enter acceptable phenotype pattern: LAT, DIAG, PS")
}

Tm = prs$env
Tm = scale(Tm)

# Get number of PCs
nPC <- ncol(vecs)

# Construct data fram to collate results
dfOut <- matrix(NA, nrow = ncol(vecs), ncol = 2)
dfOut[,1] <- seq(1,ncol(vecs))

# Compute multiple R^2 and rho(PC, FGr)
for (i in 1:nrow(dfOut)) {


  # Compute variance explained
  mod <- lm(Tm ~ vecs[,1:i])
  r2 <- cor(Tm, fitted(mod))^2
  print(r2)


  # Collect output
  dfOut[i,2] <- r2

}


# Save weights
fwrite(dfOut, out_file,row.names=F,quote=F,sep="\t", col.names = T)

