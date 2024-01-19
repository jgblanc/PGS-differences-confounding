## Select causal sites and draw effect sizes

args=commandArgs(TRUE)

if(length(args)<10){stop("Rscript draw_effects_sizes_num_causal_4PopSplit.R <frequency file> <output_file> <heritability> <alpha>
                        <test panel IDs> <popfile> <probability of flipping effect size> <test panel genotype prefix> <common snps file>")}

suppressWarnings(suppressMessages({
  library(data.table)
  library(dplyr)
  library(Matrix)
  library(pgenlibr)
  library(tidyr)
}))




# Read in all arguments
freq_file=args[1]
out_file = args[2]
h2 = as.numeric(args[3])
alpha = as.numeric(args[4])
testIDs_file = args[5]
popfile = args[6]
prob = as.numeric(args[7]) # probability effect size is positive given pC - pD is positive
nc = as.numeric(args[8])
geno_prefix = args[9]
common_snps_file = args[10]


# load variant frequency file
p = fread(freq_file)
colnames(p)=c("chr","ID","REF","ALT","ALT_FREQS","COUNT")
p=p[,c("chr","ID","ALT_FREQS")]
p[, c("CHROM", "position","ref","alt") := tstrsplit(ID, "_", fixed=TRUE)]
p = p[,c("CHROM","ID","position","ALT_FREQS")]
p$position = as.numeric(p$position)

# Read in common snps and select allele frequencies from them
cSnps <- fread(common_snps_file, header = FALSE)
p <- p %>% filter(p$ID %in% cSnps$V1)


# Function to read in genotype matrix for a set of variants
read_genos <- function(geno_prefix, betas) {

  pvar <- pgenlibr::NewPvar(paste0(geno_prefix, ".pvar"))
  d1 <- pgenlibr::NewPgen(paste0(geno_prefix, ".pgen"))
  var.ids <- betas$ID
  var.indx <- rep(0, length(var.ids))
  for (i in 1:length(var.indx)) {
    var.indx[i] <- pgenlibr::GetVariantsById(pvar,var.ids[i])
  }
  X <- ReadList(d1,var.indx, meanimpute=F)
  colnames(X) <- var.ids

  return(X)
}

# Read in population ID info
pop <- fread(popfile, header = F)
colnames(pop) <- c("FID", "IID", "POP")

# Read in test fam file
ids <- fread(testIDs_file)
fam <- fread(paste0(geno_prefix, ".psam"))

# Get only test individuals
test_inds <- inner_join(ids, pop)


# Get number of individuals in each population
pops <- unique(test_inds$POP)
n1IDs <- test_inds %>% filter(POP == pops[1]) %>% select(FID)
n1IDs <- n1IDs$FID
n2IDs <- test_inds %>% filter(POP == pops[2]) %>% select(FID)
n2IDs <- n2IDs$FID


# Select casual variants and effect sizes
if (nc == 0) { # Check if there are no causal variants and if so add a dummy one with effect size 0
  causal.variants <- sample_n(p, 1)
  causal.variants$beta <- 0
} else { # If not sample effect sizes

  # Sample nc causal variants
  causal.variants <- sample_n(p, nc)
  print(paste0("The number of variants is ", nrow(causal.variants)))

  # Now generate the effect sizes from these variants calculate the independent component of variance required
  sigma2_l = h2 / sum( sapply( causal.variants$ALT_FREQS,function(x){
    beta= ( 2*x*(1-x)) ^ (1-alpha)
    return(beta)
  }))
  print(paste0("Sigma2_1 is ", sigma2_l))

  #sample maf-dependent effects using the model above
  causal.variants$beta = sapply( causal.variants$ALT_FREQS , function(x){
    beta = rnorm( 1 , mean = 0, sd = sqrt(sigma2_l * (2*x*(1-x))^-alpha ))
  })
  causal.variants$beta <- abs(causal.variants$beta)

  #let's calculate sigma2_g to confirm that the total genetic variance is indeed h2
  sigma2_g = sum( mapply(function(b,p){ b^2* 2*p*(1-p) }, causal.variants$beta, causal.variants$ALT_FREQS))
  print(paste0("Sigma2_g is ", sigma2_g))

  # Read in genotype matrix for causal variants
  G <- read_genos(geno_prefix, causal.variants[,"ID"])
  rownames(G) <- fam$IID
  G1 <- G[rownames(G) %in% n1IDs, ]
  G2 <-  G[rownames(G) %in% n2IDs, ]

  # Calculate population specific allele frequeny
  p1 <- colMeans(G1)/2
  p2 <- colMeans(G2)/2

  # Get allele frequency difference
  diff <- p1 - p2

  # Create correlation between effect size and pop ID
  for (i in 1:nrow(causal.variants)){
    b <- causal.variants[i,"beta"]
    if (diff[i] >= 0) {
      causal.variants[i,"beta"] <- sample(c(-1, 1),1, prob = c((1-prob), prob)) * b
    } else {
      causal.variants[i,"beta"] <- sample(c(1, -1),1, prob = c((1-prob), prob)) * b
    }
  }

  # Print probability of positive beta
  indx_greater = which(diff > 0)
  indx_smaller = which(diff < 0)
  print(sum(causal.variants[indx_greater,]$beta > 0)/sum(diff > 0))
  print(sum(causal.variants[indx_smaller,]$beta < 0)/sum(diff < 0))



}


#save the effect sizes to file and use plink2 to generate PRS
fwrite(causal.variants%>%
         mutate(ALT = "T")%>%
         select(ID,ALT,beta),
           out_file,
       row.names=F,col.names=F,quote=F,sep="\t")


