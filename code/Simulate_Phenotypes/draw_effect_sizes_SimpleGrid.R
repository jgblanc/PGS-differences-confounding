## Select causal sites and draw effect sizes

args=commandArgs(TRUE)

if(length(args)<7){stop("Rscript draw_effects_sizes_SimpleGrid.R <frequency file> <output_file> <heritability> <alpha>
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
nc = as.numeric(args[6])
common_snps_file = args[7]


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


}


#save the effect sizes to file and use plink2 to generate PRS
fwrite(causal.variants%>%
         mutate(ALT = "T")%>%
         select(ID,ALT,beta),
           out_file,
       row.names=F,col.names=F,quote=F,sep="\t")


