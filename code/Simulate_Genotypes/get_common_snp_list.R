# Get list of SNPs that are at 1% MAF or greater in both the Test and GWAS panel
args=commandArgs(TRUE)

library(data.table)
library(dplyr)

# Assign all arguments
tpFile = args[1]
gpFile = args[2]
nSNP = as.numeric(args[3])
outfile = args[4]

# Test panel.freq file
tp <- fread(tpFile, header = T)

# GWAS panel .freq
gp <- fread(gpFile, header = T)

# Subset variants
t <- subset(tp, tp$ALT_FREQS > 0.01 & tp$ALT_FREQS < 0.99)
g <- subset(gp, gp$ALT_FREQS > 0.01 & gp$ALT_FREQS < 0.99)
print(paste0("The Test SNP number is ", nrow(t)))
print(paste0("The GWAS SNP number is ", nrow(g)))


# Merge to include only variants passing frequency filter in both panels
#dat <- inner_join(t,g, by=c("#CHROM","ID"))

# Randomly sample only nSNP number of SNPs
#dat <- dat %>% sample_n(nSNP)

#write.table(dat$ID, outfile, quote = F, col.names = F, row.names = F)
