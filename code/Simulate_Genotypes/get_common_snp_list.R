# Get list of SNPs that are at 5% MAF or greater in both the Test and GWAS panel
library(data.table)

# Test panel.freq file
tp <- fread(snakemake@input[[1]], header = T)

# GWAS panel .freq
gp <- fread(snakemake@input[[2]], header = T)

# Subset variants
t <- subset(tp, tp$ALT_FREQS > 0.05 & tp$ALT_FREQS < 0.95)
g <- subset(gp, gp$ALT_FREQS > 0.05 & gp$ALT_FREQS < 0.95)

# Merge to include only variants passing frequency filter in both panels
dat <- merge(t,g, by=c("#CHROM","ID"))

write.table(dat$ID, snakemake@output[[1]], quote = F, col.names = F, row.names = F)
