# Clear environment
rm(list = ls())

# Load necessary libraries
library(SNPRelate)
library(tidyverse)
library(gdsfmt)
library(RColorBrewer)

# Set run date
run_date <- format(Sys.Date(), "%d%b%Y")

# File paths
vcf_file <- "Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.141ind.vcf"
gds_file <- "Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.141ind.gds"
species_file <- "Tamias_141PCA_SpeciesList_noMinimus.txt"

# Convert VCF to GDS (if not already done)
snpgdsVCF2GDS(vcf.fn = vcf_file, out.fn = gds_file, method = "biallelic.only")

# Open GDS file
genofile <- openfn.gds(gds_file)

# Run PCA
pca <- snpgdsPCA(genofile, autosome.only = FALSE, num.thread = parallel::detectCores())
pc.percent <- round(pca$varprop * 100, 2)
print(head(pc.percent))

# Read sample info
pop_raw <- scan(species_file, what = character())
indiv <- sapply(strsplit(pop_raw, ","), `[[`, 2)

# Get sample IDs
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Create PCA table
tab <- data.frame(
  sample.id = pca$sample.id,
  pop = factor(indiv)[match(pca$sample.id, samp.id)],
  EV1 = pca$eigenvect[, 1],
  EV2 = pca$eigenvect[, 2],
  stringsAsFactors = FALSE)

# Axis labels with explained variance
xlab_text <- paste0("PC 1 (", pc.percent[1], "%)")
ylab_text <- paste0("PC 2 (", pc.percent[2], "%)")

# Use RColorBrewer palette
palette_name <- "Set1"
n_pop <- length(levels(tab$pop))
plot_colors <- brewer.pal(min(n_pop, brewer.pal.info[palette_name, "maxcolors"]), palette_name)

# Save to PDF
pdf(file = paste0("dddRAD_131Indv.RefMap_FILTERED.SNPRelatePCA_IC_", run_date, ".pdf"),
    width = 6, height = 6)
plot(tab$EV1, tab$EV2, col = plot_colors[as.integer(tab$pop)], pch = 20, cex = 2.5,
     xlab = xlab_text, ylab = ylab_text, main = "PCA Plot")
# text(tab$EV1, tab$EV2, labels = tab$sample.id, pos = 3, cex = 0.6)
legend("topright", legend = levels(tab$pop), pch = 19, col = plot_colors)
dev.off()

# Close GDS
snpgdsClose(genofile)
