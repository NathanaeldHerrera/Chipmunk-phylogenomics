# Load necessary libraries
library("gdsfmt")      # For working with GDS files
library("SNPRelate")   # For SNP analysis

# Convert a VCF file to GDS format
snpgdsVCF2GDS("Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01.172ind.vcf", 
              "Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01.172ind.gds")

# Load in the GDS file
genofile <- openfn.gds("Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01.172ind.gds")

# Perform PCA on the SNP data
pca <- snpgdsPCA(autosome.only=FALSE, genofile, snp.id=NULL, num.thread=16)

# Calculate the percentage of variance explained by each principal component
pc.percent <- pca$varprop * 100
head(round(pc.percent, 2))  # Display the first few percentages

# Load population codes from a text file
pop_code <- scan("Tamias_158PCA_SpeciesList.txt", what=character())

# Retrieve sample IDs from the GDS file
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))

# Extract individual names from the population codes
indiv <- sapply(strsplit(pop_code, split=","), "[[", 2)
head(cbind(samp.id, indiv))  # Display the first few sample IDs and individual names

# Create a data frame to hold PCA results and population information
tab <- data.frame(sample.id = pca$sample.id,
                  pop = factor(indiv)[match(pca$sample.id, samp.id)],
                  EV1 = pca$eigenvect[, 1],  # First eigenvector
                  EV2 = pca$eigenvect[, 2],  # Second eigenvector
                  stringsAsFactors = FALSE)
head(tab)  # Display the first few rows of the data frame

# Plot the PCA results
plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), pch = 19, 
     xlab="PC 1", ylab="PC 2")
legend("topright", legend=levels(tab$pop), pch=19, col=1:6)

# Save the plot to a PDF file
pdf("dddRAD_172indv.RefMap_FILTERED.SNPRelatePCA.pdf", width = 6, height = 6)
plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), pch = 20, 
     xlab="PC 1", ylab="PC 2")  # Make the plot
legend("bottomleft", legend=levels(tab$pop), pch=19, col=1:4)
dev.off()  # Close the PDF device
