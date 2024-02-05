

# Define some env variables for ease of use
IN_VCF=./populations.snps.vcf
OUT=./vcf_stats/tamias_178Ind_pop_raw

# calculate allele frequency for each variant. 
# The --freq2 just outputs the frequencies without information about the alleles, 
# --freq would return their identity. We need to add max-alleles 2 to exclude sites 
# that have more than two alleles.
vcftools --gzvcf $IN_VCF --freq2 --out $OUT --max-alleles 2

# Calculate mean depth per individual
vcftools --gzvcf $IN_VCF --depth --out $OUT

# Calculate mean depth per site
vcftools --gzvcf $IN_VCF --site-mean-depth --out $OUT

#Calculate proportion of missing data per individual
vcftools --gzvcf $IN_VCF --missing-indv --out $OUT

# Calculate proportion of missing data per site
vcftools --gzvcf $IN_VCF --missing-site --out $OUT

# Calculate heterozygosity and inbreeding coefficient per individual
# Computing heterozygosity and the inbreeding coefficient (F) for each individual can quickly 
# highlight outlier individuals that are e.g. inbred (strongly negative F), suffer from high 
# sequencing error problems or contamination with DNA from another individual leading to inflated 
# heterozygosity (high F), or PCR duplicates or low read depth leading to allelic dropout and 
# thus underestimated heterozygosity (stongly negative F). However, note that here we assume 
# Hardy-Weinberg equilibrium. If the individuals are not sampled from the same population, 
# the expected heterozygosity will be overestimated due to the Wahlund-effect. It may still be 
# worth to compute heterozygosities even if the samples are from more than one population to 
# check if any of the individuals stands out which could indicate problems.

vcftools --gzvcf $IN_VCF --het --out $OUT

#PLOT stats in R

# load necessary packages for plotting
library(tidyverse)
library(ggplot2)

# Run date for output plots
run_date <- format(Sys.Date(), "%d%b%Y")

# Variant mean Depth
var_depth <- read_delim("./tamias_178Ind_pop.ldepth.mean", delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

# Summarize mean Depth 
summary(var_depth$mean_depth)
a + theme_light() + xlim(0, 100)

# plot right here
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()           

# Output to pdf
pdf(file=paste("./vcf_metric_plots/Rplot_meanDP_",run_date,".pdf",sep=""),width=8,height=8)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() 
dev.off()

# Variant missingness
var_miss <- read_delim("./tamias_178Ind_pop.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

# Summarize variant missingness
summary(var_miss$fmiss)

# Plot for right here
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

# Output to pdf
pdf(file=paste("./vcf_metric_plots/Rplot_VariantMiss_",run_date,".pdf",sep=""),width=8,height=8)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() 
dev.off()

# Minor allele frequency
var_freq <- read_delim("./tamias_178Ind_pop.frq", delim = "\t",
                       col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)
# find minor allele frequency
var_freq$maf <- var_freq %>% select(a1, a2) %>% apply(1, function(z) min(z))

# Summarize MAF
summary(var_freq$maf)

# Plot for right here
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

# Output to pdf
pdf(file=paste("./vcf_metric_plots/Rplot_MAF_",run_date,".pdf",sep=""),width=8,height=8)
a <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() 
dev.off()

### Individual based statistics
#Mean depth per Individual
ind_depth <- read_delim("./tamias_178Ind_pop.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1)

a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

# Proportion of missing data per individual
ind_miss  <- read_delim("./tamias_178Ind_pop.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()


