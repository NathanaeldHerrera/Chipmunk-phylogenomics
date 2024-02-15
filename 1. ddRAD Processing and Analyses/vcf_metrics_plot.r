# load necessary packages for plotting
library(tidyverse)
library(ggplot2)

# Run date for output plots
run_date <- format(Sys.Date(), "%d%b%Y")

# Mean depth per Individual
ind_depth <- read_delim("./tamias_178Ind_pop.idepth", delim = "\t",
                        col_names = c("ind", "nsites", "depth"), skip = 1, 
                        show_col_types = FALSE)

# Plot for right here
a <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

# Output to pdf
pdf(file=paste("./vcf_metric_plots/Rplot_meanDP_",run_date,".pdf",sep=""),width=8,height=8)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() 
dev.off()

# Proportion of missing data per individual
ind_miss  <- read_delim("./tamias_178Ind_pop.imiss", delim = "\t",
                        col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), 
                        skip = 1, show_col_types = FALSE)

# Plot for right here
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

# Output to pdf
pdf(file=paste("./vcf_metric_plots/Rplot_indMiss_",run_date,".pdf",sep=""),width=8,height=8)
a <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() 
dev.off()

# Variant mean depth
var_depth <- read_delim("./tamias_178Ind_pop.ldepth.mean", delim = "\t",
           col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1, 
           show_col_types = FALSE)

# Summarize mean Depth 
summary(var_depth$mean_depth)
a + theme_light() + xlim(0, 100)

# plot right here
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()           

# Output to pdf
pdf(file=paste("./vcf_metric_plots/Rplot_Var_mean_DP_",run_date,".pdf",sep=""),width=8,height=8)
a <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() 
dev.off()

# Per site variant missingness
var_miss <- read_delim("./tamias_178Ind_pop.lmiss", delim = "\t",
                       col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), 
                       skip = 1, show_col_types = FALSE)

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

# Per site (Phred based) quality

var_qual <- read_delim("./tamias_178Ind_pop.lqual", delim = "\t", 
            col_names = c("chr", "pos", "qual"), skip = 1, 
            show_col_types = FALSE)

# Plot for right here 
a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light()

# Output to pdf
pdf(file=paste("./vcf_metric_plots/Rplot_VariantQual_",run_date,".pdf",sep=""),width=8,height=8)
a <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)
a + theme_light() 
dev.off()

