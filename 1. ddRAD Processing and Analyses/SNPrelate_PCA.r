library("gdsfmt")
library("SNPRelate")

snpgdsVCF2GDS("Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01.172ind.vcf", "Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01.172ind.gds")
#snpgdsSummary("JG10_MultiSmpl_MultiSampleFiltered2.recode.gds")
genofile <- openfn.gds("Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01.172ind.gds")
pca <- snpgdsPCA(autosome.only=FALSE, genofile, snp.id=NULL, num.thread=16)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

pop_code <- scan("Tamias_158PCA_SpeciesList.txt", what=character())
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
#sapply(strsplit(pop_code, split=","), "[[", 2)
indiv <- sapply(strsplit(pop_code, split=","), "[[", 2)
head(cbind(samp.id, indiv))

tab <- data.frame(sample.id = pca$sample.id,
   pop = factor(indiv)[match(pca$sample.id, samp.id)],
   EV1 = pca$eigenvect[,1], # the first eigenvector
   EV2 = pca$eigenvect[,2], # the second eigenvector
   stringsAsFactors = FALSE)
head(tab)pop_code <- scan("Tamias_158PCA_SpeciesList.txt", what=character())
samp.id <- read.gdsn(index.gdsn(genofile, "sample.id"))
#sapply(strsplit(pop_code, split=","), "[[", 2)
indiv <- sapply(strsplit(pop_code, split=","), "[[", 2)
head(cbind(samp.id, indiv))

tab <- data.frame(sample.id = pca$sample.id,
    pop = factor(indiv)[match(pca$sample.id, samp.id)],
    EV1 = pca$eigenvect[,1], # the first eigenvector
    EV2 = pca$eigenvect[,2], # the second eigenvector
    stringsAsFactors = FALSE)
head(tab)

plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), pch = 19, xlab="PC 1", ylab="PC 2")
legend("topright", legend=levels(tab$pop), pch=19, col=1:6)


pdf("dddRAD_172indv.RefMap_FILTERED.SNPRelatePCA.pdf", width = 6, height = 6)
plot(tab$EV1, tab$EV2, col=as.integer(tab$pop), pch = 20, xlab="PC 1", ylab="PC 2") # Make plot
legend("bottomleft", legend=levels(tab$pop), pch=19, col=1:4)
dev.off()
