## ddRAD data processing using [Stacks v2.55](http://catchenlab.life.illinois.edu/stacks/)

I used process_radtags first to process our raw ddRAD-seq data (PE sequencing).
Libraries are demultiplexed as below:
I have a barcode file for each library formatted with 2 columns (no headers) as: 

|          |           |
|----------|-----------|
|barcode_1 | sampleID_1|
|barcode_2 | sampleID_2|
|...       | ...       |

For our PE sequence data I used the following command to demultiplex and pre-process ddRADs
```
for i in *;
do
	process_radtags -P -p ./raw/"$i" -b ./barcodes/"$i"_barcodes.txt -o ./samples/"$i" -q -r --inline_index --barcode-dist-2 2 --renz_1 sbfI --renz_2 mspI
done
```
PCR duplicates were then removed using a custom Python script from [Peterson et al. 2012](https://pubmed.ncbi.nlm.nih.gov/22675423/)
- Requires Python 2 
```
for i in *.1.fq.gz ;
do
   name1=$(echo $i | cut -d '.' -f 1);
   ddRad_pcr_deduper.py -1 "$name1".1.fq.gz -2 "$name1".2.fq.gz -o dedupped_reads
done
```
## Raw read mapping 
Our final raw reads were mapped to the T. minimus V1 reference using [BWA-mem](https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=true)
```
for i in *1.fq.gz;
do
   name1=$(echo $i | cut -d '.' -f 1);
if [ ! -f "$name1"_PE.bam ]
 then
   echo "running"
   bwa mem -M -t 6 ./minimusNDH064_pseudohap.fasta "$name1".1.fq.gz "$name1".2.fq.gz | samtools view -Sb - > "$name1"_PE.bam 
 fi
 done
```
Add RG to bams: NOTE- fields will change depending on sequencing run.
```
for i in *_PE.bam;
do
   name1=${i%_*};
   name2=$(echo $name1 | cut -d '/' -f 2);
   
if [ ! -f "$name1"_PE.addRG.bam ]
 then
   echo "running"
   picard AddOrReplaceReadGroups -I "$name1"_PE.bam -O "$name1"_PE.addRG.bam -SO coordinate -RGID E00558 -LB tamias_ddRAD -PL illumina -PU misc -SM $name2 -VALIDATION_STRINGENCY LENIENT
 fi
 done
```
Next, we will sort and index our final bam files.
```
for i in *_PE.addRG.bam;
do
   name1=${name1%_*}; 
   samtools sort "$name1"_PE.addRG.bam > "$name1"_final.bam
done

for i in *_final.bam;
do
   name1=${i%_*}; 
   samtools index "$name1"_final.bam
done
```
## Genotype sites using gstacks
We provide a population map file (tamias_popmap.txt) defining our metapopulations for genotyping. It is formatted as:

|        |       |
|--------|-------|
|indv_01 | pop_01|
|indv_02 | pop_01|
|indv_03 | pop_02|
|indv_04 | pop_02|
|indv_05 | pop_03|
|...     | ...   |

```
gstacks -I ./final_aligned_bams -O ./stacks/gstacks_2 -M ./tamias_popmap.txt -t 24
```
## Here, we use populations to generate various call files (VCF,etc) for population genetic analyses
```
populations \
   -P gstacks_RefMap \
   -O populations_RefMap \
   -M tamias_popmap5sp.txt \
   -p 5 \
   -r 0.5 \
   -t 12 \
   --hwe \
   --smooth \
   --bootstrap \
   --fstats \
   --vcf \
   --phylip-var
```
## Before moving on I want to assess some QC metrics such as:
  - Mean depth per Individual
  - Proportion of missing data per individual
  - Variant mean Depth
  - Variant missingness

First, I generate the metric files using VCFTools:
```
# Calculate mean depth per individual
vcftools --gzvcf populations.snps.vcf --depth --out ./vcf_stats/tamias_178Ind_pop_raw

#Calculate proportion of missing data per individual
vcftools --gzvcf populations.snps.vcf --missing-indv --out ./vcf_stats/tamias_178Ind_pop_raw

# Calculate mean depth per site
vcftools --gzvcf populations.snps.vcf --site-mean-depth --out ./vcf_stats/tamias_178Ind_pop_raw

# Calculate proportion of missing data per site
vcftools --gzvcf populations.snps.vcf --missing-site --out ./vcf_stats/tamias_178Ind_pop_raw
```
## Next we can use R to plot the results to detect outliers and set DP/ missingness cutoffs using the emperical distributions.
## Phylogenetic Analysis
We first will apply some basic filtering to clean up our final alignment:
```
vcftools --vcf populations.snps.vcf --minDP 5 --minGQ 20 --max-missing 0.50 --max-alleles 2 --recode --recode-INFO-all --out iqTree.snps.filtered
```
After filtering, we kept 189 out of 189 Individuals and 310,361 out of a possible 318,618 Sites

We also want to exclude low coverage samples using a > 50% missingness cutoff.
```
vcftools --vcf iqTree.snps.filtered.recode.vcf --missing-indv --out tamias_178Ind_pop
# This creates a list of those indv. that have more than 50% missing data
mawk '$5 > 0.5' tamias_189Ind_pop.imiss | cut -f1 > lowDP.indv

vcftools --vcf iqTree.snps.filtered.recode.vcf --remove lowDP.indv --recode --recode-INFO-all --out iqTree.snps.filtered.172ind_Final
# Excludes: 17 individuals (172 out of 189)
```
Convert VCF file to a phylip file
```
vcf2phylip.py -i iqTree.snps.filtered.172ind_Final.recode.vcf -n
```
## IQ-tree ML tree estimation for the concatenated ddRAD data

Here are the final stats for this particular analysis:
172 individuals, alignment length: 358,836 bp 

```
iqtree2 -s iqTree.snps.filtered.172ind_Final.recode.min4.phy -m TEST -nt AUTO -pre tamias_172ind_RefMap_iqTree -bb 1000 -alrt 1000
```
-m Test: will test for a model of evolution
-nt AUTO: will determine best number of threads for run.
Best num threads: 11
Best model: TVM+F+R5

## Calculate site concordance factors using IQ-Tree
sCF is an additional metric we can use to assess branch support. It is the percentage of phylogentically informative sites that supports
a branch in the reference tree. 
```
iqtree2 -t tamias_172ind_RefMap_iqTree.contree -s iqTree.snps.filtered.172ind_Final.recode.min4.phy --scf 1000 --prefix tamias_172ind_RefMap_iqTree.sCF
```
_____________________________________________________________________________________________________________________________________________
### Analysis of Admixture
First, we will randomly subsample one SNP per locus
```
populations \
   -P gstacks_RefMap \
   -O populations_randomSNP \
   -M tamias_popmap5sp.txt \
   -p 5 \
   -r 0.85 \
   -t 12 \
   --vcf \
   --write-random-snp
```
## Admixture
Admixture requires a PLINK .bed file and the estimated number of ancestral populations.
We can convert our VCF file to the proper format using:
```
plink --vcf raw.lm.recode.vcf --aec --make-bed --out admix_filtered.SNPs
```
Now we can run ADMIXTURE
```
for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15 ;
do admixture --cv admix_filtered.SNPs.bed $K | tee log${K}.out; done
```
To find best K:
```
grep -h CV log*.out
```
## PCA
For the PCA, I applied (as above) some general filtering using vcftools to remove individuals and sites with a high proportion of missing data. This resulted in a final dataset of 
158 individuals and 218,700 sits

Export VCF to Plink format, filter for MAF
```
bcftools view -H populations.snps.lmDP5g85.recode.vcf | cut -f 1 | uniq | awk '{print $0"\tscf"$0}' > Tamias_MultiSmpl_Filtered.recode.vcf.chrom-map.txt
vcftools --gzvcf populations.snps.lmDP5g85.recode.vcf  --plink --chrom-map Tamias_MultiSmpl_Filtered.recode.vcf.chrom-map.txt --out Tamias_MultiSmpl_Filtered_plink

less Tamias_MultiSmpl_Filtered_plink.ped | cut -f 1,2 > Tamias_updateID.158.txt
plink --file Tamias_MultiSmpl_Filtered_plink --aec --allow-no-sex --update-ids Tamias_updateID.158.txt --make-bed --out Tamias_MultiSmpl_Filtered_plink.newID

### Filter using plink on MAF (only for within population level analyses)
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID --aec --maf 0.01 --out Tamias_MultiSmpl_Filtered_plink.newID.maf.01 --make-bed

### convert back into a new .ped file
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID.maf.01 --aec --allow-no-sex --recode --out Tamias_MultiSmpl_Filtered_plink.newID.maf.01
```
Linkage disequilibrium
PLINK has several options to estimate linkage disequilibrium (LD), the non-independent segregation of two loci, and filter variants according to certain distances and thresholds.
For biallelic markers, one of the most commonly used measures for LD is the squared coefficient of correlation (r2) 
(Hill & Robertson, 1968). It ranges between 0 and 1, where 0 indicates no correlation (= no linkage) and 1 perfect correlation (= perfect linkage disequilibrium).
```
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID --aec --indep-pairwise 1 kb 1 0.8 --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8 &> /dev/null
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID --aec --exclude Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.prune.out --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8 --make-bed
# convert back into a new .ped file 
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID.maf.01 --aec --indep-pairwise 1 kb 1 0.8 --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01 &> /dev/null
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID.maf.01 --aec --exclude Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01.prune.out --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01 --make-bed
```
Convert back to vcf file for final PCA analysis
```
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01 --aec --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01.172ind --recode vcf
```
For PCA anaysis, we are using the R package SNPRelate: See SNPRelate_PCA.r

## Species Tree Analysis
We can analyze our ddRAD data using SVDQuartets which estimates a species tree from unlinked SNP data. For SVDQ, we are going to use the same filtered VCF file from IQ-Tree but we want unlinked snps so we will add an extra thinning step
To get unlinked snps, we will use bcftools:
-w prune distance, -n means keep 1 snp per window
-N rand means select one snp at random in window
```
bcftools +prune -w 100bp -n 1 -N rand -o svdQ_SNP.filtered.172ind.100bp.thin.vcf.gz iqTree_SNP.filtered.172ind.recode.vcf.gz
```
count SNPS in final file
```
bcftools view -H svdQ_SNP.filtered.172ind.100bp.thin.vcf.gz | wc -l
```
This gives us 47,537 SNPs

Safety check, look at missingness per indv.
```
vcftools --gzvcf svdQ_SNP.filtered.172ind.100bp.thin.vcf.gz --missing-indv --out lowDP.indv
```
Convert to nexus for use in PAUP*
```
vcf2phylip.py -i svdQ_SNP.filtered.172ind.100bp.thin.vcf.gz -n
```

## SVDquartets
SVDquartets is described in [Chifman & Kubatko 2014](https://academic.oup.com/bioinformatics/article/30/23/3317/206559?login=true).
SVDQ is run through [PAUP*](https://paup.phylosolutions.com/)
First, we will construct our Paup blocks for two analyses. One will be a lineage tree and a second where we define species groupings and etimate a species tree.
initiat paup:
```
paup4a168_ubuntu64
```
# call file
Execute run_svdquartets_ddRads_lineagetree.sh

#Lineage Tree paup block: run_svdquartets_ddRads_lineagetree.sh
begin paup;
        log start file=svdQ.172ind.lineageTree.log replace;
        execute ./svdQ_SNP.filtered.172ind.100bp.thin.min4.nexus;
   outgroup minimusMSB268068 minimusNDH041 minimusNDH045 minimusNDH063 minimusNDH064 minimusNDH072 minimusNDH084 minimusZM10995 minimusZM11027 minimusZM11417 minimusZM12230 minimusZM12232 minimusZM12233 minimusZM12277 minimusZM13061 minimusZM16544 minimusZM16545;
        SVDquartets evalQuartets=all bootstrap=yes nreps=1000 nthreads=24 mrpFile=svdQ.172ind.lineageTree_qall_b1000; 
        savetrees file=svdQ.172ind.lineageTree.b1000.tre format=Newick root=yes supportValues=nodeLabels;
end;

# For the species tree, we first need to create taxa partitions (Tamias_species). I created a separate txt file called svdQ_totalSet_lineage_partitions.txt
# Here is what it looks like:
begin sets;
    taxpartition Tamias_species =
      amoenus: amoenusFMNH126103 amoenusJMG008 amoenusJMG058 amoenusJMS220 amoenusJRD008 amoenusJRD023 amoenusJRD028 amoenusJRD158 amoenusJRD161 amoenuslutiJMS297 amoenusMSB224396 amoenusMSB224913 amoenusMSB225555 amoenusMSB225676 amoenusMSB227699 amoenusMSB227795 amoenusMSB227882 amoenusMSB228278 amoenusMSB230567 amoenusMSB268064 amoenusMSB268065 amoenusMSB268069 amoenusMSB269634 amoenusMSB269635 amoenusMSB269636 amoenusMSB269638 amoenusMSB269640 amoenusMSB269856 amoenusMSB269857 amoenusMSB269858 amoenusMSB269867 amoenusMSB269868 amoenusMSB269869 amoenusMSB269870 amoenusMSB269871 amoenusMSB269882 amoenusMSB269957 amoenusMSB269966 amoenusMSB274470 amoenusMSB274506 amoenusMSB292133 amoenusMSB292142 amoenusMSB292143 amoenusMSB292145 amoenusNEOR184 amoenusRBCM19566 amoenusRBCM20037 amoenusZM12231 amoenusZM12446 amoenusZM13088 amoenusZM13089 amoenusZM16576 amoenusZM16577 amoenusZM16578 amoenusZM16579 amoenusZM16593 amoenusZM16594 amoenusZM16599,
      cratericus_A:  cratericusJRD111 cratericusJRD121 cratericusJRD122 cratericusJRD123 cratericusMSB197856 cratericusMSB264026 cratericusMSB269928 cratericusNDH002 cratericusNDH003 cratericusNDH007 cratericusNDH009 cratericusNDH011 cratericusNDH013 cratericusNDH014 cratericusNDH016 cratericusNDH018 cratericusNDH019 cratericusNDH020 cratericusNDH021 cratericusNDH022 cratericusNDH023 cratericusNDH024 cratericusNDH025 cratericusNDH026 cratericusNDH029 cratericusNDH030 cratericusNDH031 cratericusNDH032 cratericusNDH033 cratericusNDH034 cratericusNDH035 cratericusNDH036 cratericusNDH037 cratericusNDH038 cratericusZM13090 cratericusZM13093 cratericusZM13094 cratericusZM13097 cratericusZM13098 cratericusZM13100 cratericusZM13101 cratericusZM13102 cratericusZM13103 cratericusZM13105 cratericusZM13693 cratericusZM13694 cratericusZM13696 cratericusZM13740 cratericusZM13741 cratericusZM13742 cratericusZM13754 cratericusZM13755 cratericusZM13756 cratericusZM16546 cratericusZM16547 cratericusZM16548 cratericusZM16549 cratericusZM16555 cratericusZM16556 cratericusZM16580 cratericusZM16581 cratericusZM16582 cratericusZM16584 cratericusZM16585 cratericusZM16587 cratericusZM16588 cratericusZM16600 cratericusZM16601 cratericusZM16602, 
      cratericus_B:  cratericusJRD109 cratericusJRD110 cratericusJRD112 amoenusMSB274473 amoenusMSB274502 amoenusMSB274505,
      ruficaudus: ruficaudusJMG016 ruficaudusJMG134 ruficaudusJMS122 ruficaudusJMS201 ruficaudusJMS203 ruficaudusJRD142 ruficaudusZM16595 ruficaudusZM16596 ruficaudusZM16597 ruficaudusZM16598,
      minimus: minimusMSB268068 minimusNDH041 minimusNDH045 minimusNDH063 minimusNDH064 minimusNDH072 minimusNDH084 minimusZM10995 minimusZM11027 minimusZM11417 minimusZM12230 minimusZM12232 minimusZM12233 minimusZM12277 minimusZM13061 minimusZM16544 minimusZM16545,
      ;
END;

# Finally, cat the nexus dataset file with partionsFile.txt to create the final file which will be run.

#Species Tree paup block: run_svdquartets_ddRads_speciestree.sh
Execute run_svdquartets_ddRads_speciestree.sh
### Paup block run:
begin paup;
   log start file=svdQ.172ind.speciesTree.log replace;
   execute ./svdQ_SNP.filtered.172ind.100bp.thin.SpeciesSet.nexus;
   outgroup minimusMSB268068 minimusNDH041 minimusNDH045 minimusNDH063 minimusNDH064 minimusNDH072 minimusNDH084 minimusZM10995 minimusZM11027 minimusZM11417 minimusZM12230 minimusZM12232 minimusZM12233 minimusZM12277 minimusZM13061 minimusZM16544 minimusZM16545;
   SVDquartets evalQuartets=all bootstrap=yes nreps=1000 taxpartition=Tamias_species nthreads=8 mrpFile=svdQ.172ind.SpeciesTree_qall_b1000;
   savetrees file=svdQ.172ind..SpeciesTree.b1000.tre format=Newick root=yes supportValues=nodeLabels;
end;
________________________________________________________________________________________________________________________________________
#### HyDe ####
# We used HyDe to detect hybridization under the invariants framework of SVDquartets.
# HyDe allows us to estimate γ, which is the parental contribution in a putatively hybrid genome, 
# where a value of 0.5 would indicate a 50:50 genomic contribution from each parent. 

vcf2phylip.py -i svdQ_SNP.filtered.172ind.100bp.thin.vcf -n
run_hyde.py -i HyDe_SNP.filtered.172ind.100bp.thin.phy -m hyde_popMap_5sp.txt -o minimus -n 172 -t 5 -s 47537
________________________________________________________________________________________________________________________________________
