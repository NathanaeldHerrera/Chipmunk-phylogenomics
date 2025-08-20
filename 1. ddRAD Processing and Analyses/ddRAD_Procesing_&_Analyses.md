# ddRAD Data Processing & Population Genomics Pipeline ([Stacks v2.55](http://catchenlab.life.illinois.edu/stacks/))
This repository outlines an end‑to‑end pipeline for demultiplexing, mapping, genotyping, and downstream analyses of ddRAD‑seq data using Stacks v2.55 and common population‑genomics tools.

We process our raw ddRAD-seq data (PE sequencing) using process_radtags.
Libraries are demultiplexed as below:
I have a barcode file for each library formatted with 2 columns (no headers) as: 

|          |           |
|----------|-----------|
|barcode_1 | sampleID_1|
|barcode_2 | sampleID_2|
|...       | ...       |

For our PE sequence data I used the following command to demultiplex and pre-process ddRADs
```
   # Loop through all items in the current directory.
for i in *; do
    # Run the process_radtags command for each item.
    process_radtags \
        -P \  # Use paired-end processing.
        -p ./raw/"$i" \  # Input path for raw data files, using the current item as the directory name.
        -b ./barcodes/"$i"_barcodes.txt \  # Barcode file corresponding to the current item.
        -o ./samples/"$i" \  # Output directory for processed samples, named after the current item.
        -q \  # Enable quiet mode (suppress output).
        -r \  # Remove uncalled bases.
        --inline_index \  # Use inline index for barcoding.
        --barcode-dist-2 2 \  # Allow two mismatches in barcode comparison.
        --renz_1 sbfI \  # First restriction enzyme (sbfI).
        --renz_2 mspI  # Second restriction enzyme (mspI).
done
```
PCR duplicates were then removed using a custom Python script ([ddRad_pcr_deduper.py](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/1.%20ddRAD%20Processing%20and%20Analyses/ddRad_pcr_deduper.py)) from [Peterson et al. 2012](https://pubmed.ncbi.nlm.nih.gov/22675423/)
- Requires Python 2
- NOTE: The python script ddRad_pcr_deduper.py was written and provided to me for this project by [Dr. Jesse Webber](https://weberlab.integrativebiology.wisc.edu/). I made no modifications to the original script.
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
# Loop through all files that match the pattern *1.fq.gz in the current directory.
for i in *1.fq.gz; do
    # Extract the base name by removing the file extension.
    name1=$(echo $i | cut -d '.' -f 1)
    
    # Check if the corresponding BAM file does not already exist.
    if [ ! -f "$name1"_PE.bam ]; then
        echo "running"  # Print a message indicating the process is starting.

        # Run BWA to align reads, specifying the reference genome and input FASTQ files.
        bwa mem -M -t 6 ./minimusNDH064_pseudohap.fasta "$name1".1.fq.gz "$name1".2.fq.gz | \
        samtools view -Sb - > "$name1"_PE.bam  # Convert the output to BAM format.
    fi
done
```
Use Picard to add or replace read groups in the BAM file.
- fields will change depending on sequencing run.
```
for i in *_PE.bam; do
    name1=${i%_*}
    name2=$(echo $name1 | cut -d '/' -f 2)
    if [ ! -f "$name1"_PE.addRG.bam ]; then
        echo "running"
        picard AddOrReplaceReadGroups \
            -I "$name1"_PE.bam \  # Input BAM file.
            -O "$name1"_PE.addRG.bam \  # Output BAM file with added read groups.
            -SO coordinate \  # Sort order for the output file.
            -RGID E00558 \  # Read group ID.
            -LB tamias_ddRAD \  # Library name.
            -PL illumina \  # Platform (e.g., Illumina).
            -PU misc \  # Platform unit.
            -SM $name2 \  # Sample name (derived from the file path).
            -VALIDATION_STRINGENCY LENIENT  # Validation setting for processing.
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
## Run gstacks to assemble stacks of RAD tags from aligned BAM files.
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
gstacks \
    -I ./final_aligned_bams \  # Input directory containing the aligned BAM files.
    -O ./stacks/gstacks_2 \  # Output directory where gstacks will save the results.
    -M ./tamias_popmap.txt \  # Population map file that specifies sample and population assignments.
    -t 24  # Number of threads to use for processing.
```
## Population genetic analyses
Run populations to summarize genetic data from stacks created by gstacks.
```
populations \
   -P gstacks_RefMap \  # Input directory containing the output from gstacks.
   -O populations_RefMap \  # Output directory where results will be saved.
   -M tamias_popmap5sp.txt \  # Population map file that specifies sample and population assignments.
   -p 5 \  # Minimum number of individuals required to call a locus.
   -r 0.5 \  # Minimum frequency threshold for a locus to be included (50% in this case).
   -t 12 \  # Number of threads to use for processing.
   --hwe \  # Include Hardy-Weinberg Equilibrium tests in the output.
   --smooth \  # Apply smoothing to allele frequency estimates.
   --bootstrap \  # Perform bootstrap analysis to assess the robustness of estimates.
   --fstats \  # Calculate F-statistics.
   --vcf \  # Generate a VCF file with the results.
   --phylip-var  # Output a PHYLIP file with variance information.
```
## Before moving on I want to assess some QC metrics
such as:
  - Mean depth per Individual
  - Proportion of missing data per individual
  - Variant mean Depth
  - Variant missingness

First, I generate the metric files using VCFTools:
```
# Calculate depth statistics for the VCF file .
vcftools --gzvcf populations.snps.vcf.gz \
    --depth \
    --out ./vcf_stats/tamias_178Ind_pop_raw

# Calculate the proportion of missing data per individual.
vcftools --gzvcf populations.snps.vcf.gz \
    --missing-indv \
    --out ./vcf_stats/tamias_178Ind_pop_raw

# Calculate the mean depth per site.
vcftools --gzvcf populations.snps.vcf.gz \
    --site-mean-depth \
    --out ./vcf_stats/tamias_178Ind_pop_raw

# Calculate the proportion of missing data per site.
vcftools --gzvcf populations.snps.vcf.gz \
    --missing-site \
    --out ./vcf_stats/tamias_178Ind_pop_raw

# Calculate per-site (Phred-based) quality.
vcftools --gzvcf populations.snps.vcf.gz \
    --site-quality \
    --out ./vcf_stats/tamias_178Ind_pop_raw
```
Next we can use R to plot the results to detect outliers and set DP/ missingness cutoffs using the emperical distributions.

See [vcf_QC_metrics_plot.r](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/1.%20ddRAD%20Processing%20and%20Analyses/vcf_metrics_plot.r)

## Phylogenetic Analysis
We first will apply some basic filtering to clean up our final alignment:
```
# Filter a VCF file based on specific criteria and generate a new filtered VCF file.
vcftools --vcf populations.snps.vcf \  # Input VCF file containing SNPs.
    --minDP 5 \  # Minimum depth of coverage required for a variant to be included.
    --minGQ 20 \  # Minimum genotype quality required for a variant to be included.
    --max-missing 0.50 \  # Maximum proportion of missing data allowed for a site.
    --max-alleles 2 \  # Maximum number of alleles allowed at any site (biallelic variants only).
    --recode \  # Recode the VCF file to include only the variants that pass the filters.
    --recode-INFO-all \  # Retain all INFO fields in the output VCF.
    --out iqTree.snps.filtered  # Output file prefix for the filtered VCF.
```
Exclude low coverage samples using a > 50% missingness cutoff.
```
# Calculate the proportion of missing data per individual and output the results.
vcftools --vcf iqTree.snps.filtered.recode.vcf \
    --missing-indv \
    --out tamias_178Ind_pop

# Filter the results to create a list of individuals with more than 50% missing data.
mawk '$5 > 0.5' tamias_178Ind_pop.imiss | cut -f1 > lowDP.indv

# Remove individuals with high missing data from the VCF file and create a new filtered VCF.
vcftools --vcf iqTree.snps.filtered.recode.vcf --missing-indv --out tamias_178Ind_pop
# Excludes: 17 individuals (158 out of 178)
```
Convert VCF file to a phylip file
```
vcf2phylip.py -i iqTree.snps.filtered.158ind_Final.recode.vcf -n
```
## IQ-tree ML tree estimation for the concatenated ddRAD data

Here are the final stats for this particular analysis:
172 individuals, alignment length: 358,836 bp 

```
# Run IQ-tree
# 158 individuals, alignment length: 574,668 bp 
iqtree2 -s populations.snps.lowDP.trim_75.phylip -m TEST -nt AUTO -pre tamias_158ind.snps.lowDP.trim_75_iqTree -bb 1000 -alrt 1000 
```
First run we can assess best model and number of threads:
Best num threads: 7
Best model: TVM+F+G4

## Calculate site concordance factors using IQ-Tree
sCF is an additional metric we can use to assess branch support. It is the percentage of phylogentically informative sites that supports
a branch in the reference tree. 
```
iqtree2 -t tamias_158ind.snps.lowDP.trim_75_iqTree.contree -s populations.snps.lowDP.trim_75.phylip --scf 1000 --prefix tamias_158ind.snps.lowDP.trim_75_iqTree

```
## Analysis of Admixture
First, we will randomly subsample one SNP per locus
```
   populations \
   -P /home/nh253642e/ch2_wgs/ddRad_new_analyses/samples/stacks/gstacks_Final_RefMap \
   -O new_run_042125/c_Admixture \
   -M new_run_042125/tamias_popmap_5sp_042124.txt \
   -p 5 \
   -r 0.5 \
   --min-maf 0.01 \
   -t 12 \
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
# Loop through a range of values for K (the number of populations).
for K in 2 3 4 5 6 7 8 9 10 11 12 13 14 15; do
    # Run admixture analysis for the current value of K, with cross-validation.
    admixture --cv admix_filtered.SNPs.bed $K | tee log${K}.out
done
```
To find best K:
```
grep -h CV log*.out
```
## PCA
For the PCA, I applied (as above) some general filtering using vcftools to remove individuals and sites with a high proportion of missing data.

Export VCF to Plink format, filter for MAF
```
# Remove those indiviudals that sequenced poorly using the .imiss output.
# My data has most individuals showing less then .25 missing data (yay). 
# Remove individuals with more than 50% missing data:
vcftools --vcf populations.snps.vcf --missing-indv --out out

# THis creates a list of those indv. that have more than 50% missing data
mawk '$5 > 0.25' out.imiss | cut -f1 > lowDP.indv

# Remove those individuals:
vcftools --vcf populations.snps.vcf --remove lowDP.indv --recode --recode-INFO-all --out populations.snps_lowDP_nomini
# Excludes: 17 individuals

### Look at PCA without minimus
vcftools --vcf populations.snps.vcf --remove minimus.indv --recode --recode-INFO-all --out raw.lm_noMINI

### Apply some general filtering:
vcftools --vcf populations.snps_lowDP.vcf --minDP 5 --max-missing 0.85 --minGQ 20 --max-alleles 2 --recode --recode-INFO-all --out populations.snps.lmDP3g85
### Export VCF to Plink format
#Total Set FOR PCA
bgzip -c populations.snps_lowDP_131ind.vcf > populations.snps_lowDP_131ind.vcf.gz
tabix -p vcf populations.snps_lowDP_131ind.vcf.gz

bcftools view -H populations.snps_lowDP_131ind.vcf | cut -f 1 | uniq | awk '{print $0"\tscf"$0}' > Tamias_MultiSmpl_Filtered.recode.vcf.chrom-map.txt

vcftools --vcf populations.snps_lowDP_131ind.vcf  --plink --chrom-map Tamias_MultiSmpl_Filtered.recode.vcf.chrom-map.txt --out Tamias_MultiSmpl_Filtered_plink

less Tamias_MultiSmpl_Filtered_plink.ped | cut -f 1,2 > Tamias_updateID.131.txt
plink --file Tamias_MultiSmpl_Filtered_plink --aec --allow-no-sex --update-ids Tamias_updateID.131.txt --make-bed --out Tamias_MultiSmpl_Filtered_plink.newID
```

# Linkage disequilibrium
PLINK has several options to estimate linkage disequilibrium (LD), the non-independent segregation of two loci, 
and filter variants according to certain distances and thresholds.
For biallelic markers, one of the most commonly used measures for LD is the squared coefficient of correlation (r2) 
(Hill & Robertson, 1968). It ranges between 0 and 1, where 0 indicates no correlation (= no linkage) and 1 perfect correlation (= perfect linkage disequilibrium).
```
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID --aec --indep-pairwise 1 kb 1 0.8 --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8 &> /dev/null
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID --aec --exclude Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.prune.out --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8 --make-bed

### convert back into a new .ped file
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID.maf.01 --aec --indep-pairwise 1 kb 1 0.8 --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01 &> /dev/null
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID.maf.01 --aec --exclude Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01.prune.out --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.maf01 --make-bed
```
#Export to VCF for PCA analysis
```
plink --bfile Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8 --aec --out Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.131ind --recode vcf

```
For PCA anaysis, we are using the R package SNPRelate: See [SNPRelate_PCA.r](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/1.%20ddRAD%20Processing%20and%20Analyses/SNPrelate_PCA.r)

## Species Tree Analysis
We can analyze our ddRAD data using SVDQuartets which estimates a species tree from unlinked SNP data. ### For SVDQ, we are going to use the same filtered VCF file from above
Safety check, look at missingness per indv.
```
vcftools --gzvcf populations.snps.vcf --missing-indv --out lowDP.indv
```
Convert to nexus
```
/home/nh253642e/software/vcf2phylip/vcf2phylip.py -i populations.snps.vcf -n
grep -rl _final ./populations.snps.vcf | xargs sed -i 's/_final//g'
cut -d ' ' -f 1 populations.snps.min4.phy
```
## SVDquartets
SVDquartets is described in [Chifman & Kubatko 2014](https://academic.oup.com/bioinformatics/article/30/23/3317/206559?login=true).
SVDQ is run through [PAUP*](https://paup.phylosolutions.com/)
First, we will construct our Paup blocks for two analyses. One will be a lineage tree and a second where we define species groupings and etimate a species tree.

Lineage Tree paup block: run_svdquartets_ddRads_lineagetree.sh
```
begin paup;
        log start file=svdQ.158ind.lineageTree.log replace;
        execute ../populations.snps.min4.nexus;
   outgroup minimusMSB268068 minimusNDH041 minimusNDH045 minimusNDH063 minimusNDH064 minimusNDH072 minimusNDH084 minimusZM10995 minimusZM11027 minimusZM114>
        SVDquartets evalQuartets=all bootstrap=yes nreps=1000 nthreads=24 mrpFile=svdQ.158ind.lineageTree_qall_b1000; 
        savetrees file=svdQ.158ind.lineageTree.b1000.tre format=Newick root=yes supportValues=nodeLabels;
end;
```
begin paup:
```
paup4a168_ubuntu64
```
Within Paup* we can execute our analyses
```
Execute run_svdquartets_ddRads_lineagetree.sh
```
For the species tree, we first need to create taxa partition. This taxa partiotion is then added to the final nexus file at the end for the species tree analysis
```
begin sets;
    taxpartition Tamias_species =
        amoenus:        amoenusJMG058 amoenusJMS220 amoenusJRD008 amoenusJRD023 amoenusJRD158 amoenusJRD161 amoenuslutiJMS297 amoenusMSB224396 amoenusMSB224913 amoenusMSB227699 amoenusMSB227795 amoenusMSB227882 amoenusMSB228278 amoenusMSB230567 amoenusMSB268064 amoenusMSB268065 amoenusMSB268069 amoenusMSB269867 amoenusMSB269868 amoenusMSB269869 amoenusMSB269870 amoenusMSB269871 amoenusMSB269882 amoenusMSB269957 amoenusMSB269966 amoenusMSB292133 amoenusMSB292142 amoenusMSB292143 amoenusMSB292145 amoenusNEOR184 amoenusRBCM19566 amoenusRBCM20037 amoenusZM12231 amoenusZM13088 amoenusZM13089 amoenusZM16576 amoenusZM16577 amoenusZM16578 amoenusZM16579 amoenusZM16593 amoenusZM16594 amoenusZM16599 cratericusMSB269928 cratericusZM13093 cratericusZM13101 cratericusZM13102 cratericusZM13103 cratericusZM16601,
        cratericus_A:   amoenusJMG008 amoenusJRD028 amoenusMSB225555 amoenusMSB225676 amoenusMSB269634 amoenusMSB269635 amoenusMSB269636 amoenusMSB269638 amoenusMSB269640 amoenusMSB269856 amoenusMSB269857 amoenusMSB269858 amoenusMSB274470 amoenusMSB274506 amoenusZM12446 cratericusJRD111 cratericusJRD121 cratericusJRD122 cratericusJRD123 cratericusMSB264026 cratericusNDH002 cratericusNDH003 cratericusNDH007 cratericusNDH009 cratericusNDH011 cratericusNDH014 cratericusNDH018 cratericusNDH019 cratericusNDH020 cratericusNDH021 cratericusNDH022 cratericusNDH023 cratericusNDH024 cratericusNDH025 cratericusNDH026 cratericusNDH029 cratericusNDH030 cratericusNDH031 cratericusNDH032 cratericusNDH033 cratericusNDH034 cratericusNDH035 cratericusNDH036 cratericusNDH037 cratericusNDH038 cratericusZM13090 cratericusZM13094 cratericusZM13097 cratericusZM13098 cratericusZM13100 cratericusZM13105 cratericusZM13693 cratericusZM13694 cratericusZM13696 cratericusZM13740 cratericusZM13741 cratericusZM13742 cratericusZM13754 cratericusZM13755 cratericusZM13756 cratericusZM16546 cratericusZM16547 cratericusZM16548 cratericusZM16549 cratericusZM16555 cratericusZM16556 cratericusZM16580 cratericusZM16581 cratericusZM16582 cratericusZM16584 cratericusZM16585 cratericusZM16587 cratericusZM16588 cratericusZM16600 cratericusZM16602,
        cratericus_B:   amoenusMSB274473 amoenusMSB274502 amoenusMSB274505 cratericusJRD109 cratericusJRD110 cratericusJRD112 cratericusNDH013,
        ruficaudus:     ruficaudusJMG016 ruficaudusJMG134 ruficaudusJMS122 ruficaudusJMS201 ruficaudusJMS203 ruficaudusJRD142 ruficaudusZM16595 ruficaudusZM16596 ruficaudusZM16597 ruficaudusZM16598,
        minimus:        cratericusMSB197856 minimusMSB268068 minimusNDH041 minimusNDH045 minimusNDH063 minimusNDH064 minimusNDH072 minimusNDH084 minimusZM10995 minimusZM11027 minimusZM11417 minimusZM12230 minimusZM12232 minimusZM12233 minimusZM12277 minimusZM13061 minimusZM16544 minimusZM16545,
        ;
END;
```
Species Tree paup block: run_svdquartets_ddRads_speciestree.sh
Execute run_svdquartets_ddRads_speciestree.sh

```
Execute run_svdquartets_ddRads_lineagetree.sh
```
here is what our paupBlock looks like:
```
begin paup;
   log start file=svdQ.158ind.speciesTree.log replace;
   execute ./populations.snps.min4_Species_partitions.nexus;
   outgroup minimusMSB268068 minimusNDH041 minimusNDH045 minimusNDH063 minimusNDH064 minimusNDH072 minimusNDH084 minimusZM10995 minimusZM11027 minimusZM114>
   SVDquartets evalQuartets=all bootstrap=yes nreps=1000 taxpartition=Tamias_species nthreads=12 mrpFile=svdQ.158ind.SpeciesTree_qall_b1000;
   savetrees file=svdQ.158ind.SpeciesTree.b1000.tre format=Newick root=yes supportValues=nodeLabels;
end;
```
## Caster
Caster infers a species trees from site patterns using a statistical framework that is consistent with the multispecies coalescent model. We use both caster-sites (uses single-site patterns) and caster-pair (uses pairs of sites). See [Zhang et al. 2025](https://www.science.org/doi/10.1126/science.adk9688) for a full description of the method.
```
# Caster-site
/opt/ASTER/bin/caster-site -t 24 -r 16 -s 16 -i ../populations.snps.min4.phy -o caster-sites.158ind.lineageTree.tre 2> caster-sites_158ind.lineageTree.log
# Caster-pair
/opt/ASTER/bin/caster-site -t 24 -r 16 -s 16 -i ../populations.snps.min4.phy -o caster-pair.158ind.lineageTree.tre 2> caster-pair_158ind.lineageTree.log
```
## population level D-statistic using [Genomics General](https://github.com/simonhmartin/tutorials/blob/master/ABBA_BABA_whole_genome/README.md)

Convert to genotype file
```
/home/nh253642e/software/genomics_general/VCF_processing/parseVCF.py -i Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.150ind.vcf.gz --skipIndels | bgzip > Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.150ind.geno.gz
```
# calculate genome wide allele frequencies
we need to first determine the frequency of the derived allele in each populaton at each polymorphic site in the genome. 
The input file has already been filtered to contain only bi-allelic sites. We use the provided script (freq.py) to calculate alellic frequencies.
The script requires we define populatuions. In our case, we are setting each cratericus A pop as an independent pop. 
By setting --target derived we obtain the frquency of the derived allele in each population at each site. 
This is based on using the final population specified (T. miniums) as the outgroup. 
Sites at which this population is not fixed for the ancestral state are discarded.

# craA populations 
```
/home/nh253642e/software/genomics_general/freq.py -g Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.150ind.geno.gz \
-p cra_1 -p cra_2 -p cra_3 -p cra_4 -p cra_5 \
-p cra_6 -p cra_7 -p cra_8 -p cra_9 -p cra_10 \
-p cra_11 -p cra_12 -p cra_13 -p cra_14 \
-p cra_15 -p cra_16 -p cra_17 -p cra_18 -p cra_19 \
-p cra_20 -p cra_B -p amoenus -p ruficaudus -p minimus \
--popsFile tamias_150ind_RADs_craA_pops.txt --target derived \
-t 24 -o Tamias_MultiSmpl_Filtered_plink.newID.ldkb1r0.8.150ind_cra_A_populations.tsv.gz
```
Finally, we can use the Rscript [ddrad_dstat.R](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/1.%20ddRAD%20Processing%20and%20Analyses/ddrad_dstats.R) to calculate site pattern counts, D-stat and calculate Z scores. 

#### Go back to [main page](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/tree/main)
