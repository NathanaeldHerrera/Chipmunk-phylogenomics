## Raw read processing using Nexflow workflow
# Workflow Overview
This pipeline automates read cleaning, alignment, QC, quantification, and reporting.
There are two files, the workflow and a config file for job submission via slurm.
[clean-map-qc-dna.nf]()
[clean-map-qc-dna.config]()

Steps:
Fastp – Adapter trimming and quality filtering
BWA-mem + Picard – Paired-end read alignment and read group tagging
samtools index
Qualimap – BAM-level QC
MultiQC – Summary reports (Fastp + Qualimap)

Inputs:

Paired-end FASTQ files (with _1.fq.gz/_2.fq.gz naming)
Reference genome index (bwa-mem, samtools faidx, picard)
Command Example:

nextflow run clean-map-qc-dna.nf -c clean-map-qc-dna.config -with-trace
Key Outputs:

results/fastp_cleaned_reads/: Trimmed reads
results/bams/: Sorted and indexed BAM files
results/qualimap/: Alignment quality reports

## Indel realignment
We are going to use UnifiedGenotyper from the GATK. It runs off of a java version that does not play well with updated java version required by Nextflow so we carry this step out independently. 
```
conda activate gatk3.8-env

for i in *deduped.bam;
	name1=${i%.*};
	name2=${name1%_*};
	if [ ! -f $name2.realignment_targets.list ];              
	then
		echo "running"
		# Create a list of intervals to realign
		gatk -T RealignerTargetCreator -R $REF -I "$name2"_deduped.bam -o "$name2".realignment_targets.list
	fi
	if [ ! -f "$name2"_realigned.bam ];                      
	then
		echo "running"
		# Perform indel realignment
		gatk -T IndelRealigner -Xmx50g -R $REF -I "$name2"_deduped.bam -targetIntervals $name2.realignment_targets.list -o "$name2"_realigned.bam TMP_DIR=$TMP
	fi
done
```
#### Go back to [main page](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/tree/main)
