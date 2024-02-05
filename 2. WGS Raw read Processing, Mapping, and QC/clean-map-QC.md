## Process raw reads with [FastP](https://github.com/OpenGene/fastp)
```
for i in *1.fastq.gz;
do
  name1=${i%_*};
  name2=$(echo $name1 | cut -d '_' -f 1);
  fastp -i "$name2"_R1_001.fastq.gz -I "$name2"_R2_001.fastq.gz -m --merged_out "$name2"_fastp_ME.fastq.gz --out1 "$name2"_fastp_R1.fastq.gz --out2 "$name2"_fastp_R2.fastq.gz --unpaired1 "$name2"_fastp_SE.fastq.gz --unpaired2 "$name2"_fastp_SE.fastq.gz --detect_adapter_for_pe --cut_front --cut_front_window_size 5 --cut_front_mean_quality 20 -l 25 -j "$name2"_fastp.json -h "$name2"_fastp.html 2> "$name2".log
done
````
 Create single sample directories for all cleaned reads
```
for i in *_fastp_PE.fastq.gz; 
do echo $i; 
barcode=$(echo $i | cut -d '_' -f1); 
mkdir $barcode; 
mv *$barcode* $barcode; 
done
```
## FastP QC using MultiQC
We can use multiQC to consolodate the FastP log files into a nice HTML to assess quality metrics of our cleaned read data
```
multiqc .
```

## Read mapping using [BWA-MEM](https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=true)
This single command will map the PE, SE, ME reads to the T. minimus reference; sort and merge; add read groups; and finally 
mark duplicates with picard tools.
- All samples processed here are from a single sequencing effort. IF you have multiple seq. efforts,
you may want to add read groups as an independent step for each sequencing effort so as to maintain proper RG IDs.  
```
for i in */*R1.fastq.gz;
do
        name1=${i%.*};
        name2=$(echo $name1 | cut -d '_' -f 2);
        name3=$(echo $name2 | cut -d '/' -f 2); # This is only used for add RG (removes path from sample_ID/sample_ID)
# PE reads
if [ ! -f "$name2"_PE.bam ]
then
        bwa mem -M -t 6 $REF "$name2"_fastp_R1.fastq.gz "$name2"_fastp_R2.fastq.gz | samtools view -Sb - > "$name2"_PE.bam
fi
# For PE merged reads
if [ ! -f "$name2"_ME.bam ]
then
        bwa mem -M -t 6 $REF "$name2"_fastp_ME.fastq.gz | samtools view -Sb - > "$name2"_ME.bam
fi
# For SE reads
if [ ! -f "$name2"_SE.bam ]
then
        bwa mem -M -t 6 $REF "$name2"_fastp_SE.fastq.gz | samtools view -Sb - > "$name2"_SE.bam
fi
# Merging and sort the bams
if [ ! -f "$name2"_merge.bam ]
then
        samtools merge "$name2"_merge.bam "$name2"_PE.bam "$name3"_ME.bam "$name2"_SE.bam
fi
# Add readgroups
if [ ! -f "$name2"_merge_addRG.bam ]
then
        picard -Xmx50g AddOrReplaceReadGroups I="$name2"_merge.bam O="$name2"_merge_addRG.bam SO=coordinate RGID=A00738 LB=ch2_WGS PL=illumina PU=misc SM=$name3 TMP_DIR=$TMP
fi
# Mark duplicates
if [ ! -f "$name2"_deduped.bam ]
then
        picard -Xmx50g MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 I="$name2"_merge_addRG.bam O="$name2"_deduped.bam METRICS_FILE="$name2"_dedupe_metrics.txt VALIDATION_STRINGENCY=LENIENT TMP_DIR=$TMP
fi
done
```
## Index final bams
```
for i in */*_deduped.bam;
do
	name1=${i%.*};
	name2=${name1%_*}; 
	samtools index "$name2"_deduped.bam
done
```
## QualiMAP bam QC to assess mapping rate and coverage
```
for i in */*deduped.bam;
do
        name1=${i%_*};
        qualimap bamqc -nt 40 --java-mem-size=200G -sd -bam "$name1"_deduped.bam -outformat HTML -outdir $OUT -outfile "$name1"_bamQC_report
done
```
As before, we can consolodate our Qualimap reports with multiqc
```
multiqc .
```
## Variant discovery using GATK version 3.8 
We are now ready to move on to variant calling using the GATK (v 3.8) pipeline.

First, we perform indel realignment
```
for i in */*deduped.bam;
do
	name1=${i%.*};
	name2=${name1%_*};

if [ ! -f $name2.realignment_targets.list ]
 then
 	echo "running"
	gatk3 -T RealignerTargetCreator -R $REF -I "$name2"_deduped.bam -o "$name2".realignment_targets.list
fi
if [ ! -f "$name2"_realigned.bam ]
then
	echo "running"
	gatk3 -T IndelRealigner -Xmx50g -R $REF -I "$name2"_deduped.bam -targetIntervals $name2.realignment_targets.list -o "$name2"_realigned.bam TMP_DIR=$TMP
fi
done
```

## We can now Remove intermediate files to maintain a small footprint
```
rm *_PE.ba* *_SE.ba* *_Merged.ba* *_addRG.ba* *deduped.ba* *AllCalls.vcf* *.bed *_consensus.fa
```
