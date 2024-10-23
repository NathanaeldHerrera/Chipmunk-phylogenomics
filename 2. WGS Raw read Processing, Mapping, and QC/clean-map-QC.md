## Process raw reads with [FastP](https://github.com/OpenGene/fastp)
```
# Loop through all files matching the pattern '*1.fastq.gz'
for i in *1.fastq.gz; 
do
  # Extract the base name by removing the trailing '_1.fastq.gz' from the filename
  name1=${i%_*};

  # Extract the sample name (first part before the underscore) from the base name
  name2=$(echo $name1 | cut -d '_' -f 1);

  # Run fastp to process the paired-end reads
  fastp \
    -i "$name2"_R1_001.fastq.gz \  # Input R1 file
    -I "$name2"_R2_001.fastq.gz \  # Input R2 file
    -m \                             # Enable automatic adapter detection and trimming
    --merged_out "$name2"_fastp_ME.fastq.gz \  # Output merged reads
    --out1 "$name2"_fastp_R1.fastq.gz \        # Output trimmed R1 reads
    --out2 "$name2"_fastp_R2.fastq.gz \        # Output trimmed R2 reads
    --unpaired1 "$name2"_fastp_SE.fastq.gz \   # Output unpaired R1 reads
    --unpaired2 "$name2"_fastp_SE.fastq.gz \   # Output unpaired R2 reads
    --detect_adapter_for_pe \                    # Automatically detect adapters for paired-end data
    --cut_front \                                # Enable trimming from the front of reads
    --cut_front_window_size 5 \                  # Size of the window for the front cut
    --cut_front_mean_quality 20 \                # Mean quality threshold for cutting
    -l 25 \                                      # Minimum length for reads after trimming
    -j "$name2"_fastp.json \                     # Output JSON report
    -h "$name2"_fastp.html \                     # Output HTML report
    2> "$name2".log                              # Redirect error messages to a log file
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
We can use multiQC to aggregate and visualize (as an HTML) the FastP log files to assess quality metrics of our cleaned read data
```
multiqc .
```

## Read mapping using [BWA-MEM](https://academic.oup.com/bioinformatics/article/25/14/1754/225615?login=true)
This single command will map the PE, SE, ME reads to the T. minimus reference; sort and merge; add read groups; and finally 
mark duplicates with picard tools.
- All samples processed here are from a single sequencing effort. IF you have multiple seq. efforts,
you may want to add read groups as an independent step for each sequencing effort so as to maintain proper RG IDs.  
```
for i in */*R1.fastq.gz; do
    name1=${i%.*};
    name2=$(echo $name1 | cut -d '_' -f 2);  # Extract sample name
    name3=$(echo $name2 | cut -d '/' -f 2);  #removes path

    # Step 1: Align paired-end (PE) reads
    if [ ! -f "$name2"_PE.bam ]; then
        # Align R1 and R2 using BWA MEM and convert to BAM using Samtools
        bwa mem -M -t 6 $REF "$name2"_fastp_R1.fastq.gz "$name2"_fastp_R2.fastq.gz | samtools view -Sb - > "$name2"_PE.bam
    fi
    # Step 2: Align merged reads
    if [ ! -f "$name2"_ME.bam ]; then
        # Align merged reads using BWA MEM
        bwa mem -M -t 6 $REF "$name2"_fastp_ME.fastq.gz | samtools view -Sb - > "$name2"_ME.bam
    fi
    # Step 3: Align single-end (SE) reads
    if [ ! -f "$name2"_SE.bam ]; then
        # Align single-end reads using BWA MEM
        bwa mem -M -t 6 $REF "$name2"_fastp_SE.fastq.gz | samtools view -Sb - > "$name2"_SE.bam
    fi
    # Step 4: Merge BAM files (PE, ME, SE)
    if [ ! -f "$name2"_merge.bam ]; then
        # Merge the different BAM files into a single BAM file
        samtools merge "$name2"_merge.bam "$name2"_PE.bam "$name3"_ME.bam "$name2"_SE.bam
    fi
    # Step 5: Add read groups to the merged BAM file
    if [ ! -f "$name2"_merge_addRG.bam ]; then
        picard -Xmx50g AddOrReplaceReadGroups \        # Using Picard tools with a maximum heap size of 50 GB to handle large files (this flag depends on HPC set-up and may not be required for some instances)
            I="$name2"_merge.bam \              # Specifies the input BAM file to be processed (merged BAM from earlier steps)
            O="$name2"_merge_addRG.bam \        # Defines the output BAM file name, which will include read group information
            SO=coordinate \                       # Sets the sorting order of the output BAM to 'coordinate' (sorted by genomic position)
            RGID=A00738 \                        # Assigns a unique read group ID
            LB=ch2_WGS \                         # Specifies the library name
            PL=illumina \                        # Indicates the platform used for sequencing
            PU=misc \                            # Sets the platform unit to 'misc', which may indicate a mixed or non-specific run
            SM=$name3 \                         # Specifies the sample name using the variable $name3, which holds the sample identifier
            TMP_DIR=$TMP                         # Designates a temporary directory for any temporary files generated during processing

    fi
    # Step 6: Mark duplicates in the merged BAM file
    if [ ! -f "$name2"_deduped.bam ]; then
        picard -Xmx50g MarkDuplicates \ 
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \  # Sets the maximum number of file handles for read ends mapping to optimize resource usage
            I="$name2"_merge_addRG.bam \
            O="$name2"_deduped.bam \
            METRICS_FILE="$name2"_dedupe_metrics.txt \
            VALIDATION_STRINGENCY=LENIENT \             # Allows for some leniency in validation errors to avoid halting the process for minor issues
            TMP_DIR=$TMP
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
        qualimap bamqc -nt 40 \
            --java-mem-size=200G \
            -sd \
            -bam "$name1"_deduped.bam \
            -outformat HTML \
            -outdir $OUT \
            -outfile "$name1"_bamQC_report
done
```
As before, we can consolodate our Qualimap reports with multiqc
```
multiqc .
```
## Indel realignment
This will result in our final mapped bams

```
for i in */*deduped.bam;
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

## We can now Remove intermediate bam files to maintain a small(er) footprint
```
rm *_PE.ba* *_SE.ba* *_ME.ba* *_addRG.ba* *deduped.ba*
```

#### Go back to [main page](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/tree/main)
