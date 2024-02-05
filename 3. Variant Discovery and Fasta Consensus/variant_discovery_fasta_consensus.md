## Variant discovery using GATK version 3.8 
We are now ready to move on to variant calling using the GATK (v 3.8) pipeline.
This pipeline will generate generate a single sample VCF where we emit all sites (that is, it will make a call at every position in our reference regardless if it is homozygous REF or a no call ./.)
using UnifiedGenotyper (GATK v3.8)

```
for i in */*realigned.bam;
do
	name1=${i%.*};
	name2=${name1%_*};
if [ ! -f "$name2"_AllCalls.vcf ]
then
	echo "running"
	gatk3 -Xmx100g -T UnifiedGenotyper -nct 20 -R $REF -I "$name2"_realigned.bam --genotyping_mode DISCOVERY --output_mode EMIT_ALL_SITES -stand_call_conf 30 -o "$name2"_AllCalls.vcf
fi
done
```
## Soft variant filtering using GATK Variantfiltration
This gives us our final filtered single sample VCF files.
Filters applied are
- Quality by depth (QD)
- RMS Mapping quality (MQ)
- Per site depth (DP)

See here for information on [filters](https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variantsjava)  
``` 
for i in */*_AllCalls.vcf;
do
        name1=${i%.*};
        name2=${name1%_*};

if [ ! -f "$name2"_AllCalls.filtered.vcf ]
then
        echo "running"
        gatk3 -Xmx100g -T VariantFiltration -R $REF -V "$name2"_AllCalls.vcf --filterExpression "QD < 2.0" --filterName "LowQD" --filterExpression "MQ < 10.0" --filterName "LowMQ" --filterExpression "DP < 5" --filterName "LowDP" -l ERROR -o "$name2"_AllCalls.filtered.vcf
fi
done
```
## FASTA consensus
We use FastaAlternateReferenceMaker to create a single sample FASTA consensus file. This will inject all calls into the reference 
WITH IUPAC ambiguity codes (-IUPAC flag).
```
for i in */*_AllCalls.filtered.vcf.gz;
do
        name1=${i%.*}; 
        name2=${name1%_*};
        name3=$(echo $name2 | cut -d '/' -f 2);
         
if [ ! -f "$name2"_consensus.fa ]
then
        echo "running"
        gatk3 -Xmx50g -T FastaAlternateReferenceMaker -R $REF -V "$name2"_AllCalls.filtered.vcf.gz -o "$name2"_consensus.fa -IUPAC "$name3"
fi
done
```

## Mask out no call positions and filtered sites
This call is specific to UnifedGenotpyer output
It is IMPORTANT to make sure each grep is specified to the actual filter names specified above.

No calls
```
for i in */*_AllCalls.filtered.vcf.gz;
do
        name1=${i%.*};
        name2=${name1%_*};

if [ ! -f "$name2"_nocalls.combined.bed ]
then
        echo "running"
        zgrep "\./\." "$name2"_AllCalls.filtered.vcf.gz | awk '{{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > "$name2"_nocalls.combined.bed
fi
done
```
Low QD
```
for i in */*_AllCalls.filtered.vcf.gz;
do
        name1=${i%.*};
        name2=${name1%_*};

if [ ! -f "$name2"_LowQD.combined.bed ]
then
        echo "running"
        zgrep "LowQD" "$name2"_AllCalls.filtered.vcf.gz | awk '{{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > "$name2"_LowQD.combined.bed
fi
done
```
Low MQ
```
for i in */*_AllCalls.filtered.vcf.gz;
do
        name1=${i%.*};
        name2=${name1%_*};

if [ ! -f "$name2"_LowMQ.combined.bed ]
then
        echo "running"
        zgrep "LowMQ" "$name2"_AllCalls.filtered.vcf.gz | awk '{{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > "$name2"_LowMQ.combined.bed
fi
done
```
Low DP
```
for i in */*_AllCalls.filtered.vcf.gz;
do
        name1=${i%.*};
        name2=${name1%_*};

if [ ! -f "$name2"_LowDP.combined.bed ]
then
        echo "running"
        zgrep "LowDP" "$name2"_AllCalls.filtered.vcf.gz | awk '{{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > "$name2"_LowDP.combined.bed
fi
done
```
Low Qual
```
for i in */*_AllCalls.filtered.vcf.gz;
do
        name1=${i%.*};
        name2=${name1%_*};

if [ ! -f "$name2"_LowQual.combined.bed ]
then
        echo "running"
        zgrep "LowQual" "$name2"_AllCalls.filtered.vcf.gz | awk '{{OFS="\t"; if ($0 !~ /\#/); print $1, $2-1, $2}}' | bedtools merge -i - > "$name2"_LowQual.combined.bed
fi
done
```
Now we can combine all of the above failed sites .bed files into a single all positions to mask bed file
```
for i in */*_AllCalls.filtered.vcf.gz;
do
        name1=${i%.*};
        name2=${name1%_*};

if [ ! -f "$name2"_all.combined.bed ]
then
        echo "running"
        cat "$name2"_nocalls.combined.bed "$name2"_LowQD.combined.bed "$name2"_LowMQ.combined.bed "$name2"_LowDP.combined.bed "$name2"_LowQual.combined.bed  > "$name2"_all.combined.bed
fi
done
```
Sort and merge final records. This creates 'all_positions_to_mask.bed' which is what we use to mask our final FASTA files. It does this by converting failed sites/ no calls 
to an 'N' in the fasta file.
```
for i in */*_AllCalls.filtered.vcf.gz;
do
        name1=${i%.*};
        name2=${name1%_*};

if [ ! -f "$name2"_all_positions_to_mask.bed ]
then
        echo "running"
        bedtools sort -i "$name2"_all.combined.bed | bedtools merge -i - > "$name2"_all_positions_to_mask.bed
fi
done
```
Unfortunaltey FastaAlternateReferenceMaker modifies the refence scaffold names to a simplified sequential order. We can use [rename.py](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/3.%20Variant%20Discovery%20and%20Fasta%20Consensus/rename.py) to make sure reference genome scaffold
names match single sample FASTA scaffold names.
```
for i in */*_consensus.fa ;
do
        name1=${i%.*};
        name2=${name1%_*};
        python rename.py $REF "$name2"_consensus.fa
done    
```
Finally, we mask out the S@*t positions from the final consensus.fa file
```
for i in */*_consensus.fa;
do
        name1=${i%.*};
        name2=${name1%_*};
        bedtools maskfasta -fi "$name2"_consensus.fa -fo "$name2"_consensus.masked.fa -bed "$name2"_all_positions_to_mask.bed
done
```
As a check we can count the number of nucleotides in each fasta file. They should be identical.
```
for i in *_consensus.masked.fa ;
do
        name1=${i%.*};
        name2=${name1%_*};
        grep -v ">" "$name2"_consensus.masked.fa | wc | awk '{print $3-$1}'
done
```
