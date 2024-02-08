## Multi sequence alignments 
Here, we will create mutli sequence alignments (MSA) for phylogenetic analyses. There are two primary parts to MSA creation
1) MSA for Maximum Likelihood analysis of the concatenated data set. 
2) MSAs for gene tree analyses. 

We will also treat the X-linked scaffolds and autosomes seperately. For our referenece, X-linked scaffolds are scf 102, scf 163, and scf 370.

To begin, we are going to create per scaffold MSAs. We can do this quite easily because we have single sample FASTA files that all share the exact same coordinate system as the original reference (see previous section Variant Discovery and FASTA Consensus). Therefore, we can just create a per scaffold alignment using a custom R script that will take each FASTA file and split out each scaffold and then cat to a new scaffold MSA for each individual. 

I added all of the single sample FASTA to a single directory and then call the R script [BioConduct_Contig.aln.R](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/4.%20Multi%20Sequence%20Alignments/BioConduct_Contig.aln.R)
Requires R 3.4
```
Rscript BioConduct_Contig.aln.R
```
After running this, we now have scaffold level MSAs that do not require any downstream alignment tools because the coordinate system is retained across all samples. 

For each scaffold we use msa_split from the [phast](https://academic.oup.com/bib/article/12/1/41/244593?login=true) software package, to split each scaffold into 50 kb "gene tree" windows (we also analysed 100 kb windows but found no difference in results).

```
for i in *.fasta;
do
        name1=${i%.*}; 
        msa_split -i FASTA -o FASTA -r "$name1"_50kb --windows 50000,0 "$name1".fasta
done
```
Next we will us [AMAS](https://github.com/marekborowiec/AMAS) to trim each 50 kb window to include only sites with 75% or greater sample coverage. 
```
for i in *.fa;
do
        name1=${i%.*};
        name2=${name1%.*};
        AMAS.py trim -t 0.75 -f fasta -d dna -u fasta -i "$name2".fasta -o "$name2"_trim.75.fa
done
```
At this point, I want to exclude any 50 kb window that has fewer than 25 kb sites after trimming. 
```
for i in *.fa ;
do
    name1=${i%_*};
bioawk -c fastx 'length($seq) >99999 {print ">"$name"\n"$seq}' "$name1"_0.75trim.fa > "$name1"_0.75trim.fasta
done
```

This results in ready to analzye gene tree MSAs using IQ-Tree (see phylogenetic analyses section). Before we move on, we also want to create concatenated sequence alignments for the autosomal scaffolds and X-linked scaffolds for ML concat analysis. This is easy to do using AMAS:

For autosomal:
```
AMAS.py concat -f fasta -d dna -i *.fasta --part-format raxml -u phylip -t tamias_12ind_auto_trim75.phylip
```
For X-linked scaffolds:
```
AMAS.py concat -f fasta -d dna -i *.fasta --part-format raxml -u phylip -t phyllotis_45ind_X-linked_trim75.phylip
```
Now, we have gene tree alignments for gene tree analyses and species delimitation as well as a concatenated alignments for autosomes and X-linked scaffolds for ML analyses.
