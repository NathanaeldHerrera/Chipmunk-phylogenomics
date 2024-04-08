## Gene Tree Alignments and Analysis
This pipeline is designed to create "gene tree" alignments from the Tamias WGS data where coordinates will be preserved across the genome.
We are using a 50 kb stepped interval across each scaffold (for those that are 10 mb or larger).

1) MSA for Maximum Likelihood analysis of the concatenated data set. 
2) MSAs for gene tree analyses. 

We will also treat the X-linked scaffolds and autosomes seperately. For our referenece, X-linked scaffolds are scf 102, scf 163, and scf 370.

To begin, we are going to create per scaffold MSAs. We can do this quite easily because we have single sample FASTA files that all share the exact same coordinate system as the original reference (see previous section Variant Discovery and FASTA Consensus). Therefore, we can just create a per scaffold alignment using a custom R script that will take each FASTA file and split out each scaffold and then cat to a new scaffold MSA for each individual. 

I added all of the single sample FASTA to a single directory and then call the R script [BioConduct_Contig.aln.R](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/4.%20Multi%20Sequence%20Alignments/BioConduct_Contig.aln.R)
(Requires R 3.4).
```
Rscript BioConduct_Contig.aln.R
```
We now have scaffold level MSAs that do not require any downstream alignment tools because the coordinate system is retained across all samples. 
