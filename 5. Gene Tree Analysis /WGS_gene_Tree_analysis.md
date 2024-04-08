## Gene Tree Alignments and Analysis
This pipeline is designed to create gene tree alignments from the Tamias WGS data where coordinates will be preserved across the genome.
We are using a 50 kb stepped intervals across each scaffold (for those that are 10 mb or larger).

There are a few ways to do this. More recently, I have incorportated the Tree House Explorer (THEx) pipleine created by [Sam Harris](http://www.eutherialab.org/treehouseexplorer/) of the [Murphy lab](http://www.eutherialab.org/) into my workflow as it is well documented and can streamline many elements. I will also include the pipeline elements that do not include THEx for the analysis of gene trees using IQ-Tree and {GNU parallel](https://www.gnu.org/software/parallel/).

For our referenece, X-linked scaffolds are scf 102, scf 163, and scf 370. For all analayses, we treat the X-linked scaffolds seperately from the autosomes. 

First, we will generate 50 kb stepped windows across each scaffold. This will retain coordinate possition. For THEx, All of these steps can be incorporated into a single config file. For simplicity and claritry, I show each step here as an individual unit. 

The directory: multialignmentdir contains each of the scaffold level MSAs generated in the [previous section]()

```
thexb --minifastas \
  -i ./multiAlignmentDir/ \
  -o ./outdir \
  --window_size 50kb
```
Now, we can use [TrimAL](), to trim our 50 kb alignments. I am using a gap threashold of 0.75 and a minimum sequence length of 25 kb and specify to drop filtered windows. 
```
thexb --trimal \
  -i ./outdir/windowed_fastas \
  -o ./outdir \
  --trimal_min_seq_len 25000bp \
  --trimal_gap_threshold 0.75
```
Next, we can run IQ-Tree on our individual gene trees.

```
thexb --iqtree \
  -i ./outdir/trimal_filtered_windows/ \
  -o ./outdir/ \
  --iqtree_model GTR*H4 \
  --iqtree_bootstrap 1000 \
  --tv_file_name TreeViewer_input_tamias_file.xlsx
```
NOTE: This is how I analyzed gene tree alignments before THEx was available.
Depending on compuational resources, you can also parallelize IQ-Tree using GNU parallel as shown here:
```
## Create a call list of gene trees
ls *.fa > tamias_50kb_contigs.txt
##IQ-Tree paralleliation
parallel -a tamias_50kb_contigs.txt -j 120 iqtree2 -s {} -m TEST -nt 1 -pre {} -bb 1000 -alrt 1000
```




