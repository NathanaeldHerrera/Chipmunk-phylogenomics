# Gene Tree Alignments and Analysis
This pipeline is designed to create gene tree alignments from the Tamias WGS data where coordinates will be preserved across the genome.
We are using 100 kb stepped intervals across each scaffold (for those that are 10 mb or larger).

There are a few ways to do this. More recently, I have incorportated the Tree House Explorer (THEx) pipleine created by [Sam Harris](http://www.eutherialab.org/treehouseexplorer/) of the [Murphy lab](http://www.eutherialab.org/) into my workflow as it is well documented and can streamline many elements.

For our referenece, X-linked scaffolds are scf 102, scf 163, and scf 370 (excluded to being < 10 Mb). For all analayses, we treat the X-linked scaffolds seperately from the autosomes. 

First, we will generate 50 kb stepped windows across each scaffold. This will retain coordinate possition. For THEx, All of these steps can be incorporated into a single config file. For simplicity and claritry, I show each step here as an individual unit. 

The directory: multialignmentdir contains each of the scaffold level MSAs generated in the [previous section](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/4.%20Multi%20Sequence%20Alignments/MSA_generation.md)

```
thexb --minifastas \
  -i ./multiAlignmentDir/ \
  -o ./outdir \
  --window_size 100kb
```
Now, we can use [TrimAL](https://vicfero.github.io/trimal/), to trim our 100 kb alignments. I am using a gap threashold of 0.75, a minimum sequence length of 25 kb, and specify to drop filtered windows (This will create an empty file that will maintiain coordinate across the geneome, regardless if there is alignment or not). 
```
thexb --trimal \
  -i ./outdir/windowed_fastas \
  -o ./outdir \
  --trimal_min_seq_len 25000bp \
  --trimal_gap_threshold 0.75
```
Next, we can run [IQ-Tree](http://www.iqtree.org/doc/) using the GTR + GHOST model of rate heterogeneity (see methods section of manuscript for an explaination of model selection) on our individual gene trees.

```
thexb --iqtree \
  -i ./outdir/trimal_filtered_windows/ \
  -o ./outdir/ \
  --iqtree_model GTR*H4 \
  --iqtree_bootstrap 1000 \
  --tv_file_name TreeViewer_input_tamias_file.xlsx
```
This will yeild two outputs. First, we will have ML gene trees that we can use for downstream species tree delimiation approaches (see: section 6. [Species Delimiation](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/6.%20Species%20Delimitation%20and%20Network%20Analyses/species_delimiation_%26_network_analyses.md)) and second, a .xlsx file of each coordinate specific locus tree and the topology. THEx can then be used to rank the topologies and plot them to look at the distribution of locus trees across the geneome. 

I outline the species delimitation pipeline is [section 6](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/6.%20Species%20Delimitation%20and%20Network%20Analyses/species_delimiation_%26_network_analyses.md). Here, we will continue with the gene tree analysis in THEx. 

We now want to assign rank order (the proportion each locus tree topology is represented across our genome) using THEx topobinner.
```
thexb --topobinner \
  -i ./outdir/TreeViewer_input_tamias_file.xlsx \
  -o ./outdir/
```
As per THEx documentation: "Topobinner is a tool used to organize and label identical tree topologies based on RF-distance. Trees are treated as equal topologies when their RF-distance is equal to 0. Trees are binned and then labeled based on their whole genome frequency, meaning the most frequent topology in the genome will be labeled Tree1, then Tree2 for the second most frequent topology, etc. Once this step is completed, an updated Tree Viewer file with binned topologies will be generated and placed in the defined output directory."

This gives us our final output .xlsx file for plotting in THEx. In conjunction to the ranked locus tree file, you must provide a bed file of the genome. 
