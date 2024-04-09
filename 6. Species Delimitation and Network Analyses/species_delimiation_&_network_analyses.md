# Speices Tree and Network Analyses
This pipeline outlines the three approaches I used to explore fine-scale patterns of gene tree-species tree discordance under a multispecies-coalescent model.
For ASTRAL and MP-EST, I used the gene trees where coordinate is retained from the previous section: [gene tree analysis](). Again, we analyze both atuosomal and X-linked scaffolds seperately.

## ASTRAL
Assuming sets of independent and accurately estimated gene trees, ASTRAL estimates an unrooted species tree by finding the species tree that has the maximum number of shared induced quartet trees with the given set of gene trees. See [Zhang et al. 2018](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2129-y) for a description of the theory behind ASTRAL.

For ASTRAL, we want unrooted gene trees.

To unroot gene trees we can use the R package [Ape](https://cran.r-project.org/web/packages/ape/) as implemented in the R script [ape_unroot_geneTrees.R](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/6.%20Species%20Delimitation%20and%20Network%20Analyses/ape_unroot_geneTrees.R). To run this, call script in the directory containing all .contree gene tree files from IQ-Tree.

Create a list of all gene trees for X-linked and autos:
```
find . -name "*contree" -exec cat '{}' ';' > tamias_50KB_conTrees_autosomes.tre
find . -name "*contree" -exec cat '{}' ';' > tamias_50KB_conTrees_x-linked.tre
```

We also want to collapse branches with excessively low bootstrap support. For this, I used [newick_utils](https://github.com/tjunier/newick_utils) to collapse branches with a BP < 10.  
```
nw_ed tamias_50KB_conTrees_autosomes.tre 'i & b<=10' o > tamias_50KB_conTrees_autosomes_BS10.tre
nw_ed tamias_50KB_conTrees_x-linked.tre 'i & b<=10' o > tamias_50KB_conTrees_x-linked_BS10.tre
```
ASTRAL output options:
Alternative posteriors only (-t 4): When this option is used, ASTRAL outputs three local posterior 
probabilities: one for the main topology, and one for each of the two alternatives (RS|LO and RO|LS, in that order). 
The posterior of the three topologies adds up to one. This is because of our locality assumption, 
which basically asserts that we assume the four groups around the branch (L, R, S, and O) are each 
correct and therefore, there are only three possible alternatives.

Full annotation (-t 2): When you use this option, for each branch you get a lot of different measurements:
q1,q2,q3: these three values show quartet support (as defined in the description of -t 1) for the main topology, 
the first alternative, and the second alternative, respectively.
f1, f2, f3: these are the total number of quartet trees in all the gene trees that support the main topology, the first alternative, and the second alternative, respectively.
pp1, pp2, pp3: these are the local posterior probabilities (as defined in the description of -t 4) for the main topology, the first alternative, and the second alternative, respectively.
QC: this is the total number of quartets defined around each branch.
EN: the effective number of genes for the branch.
### Autosomes 
Lineage Tree:
```
astral -t 4 -i tamias_50KB_conTrees_autosomes_BS10.tre -o tamias_50KB_autosomes_BS10_ASTRAL_LineageTree.tre | tee tamias_50KB_autosomes_BS10_ASTRAL_LineageTree.log
```
Species Tree:
```
astral -t 4 -i tamias_50KB_conTrees_autosomes_BS10.tre -a 12ind_species_list.txt -o tamias_50KB_autosomes_BS10_ASTRAL_SpeciesTree.tre | tee tamias_50KB_autosomes_BS10_ASTRAL_SpeciesTree.log
```
For the species, we have multiple individuals from the same species. In this case, we are going to force the species to be monophyletic by supplying a mapping file using the -a option. This mapping file is formatted as:

species_name1:individual_1,individual_2,...
species_name2:individual_1,individual_2,...

### X-linked scaffolds 
Lineage Tree:
```
astral -t 4 -i tamias_50KB_conTrees_x-linked_BS10.tre -o tamias_50KB_x-linked_BS10_ASTRAL_LineageTree.tre | tee tamias_50KB_x-linked_BS10_ASTRAL_LineageTree.log
```
Species Tree:
```
astral -t 4 -i tamias_50KB_conTrees_x-linked_BS10.tre -a 12ind_species_list.txt -o tamias_50KB_x-linked_BS10_ASTRAL_SpeciesTree.tre | tee tamias_50KB_x-linked_BS10_ASTRAL_SpeciesTree.log
```

## MP-EST
MP-EST uses the frequency of triplets of taxa among the gene trees to estimate the topology and branch lengths of the overall species tree. See [Lue et al. 2010]() for the theory behind MP-EST.

MP-EST takes in gene tree topologies without branch support. To do so, we will use the R script [ape_strip_BS.R](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/blob/main/6.%20Species%20Delimitation%20and%20Network%20Analyses/ape_strip_BS.R) to remove the branch support from the IQ-Tree gene trees. 

MP-EST takes ROOTED gene trees. NOTE, you can specify whether or not to root trees in IQ-Tree. IF you have unrooted topologies, you do not need to unroot with R script. IF you need to re-root gene trees, this can be done easily with newick_utils:
```
for i in *.contree_nBS.tree;
do
        name1=${i%.*};
        name2=${name1%_*};
        name3=${name2%.*};
        name4=${name3%.*};
nw_reroot -l "$name4".fa.contree_nBS.tree JPJ2177 JPJ412 > "$name4"_nBS_rooted.tre
done
```
create a list of gene trees for MP-EST
```
find . -name "*tre" -exec cat '{}' ';' > tamis_mpest_conTrees_autosomes_nBS_rooted.tre
find . -name "*tre" -exec cat '{}' ';' > tamis_mpest_conTrees_x-linked_nBS_rooted.tre
```
Run using v2.1
```
./mpest tamias_control_speciesTree
```

## SVDquartets
SVDquartets assesses support for quartets of taxa directly from site-pattern frequencies of variable sites only. This approach differs from summary methods, such as ASTRAL and MP-EST, because it does not independently estimate gene trees, avoiding the issue of gene-tree estimation error. For the theory behind SVDQ see [Chifman and Kubatko 2014](https://academic.oup.com/bioinformatics/article/30/23/3317/206559?login=true).

#### Go back to [main page](https://github.com/NathanaeldHerrera/Chipmunk-phylogenomics/tree/main)
