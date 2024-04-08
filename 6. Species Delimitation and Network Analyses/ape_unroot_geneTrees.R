#!/usr/bin/env Rscript
library(ape)

tree_list<-list.files(pattern="*contree") #rooted gene trees from IQ-Tree
 
for (i in 1:length(tree_list)){
    tr_i<-read.tree(tree_list[i])
    tr_unroot_i<-unroot(tr_i)
    outtree2=paste(tree_list[i],"_unrooted.tree",sep="")
    write.tree(tr_unroot_i,file=outtree2)
}
