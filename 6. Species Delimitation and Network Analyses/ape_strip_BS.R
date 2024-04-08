library(ape)
tree_list<-list.files(pattern="*contree")
 
for (i in 1:length(tree_list)){
    tr_i<-read.tree(tree_list[i])
    tr_i$node.label <- NULL
    outtree2=paste(tree_list[i],"_nBS.tree",sep="")
    write.tree(tr_i,file=outtree2)
}
