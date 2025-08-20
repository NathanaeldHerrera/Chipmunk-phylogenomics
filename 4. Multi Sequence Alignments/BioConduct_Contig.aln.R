library(Biostrings)
library(seqinr)
#rm(list=ls())

fastafiles <- list.files(pattern = ".fa") # this will pull any file with a .fa extension
fastafiles
### This will iterate through your fasta files and it will pull out each contig from each species file and append them to a 
### contig specific file. 
for (fasta in fastafiles) {
  sample <- readDNAStringSet(fasta)
  name <- strsplit(fasta, "_")[[1]][1]
    for(i in names(sample)){
#      print(sample[i])
#      print(name) 
      a <- sample[i]
#      print(a) 
      names(a) <- paste(name)
      writeXStringSet(a, filepath = paste(i, ".fasta", sep = ''), append = TRUE, format="fasta")
  }
}


# writeXStringSet(a, filepath = paste("Contig_",i, ".fasta", sep = ''), append = TRUE, format="fasta")
