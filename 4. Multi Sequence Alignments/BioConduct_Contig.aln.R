# Load necessary libraries
library(Biostrings)
library(seqinr)

# List all FASTA files in the current directory with a .fa extension
fastafiles <- list.files(pattern = ".fa") 
print(fastafiles)  # Display the list of FASTA files

# Iterate through each FASTA file to process the sequences
for (fasta in fastafiles) {  
  sample <- readDNAStringSet(fasta) # Read the DNA sequences from the FASTA file
  name <- strsplit(fasta, "_")[[1]][1] # Extract the sample name from the file name
  
  # Loop through each contig (sequence) in the sample
  for (i in names(sample)) {
    a <- sample[i] # Extract the specific contig
    names(a) <- paste(name) # Set the name of the contig
    writeXStringSet(a, filepath = paste(i, ".fasta", sep = ''), append = TRUE, format = "fasta") # Write the contig to a new FASTA file, appending if it already exists
  }
}
