# Clear workspace
rm(list=ls())

# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)
library(Biostrings)  # for reading FASTA files

# set seed
set.seed(123)

# Step 1: Read in the aligned sequences and metadata
fasta_file <- "pData/sars2_subsampled_unclassfil_filtered_refaligned_noN_trimmed.fasta"

# Read the aligned sequences
fasta <- readDNAStringSet(fasta_file)
keep_indices = c()
for (i in 1:length(fasta)) {
  seq <- fasta[[i]]  # Extract the sequence
  seq_char <- as.character(seq)  # Convert the DNAString to a character vector
  
  # Compute the number of 'N's
  prop_n <- sum(str_count(seq_char, "-")) / nchar(seq_char)
  
  # check if the time has two "-" in the header
  # if so, remove the sequence
  header = strsplit(names(fasta)[[i]], split="\\|")[[1]]
  # get teh final one
  tmp = header[length(header)]
  # check if it has two "-"
  if (sum(str_count(tmp, "-")) != 2) {
    next
  }
  
  if (prop_n <= 0.001) {
    keep_indices <- c(keep_indices, i)  # Keep track of the indices of valid sequences
  }
}

# randomly pick 5k of the keep indices
keep_indices <- sample(keep_indices, 5000)

# Step 5: Write the new FASTA file with the keep_indices only
valid_seqs <- fasta[keep_indices]

output_file <- "fasta/sars_all.fasta"
writeXStringSet(valid_seqs, output_file)

# Output confirmation
cat("New FASTA file has been written:", output_file, "\n")



# Step 1: Read in the aligned sequences and metadata
fasta_file <- "pData/mpox_unclassfil_filtered_refaligned_noN_trimmed.fasta"

# Read the aligned sequences
fasta <- readDNAStringSet(fasta_file)
keep_indices = c()
for (i in 1:length(fasta)) {
  seq <- fasta[[i]]  # Extract the sequence
  seq_char <- as.character(seq)  # Convert the DNAString to a character vector
  
  # Compute the number of 'N's
  prop_n <- sum(str_count(seq_char, "-")) / nchar(seq_char)
  
  # check if the time has two "-" in the header
  # if so, remove the sequence
  header = strsplit(names(fasta)[[i]], split="\\|")[[1]]
  # get teh final one
  tmp = header[length(header)]
  # check if it has two "-"
  if (sum(str_count(tmp, "-")) != 2) {
    next
  }
  
  if (prop_n <= 0.001) {
    keep_indices <- c(keep_indices, i)  # Keep track of the indices of valid sequences
  }
}

# Step 5: Write the new FASTA file with the keep_indices only
valid_seqs <- fasta[keep_indices]

output_file <- "fasta/mpxv_all.fasta"
writeXStringSet(valid_seqs, output_file)

# Output confirmation
cat("New FASTA file has been written:", output_file, "\n")
