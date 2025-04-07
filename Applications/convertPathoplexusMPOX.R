# Clear workspace
rm(list=ls())

# Load necessary libraries
library(readr)
library(dplyr)
library(stringr)
library(Biostrings)  # for reading FASTA files

# Step 1: Read in the aligned sequences and metadata
fasta_file <- "pathoplexusData/mpox_aligned-nuc_2025-02-09T1921.fasta"
metadata_file <- "pathoplexusData/mpox_metadata_2025-02-09T1920.tsv"

# Read the aligned sequences
fasta <- readDNAStringSet(fasta_file)

# Read the metadata
metadata <- read_tsv(metadata_file)

# Step 2: Filter metadata to remove entries where the sampling time is not in 'yyyy-mm-dd' format
metadata <- metadata %>%
  filter(str_detect(sampleCollectionDate, "^\\d{4}-\\d{2}-\\d{2}$"))

# Step 3: Remove sequences with more than 10% N's using a for loop
valid_seqs <- list()  # Initialize an empty list to store valid sequences
valid_headers <- character()  # Initialize an empty character vector to store the new headers

for (i in 1:length(fasta)) {
  seq <- fasta[[i]]  # Extract the sequence
  seq_char <- as.character(seq)  # Convert the DNAString to a character vector
  
  # Compute the number of 'N's
  prop_n <- sum(str_count(seq_char, "N")) / nchar(seq_char)
  
  if (prop_n <= 0.10) {
    # Get the corresponding sampling time from the metadata (match by sequence name)
    seq_name = strsplit(names(fasta)[i], "\\.")[[1]][[1]]
    
    sampling_time <- metadata$sampleCollectionDate[match(seq_name, metadata$accession)]
    
    # Check if sampling time is not NA and format the header as accession|sampling_time
    if (!is.na(sampling_time)) {
      valid_seqs[[length(valid_seqs) + 1]] <- seq  # Keep the valid sequence
      valid_headers <- c(valid_headers, paste(seq_name, sampling_time, sep = "|"))
    }
  }
}

# Step 4: Create a new DNAStringSet with the new headers
valid_seqs <- DNAStringSet(valid_seqs)  # Convert the list of valid sequences back to DNAStringSet
names(valid_seqs) <- valid_headers  # Assign the new headers to the valid sequences

# Step 5: Write the new FASTA file with the updated headers
output_file <- "fasta/mpxv_all.fasta"
writeXStringSet(valid_seqs, output_file)

# Output confirmation
cat("New FASTA file has been written:", output_file, "\n")
