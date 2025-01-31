
# Clear workspace
rm(list=ls())

burnin = 0.0911

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# Load required packages
library(Biostrings)  # for reading FASTA files
library(dplyr)       # for handling random sampling

# remove all files in xmls
unlink("xmls", recursive = TRUE)
# make a new folder called xmls
dir.create("xmls")

# Define file directories and parameters
fasta_dir <- "fasta/"
xml_dir <- "xmls/"
template_dir <- "Template/"
percentages <- c(0.1, 0.5, 1)
template <- c('Default', 'Targeted', 'Intervals', 'Topology')
cl <- c(10, 20, 50, 100) * 10^5
ll <- c(10, 20, 50, 100) * 10^2

# Get list of FASTA files
fastafiles <- list.files(fasta_dir, pattern = "*.fasta", full.names = TRUE)

for (fasta_file in fastafiles) {
  # Read the FASTA file
  fasta <- readDNAStringSet(fasta_file)
  fname <- tools::file_path_sans_ext(basename(fasta_file))
  
  for (perc in percentages) {
    # Randomly sample sequences
    indices <- sample(1:length(fasta), size = round(length(fasta) * perc))
    
    for (tmpl in template) {
      # Read template file and create output XML file
      template_file <- file.path(template_dir, paste0(tmpl, ".xml"))
      output_file <- file.path(xml_dir, paste0(fname, "_", tmpl, "_", perc, ".xml"))
      
      template_lines <- readLines(template_file)
      output_lines <- character(0)  # Initialize empty vector to store output lines
      
      for (line in template_lines) {
        if (grepl("insert_sequences", line)) {
          for (idx in indices) {
            seq <- fasta[idx]
            seq_id <- names(seq)
            seq_value <- as.character(seq)
            output_lines <- c(output_lines, paste0("\t<sequence id=\"seq_", seq_id, "\" spec=\"Sequence\" taxon=\"", seq_id, "\" totalcount=\"4\" value=\"", seq_value, "\"/>\n"))
          }
        } else if (grepl("insert_chain_length", line)) {
          output_lines <- c(output_lines, sub("insert_chain_length", as.character(cl[which(percentages == perc)]), line))
        } else if (grepl("insert_logEvery", line)) {
          output_lines <- c(output_lines, sub("insert_logEvery", as.character(ll[which(percentages == perc)]), line))
        } else {
          output_lines <- c(output_lines, line)
        }
      }
      
      # Write output lines to XML file
      writeLines(output_lines, output_file)
    }
  }
}
