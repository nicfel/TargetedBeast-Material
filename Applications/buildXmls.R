
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

template <- c('Default', 'Targeted', 'Intervals')

# leafs = c("infB"=1978, "h3n2"=5000, "wnv"=2664, "mpxv"=1579)
n_samples <- round(exp(seq(log(50),log(200), length.out = 5)))
cl = n_samples  * 10^4*5
ll = round(n_samples * 10^1/4)

# Get list of FASTA files
fastafiles <- list.files(fasta_dir, pattern = "*.fasta", full.names = TRUE)

# remove files that contains h3n2_HA mpxv_all sars_all
# fastafiles <- fastafiles[!grepl("h3n2_HA", fastafiles)]
# fastafiles <- fastafiles[!grepl("infB", fastafiles)]
# fastafiles <- fastafiles[!grepl("sars_all", fastafiles)]
# fastafiles <- fastafiles[!grepl("wnv", fastafiles)]

for (fasta_file in fastafiles) {
  # Read the FASTA file
  fasta <- readDNAStringSet(fasta_file)
  fname <- tools::file_path_sans_ext(basename(fasta_file))
  
  for (perc in n_samples) {
    # Randomly sample sequences
    indices <- sample(1:length(fasta), size = round(perc))
    
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
          output_lines <- c(output_lines, sub("insert_chain_length", as.character(cl[which(n_samples == perc)]), line))
        } else if (grepl("insert_logEvery", line)) {
          output_lines <- c(output_lines, sub("insert_logEvery", as.character(ll[which(n_samples == perc)]), line))
        }else if (grepl("weight=\"0.05\"", line)){
          output_lines <- c(output_lines, sub("weight=\"0.05\"", "weight=\"0.25\"", line))   
        } else {
          output_lines <- c(output_lines, line)
        }
      }
      
      # Write output lines to XML file
      writeLines(output_lines, output_file)
    }
  }
}
