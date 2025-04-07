
# Clear workspace
rm(list=ls())

burnin = 0.0911

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)


# Load required packages
library(Biostrings)  # for reading FASTA files
library(dplyr)       # for handling random sampling

# Define file directories and parameters
fasta_dir <- "fasta/"
xml_dir <- "xmls_long/"
template_dir <- "Template/"

template <- c('Default', 'Intervals', 'Targeted')


percentages <- c(1000)
cl <- c(10, 10) * 10^5 * 10
ll <- c(1,1) * 10^3 * 5

# Get list of FASTA files
fastafiles <- list.files(fasta_dir, pattern = "*.fasta", full.names = TRUE)

# remove h3n2_HA
# fastafiles <- fastafiles[!grepl("h3n2_HA", fastafiles)]
# fastafiles <- fastafiles[!grepl("infB", fastafiles)]
# fastafiles <- fastafiles[!grepl("sars_all", fastafiles)]
# fastafiles <- fastafiles[!grepl("wnv", fastafiles)]

set.seed(24354)  # Set seed for reproducibility


for (fasta_file in fastafiles) {
  # Read the FASTA file
  fasta <- readDNAStringSet(fasta_file)
  fname <- tools::file_path_sans_ext(basename(fasta_file))
  
  # # build a iqtree with the following commands
  # # system(['/opt/homebrew/bin/iqtree2 -nt 11 -s  raxml/'  Viruses{v} '_' num2str(y) '_' Segments{s} '_raxml.fasta -m GTR --prefix raxml/'  Viruses{v} '_' num2str(y) '_' Segments{s} '']);
  # system(paste0("/opt/homebrew/bin/iqtree2 -nt 11 -s ", fasta_file, " -m GTR --prefix tree/", fname))
  # ds
  
  for (perc in percentages) {
    # Randomly sample sequences
    total_samples = length(fasta)
    
    
    indices <- sample(1:length(fasta), size = min(total_samples, perc))
    
    for (tmpl in template) {
      # Read template file and create output XML file
      template_file <- file.path(template_dir, paste0(tmpl, ".xml"))
      output_file <- file.path(xml_dir, paste0(fname, "_", tmpl, "_", length(indices), "_rep0.xml"))
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
          # line = '\t<run id="mcmc" spec="coupledMCMC.CoupledMCMC" chains="4" deltaTemperature="0.0" optimise="false" resampleEvery="1000" chainLength="4e+07">\n'
          line = sub("insert_chain_length", as.character(cl[which(percentages == perc)]), line)
          output_lines <- c(output_lines, line)
        } else if (grepl("insert_logEvery", line)) {
          output_lines <- c(output_lines, sub("insert_logEvery", as.character(ll[which(percentages == perc)]), line))
        } else {
          output_lines <- c(output_lines, line)
        }
      }
      # Write output lines to XML file
      writeLines(output_lines, output_file)
      # make two more copies of the xml by replacing rep0 with rep1 and rep2
      # use copyfile
      system(paste0("cp ", output_file, " ", sub("_rep0.xml", "_rep1.xml", output_file)))
      system(paste0("cp ", output_file, " ", sub("_rep0.xml", "_rep2.xml", output_file)))
      
    }
    
  }
  

}
