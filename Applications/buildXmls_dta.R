
# Clear workspace
rm(list=ls())

burnin = 0.0911

# Set the directory to the directory of the file
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
df_regions <- data.frame(
  location = c(
    "New_York","Tennessee","Michigan","NovaScotia","South_Carolina",
    "Georgia","Massachusetts","California","Utah","Florida",
    "New_York_City","New_Hampshire","Colorado","Texas","Idaho",
    "Ohio","Minnesota","South_Dakota","Hawaii","British_Columbia",
    "Manitoba","Delaware",NA,"Baltimore","Missouri",
    "Hawaii-TAMC","Nevada","Rhode_Island","North_Carolina","Costa_Rica",
    "Honduras","Virginia","Washington","Pennsylvania","Los_Angeles_County",
    "Maine","Montana","Iowa","Illinois","North_Dakota",
    "New_Jersey","West_Virginia","Wisconsin","District_Of_Columbia","Arizona",
    "Kansas","Mexico_City","MIXCO","Oregon","Maryland",
    "New_Mexico","Vermont","Nebraska","Indiana","Human",
    "Alabama","Saskatchewan","Guatemala","Kentucky","Connecticut",
    "SouthAustralia","Wyoming","Nunavut","Trinidad","USA",
    "Mexico","Canada","Louisiana","Alberta","Arkansas",
    "Ontario","Oklahoma","Panama","Barbados","HONDURAS",
    "Mississippi","DE-DHSS-2","SALAMA","Manitoba744683","North_America",
    "Quebec","New_Brunswick","DE","Puerto_Rico","Jamaica",
    "GUATEMALA","Sanarate","San_Agustin_Acasaguastlan","San_Cristobal_Acasaguastlan","Mixco",
    "San_Antonio_La_Paz","CHINAUTLA","SACATEPEQUEZ","RIO_HONDO","SAN_MARCOS_",
    "San_Pedro_Ayampuc","VILLA_CANALES","Manitoba744748","Manitoba744744","Manitoba745430",
    "Nicaragua","DistrictOfColumbia","District_of_Columbia","St._Vincent_And_Grenadines","El_Salvador",
    "Alaska","FORT-DE-FRANCE","Dominican_Republic","CALIFORNIA","PA",
    "Belize","Martinique","SANARATE","Chisec","Coban",
    "AL","CHIQUIMULA","Malacatan","Chiquimulilla","Unknown",
    "Chiquimula","Canberra","DISTRICT_OF_COLUMBIA","Nova_Scotia","San_Jose_Pinula",
    "San_Marcos","Tecpan","SAN_MARCOS","Panzos","San_Jose_La_Arada",
    "SUCHITEPEQUEZ","San_Juan_Sacatepequez","Chinautla","PATZUN","SAN_PEDRO_SACATEPEQUEZ",
    "ZACAPA","Palencia","CHICHICASTENANGO","PAJAPITA","SAN_ANTONIO_LA_PAZ",
    "Villa_Nueva","Santa_Catarina_Pinula","Ipala","Cuilapa","Manitoba739850",
    "Manitoba716522","Newfoundland","Manitoba720336","St_Kitts","Guadeloupe",
    "POINTE-A-PITRE","TECPAN","_Pennsylvania","Bermuda","Dominica",
    "Anguilla"
  ),
  final_region = c(
    #  1) New_York
    "USA_East",
    #  2) Tennessee
    "USA_Central",
    #  3) Michigan
    "USA_Central",
    #  4) NovaScotia
    "Canada",
    #  5) South_Carolina
    "USA_East",
    #  6) Georgia
    "USA_East",
    #  7) Massachusetts
    "USA_East",
    #  8) California
    "USA_West",
    #  9) Utah
    "USA_West",
    # 10) Florida
    "USA_East",
    # 11) New_York_City (treated like New_York)
    "USA_East",
    # 12) New_Hampshire
    "USA_East",
    # 13) Colorado
    "USA_West",
    # 14) Texas
    "USA_Central",
    # 15) Idaho
    "USA_West",
    # 16) Ohio
    "USA_Central",
    # 17) Minnesota
    "USA_Central",
    # 18) South_Dakota
    "USA_Central",
    # 19) Hawaii
    "USA_West",
    # 20) British_Columbia
    "Canada",
    # 21) Manitoba
    "Canada",
    # 22) Delaware
    "USA_East",
    # 23) NA
    NA,
    # 24) Baltimore (Maryland)
    "USA_East",
    # 25) Missouri
    "USA_Central",
    # 26) Hawaii-TAMC (still Hawaii)
    "USA_West",
    # 27) Nevada
    "USA_West",
    # 28) Rhode_Island
    "USA_East",
    # 29) North_Carolina
    "USA_East",
    # 30) Costa_Rica
    "Central_America",
    # 31) Honduras
    "Central_America",
    # 32) Virginia
    "USA_East",
    # 33) Washington
    "USA_West",
    # 34) Pennsylvania
    "USA_East",
    # 35) Los_Angeles_County (California)
    "USA_West",
    # 36) Maine
    "USA_East",
    # 37) Montana
    "USA_West",
    # 38) Iowa
    "USA_Central",
    # 39) Illinois
    "USA_Central",
    # 40) North_Dakota
    "USA_Central",
    # 41) New_Jersey
    "USA_East",
    # 42) West_Virginia
    "USA_East",
    # 43) Wisconsin
    "USA_Central",
    # 44) District_Of_Columbia
    "USA_East",
    # 45) Arizona
    "USA_West",
    # 46) Kansas
    "USA_Central",
    # 47) Mexico_City
    "Mexico",
    # 48) MIXCO (Guatemala)
    "Central_America",
    # 49) Oregon
    "USA_West",
    # 50) Maryland
    "USA_East",
    # 51) New_Mexico
    "USA_West",
    # 52) Vermont
    "USA_East",
    # 53) Nebraska
    "USA_Central",
    # 54) Indiana
    "USA_Central",
    # 55) Human
    NA,
    # 56) Alabama
    "USA_Central",
    # 57) Saskatchewan
    "Canada",
    # 58) Guatemala
    "Central_America",
    # 59) Kentucky
    "USA_Central",
    # 60) Connecticut
    "USA_East",
    # 61) SouthAustralia
    NA,
    # 62) Wyoming
    "USA_West",
    # 63) Nunavut
    "Canada",
    # 64) Trinidad
    "Caribbean",
    # 65) USA
    "USA",
    # 66) Mexico
    "Mexico",
    # 67) Canada
    "Canada",
    # 68) Louisiana
    "USA_Central",
    # 69) Alberta
    "Canada",
    # 70) Arkansas
    "USA_Central",
    # 71) Ontario
    "Canada",
    # 72) Oklahoma
    "USA_Central",
    # 73) Panama
    "Central_America",
    # 74) Barbados
    "Caribbean",
    # 75) HONDURAS
    "Central_America",
    # 76) Mississippi
    "USA_Central",
    # 77) DE-DHSS-2
    "USA_East",
    # 78) SALAMA (Guatemala)
    "Central_America",
    # 79) Manitoba744683
    "Canada",
    # 80) North_America (treated as USA per keywords)
    "USA",
    # 81) Quebec
    "Canada",
    # 82) New_Brunswick
    "Canada",
    # 83) DE (Delaware?)
    "USA_East",
    # 84) Puerto_Rico
    "Caribbean",
    # 85) Jamaica
    "Caribbean",
    # 86) GUATEMALA
    "Central_America",
    # 87) Sanarate (Guatemala)
    "Central_America",
    # 88) San_Agustin_Acasaguastlan (Guatemala)
    "Central_America",
    # 89) San_Cristobal_Acasaguastlan (Guatemala)
    "Central_America",
    # 90) Mixco (Guatemala)
    "Central_America",
    # 91) San_Antonio_La_Paz (Guatemala)
    "Central_America",
    # 92) CHINAUTLA (Guatemala)
    "Central_America",
    # 93) SACATEPEQUEZ (Guatemala)
    "Central_America",
    # 94) RIO_HONDO (Guatemala)
    "Central_America",
    # 95) SAN_MARCOS_ (Guatemala)
    "Central_America",
    # 96) San_Pedro_Ayampuc (Guatemala)
    "Central_America",
    # 97) VILLA_CANALES (Guatemala)
    "Central_America",
    # 98) Manitoba744748
    "Canada",
    # 99) Manitoba744744
    "Canada",
    #100) Manitoba745430
    "Canada",
    #101) Nicaragua
    "Central_America",
    #102) DistrictOfColumbia
    "USA_East",
    #103) District_of_Columbia
    "USA_East",
    #104) St._Vincent_And_Grenadines
    "Caribbean",
    #105) El_Salvador
    "Central_America",
    #106) Alaska
    "USA_West",
    #107) FORT-DE-FRANCE (Martinique region)
    "Caribbean",
    #108) Dominican_Republic
    "Caribbean",
    #109) CALIFORNIA
    "USA_West",
    #110) PA (Pennsylvania)
    "USA_East",
    #111) Belize
    "Central_America",
    #112) Martinique
    "Caribbean",
    #113) SANARATE (Guatemala)
    "Central_America",
    #114) Chisec (Guatemala)
    "Central_America",
    #115) Coban (Guatemala)
    "Central_America",
    #116) AL (ambiguous - not matched to 'alabama' by exact code => NA)
    NA,
    #117) CHIQUIMULA (Guatemala)
    "Central_America",
    #118) Malacatan (Guatemala)
    "Central_America",
    #119) Chiquimulilla (Guatemala)
    "Central_America",
    #120) Unknown
    NA,
    #121) Chiquimula (Guatemala)
    "Central_America",
    #122) Canberra (Australia)
    NA,
    #123) DISTRICT_OF_COLUMBIA
    "USA_East",
    #124) Nova_Scotia
    "Canada",
    #125) San_Jose_Pinula (Guatemala)
    "Central_America",
    #126) San_Marcos (Guatemala)
    "Central_America",
    #127) Tecpan (Guatemala)
    "Central_America",
    #128) SAN_MARCOS (Guatemala)
    "Central_America",
    #129) Panzos (Guatemala)
    "Central_America",
    #130) San_Jose_La_Arada (Guatemala)
    "Central_America",
    #131) SUCHITEPEQUEZ (Guatemala)
    "Central_America",
    #132) San_Juan_Sacatepequez (Guatemala)
    "Central_America",
    #133) Chinautla (Guatemala)
    "Central_America",
    #134) PATZUN (Guatemala)
    "Central_America",
    #135) SAN_PEDRO_SACATEPEQUEZ (Guatemala)
    "Central_America",
    #136) ZACAPA (Guatemala)
    "Central_America",
    #137) Palencia (Guatemala)
    "Central_America",
    #138) CHICHICASTENANGO (Guatemala)
    "Central_America",
    #139) PAJAPITA (Guatemala)
    "Central_America",
    #140) SAN_ANTONIO_LA_PAZ (Guatemala)
    "Central_America",
    #141) Villa_Nueva (Guatemala)
    "Central_America",
    #142) Santa_Catarina_Pinula (Guatemala)
    "Central_America",
    #143) Ipala (Guatemala)
    "Central_America",
    #144) Cuilapa (Guatemala)
    "Central_America",
    #145) Manitoba739850
    "Canada",
    #146) Manitoba716522
    "Canada",
    #147) Newfoundland
    "Canada",
    #148) Manitoba720336
    "Canada",
    #149) St_Kitts
    "Caribbean",
    #150) Guadeloupe
    "Caribbean",
    #151) POINTE-A-PITRE
    "Caribbean",
    #152) TECPAN (Guatemala)
    "Central_America",
    #153) _Pennsylvania
    "USA_East",
    #154) Bermuda
    "Caribbean",
    #155) Dominica
    "Caribbean",
    #156) Anguilla
    "Caribbean"
  ),
  stringsAsFactors = FALSE
)

# set the ones with USA as location to NA
df_regions$final_region[df_regions$location == "USA"] <- NA
# same for North_America
df_regions$final_region[df_regions$location == "North_America"] <- NA


# Load required packages
library(Biostrings)  # for reading FASTA files
library(dplyr)       # for handling random sampling

# Define file directories and parameters
fasta_dir <- "fasta/"
xml_dir <- "xmls_dta/"
template_dir <- "Template/"

template <- c('Default_dta', 'Targeted_dta')


percentages <- c(10000)
cl <- c(10, 10) * 10^5 * 10
ll <- c(1,1) * 10^3 * 50

# Get list of FASTA files
fastafiles <- list.files(fasta_dir, pattern = "*.fasta", full.names = TRUE)

# remove h3n2_HA
fastafiles <- fastafiles[!grepl("h3n2_HA", fastafiles)]
fastafiles <- fastafiles[!grepl("infB", fastafiles)]
fastafiles <- fastafiles[!grepl("sars_all", fastafiles)]
fastafiles <- fastafiles[!grepl("wnv", fastafiles)]
fastafiles <- fastafiles[!grepl("mpxv", fastafiles)]

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
    # split each name by / and take the second group as the location
    locs <- sapply(names(fasta), function(x) unlist(strsplit(x, "/"))[2])
    # for each locs, get the corresponding final_region
    locs <- df_regions$final_region[match(locs, df_regions$location)]
    
    
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
        } else if (grepl("insert_traits", line)){
          for (idx in indices) {
            seq <- fasta[idx]
            seq_id <- names(seq)
            loc=locs[idx]
            # if loc is NA, then make it a ?
            if(is.na(loc)){
              loc = "?"
            }
            if (idx == indices[length(indices)]) {
              output_lines <- c(output_lines, paste0("\t\t\t\t\t\t\t", seq_id, "=", loc))
            } else {
              output_lines <- c(output_lines, paste0("\t\t\t\t\t\t\t", seq_id, "=", loc, ","))
            }
          }
        } else if (grepl("insert_code_map", line)) {
          unique_locs = unique(locs)
          # make code map of typ codeMap="chi=0,dog=1,duc=2,kid=3,pig=4,? = 0 1 2 3 4"
          cmap = ""
          # remove na from unique_locs
          unique_locs = unique_locs[!is.na(unique_locs)]
          for (loc in unique_locs) {
            cmap = paste0(cmap, loc, "=", which(unique_locs == loc)-1, ",")
          }
          # set ?=0 1 2 3...
          str = seq(0, length(unique_locs)-1, 1)
          cmap = paste0(cmap, "?=", paste0(str, collapse=" "), collapse=",")
          tmp_line = sub("insert_code_map", cmap, line)
          # sub insert_no_states with length of unique_locs
          output_lines <- c(output_lines, sub("insert_no_traits", as.character(length(unique_locs)), tmp_line))

        } else if (grepl("insert_no_traits", line)) {
          # replace insert_freqs with 1/length(unique_locs) such that they add to 1
          freq_vals = rep(1/length(unique_locs), length(unique_locs))
          # make sure that the sum of freq_vals is 1
          freq_vals[length(freq_vals)] = 1 - sum(freq_vals[-length(freq_vals)])
          line_tmp = sub("insert_freqs", paste0(freq_vals, collapse=" "), line)
          output_lines <- c(output_lines, sub("insert_no_traits", as.character(length(unique_locs)), line_tmp))
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
