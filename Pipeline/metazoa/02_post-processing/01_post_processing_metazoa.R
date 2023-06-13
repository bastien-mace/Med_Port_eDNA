# This script reads eDNA metabarcoding data from the SWARM clustering pipeline 
# It cleans and formats data correctly for the metazoa marker

# Packages
library(tidyverse)
library(dplyr)
library(seqinr)

# Source functions
source("00_functions_metazoa.R")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 1 - Assemble & Clean
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

# Clean index-hoping and tag-jump

# Path to the swarm outputs
swarm_path <- "../01_swarm/"
  
# Metadata file 
metadata <- read.table(paste0(swarm_path, "all_samples_metazoa.csv"),
                       sep=";", stringsAsFactors = F, header = F)
colnames(metadata) <- c("plate", "run", "sample_name", "project", "marker", "batch")

# Check for no NA values on the run
verif_metadata <- !is.na(metadata$run)
if(length(verif_metadata[verif_metadata == FALSE]) > 0 ) stop(print("Error: Some samples do not have metadata fields"))

# Open files 
files <- list.files(swarm_path, pattern = "(.*)metazoa(.*)tsv")

# ----- # Project files
project_table <- read.table(paste0(swarm_path, grep("Port(.*)table", files, value = T)),
                            sep="\t", stringsAsFactors = F, header = T)
  
project_taxo <- read.table(paste0(swarm_path, grep("Port(.*)ecotag", files, value = T)),
                           sep="\t", stringsAsFactors = F, header = T)
# Assemble
project_data <- assemble_data(table_otu = project_table, taxo_otu = project_taxo) %>%
  left_join(., metadata)

# ----- # Other files
other_table <- read.table(paste0(swarm_path, grep("Other(.*)table", files, value = T)),
                          sep="\t", stringsAsFactors = F, header = T)
other_taxo <- read.table(paste0(swarm_path, grep("Other(.*)ecotag", files, value = T)),
                         sep="\t", stringsAsFactors = F, header = T)
  # Assemble
other_data <- assemble_data(table_otu = other_table, taxo_otu = other_taxo) %>%
  left_join(., metadata)
  
# ----- # Blank files 
# Apply the code using the blanks only if those are present 
if(length(grep("Blank(.*)", files, value = T)) == 0){message(paste0("There is no Blank files for the metazoa data"))}

if(length(grep("Blank(.*)", files, value = T)) != 0){

  blank_table <- read.table(paste0(swarm_path, grep("Blank(.*)table", files, value = T)),
    sep="\t", stringsAsFactors = F, header = T)
  blank_taxo <- read.table(paste0(swarm_path, grep("Blank(.*)ecotag", files, value = T)),
    sep="\t", stringsAsFactors = F, header = T)
  
  # Assemble
  blank_data <- assemble_data(table_otu = blank_table, taxo_otu = blank_taxo) %>%
    left_join(., metadata)
  
  # Clean index-hoping
  project_cleaned <- clean_index_hoping(file_edna = rbind(project_data, other_data), 
                                        file_blank = blank_data)
  
  project_data <- project_cleaned[[1]]
  threshold <- project_cleaned[[3]]
  message(paste0("There is ", dim(project_cleaned[[2]])[1], " lost MOTUs during index-hopping cleaning"))
  message(paste0("\nThe Blank threshold for the ", threshold$batch, " batch in the metazoa project is ", threshold$threshold_blank))
}

# Clean tag-jump
project_data_tag <- clean_tag_jump(file_edna = project_data, file_other = other_data)

# Store in list 
list_step1 <- project_data_tag[[1]]

message(paste0("There are ", dim(project_data_tag[[2]])[1], " lost MOTUs during tag-jump cleaning"))

save(list_step1, file = "metazoa_workspace.Rdata")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 2: Complete taxonomy
# load("metazoa_workspace.Rdata")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

# Remove null elements from list
list_step1 <- list_step1[!sapply(list_step1, is.null)]

# Clean column names 
columns_to_remove <- c("amplicon", "OTU", "id", "total", "count",  "cloud", "length", "abundance", "spread", "identity", "taxonomy", "references",
                       "taxid_by_db.db_ncbi_mito", "best_match.db_ncbi_mito", "match_count.db_ncbi_mito", "species_list.db_ncbi_mito",
                       "plate", "run", "project", "marker", "batch", "somme_tot", "threshold",
                       "order", "order_name", "family", "family_name", "genus", "genus_name", "species", "species_name")

# Add classification
file_taxo <- as.data.frame(list_step1) %>%
  select(-one_of(columns_to_remove)) %>%
  add_classification()

list_step2 <- file_taxo[[1]]
archive_classif_ncbi <- file_taxo[[2]]

# Save data and the new archive file 

save(archive_classif_ncbi, file = "metazoa_archive.Rdata")
save(list_step1, list_step2, file = "metazoa_workspace.Rdata")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 3 - Data filtering
# load("metazoa_workspace.Rdata")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

# Clean with default values
data_cleaned <- clean_motus(list_step2, min_PCR = 0) # turn min_PCR to 0 if replicates were pooled

list_step3 <- data_cleaned[[1]]

# Extract data and species present in only one PCR and which are discarded (useless if min_PCR = 0)
if (is.null(data_cleaned[[2]]) == FALSE){
  list_step3_0PCR <- data_cleaned[[2]]
  species_0_PCR <- list_step3_0PCR %>%
    bind_rows() %>%
    filter(kingdom_name %in% c("Metazoa")) %>%
    distinct(sequence,best_identity_database, scientific_name)
  message(paste0("There is ", length(unique(species_0_PCR$scientific_name)), " lost animal taxa during the PCR cleaning step"))
}

save(list_step1, list_step2, list_step3, file = "metazoa_workspace.Rdata")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 4 - LULU post-clustering
# load("metazoa_workspace.Rdata")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

# Format data
list_step3$definition <- gsub(" ", "", list_step3$definition)

# Create necessary files
otu_seq <- list_step3 %>%
  distinct(definition, sequence) %>%
  mutate(sequence = toupper(sequence))

# Convert and write to fasta
write.fasta(as.list(otu_seq$sequence), names = otu_seq$definition, file.out = "OTU_sequences_metazoa2.fasta", open = "w")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #
# Create the match list necessary for lulu in a bash terminal (UNIX OS system) using the blastn tool

# ----- Launch blastn in the terminal
# $ makeblastdb -in OTU_sequences_metazoa.fasta -parse_seqids -dbtype nucl

# ----- Create the matchlist
# $ blastn -db OTU_sequences_metazoa.fasta -outfmt '6 qseqid sseqid pident' -out  match_list_metazoa.txt -task megablast -qcov_hsp_perc 80 -perc_identity 84 -query OTU_sequences_metazoa.fasta

# Transfer the matchlist into the lulu folder
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- #

# Apply lulu
list_lulu <- apply_lulu(list_step3, path_lulu = paste0(getwd(), "/"))

# MOTUs to keep
lulu_motu_keep <- list_lulu$curated_otus

# Filter out discarded MOTUs
list_step4 <- list_step3 %>%
  filter(definition %in% lulu_motu_keep)

# Compare the number of MOTUs before  and after LULU
count_motu_before_lulu <- length(unique(list_step3$definition))
count_motu_after_lulu <- length(unique(list_step4$definition))
message(paste0("There is ", count_motu_before_lulu - count_motu_after_lulu, " lost MOTUs after lulu curation"))

save(list_step1, list_step2, list_step3, list_step4, file = "metazoa_workspace.Rdata")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 5 - Last cleaning
# load("metazoa_workspace.Rdata")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

# Samples we want to keep
samples_to_keep <- unique(grep(
"SPY213525|SPY213527|SPY202505|SPY202519|SPY202520|SPY213512|SPY213513|SPY213514|SPY213515|SPY213516|SPY213517|SPY213518|SPY213519|SPY213520|SPY213521|SPY213522|SPY213523|SPY221739|SPY221740|SPY221741|SPY221742|SPY221743|SPY221744|SPY221745|SPY221746|SPY221747|SPY221748|SPY221749|SPY221750|SPY221751|SPY221753|SPY224119|SPY224120|SPY224121|SPY224122|SPY224123|SPY224124|SPY224125|SPY224126|SPY224127|SPY224128|SPY224129|SPY224130",
list_step4$sample_name, value = T))

list_step5 <- list_step4 %>%
  select(-c("n_PCR")) %>%
  filter(sample_name %in% samples_to_keep)

save(list_step1, list_step2, list_step3, list_step4, list_step5, file = "metazoa_workspace.Rdata")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Step 6 - Create a Jaccard matrix
# load("metazoa_workspace.Rdata")
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 

# Unnecessary columns to remove
columns_to_remove <- c("sequence",  "count_reads", "best_identity_database", "taxid",
                       "rank", "kingdom", "phylum", "class", "order", "family", "genus", "species",
                       "kingdom_name", "phylum_name", "class_name", "order_name", "family_name", "genus_name", "species_name")

# Build the table of MOTUs occurences among samples 
occurences_table <- list_step5 %>%
  spread(definition, scientific_name) %>%
  select(-one_of(columns_to_remove))

samples <- unique(occurences_table$sample_name)

list_step6 <- data.frame()

# Build the sample composition table (= Jaccard matrix)
for (i in 1:length(samples)){
  sublist <- filter(occurences_table, sample_name == samples[i]) %>% 
    subset(select = -c(sample_name)) %>%
    replace(!is.na(.), "1") %>%  # Non NA values are replaced by 1
    replace(is.na(.), "0") %>%   # NA values are replaced by 0
    lapply(as.numeric) %>%
    as.data.frame(check.names = F) 
  
  sublist_sum <- sublist %>%
    colSums() %>% # We sum for the whole sample
    t() %>%
    as.data.frame(check.names = F)
  
  rownames(sublist_sum) <- samples[i] # We give the sample names to the rownames
  list_step6 <- rbind(list_step6, sublist_sum)
}

# Taxonomy corresponding to MOTUs
motu_taxo <- list_step5 %>%
  select(-c(sample_name)) %>%
  unique()

unique_motu_taxo <- motu_taxo[!duplicated(motu_taxo$definition), ]

save(list_step1, list_step2, list_step3, list_step4, list_step5, list_step6, file = "metazoa_workspace.Rdata")

# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 
# Write output file
write.table(list_step6, file = "metazoa_jaccard.tsv", col.names = NA)
write.table(unique_motu_taxo, file = "metazoa_motu_taxo.tsv", col.names = NA)
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- # 