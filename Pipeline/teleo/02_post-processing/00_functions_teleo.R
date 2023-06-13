library(rentrez)
set_entrez_key("e0a14b0fc98a9389890653925881e0f04107") # My API key
Sys.getenv("ENTREZ_KEY") 

"%ni%" <- Negate("%in%")

# ---------------------------------------------------------------------------------------------------------------- # 

# Function to assemble MOTU table and taxonomy

assemble_data <- function(table_otu,
                          taxo_otu){
  
  # Normalize tables
  taxo_otu$amplicon <- str_trim(as.character(taxo_otu$definition))
  taxo_otu$sequence <- toupper(taxo_otu$sequence)                   
  
  # Checks that all MOTUs from the table are in the taxo csv (The opposite is not true due to sequence loss while cleaning rapidrun data)
  verif_motu <- table_otu$amplicon %in% taxo_otu$amplicon
  if(length(verif_motu[verif_motu == FALSE]) > 0) stop(message("Error: Some MOTU are in table but not in taxonomy file"))
  
  # Merge data
  table_taxo <- merge(taxo_otu, table_otu, by = c("amplicon", "sequence"))
  
  # Extract all sample names to prepare formatting
  samples <- table_taxo %>%
    select(starts_with("SPY"), starts_with("CNEG"), starts_with("Other"), starts_with("Blank")) %>%
    colnames()
  
  # Format to have each sequence of the samples in a single table
  merged <- table_taxo %>%
    gather(key = "sample_name", value = "count_reads", all_of(samples)) %>%
    filter(count_reads > 0)
  
  # Standardise and rename the % ID column
  database_col_name <- grep("best_identity", colnames(merged), value = T)
  merged[ ,"best_identity_database"] <- merged[[database_col_name]]
  merged[[database_col_name]] <- NULL
  
  return(merged)
}

# ---------------------------------------------------------------------------------------------------------------- # 

# Function to clean tag-jump

clean_tag_jump <- function(file_edna, file_other,
                           tag_jump_value = 0.001){
  
  # Add project type
  file_edna$clean <- "Project"
  file_other$clean <- "Other"
  
  file_all <- rbind(file_edna, file_other)
  
  # Cleaning
  file_edna_clean <- file_all %>%
    group_by(run, amplicon) %>%
    mutate(somme_tot = sum(count_reads), 
           threshold = tag_jump_value * somme_tot) %>%
    ungroup() %>%
    filter(count_reads > threshold) %>%
    filter(clean == "Project") %>%
    select(-clean)
  
  # Discarded MOTUs
  file_all_discarded <- file_all %>%
    group_by(run, amplicon) %>%
    mutate(somme_tot = sum(count_reads), 
           threshold = tag_jump_value * somme_tot) %>%
    ungroup() %>%
    filter(count_reads < threshold) %>%
    filter(clean == "Project") %>%
    select(-clean)
  
  return(list(file_edna_clean, 
              file_all_discarded))
}

# ---------------------------------------------------------------------------------------------------------------- # 

# Function to clean index-hoping

clean_index_hoping <- function(file_edna,
                               file_blank){
  
  # Check if 'batch' field in metadata exists, if not make a dummy considering all runs 
  if(is.null(metadata$batch)){
    file_edna[, "batch"] <- "A"
    file_blank[, "batch"] <- "A"
  }

  # Add project type
  file_edna$clean <- "other_or_project"
  file_blank$clean <- 'blank'

  all <- rbind(file_edna, file_blank)

  # Calculate intra-library tag-jump
  tag_jump_value <- 0.001

  table_counts <- all %>%
    # TAG JUMP
    group_by(run, sequence) %>%
    mutate(somme_tot = sum(count_reads),
           threshold = tag_jump_value * somme_tot) %>%
    ungroup() %>%
    filter(count_reads > threshold) %>%
    # TAG JUMP OVER
    group_by(batch, plate, sequence) %>%
    summarise(n_reads_tots = sum(count_reads),
              n_reads_other_project = sum(count_reads[clean == "other_or_project"]),
              n_reads_blanks = sum(count_reads[clean == "blank"])) %>%
    mutate(threshold_blank = n_reads_blanks/n_reads_tots) %>%
    filter(n_reads_blanks == 0 | n_reads_blanks > 10) %>%
    ungroup() %>%
    left_join(x = ., y = all %>% distinct(sequence, order_name), by = 'sequence')

  # Retain the maximal threshold for each batch
  threshold_blank_df <- table_counts %>%
    filter(threshold_blank < 1) %>% # remove sequences with only blank sequenced
    group_by(batch) %>%
    summarise(threshold_blank = max(threshold_blank))
  
  # Warning if threshold is too high 
  if(max(threshold_blank_df$threshold_blank) > 0.01) stop(message("Error: The cleaning threshold using blanks seems too high"))
  
  # Clean the project_data
  project_counts <- all %>%
    # Joint
    left_join(., threshold_blank_df) %>%
    group_by(batch, plate, sequence) %>%
    mutate(somme_tot_plate_batch = sum(count_reads), 
           threshold_plate_batch = threshold_blank * somme_tot_plate_batch,
           discard = ifelse(count_reads < threshold_plate_batch, "yes", "no")) %>%         
    # Clean for output
    ungroup() %>% 
    filter(project == "Port") %>% 
    select(-clean)
  
  # Filter
  project_clean <- project_counts %>%
    filter(discard == "no") %>%
    select(colnames(project_data))
  project_discarded <- project_counts %>%
    filter(discard == "yes")
  threshold_blank_df$threshold_blank <- round(threshold_blank_df$threshold_blank, 5)
  
  return(list(project_clean,
              project_discarded,
              threshold_blank_df))
}

# ---------------------------------------------------------------------------------------------------------------- # 

# Function to add complete classification, from superkingdom to species

# Create an archive to avoid running a long time ! 
# This outputs a list, the [[1]] is the file, the [[2]] is the updated archive file 

add_classification <- function(edna_file,
                     archive_file = NULL){
  require(taxize) 
  
  # Control if archive file is null or not
  if(is.null(archive_file)){
    archive_file <- setNames(data.frame(matrix(ncol = 11, nrow = 0)), c("class", "class_name",
                                                                        "order", "order_name",
                                                                        "family", "family_name",
                                                                        "genus", "genus_name",
                                                                        "species", "species_name",
                                                                        "scientific_name"))
  }
  
  # All names
  taxa_list <- unique(edna_file$scientific_name)
  
  # Are they in archive file? 
  taxa_list_present <- taxa_list[taxa_list %in% archive_file$scientific_name]
  taxa_list_query <- taxa_list[taxa_list %ni% archive_file$scientific_name]
  
  # Get the taxonomy
  list_classification <- classification(taxa_list_query, db = "ncbi", rows = 1)
  
  # Transform into data frame
  list_classification <- map_df(list_classification, ~ as.data.frame(.x), .id = "initial_name")
  taxa_list <- unique(list_classification$initial_name)
  
  # Add complete classification
  for (i in taxa_list){
    taxon_classification_complete <- list_classification %>%
      filter(initial_name == i)

    # Add class
    taxon_class <- taxon_classification_complete %>%
      filter(rank == "class") %>%
      rename(class_name = name) %>%
      rename(class = id) %>%
      select(class, class_name)
    
    # Add order
    taxon_order <- taxon_classification_complete %>%
      filter(rank == "order") %>%
      rename(order_name = name) %>%
      rename(order = id) %>%
      select(order, order_name)
    
    # Add family
    taxon_family <- taxon_classification_complete %>%
      filter(rank == "family") %>%
      rename(family_name = name) %>%
      rename(family = id) %>%
      select(family, family_name)
    
    # Add genus
    taxon_genus <- taxon_classification_complete %>%
      filter(rank == "genus") %>%
      rename(genus_name = name) %>%
      rename(genus = id) %>%
      select(genus, genus_name)
  
    # Add species
    taxon_species <- taxon_classification_complete %>%
      filter(rank == "species") %>%
      rename(species_name = name) %>%
      rename(species = id) %>%
      select(species, species_name)
    
    # Classification of the taxon (NA when rank is absent)
    taxon_classification <- do.call(cbind, lapply(c(taxon_class,
                                                    taxon_order,
                                                    taxon_family,
                                                    taxon_genus,
                                                    taxon_species,
                                                    scientific_name = i), function(x) if (length(x) == 0)  NA else x))
    
    # Add to archive
    archive_file <- rbind(archive_file, taxon_classification)
  }
  
  # Add the complete classification to the dataframe
  edna_with_classif <- edna_file %>%
    left_join(., archive_file, by=c("scientific_name"))
  
  return(list(edna_with_classif,
              archive_file))
}

# ---------------------------------------------------------------------------------------------------------------- # 

# Function to clean MOTUs based un user-input thresholds 

clean_motus <- function(edna_file,
                        min_reads = 10,
                        min_PCR = 1,
                        remove_PCR_blanks = TRUE,
                        remove_chimera = TRUE,
                        remove_not_classif_taxize = TRUE){

  # Fish classes kept 
  fish_class <- c("Actinopteri", "Chondrichthyes")
  
  # Verification that the class column exists
  if(remove_not_classif_taxize & !('class_name' %in% colnames(edna_file)) ) stop(message("Error: The column 'class_name' does not exist. Please use the 'add_classif_name' function before or set the 'remove_not_classif_taxize' option to FALSE."))
  
  # Isolate the PCR blanks for the filter
  amplicon_control <- edna_file %>%
    # Base filters: > 10 reads & Chimeras removal
    filter(count_reads > min_reads) %>% 
    `if`(remove_chimera, filter(., chimera == "N"), .) %>%
    # Keep only the control samples
    filter(!str_detect(sample_name, "SPY"))
  
  # Check condition of non NULL for amplicon_control 
  if(nrow(amplicon_control) == 0){
    amplicon_control <- data.frame(NA)
    amplicon_control$amplicon <- NA
  } else {
    
    # Print the count of MOTUs in blank
    message(paste0("\nThere is ", length(unique(amplicon_control$sequence)), " MOTUs in the blanks"))
  }
  
  # MOTUs cleaning
  edna_file_filtered <- edna_file %>%
    # Remove controls
    filter(str_detect(sample_name, "SPY")) %>% 
    # Read counts
    filter(count_reads > min_reads) %>% 
    # Chimeras
    `if`(remove_chimera, filter(., chimera == "N"), .) %>%
    select(-chimera) %>%
    # Remove MOTUs in blanks
    `if`(remove_PCR_blanks, filter(., !sequence %in% amplicon_control$sequence), .) %>%
    # Remove non fish taxa
    `if`(remove_not_classif_taxize, filter(., class_name %in% fish_class), .) %>%
    # PCR
    group_by(sequence) %>%
    mutate(n_PCR = n_distinct(sample_name)) %>%
    ungroup() %>%
    filter(n_PCR > min_PCR) %>%
    as.data.frame()
  
  # Isolate data present in only one PCR
  if (min_PCR > 0){
    edna_min_PCR <- edna_file %>%
      # Remove controls
      filter(str_detect(sample_name, "SPY")) %>% 
      # Read counts
      filter(count_reads > min_reads) %>% 
      # Chimeras
      `if`(remove_chimera, filter(., chimera == "N"), .) %>%
      select(-chimera) %>%
      # Remove MOTUs in blank
      `if`(remove_PCR_blanks, filter(., !sequence %in% amplicon_control$sequence), .) %>%
      # Remove non fish taxa
      `if`(remove_not_classif_taxize, filter(., class_name %in% fish_class), .) %>%
      # PCR
      group_by(sequence) %>%
      mutate(n_PCR = n_distinct(sample_name)) %>%
      ungroup() %>%
      filter(n_PCR == min_PCR) %>%
      as.data.frame()
  } else {
    edna_min_PCR <- NULL
  }
  
  return(list(edna_file_filtered,
              edna_min_PCR))
}

# ---------------------------------------------------------------------------------------------------------------- # 

# Function to clean MOTUs based on LULU
  
apply_lulu <- function(edna_file,
                       path_lulu,
                       match_lulu = 84,
                       co_occurence_lulu = 0.95){
  require(lulu)

  # Format for LULU function
  otutab <- edna_file %>%
    select(definition, count_reads, sample_name) %>%
    spread(sample_name, count_reads) %>%
    replace(., is.na(.), "0") %>%
    as.data.frame()
  
  # Put the OTU name in row name
  otutab2 <- otutab
  otutab$definition <- NULL
  otutab <- otutab %>%
    mutate_if(is.character, as.numeric)
  rownames(otutab) <- otutab2$definition
  
  # Import matchlist
  matchlist <- read.table(paste(path_lulu, "match_list_teleo.txt", sep=""), header = F, as.is = T, stringsAsFactors=F)
  
  # Run LULU 
  lulu_edna <- lulu(otutab, matchlist, minimum_match = match_lulu, minimum_relative_cooccurence = co_occurence_lulu)
  
  # Clean the logs 
  system("rm lulu.log*")
  
  return(lulu_edna)
}

# ---------------------------------------------------------------------------------------------------------------- # 