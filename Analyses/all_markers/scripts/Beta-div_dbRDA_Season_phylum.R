# Load the required packages
library(tidyverse)
library(reshape2)
library(vegan)
library(ggalt)
library(ggrepel)


# Extract MOTU phylum names
motu_phylum <- read.csv2("../data/all_motu_taxo.csv") %>%
  select(c("phylum_name")) %>%
  na.omit() %>%
  unique()
motu_phylum

# Extract MOTUs corresponding to these phylum
motu_table <- read.csv2("../data/all_motu_taxo.csv") %>%
  filter(phylum_name %in% motu_phylum$phylum_name) %>%
  select(c("definition", "phylum_name"))

# Load metadata
metadata <- read.csv2("../../00_metadata/metadata_pooled.csv") %>%
  filter(Habitat %in% c("Marina_Summer", "Marina_Autumn"))
  
# Extract where they occur
occurence_motu <- read.csv2("../data/all_jaccard_pooled.csv") %>%
  inner_join(., metadata) %>%
  select(-c("Habitat", "Season", "Longitude", "Latitude")) %>%
  column_to_rownames("Site") %>%
  select(motu_table$definition) %>%
  t(.) %>%
  as.data.frame() %>%
  rownames_to_column("definition") %>%
  inner_join(., motu_table) %>%
  select(-definition)


# Build a matrix with the number of MOTUs belonging to each phylum in each site
occurence_phylum <- data.frame(matrix(ncol = length(colnames(occurence_motu)) - 1, nrow = length(motu_phylum$phylum_name)))
colnames(occurence_phylum) <- c(colnames(occurence_motu)[1:length(colnames(occurence_motu)) - 1])
rownames(occurence_phylum) <- motu_phylum$phylum_name

## Calculate the abundance of each phylum in each site
for (i in motu_phylum$phylum_name){
  occurence_phylum_i <- occurence_motu %>%
    filter(phylum_name == i)
  occurence_phylum[i, ] <- colSums(occurence_phylum_i[, -c(dim(occurence_phylum_i)[2])])
}


# Prepare data for dbRDA
biodiv_table <- occurence_phylum %>%
  t(.) %>%
  as.data.frame() %>%
  rownames_to_column("Site") %>%
  inner_join(., metadata)

# Define dbRDA variables
## Select the Ys
data_dbRDA_phylum <- biodiv_table %>%
  select(-c("Site", "Habitat", "Season", "Longitude", "Latitude"))

## Select the Xs
meta_dbRDA_phylum <- biodiv_table %>%
  select(c("Site", "Habitat", "Season", "Longitude", "Latitude"))

meta_dbRDA_phylum$Longitude <- as.numeric(meta_dbRDA_phylum$Longitude)
meta_dbRDA_phylum$Latitude <- as.numeric(meta_dbRDA_phylum$Latitude)

row.names(data_dbRDA_phylum) <- meta_dbRDA_phylum$Site

## Remove absent MOTUs (those only present in reserves)
data_dbRDA_phylum <- data_dbRDA_phylum[, colSums(data_dbRDA_phylum) != 0]


# Set the model
dbRDA_phylum <- capscale(data_dbRDA_phylum ~ Season + Condition(Longitude + Latitude),
                         meta_dbRDA_phylum, distance ='bray')

# Get the inertia statistics                   
dbRDA_phylum

# Test the significance of the model
anova.cca(dbRDA_phylum, permutations = 9999)

# Calculate the adjusted-R²
RsquareAdj(dbRDA_phylum)


# Extract dbRDA data
sumdbRDA <- summary(dbRDA_phylum)

# Get the scores
scrs <- scores(dbRDA_phylum, display = c("species", "sites", "bp"), scaling = 2)

## Extract site names
site_names <- biodiv_table %>%
  arrange(match(Site, rownames(scrs$sites))) %>%
  pull(Site)

## Extract season names
season_names <- biodiv_table %>%
  arrange(match(Site, rownames(scrs$sites))) %>%
  pull(Season)

## Extract site scores and add site and habitat names
sites_centroids <- data.frame(scrs$sites)
sites_centroids$SITE <- site_names
sites_centroids$SEASON <- season_names

## Remove the suffix
sites_centroids[sites_centroids == "Marina_Summer"] <- "Summer"
sites_centroids[sites_centroids == "Marina_Autumn"] <- "Autumn"

sites_centroids$SITE <- gsub("*_Summer", "", sites_centroids$SITE)
sites_centroids$SITE <- gsub("*_Autumn", "", sites_centroids$SITE)

# Plot the data
dbRDA_site_scores <- ggplot(data = sites_centroids,
                    aes(x = CAP1, y = MDS1, color = SEASON, label = SITE)) +
  geom_point(size = 3) +
  geom_text_repel(color = "black") +
  geom_encircle(aes(group = SEASON, linetype = SEASON, fill = SEASON), s_shape = 1, expand = 0, alpha = 0.5, show.legend = FALSE) +
  scale_color_manual(values = c("#FDC086", "#BEAED4")) +
  scale_fill_manual(values = c("#FDC086", "#BEAED4")) +
  geom_hline(yintercept = 0, lty = "dotted") +
  geom_vline(xintercept = 0, lty = "dotted") +
  theme_classic() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 15, colour = "black"),
        legend.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        legend.position = c(0.9, 0.1),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "white", color = "black")) +
  labs(x = paste0("Constrained (", round(sumdbRDA$cont$importance[2,1]*100, 2), "%)"),
       y = paste0("Unconstrained - Axis 1 (", round(sumdbRDA$cont$importance[2,2]*100, 2), "%)"),
       color = "Habitat")
dbRDA_site_scores

# ggsave("../plots/Beta-div_dbRDA_Season_phylum_sites.tiff", dbRDA_site_scores, dpi = 600, width = 20, height = 15, units = "cm")


## Extract species scores and add species names
species_centroids <- data.frame(scrs$species)
species_centroids$SPECIES <- rownames(species_centroids)

### Identify the 10% most contributing species
abs_CAP1 <- abs(species_centroids$CAP1)
### Calculate the threshold for the top 10% values
threshold <- quantile(abs_CAP1, 0.9)
### Extract the indices of the top 10% values
top_indices <- which(abs_CAP1 >= threshold)
### Extract the corresponding MOTUs
species_centroids_10 <- species_centroids[top_indices, ]
species_centroids_10 <- cbind(species_centroids_10, Group = gsub("_OTU.*", "", species_centroids_10$SPECIES))

## Create a table with the markers and their corresponding phyla
markers_phylum <- motu_table %>%
  rename("SPECIES" = "phylum_name")
markers_phylum$definition <- gsub("teleo.*", "teleo", markers_phylum$definition)
markers_phylum$definition <- gsub("metazoa.*", "metazoa", markers_phylum$definition)
markers_phylum$definition <- gsub("euka2.*", "euka2", markers_phylum$definition)
markers_phylum$definition <- gsub("bact2.*", "bact2", markers_phylum$definition)

### Dereplicate the data
duplicate_cols <- duplicated(markers_phylum)
markers_phylum_unique <- markers_phylum[!duplicate_cols, ]

species_centroids_10 <- inner_join(species_centroids_10, markers_phylum_unique)


# Plot the data
dbRDA_species_scores <- ggplot() + 
  geom_segment(data = species_centroids, aes(x = 0, xend = CAP1, y = 0, yend = MDS1), col = "grey",
               arrow = arrow(length = unit(0.01,"npc"))) + # all species
  geom_segment(data = species_centroids_10, aes(x = 0, xend = CAP1, y = 0, yend = MDS1), col = "black",
               arrow = arrow(length = unit(0.01,"npc"))) + # most differentiated species
  geom_hline(yintercept = 0, lty = 2, col = "grey") +
  geom_vline(xintercept = 0, lty = 2, col = "grey") +
  geom_label_repel(data = species_centroids_10, aes(x = CAP1, y = MDS1, label = SPECIES, color = definition)) +
  guides(color = guide_legend(override.aes = list(labels = c("Prokaryotes", "Eukaryotes", "Metazoans", "Fish")),
                              keywidth = 4)) +
  scale_color_manual(values = c("#E78AC3", "#66C2A5", "#FC8D62", "#8DA0CB")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA),
        axis.text = element_text(colour = "black", size = 10),
        axis.title = element_text(size = 15, colour = "black"),
        legend.text = element_blank(),
        legend.title = element_text(size = 10, colour = "black", face = "bold"),
        legend.position = c(0.915, 0.135),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "white", color = "black")) +
  labs(x = paste0("Constrained (", round(sumdbRDA$cont$importance[2,1]*100, 2), "%)"),
       y = paste0("Unconstrained - Axis 1 (", round(sumdbRDA$cont$importance[2,2]*100, 2), "%)"),
       color = "Group")
dbRDA_species_scores

# ggsave("./plots/Beta-div_dbRDA_Habitat_phylum_species.tiff", dbRDA_species_scores, dpi = 600, width = 20, height = 15, units = "cm")