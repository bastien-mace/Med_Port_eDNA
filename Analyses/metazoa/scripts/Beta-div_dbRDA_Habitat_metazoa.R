# Load the required packages
library(tidyverse)
library(vegan)
library(ggalt)
library(ggrepel)


# Load the data
## Load the Jaccard presence/absence matrix
jaccard <- read.csv2("../data/metazoa_jaccard.csv")

## Load metadata
metadata <- read.csv2("../../00_metadata/metadata.csv") %>%
  filter(Season == "Summer") # as reserves were sampled in summer only, they are
                             # compared with marinas samples collected in summer 
## Combine data
biodiv_table <- inner_join(jaccard, metadata) 


# Define dbRDA variables
## Select the Ys
data_dbRDA <- biodiv_table %>%
  select(-c("Sample", "Site", "Habitat", "Season", "Longitude", "Latitude"))

## Select the Xs
meta_dbRDA <- biodiv_table %>%
  select(c("Sample", "Site", "Habitat", "Season", "Longitude", "Latitude"))

meta_dbRDA$Longitude <- as.numeric(meta_dbRDA$Longitude)
meta_dbRDA$Latitude <- as.numeric(meta_dbRDA$Latitude)

## Remove absent MOTUs (those only present in marinas in autumn)
data_dbRDA <- data_dbRDA[, colSums(data_dbRDA) != 0]


# Set the model
dbRDA <- capscale(data_dbRDA ~ Habitat + Condition(Longitude + Latitude),
                  meta_dbRDA, distance ='jaccard')

# Get the inertia statistics                   
dbRDA

# Test the significance of the model
anova.cca(dbRDA, permutations = 9999)

# Calculate the adjusted-R²
RsquareAdj(dbRDA)


# Extract dbRDA data
sumdbRDA <- summary(dbRDA)

# Get the scores
scrs <- scores(dbRDA, display = c("species", "sites", "bp"), scaling = 2)

## Extract site names
site_names <- biodiv_table %>%
  arrange(match(Site, rownames(scrs$sites))) %>%
  pull(Site)

## Extract habitat names
habitat_names <- biodiv_table %>%
  arrange(match(Site, rownames(scrs$sites))) %>%
  pull(Habitat)

## Extract site scores and add site and habitat names
sites_centroids <- data.frame(scrs$sites)
sites_centroids$SITE <- site_names
sites_centroids$HABITAT <- habitat_names

## Remove the suffix
sites_centroids$SITE <- gsub("*.Summer", "", sites_centroids$SITE)


# Plot the data
dbRDA_site_scores <- ggplot(data = sites_centroids,
                            aes(x = CAP1, y = MDS1, color = HABITAT, label = SITE)) +
  geom_point(size = 3) +
  geom_text_repel(color = "black") +
  geom_encircle(aes(group = HABITAT, linetype = HABITAT, fill = HABITAT), s_shape = 1, expand = 0, alpha = 0.5, show.legend = FALSE) +
  scale_color_brewer(palette = "Accent", direction = -1) +
  scale_fill_brewer(palette = "Accent", direction = -1)  +
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
        legend.position = c(0.1, 0.9),
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "white", color = "black")) +
  labs(x = paste0("Constrained (", round(sumdbRDA$cont$importance[2,1]*100, 2), "%)"),
       y = paste0("Unconstrained - Axis 1 (", round(sumdbRDA$cont$importance[2,2]*100, 2), "%)"),
       color = "Habitat")
dbRDA_site_scores

# ggsave("../plots/Beta-div_dbRDA_Habitat_metazoa.tiff", dbRDA_site_scores, dpi = 600, width = 20, height = 15, units = "cm")
