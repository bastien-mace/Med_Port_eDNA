# Load the required packages
library(tidyverse)
library(vegan)
library(ade4)
library(ggalt)


# Load the data
## Load the Jaccard presence/absence matrix
jaccard <- read.csv2("../data/all_jaccard_pooled.csv")

## Load metadata
meta <- read.csv2("../../00_metadata/metadata_pooled.csv") %>%
  filter(Season == "Summer") # as reserves were sampled in summer only, they are
                             # compared with marinas samples collected in summer 
## Combine data
biodiv_table <- inner_join(jaccard, meta)


# Define dbRDA variables
## Select the Ys
data_motu <- biodiv_table %>%
  select(-c("Site", "Habitat", "Season", "Longitude", "Latitude"))

## Select the Xs
meta_motu <- biodiv_table %>%
  select(c("Site", "Habitat", "Season", "Longitude", "Latitude"))

meta_motu$Longitude <- as.numeric(meta_motu$Longitude)
meta_motu$Latitude <- as.numeric(meta_motu$Latitude)

## Remove absent MOTUs (those only present in marinas in autumn)
data_motu <- data_motu[, colSums(data_motu) != 0]


# Compute Jaccard distance between samples
motu.bray <- vegdist(data_motu, method = "jaccard")

# Compute PCoA
pcoa_motu <- dudi.pco(motu.bray, scannf = FALSE)


# Extract coordinates
pcoa_motu_2_plot <- pcoa_motu$li
pcoa_motu_2_plot$Site <- meta_motu$Site
pcoa_motu_2_plot$Habitat <- meta_motu$Habitat

# Remove the suffix
pcoa_motu_2_plot[pcoa_motu_2_plot == "Marina_Summer"] <- "Marina"
pcoa_motu_2_plot$Site <- gsub("*_Summer", "", pcoa_motu_2_plot$Site)


# Plot the data
PCoA_plot <- ggplot(pcoa_motu_2_plot,
                    aes(x = A1, y = A2, color = Habitat, label = Site))+
  geom_point(size = 3)+
  geom_text_repel(color = "black") +
  scale_color_brewer(palette = "Accent", direction = -1) +
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
  labs(x = paste0("PCoA 1 (", round(pcoa_motu$eig[1]/sum(pcoa_motu$eig)*100, 2), "%)"),
       y = paste0("PCoA 2 (", round(pcoa_motu$eig[2]/sum(pcoa_motu$eig)*100, 2), "%)"),
       color = "Habitat")
PCoA_plot

# ggsave("../plots/Beta-div_PCoA_Habitat.tiff", PCoA_plot, dpi = 600, width = 20, height = 15, units = "cm")
