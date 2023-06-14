# Load the required package
library(rstatix)

# Taxonomic level: MOTUs
## Extract the data already used for the dbRDA
source("Beta-div_dbRDA_Habitat_MOTU.R")
row.names(data_dbRDA_motu) <- meta_dbRDA_motu$Site

## Extract MOTU from marinas
marinas <- meta_dbRDA_motu %>%
  filter(Habitat == "Marina_Summer") %>%
  select(Site)

marina_motu <- subset(data_dbRDA_motu, rownames(data_dbRDA_motu) %in% marinas[, 1])

## Extract MOTU from reserves
reserves <- meta_dbRDA_motu %>%
  filter(Habitat == "Reserve") %>%
  select(Site)

reserve_motu <- subset(data_dbRDA_motu, rownames(data_dbRDA_motu) %in% reserves[, 1])

## Compute the Jaccard distances between samples
dist_marina_natural_motu <- vegdist(rbind(marina_motu, reserve_motu), method = "jaccard")
dist_marina_natural_motu

## Only keep Marina-Reserve distances
list_beta_motu <- as.numeric(dist_marina_natural_motu[c(6:9, 14:17, 21:24, 27:30, 32:35, 36:39)])


# Taxonomic level: Families
## Extract the data already used for the dbRDA
source("Beta-div_dbRDA_Habitat_family.R")
row.names(data_dbRDA_family) <- meta_dbRDA_family$Site

## Extract MOTU from marinas
marina_family <- subset(data_dbRDA_family, rownames(data_dbRDA_family) %in% marinas[, 1])

## Extract MOTU from reserves
reserve_family <- subset(data_dbRDA_family, rownames(data_dbRDA_family) %in% reserves[, 1])

## Compute the Bray-Curtis distances between samples
dist_marina_natural_family <- vegdist(rbind(marina_family, reserve_family), method = "bray")
dist_marina_natural_family

## Only keep Marina-Reserve distances
list_beta_family <- as.numeric(dist_marina_natural_family[c(6:9, 14:17, 21:24, 27:30, 32:35, 36:39)])


# Taxonomic level: Orders
## Extract the data already used for the dbRDA
source("Beta-div_dbRDA_Habitat_order.R")
row.names(data_dbRDA_order) <- meta_dbRDA_order$Site

## Extract MOTU from marinas
marina_order <- subset(data_dbRDA_order, rownames(data_dbRDA_order) %in% marinas[, 1])

## Extract MOTU from reserves
reserve_order <- subset(data_dbRDA_order, rownames(data_dbRDA_order) %in% reserves[, 1])

## Compute the Bray-Curtis distances between samples
dist_marina_natural_order <- vegdist(rbind(marina_order, reserve_order), method = "bray")
dist_marina_natural_order

## Only keep Marina-Reserve distances
list_beta_order <- as.numeric(dist_marina_natural_order[c(6:9, 14:17, 21:24, 27:30, 32:35, 36:39)])


# Taxonomic level: Classes
## Extract the data already used for the dbRDA
source("Beta-div_dbRDA_Habitat_class.R")
row.names(data_dbRDA_class) <- meta_dbRDA_class$Site

## Extract MOTU from marinas
marina_class <- subset(data_dbRDA_class, rownames(data_dbRDA_class) %in% marinas[, 1])

## Extract MOTU from reserves
reserve_class <- subset(data_dbRDA_class, rownames(data_dbRDA_class) %in% reserves[, 1])

## Compute the Bray-Curtis distances between samples
dist_marina_natural_class <- vegdist(rbind(marina_class, reserve_class), method = "bray")
dist_marina_natural_class

## Only keep Marina-Reserve distances
list_beta_class <- as.numeric(dist_marina_natural_class[c(6:9, 14:17, 21:24, 27:30, 32:35, 36:39)])


# Taxonomic level: Phyla
## Extract the data already used for the dbRDA
source("Beta-div_dbRDA_Habitat_phylum.R")
row.names(data_dbRDA_phylum) <- meta_dbRDA_phylum$Site

## Extract MOTU from marinas
marina_phylum <- subset(data_dbRDA_phylum, rownames(data_dbRDA_phylum) %in% marinas[, 1])

## Extract MOTU from reserves
reserve_phylum <- subset(data_dbRDA_phylum, rownames(data_dbRDA_phylum) %in% reserves[, 1])

## Compute the Bray-Curtis distances between samples
dist_marina_natural_phylum <- vegdist(rbind(marina_phylum, reserve_phylum), method = "bray")
dist_marina_natural_phylum

## Only keep Marina-Reserve distances
list_beta_phylum <- as.numeric(dist_marina_natural_phylum[c(6:9, 14:17, 21:24, 27:30, 32:35, 36:39)])


# Create a data frame with the dissimilarity values for each taxonomic level
matrix_taxo <- melt(data.frame(MOTU = list_beta_motu,
                               Family = list_beta_family,
                               Order = list_beta_order,
                               Class = list_beta_class,
                               Phylum = list_beta_phylum))


# Test if there are differences in pairwise dissimilarity between taxonomic levels
## Kruskal-Wallis test: Overall difference
kruskal.test(value ~ variable, data = matrix_taxo)

## Dunn test: Differences between each combinations
dunn_test(value ~ variable, data = matrix_taxo)


# Plot the data
## Add the number of MOTU/Family/Order/Class/Phylum
matrix_taxo$variable <- gsub("MOTU", "MOTU\n(1213)", matrix_taxo$variable)
matrix_taxo$variable <- gsub("Family", "Family\n(178)", matrix_taxo$variable)
matrix_taxo$variable <- gsub("Order", "Order\n(113)", matrix_taxo$variable)
matrix_taxo$variable <- gsub("Class", "Class\n(54)", matrix_taxo$variable)
matrix_taxo$variable <- gsub("Phylum", "Phylum\n(31)", matrix_taxo$variable)

matrix_taxo$variable <- factor(matrix_taxo$variable,                                   
                                levels = c("MOTU\n(1213)",
                                           "Family\n(178)",
                                           "Order\n(113)",
                                           "Class\n(54)",
                                           "Phylum\n(31)"))

beta_div_taxo <- ggplot(matrix_taxo,aes(x = variable, y = value, fill = variable)) +
  geom_violin(show.legend = F) +
  stat_summary(fun = mean, geom = "point", color = "black", size = 3,
               show.legend = F) +
  stat_summary(fun = mean, geom = "text", color = "black", vjust = 2,
              aes(label = round(..y.., digits = 2)),
              show.legend = F) +
   scale_fill_brewer(palette = "YlGn", direction = -1) +
  ylab("Dissimilarity index") +
  xlab("Taxonomic level") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_text(colour = "black", size = 15, face = "bold"))   
beta_div_taxo

# ggsave("../plots/Beta-div_taxonomic_levels.tiff", beta_div_taxo, dpi = 600, width = 25, height = 17.5, units = "cm")

