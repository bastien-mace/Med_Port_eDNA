# Load the required packages
library(tidyverse)
library(spaMM)
library(rredlist)
library(ggsignif)

# Set the 'Negate' function
"%ni%" <- Negate("%in%")

# Set the IUCN API key
Sys.setenv(IUCN_KEY = "9bb4facb6d23f48efbf424bb05c0c1ef1cf6f468393bc745d42179ac4aca5fee")
API_key <- Sys.getenv("IUCN_KEY")


# Load scientific names of species identified with 100% identity
motu_taxo <- read.csv2("../data/teleo_motu_taxo.csv") %>%
  filter(best_identity_database >= 0.99) %>%
  select(c("species_name")) %>%
  na.omit() %>%
  unique()
motu_taxo


# Check IUCN status of these species in the Mediterranean Red List
species_red_list_med <- data.frame(species_name = motu_taxo$species_name) %>% 
  mutate(iucn_query = map(species_name, rl_search, region = "mediterranean", key = API_key)) 

# Extract their categories
species_categories <- species_red_list_med %>% 
  mutate(category = map_chr(iucn_query, pluck, "result", "category", .default = NA_character_)) %>% 
  select(species_name, category) %>% 
  mutate_at("category", ~ replace_na(.,"NE")) %>%
  cbind(region = "med")


# Check global IUCN status for species unassigned in the Mediterranean Red List
species_unassigned_med <- species_categories %>%
  filter(category %in% c("NE", "DD", "NA"))

species_red_list_global <- data.frame(species_name = species_unassigned_med$species_name) %>% 
  mutate(iucn_query = map(species_name, rl_search, region = "global", key = API_key))

# Extract their categories
species_categories_global <- species_red_list_global %>% 
  mutate(category = map_chr(iucn_query, pluck, "result", "category", .default = NA_character_)) %>% 
  select(species_name, category) %>% 
  mutate_at("category", ~ replace_na(.,"NE")) %>%
  cbind(region = "global")

# Replaced unassigned status in the IUCN Mediterranean Red List by their status
# in the global IUCN Red List
for (i in species_categories_global$species_name){
  species_categories[species_categories$species_name == i, ] <- species_categories_global[species_categories_global$species_name == i, ]
}

# Extract species for which NCBI scientific name is not recognized by IUCN database
species_to_correct <- species_categories %>%
  filter(category %in% c("NE"))
species_to_correct

# NCBI scientific name         # IUCN scientific name - Status - Database
#-------------------------------------------------------------------------
# Epinephelus marginatus      -> Mycteroperca marginatus - VU - global
species_categories[species_categories$species_name == "Epinephelus marginatus", 2:3] <- c("VU", "global")
# Lepadogaster candolii       -> Lepadogaster candollei - LC - med
species_categories[species_categories$species_name == "Lepadogaster candolii", 2:3] <- c("LC", "med")
# Aphia minuta                -> Aphia minuta - NE - global
species_categories[species_categories$species_name == "Aphia minuta", 2:3] <- c("NE", "global")
# Pegusa nasuta               -> Pegusa lascaris - LC - global
species_categories[species_categories$species_name == "Pegusa nasuta", 2:3] <- c("LC", "global")
# Zosterisessor ophiocephalus -> Gobius ophiocephalus - LC - med
species_categories[species_categories$species_name == "Zosterisessor ophiocephalus", 2:3] <- c("LC", "med")
# Cololabis saira             -> Cololabis saira - NE - global
species_categories[species_categories$species_name == "Cololabis saira", 2:3] <- c("NE", "global")
#-------------------------------------------------------------------------

species_redlist <- species_categories %>%
  filter(category %in% c("VU", "EN", "CR"))
species_redlist


# Extract MOTUs corresponding to these species
motus_redlist <- read.csv2("../data/teleo_motu_taxo.csv") %>%
  filter(scientific_name %in% species_redlist$species_name) %>%
  select(c("definition", "scientific_name"))


# Import metadata
metadata <- read.csv2("../../00_metadata/metadata.csv")


# Restrict the presence/absence matrix to the Red List species
redlist_sp <- read.csv2("../data/teleo_jaccard.csv") %>%
  select("Sample", motus_redlist$definition) %>%
  left_join(., metadata) %>%
  filter(Season == "Summer") %>%
  select(-c("Season", "Longitude", "Latitude"))

redlist_sp_jaccard <- select(redlist_sp, -c("Sample", "Site", "Habitat"))

names(redlist_sp_jaccard) <- motus_redlist$scientific_name


# Count the number of Red List species found in each sample
redlist_sp_by_sample <- c()
for (i in 1:dim(redlist_sp_jaccard)[1]){
  nb_in_sample <- sum(redlist_sp_jaccard[i, ])
  redlist_sp_by_sample <- c(redlist_sp_by_sample, nb_in_sample)
}


# Create a data frame with threatened species richness by sample, with their
# corresponding site and habitat
redlist_richness <- data.frame(Sample = redlist_sp$Sample,
                                Richness = redlist_sp_by_sample,
                                Habitat = redlist_sp$Habitat,
                                Site = redlist_sp$Site) %>%
  inner_join(., metadata)

redlist_richness$Longitude <- as.numeric(redlist_richness$Longitude)
redlist_richness$Latitude <- as.numeric(redlist_richness$Latitude)


# Data exploration
## Summary: Descriptive statistics
summary(redlist_richness$Richness)

## Distribution of MOTU richness
par(mfrow = c(2, 2))
### Boxplot
boxplot(redlist_richness$Richness, col = 'blue', ylab = 'Richness')
### Cleveland plot
dotchart(redlist_richness$Richness, pch = 16, col = 'blue', xlab = 'Richness')
### Histogram
hist(redlist_richness$Richness, col = 'blue', xlab = "Richness", main = "")
### Quantile-Quantile plot
qqnorm(redlist_richness$Richness, pch = 16, col = 'blue', xlab = '')
qqline(redlist_richness$Richness, col = 'red')


# Build the LMM model (ML method)
model_ML <- fitme(Richness ~ Habitat + Matern(1|Longitude + Latitude),
                  data = redlist_richness,
                  method = "ML")

# Build the LMM model (REML method)
model_REML <- fitme(Richness ~ Habitat + Matern(1|Longitude + Latitude),
                    data = redlist_richness,
                    method = "REML")

# Test
fixedLRT(Richness ~ 1 + Matern(1|Longitude + Latitude),
         Richness ~ Habitat + Matern(1|Longitude + Latitude),
         data = redlist_richness,
         method = "ML")

# Obtain the pseudo-R²
pseudoR2(model_ML, nullform = Richness ~ 1 + Matern(1|Longitude + Latitude))


# Check normality of residuals
par(mfrow = c(1, 2))
## Histogram
hist(residuals(model_REML), col = 'blue', xlab = "Residuals", main = "Check Normality")
## Quantile-Quantile plot
qqnorm(residuals(model_REML), pch = 16, col = 'blue', xlab = '')
qqline(residuals(model_REML), col = 'red')

## Shapiro test
shapiro.test(residuals(model_REML))


# Check homogeneity of variances
par(mfrow = c(1, 2))
# Residuals vs Fitted
plot(residuals(model_REML) ~ fitted(model_REML), col = 'blue', pch = 16,
     ylab = "Residuals",
     xlab = "Fitted",)
abline(h = 0)
# Residuals vs Habitat
boxplot(residuals(model_REML) ~ redlist_richness$Habitat, 
        varwidth = TRUE,
        ylab = "Residuals",
        xlab = "Habitat")
abline(h = 0)


# Extract fitted values
predict_values <- predict(model_REML)
redlist_richness$predict <- predict_values[, 1]

# Rename the categories with the number of samples
redlist_richness[redlist_richness == "Marina_Summer"] <- "Marina\n(n = 12)"
redlist_richness[redlist_richness == "Reserve"] <- "Reserve\n(n = 8)"

# Plot the data
alpha_div_predicted <- ggplot(redlist_richness, aes(x = Habitat, y = predict, fill = Habitat)) +
  geom_violin(show.legend = F) +
  stat_summary(fun = mean, geom = "point", color = "darkred", size = 3,
               show.legend = F) +
  stat_summary(fun = mean, geom = "text", color = "darkred", vjust = -1.5,
               aes(label = round(..y.., digits = 2)),
               show.legend = F) +
  geom_signif(comparisons = list(c("Marina\n(n = 12)", "Reserve\n(n = 8)")),
              annotations = c("*"), textsize = 5, fontface = "bold") +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  ylab("Threatened species richness (LMM predicted values)") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_blank()) 
alpha_div_predicted

# ggsave("../plots/RedList_teleo_predicted.tiff", alpha_div_predicted, dpi = 600, width = 15, height = 15, units = "cm")

## Plot raw data
alpha_div <- ggplot(redlist_richness, aes(x = Habitat, y = Richness, fill = Habitat)) +
  geom_violin(show.legend = F) +
  stat_summary(fun = mean, geom = "point", color = "darkred", size = 3,
               show.legend = F) +
  stat_summary(fun = mean, geom = "text", color = "darkred", vjust = -1,
               aes(label = round(..y.., digits = 2)),
               show.legend = F) +
  scale_fill_brewer(palette = "Accent", direction = -1) +
  ylab("Threatened species richness (LMM predicted values)") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_blank()) 
alpha_div

# ggsave("../plots/RedList_teleo.tiff", alpha_div, dpi = 600, width = 15, height = 15, units = "cm")
