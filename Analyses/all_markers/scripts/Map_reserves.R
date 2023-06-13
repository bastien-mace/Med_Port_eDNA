# Load the required packages
library(tidyverse)
library(reshape2)
library(magick)
library(mapdata)


# Import MOTU richness for each marker
teleo_richness <- read.csv2("../data/teleo_richness_reserve.csv") %>%
  as.data.frame() %>%
  rename("teleo" = "Richness")

metazoa_richness <- read.csv2("../data/metazoa_richness_reserve.csv") %>%
  as.data.frame() %>%
  rename("metazoa" = "Richness")

euka2_richness <- read.csv2("../data/euka2_richness_reserve.csv") %>%
  as.data.frame() %>%
  rename("euka2" = "Richness")

bact2_richness <- read.csv2("../data/bact2_richness_reserve.csv") %>%
  as.data.frame() %>%
  rename("bact2" = "Richness")

# Merge all data
table_richness <- teleo_richness %>%
  inner_join(., metazoa_richness) %>%
  inner_join(., euka2_richness) %>%
  inner_join(., bact2_richness) %>%
  melt() %>%
  dplyr::rename("Marker" = "variable")


# Download the map of the Mediterranean Sea
med <- map_data("worldHires",  xlim = c(3, 10), ylim = c(42, 44))
## Create the map and save it
fig <- image_graph(width = 2400, height = 2200, res = 300)
ggplot() +
  geom_polygon(data = med, aes(x = long, y = lat, group = group), fill = "slategrey") +
  coord_fixed(xlim = c(3, 9.5), ylim = c(42, 44), ratio = 1.2)+
  theme(panel.background = element_blank(),
        axis.title = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
dev.off()


# Create a list of barplots
sites <- unique(table_richness$Site)
bp <- list()
for (i in 1:length(sites)) {
  site_to_plot <- table_richness[which(table_richness$Site %in% sites[i]), ]
  
  bp[[i]] <- ggplot(data = site_to_plot, aes(x = Marker, y = value, fill = interaction(Marker, Season))) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("#CBD5E8", "#FDCDAC", "#B3E2CD", "#F4CAE4", "#8DA0CB", "#FC8D62", "#66C2A5", "#E78AC3")) +
    labs(title = sites[i]) +
    ylim(0, 500) +
    theme(panel.background = element_rect(fill = "transparent"),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text.y = element_text(colour="black", size = 5),
          axis.text.x = element_text(colour="black", size = 5, angle = 20),
          axis.title.y = element_blank(),
          axis.title.x = element_blank(),
          plot.title = element_text(colour = "black", size = 10, face = "bold"))
  
}
names(bp) <- sites

# Plot them
barfig1 <- image_graph(width = 350, height = 350, res = 300, bg = 'transparent')
bp[[1]]
dev.off()
barfig2 <- image_graph(width = 350, height = 350, res = 300, bg = 'transparent')
bp[[2]]
dev.off()
barfig3 <- image_graph(width = 350, height = 350, res = 300, bg = 'transparent')
bp[[3]]
dev.off()
barfig4 <- image_graph(width = 350, height = 350, res = 300, bg = 'transparent')
bp[[4]]
dev.off()

# Add them to the map
final <- image_composite(fig, barfig2, offset = "+200+1125")
final <- image_composite(final, barfig1, offset = "+750+850")
final <- image_composite(final, barfig4, offset = "+1100+925")
final <- image_composite(final, barfig3, offset = "+1400+625")
final

# ggsave("../plots/Alpha-div_Map_reserve.tiff", image_ggplot(final), width = 2400, height = 2200, unit = "px", dpi = 600)
