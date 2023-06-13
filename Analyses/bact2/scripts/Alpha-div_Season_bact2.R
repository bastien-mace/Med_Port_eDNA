# Load the required packages
library(tidyverse)
library(spaMM)
library(ggsignif)


# Import metadata
metadata <- read.csv2("../../00_metadata/metadata.csv") %>%
  filter(Habitat %in% c("Marina_Summer", "Marina_Autumn"))


# Import MOTU richness file
richness_table <- read.csv2("../data/bact2_richness.csv") %>%
  inner_join(., metadata) # all data are merged

richness_table$Longitude <- as.numeric(richness_table$Longitude)
richness_table$Latitude <- as.numeric(richness_table$Latitude)


# Data exploration
## Summary: Descriptive statistics
summary(richness_table$Richness)

## Distribution of MOTU richness
par(mfrow = c(2, 2))
### Boxplot
boxplot(richness_table$Richness, col = 'blue', ylab = 'Richness')
### Cleveland plot
dotchart(richness_table$Richness, pch = 16, col = 'blue', xlab = 'Richness')
### Histogram
hist(richness_table$Richness, col = 'blue', xlab = "Richness", main = "")
### Quantile-Quantile plot
qqnorm(richness_table$Richness, pch = 16, col = 'blue', xlab = '')
qqline(richness_table$Richness, col = 'red')


# Build the LMM model (ML method)
model_ML <- fitme(Richness ~ Season + Matern(1|Longitude + Latitude),
                  data = richness_table,
                  method = "ML")

# Build the LMM model (REML method)
model_REML <- fitme(Richness ~ Season + Matern(1|Longitude + Latitude),
                    data = richness_table,
                    method = "REML")

# Test
fixedLRT(Richness ~ 1 + Matern(1|Longitude + Latitude),
         Richness ~ Season + Matern(1|Longitude + Latitude),
         data = richness_table,
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
# Residuals vs Season
boxplot(residuals(model_REML) ~ richness_table$Season, 
        varwidth = TRUE,
        ylab = "Residuals",
        xlab = "Season")
abline(h = 0)


# Extract fitted values
predict_values <- predict(model_REML)
richness_table$predict <- predict_values[, 1]

# Rename the categories with the number of samples
richness_table[richness_table == "Marina_Summer"] <- "Summer\n(n = 12)"
richness_table[richness_table == "Marina_Autumn"] <- "Autumn\n(n = 12)"

# Plot the data
alpha_div_predicted <- ggplot(richness_table, aes(x = Habitat, y = predict, fill = Habitat)) +
  geom_violin(show.legend = F) +
  stat_summary(fun = mean, geom = "point", color = "darkred", size = 3,
               show.legend = F) +
  stat_summary(fun = mean, geom = "text", color = "darkred", vjust = -1,
               aes(label = round(..y.., digits = 2)),
               show.legend = F) +
  geom_signif(comparisons = list(c("Autumn\n(n = 12)", "Summer\n(n = 12)")),
              annotations = c("***"), textsize = 5, fontface = "bold") +
  scale_fill_manual(values = c("#FDC086", "#BEAED4")) +
  ylab("MOTU richness (LMM predicted values)") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_blank()) 
alpha_div_predicted

# ggsave("../plots/Alpha-div_Season_bact2_predicted.tiff", alpha_div_predicted, dpi = 600, width = 15, height = 15, units = "cm")

## Plot raw data
alpha_div <- ggplot(richness_table, aes(x = Habitat, y = Richness, fill = Habitat)) +
  geom_violin(show.legend = F) +
  stat_summary(fun = mean, geom = "point", color = "darkred", size = 3,
               show.legend = F) +
  stat_summary(fun = mean, geom = "text", color = "darkred", vjust = -1,
               aes(label = round(..y.., digits = 2)),
               show.legend = F) +
  scale_fill_manual(values = c("#FDC086", "#BEAED4")) +
  ylab("MOTU richness") +
  theme_classic() +
  theme(axis.text.y = element_text(colour = "black", size = 10),
        axis.text.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15, face = "bold"),
        axis.title.x = element_blank()) 
alpha_div

# ggsave("../plots/Alpha-div_Season_bact2.tiff", alpha_div, dpi = 600, width = 15, height = 15, units = "cm")
