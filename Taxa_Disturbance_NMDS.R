# Load necessary libraries
library(dplyr)
library(vegan)
library(ggplot2)
library(tidyr)
library(ggrepel)

# Step 1: Load and preprocess data
community <- read.csv("Typology1_Macro_Taxa_Count_Complete.csv")
df <- community[, c("Genus", "Site_ID", "Count")]

# Create a matrix of Genus counts per Site_ID
df_matrix <- df %>%
  group_by(Genus, Site_ID) %>%
  summarise(Count = sum(Count)) %>%
  spread(key = Genus, value = Count, fill = 0)

# Step 2: Compute NMDS
dist_matrix <- vegdist(df_matrix[, -1], method = "bray")
nmds <- metaMDS(dist_matrix, k=2, trymax=100, distance = "bray", wascores= TRUE)
scores <- as.data.frame(scores(nmds))
scores$Site_ID = df_matrix$Site_ID

# Step 3: Import and categorize disturbance data
disturbance_data <- read.csv("Typology_1-Stream_Site_IDisturbance_Impact.csv")
quantiles_PAS <- quantile(disturbance_data$Ws_PAS_.prop., c(0.33, 0.66))
quantiles_BMP <- quantile(disturbance_data$ECO_BMP_Intensity, c(0.33, 0.66))
disturbance_data$Pasture <- cut(disturbance_data$Ws_PAS_.prop., breaks = c(-Inf, quantiles_PAS, Inf), labels = c("low", "medium", "high"))
disturbance_data$BMP <- cut(disturbance_data$ECO_BMP_Intensity, breaks = c(-Inf, quantiles_BMP, Inf), labels = c("low", "medium", "high"))

# Merge disturbance data with scores based on Site_ID
scores <- merge(scores, disturbance_data[, c("Site_ID", "Pasture", "BMP")], by = "Site_ID")
scores$Disturbance_Intensity <- paste(scores$Pasture, scores$BMP, sep="-")

# Step 4: Conduct Envfit analysis and filter significant taxa
envfit_scores <- envfit(nmds, df_matrix[, -1], permutations = 999) 
taxa_scores <- as.data.frame(scores(envfit_scores, display = "vectors"))
taxa_scores <- cbind(taxa_scores, Taxon = rownames(taxa_scores), pval = envfit_scores$vectors$pvals, r2 = envfit_scores$vectors$r)
sig_taxa_scores <- subset(taxa_scores, pval <= 0.09)
sig_taxa_scores <- sig_taxa_scores[order(sig_taxa_scores$r2, decreasing=TRUE), ]

# Step 5: Calculate convex hulls for each Disturbance_Intensity group
compute_hull <- function(df) {
  hpts <- chull(df$NMDS1, df$NMDS2)
  df[hpts, ]
}
# For Disturbance_Intensity
hulls_DI <- scores %>% 
  group_by(Disturbance_Intensity) %>% 
  do(compute_hull(.))

# For BMP
hulls_BMP <- scores %>% 
  group_by(BMP) %>% 
  do(compute_hull(.))

# For Pasture
hulls_PAS <- scores %>% 
  group_by(Pasture) %>% 
  do(compute_hull(.))

# Step 6: Visualize NMDS plots with Disturbance_Intensity and significant taxa information
# NMDS plot for Disturbance_Intensity
NMDS_plot_Dist <- ggplot() + 
  geom_polygon(data = hulls_DI, aes(x = NMDS1, y = NMDS2, fill = Disturbance_Intensity, group = Disturbance_Intensity), alpha = 0.3) +
  geom_text(data = scores, aes(x = NMDS1, y = NMDS2, label = Site_ID), size = 2.5, color = "black") +
  geom_text_repel(data = sig_taxa_scores, aes(x = NMDS1, y = NMDS2, label = Taxon), size = 3, box.padding = 0.5, point.padding = 0.5, segment.color = "grey50") + 
  scale_fill_manual(values = color_palette) +   
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 12), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2")
# NMDS plot for BMP
# colorblind-friendly BMP color palette & shape palette
BMP_palette <- c("high" = "#D78F02", "medium" = "#3B9E77", "low" = "#7540B3")

NMDS_plot_BMP <-  ggplot() + 
  geom_polygon(data = hulls_BMP, aes(x = NMDS1, y = NMDS2, fill = BMP, group = BMP), alpha = 0.4) +
  geom_point(data = scores, aes(x = NMDS1, y = NMDS2, shape = BMP), size = 1.5) + 
  geom_point(data = sig_taxa_scores, aes(x = NMDS1, y = NMDS2), shape = 1, size = 1.5) +
 geom_text_repel(data = scores, aes(x = NMDS1, y = NMDS2, label = Site_ID), size = 2.7, color = "black", box.padding = 0.3, point.padding = 0.55, segment.color = "grey50") +
geom_text_repel(data = sig_taxa_scores, aes(x = NMDS1, y = NMDS2, label = Taxon), 
                  size = 3, box.padding = 0.12, point.padding = 0.2, 
                  segment.color = "grey50", max.overlaps = 9) +  
  scale_fill_manual(values = BMP_palette) +     
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 12), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2")

# NMDS plot for Pasture
# colorblind-friendly, distinct color palette for PAS categories
PAS_palette <- c("high" = "#D95F02", "medium" = "#1B9E77", "low" = "#7570B3")
NMDS_plot_PAS <- ggplot() + 
  geom_polygon(data = hulls_PAS, aes(x = NMDS1, y = NMDS2, fill = Pasture, group = Pasture), alpha = 0.2) +
  geom_text(data = scores, aes(x = NMDS1, y = NMDS2, label = Site_ID), size = 3, color = "black") +
  geom_text_repel(data = sig_taxa_scores, aes(x = NMDS1, y = NMDS2, label = Taxon), size = 3, box.padding = 0.5, point.padding = 0.5, segment.color = "grey50") + 
  scale_fill_manual(values = PAS_palette) +   # Updated palette name
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 12), 
        axis.title.x = element_text(face = "bold", size = 12, colour = "black"), 
        legend.title = element_text(size = 12, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, linewidth = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2")
