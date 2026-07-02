########################################
### Tree growth response to competition and health condition in a 
### subtropical savanna
###
### Principal Coordinates Analysis
###
### D. Alex Bowers
### Daniel J. Johnson
########################################

# required packages
library(vegan)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(cowplot)

# needed data
trees <- read.csv("data/trees_PCoA.csv")


##### Run PCoA #####

# function to run PCoA
pcoa <- function(data, species, method) {
  # select species
  pcoa_data <- subset(data, sp  == species)
  # select FADs
  predictors_health <- pcoa_data %>% dplyr::select(broken_stem, defoliation, rot, wound, 
                                                   fire, insect, canker, leaf_damage, 
                                                   hollow_stem, overtopped, animal, 
                                                   lightning, fungi, crushed,
                                                   bark_beetles, unknown,
                                                   postfire_leaf_loss)
  
  # ensure no columns have all 0s
  predictors_health <- predictors_health %>% select_if(~ any(. != 0))
  
  # gower dissimilarity index
  g_predictors_health <- vegdist(predictors_health, method)
  
  # run PCoA
  ### k is n col - 1
  cmd_health <- cmdscale(g_predictors_health, (k = 13), eig = TRUE)
  
  # PCoA table
  eigenvalues_health <- cmd_health$eig[1:13]
  propVar_health <- eigenvalues_health/sum(eigenvalues_health)
  cumVar_health <- cumsum(propVar_health)
  PCoA_Table_health <- cbind(eigenvalues_health, propVar_health, cumVar_health)
  
  # Scree plot:
  scree <- data.frame(Index = c(1,2,3,4,5,6,7,8,9,10,11,12,13), Eigen = eigenvalues_health)
  scree_plot <- ggplot(scree, aes(x = Index, y = Eigen)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(
      strip.background = element_rect(fill = NA, linewidth = 1),
      panel.border = element_rect(linewidth = 1),
    ) +
    labs(
      x = "Index",
      y = "Eigen Value"
    )
  
  
  # Convert PCoA coordinates to a data frame
  pcoa_coords_health <- as.data.frame(cmd_health$points)
  pcoa_coords_health <- pcoa_coords_health[1:3]
  colnames(pcoa_coords_health) <- c("PC1", "PC2", "PC3")
  pcoa_coords_health <- merge(pcoa_coords_health, pcoa_data, by = "row.names")
  
  
  # get predictor scores (weighted averages)
  predictors_scores_health <- wascores(cmd_health$points[, 1:3], predictors_health)
  predictors_scores_health <- data.frame(predictors_scores_health)
  
  
  # pca_coords_health is data frame with all values
  # PCoA_Table_health shows cumulative variance for deciding which axes to keep
  # preictors_scores_scores health shows weights of each variable
  # scree plot
  # 5 and 6 for veg dist 
  return(list(pcoa_coords_health, PCoA_Table_health, predictors_scores_health, scree_plot, cmd_health, predictors_health))
}


# run pcoa for pines and oaks
pipa_pcoa <- pcoa(trees, "PIPA", "gower")
qula_pcoa <- pcoa(trees, "QULA", "gower")

# correct oak PCoA for ease of understanding
qula_pcoa[[1]]$PC1 <- qula_pcoa[[1]]$PC1 * -1
qula_pcoa[[3]]$X1 <- qula_pcoa[[3]]$X1 * -1
qula_pcoa[[1]]$PC2 <- qula_pcoa[[1]]$PC2 * -1
qula_pcoa[[3]]$X2 <- qula_pcoa[[3]]$X2 * -1

# view PCoA
pipa_pcoa[[2]]
qula_pcoa[[2]]

# data for regression
pipa_lm <- pipa_pcoa[[1]]
qula_lm <- qula_pcoa[[1]]


##### Find significant factors #####

# function to find significant variables in PCoA
# check which loadings are significant 
sig_ordination <- function(df) {
  set.seed(483729)
  envfit <- envfit(df[[5]], df[[6]])
  vectors <- as.data.frame(scores(envfit, display = "vectors"))
  vectors$r2 <- envfit$vectors$r
  vectors$pval <- envfit$vectors$pvals
  # merge as a data frame with weigted averages
  vectors <- cbind(df[[3]], vectors)
  vectors <- vectors[order(-vectors$r2), ]
  return(vectors)
}

# check to see which loadings are signigificant by fitting an 
# environmental vector to ordination model using envfit
sig_ord_pipa <- sig_ordination(pipa_pcoa)
sig_ord_pipa
sig_ord_qula <- sig_ordination(qula_pcoa)
sig_ord_qula



##### Plot the PCoA

# function to plot pcoa
plot_pcoa <- function(sp_pcoa, highlight_vars = NULL) {
  
  # separate lists
  pcoa_coords_health <- sp_pcoa[[1]]
  predictors_scores_health <- sp_pcoa[[3]]
  colnames(predictors_scores_health) <- c("PC1", "PC2", "PC3")
  
  # convert to dataframe
  predictors_scores_health <- as.data.frame(predictors_scores_health)
  
  # classify based on highlight_vars instead of thresholds
  predictors_scores_health$importance <- ifelse(
    rownames(predictors_scores_health) %in% highlight_vars, 
    "high", "low"
  )
  
  # add variables for plot
  PCa <- pcoa_coords_health[, "PC1"]
  PCb <- pcoa_coords_health[, "PC2"]
  Xa <- predictors_scores_health[, "PC1"]
  Xb <- predictors_scores_health[, "PC2"]
  pcoa_labels <- sp_pcoa[[2]][1:3, ]
  rownames(pcoa_labels) <- c("PC1", "PC2", "PC3")
  
  # plot
  ggplot() + 
    geom_hline(yintercept = 0, linetype = "solid", color = "black") +
    geom_vline(xintercept = 0, linetype = "solid", color = "black") +
    
    # plot tree points
    geom_jitter(data = pcoa_coords_health, aes(x = PCa, y = PCb), size = 3, 
                alpha = .5, width = 0.005, height = 0.005) +
    
    # plot variable arrows and labels
    geom_text_repel(data = predictors_scores_health, 
                    aes(x = PC1, y = PC2, 
                        label = rownames(predictors_scores_health), 
                        color = importance),
                    size = 8, max.overlaps = 115, box.padding = 0.55, 
                    segment.color = "#545454", segment.size = 1, 
                    segment.curvature = 0.3, segment.ncp = 5, point.padding = 0.1) +
    
    geom_point(data = predictors_scores_health, 
               aes(x = PC1, y = PC2, color = importance),
               size = 3, shape = 16) +
    
    scale_color_manual(values = c("high" = "red", "low" = "grey")) +
    theme(legend.position = "none") +
    xlim(-0.3, 0.3) +
    ylim(-0.3, 0.3) +
    labs(
      x = paste0("PCoA 1", " (Explains ", round(pcoa_labels["PC1", 2] * 100, 1), "% of variation)"), 
      y = paste0("PCoA 2", " (Explains ", round(pcoa_labels["PC2", 2] * 100, 1), "% of variation)")
    )
}

# plot pines and oaks
# pine pc1 and 2
fig1a <- plot_pcoa(pipa_pcoa,
                   c("postfire_leaf_loss", "fire", "wound", "rot")) +
  ggtitle(expression(italic("Pinus palustris")))
# oak pc1 and 2
fig1b <- plot_pcoa(qula_pcoa,
                   c("defoliation", "rot", "wound", "fire", "broken_stem", 
                     "postfire_leaf_loss", "hollow_stem", "leaf_damage")) +
  ggtitle(expression(italic("Quercus laevis")))

# theme for manuscript
theme <- theme_bw() + 
  theme(
    strip.background = element_rect(fill = NA, linewidth = 1),
    panel.border = element_rect(linewidth = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )

# figure 1
plot_grid(
  fig1a + theme,
  fig1b + theme,
  ncol = 2
)
]