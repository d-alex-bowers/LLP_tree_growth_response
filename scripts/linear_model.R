########################################
### Tree growth response to competition and health condition in a 
### subtropical savanna
###
### Multiple Linear Regression Models
###
### D. Alex Bowers
### Daniel J. Johnson
########################################


# required packages
library(vegan)
library(car)
library(sf)
library(spdep)
library(ggplot2)

# needed data
pipa_lm <- read.csv("data/pipa_lm.csv") # Pinus palustris
qula_lm <- read.csv("data/qula_lm.csv") # Quercus laevis



##### linear regression #####

# function to clean data for model including outlier tests
clean_health_comp_model <- function(pcoa_coords_health) {
  # Select variables needed for model
  lm_sp <- dplyr::select(pcoa_coords_health, StemTag, bai, PC1, PC2, PC3, 
                         idw_con, idw_het, ba_m_2019, crown_living, water)
  
  # Compute Z-scores
  lm_sp$z_scores <- scale(lm_sp$bai)
  
  # Define threshold
  lm_sp <- lm_sp[abs(lm_sp$z_scores) < 5, ]
  
  # remove NAs
  lm_sp <- na.omit(lm_sp)
  
  return(lm_sp)
}


# clean data 
pipa_clean <- clean_health_comp_model(pipa_lm)
qula_clean <- clean_health_comp_model(qula_lm)

# model to run multiple linear regression
health_comp_model <- function(pcoa_coords_health, model_ac = FALSE, ac) {
  if (model_ac) {
    # remove inf from neighbors
    pcoa_coords_health$ac_term <- ac 
    pcoa_coords_health <- pcoa_coords_health[!is.infinite(pcoa_coords_health$ac_term), ]
    # run model
    model <- lm(bai ~ scale(PC1) + scale(PC2) + scale(idw_con) +
                  scale(idw_het) + scale(crown_living) +
                  scale(ba_m_2019) + scale(water) + scale(ac_term),
                data = pcoa_coords_health) 
  } else {
    model <- lm(bai ~ scale(PC1) + scale(PC2) + scale(idw_con) +
                  scale(idw_het) + scale(crown_living) +
                  scale(ba_m_2019) + scale(water),
                data = pcoa_coords_health)
  }
  # Return model summary
  return(model)
}


# run model for both species
pipa_model <- health_comp_model(pipa_clean)
qula_model <- health_comp_model(qula_clean)
# view summary
summary(pipa_model)
summary(qula_model)

# check for multicolinerality
vif(pipa_model)
vif(qula_model)

# check leverage using cooks distance
pipa_clean$cooks <- cooks.distance(pipa_model)
qula_clean$cooks <- cooks.distance(qula_model)
# find poor values
pipa_cooks <- subset(pipa_clean, pipa_clean$cooks <= 4/nrow(pipa_clean))
qula_cooks <- subset(qula_clean, qula_clean$cooks <= 4/nrow(qula_clean))
# rerun model
pipa_cooks_model <- health_comp_model(pipa_cooks)
qula_cooks_model <- health_comp_model(qula_cooks)
# view summary
summary(pipa_cooks_model)
summary(qula_cooks_model)



##### plot the model #####

# Extract coefficients and related statistics
# pines
pipa_model_summary <- summary(pipa_cooks_model)
pipa_coef_df <- as.data.frame(pipa_model_summary$coefficients)
names(pipa_coef_df) <- c("Estimate", "Std.Error", "z.value", "p.value")
pipa_coef_df$parameter <- rownames(pipa_coef_df)
pipa_coef_df$sp <- "Pinus palustris"
# oaks
qula_model_summary <- summary(qula_cooks_model)
qula_coef_df <- as.data.frame(qula_model_summary$coefficients)
names(qula_coef_df) <- c("Estimate", "Std.Error", "z.value", "p.value")
qula_coef_df$parameter <- rownames(qula_coef_df)
qula_coef_df$sp <- "Quercus laevis"
# combine
coef_df <- rbind(pipa_coef_df, qula_coef_df)
rownames(coef_df) <- NULL

# Add a column for significance
alpha <- 0.05
coef_df$Significant <- ifelse(coef_df$p.value < alpha, 19, 1)

# remove Intercept
coef_df <- coef_df %>%
  filter(parameter != "(Intercept)")

# reorder columns
coef_df$parameter <- factor(coef_df$parameter, levels = 
                              c(
                                "scale(ac_term)",
                                "scale(water)",
                                "scale(ba_m_2019)",
                                "scale(idw_het)",
                                "scale(idw_con)",
                                "scale(crown_living)",
                                "scale(PC2)",
                                "scale(PC1)"))
# italizise speices names
coef_df$sp <- ifelse(
  coef_df$sp == "Pinus palustris",
  "italic(Pinus~palustris)",
  "italic(Quercus~laevis)"
)

# plot
ggplot(coef_df, aes(
  x = parameter,
  y = Estimate,
  ymin = Estimate - 2 * Std.Error,
  ymax = Estimate + 2 * Std.Error,
  shape = Significant
)) +
  geom_linerange() +
  geom_point(size = 4) +
  scale_shape_identity() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  coord_flip() +
  labs(
    y = "Parameter estimate",
    x = NULL
  ) +
  scale_x_discrete(labels = c(
    "scale(PC1)" = "Physical damage (PCoA1)",
    "scale(PC2)" = "Leaf loss (PCoA2)",
    "scale(idw_con)" = "Intraspecific competition",
    "scale(idw_het)" = "Interspecific competition",
    "scale(crown_living)" = "Percent crown living",
    "scale(ba_m_2019)" = "Initial basal area",
    "scale(water)" = "Depth to water table",
    "scale(ac_term)" = "Autocovariate"
  )) +
  facet_wrap(
    ~ sp,
    labeller = label_parsed
  ) + 
  theme_bw() + 
  theme(
    strip.background = element_rect(fill = NA, linewidth = 1),
    panel.border = element_rect(linewidth = 1),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 18),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )
