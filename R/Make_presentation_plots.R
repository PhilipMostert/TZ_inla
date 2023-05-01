rm(list = ls())
library(INLA)
library(INLAutils)
library(brinla)
library(tidyverse)
library(GGally)
library(raster)
library(readxl)
library(ggpubr)
library(patchwork)
library(ggthemes)
library(ztable)
library(magrittr)
library(terra)

# ------------------------------------------------------------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------------------------------------------------------------
load('')
source("source/misc_functions.R")
load('model_data/regression_data.RData')
load("results/traits_model_list.RData")

# ------------------------------------------------------------------------------------------------------------------------
# Raw data plot
# ------------------------------------------------------------------------------------------------------------------------
plot_80s <- ggplot() +
      geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group),
                   colour = "white", fill = NA) +
      geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 0 & as.data.frame(atlas_sp)$date_index == 1, ],
                 aes(x = x, y = y, col = "Absence"), alpha = 0.8, pch = 15, cex = 2) +
      geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 1 & as.data.frame(atlas_sp)$date_index == 1, ],
                 aes(x = x, y = y, col = "Presence"), alpha = 0.8, pch = 15, cex = 3) +
      geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 0 & as.data.frame(ebird_sp)$date_index == 1, ],
                 aes(x = x, y = y, fill = "Absence"), pch = 21, alpha = 0.5) +
      geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 1 & as.data.frame(ebird_sp)$date_index == 1, ],
                 aes(x = x, y = y, fill = "Presence"), pch = 23, cex = 2) +
      theme_void() +
      coord_equal() +
      scale_fill_discrete(name = 'eBird') +
      scale_color_discrete(name = 'Atlas') +
      ggtitle("1980-1999") +
      theme(plot.title = element_text(colour = "white", hjust = 0.5),
            legend.text=element_text(colour = "white", size=13),
            legend.title=element_text(colour = "white", size=13))

plot_2000s <- ggplot() +
      geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group),
                   colour = "white", fill = NA) +
      geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 0 & as.data.frame(atlas_sp)$date_index == 2, ],
                 aes(x = x, y = y, col = "Absence"), alpha = 0.8, pch = 15, cex = 2) +
      geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 1 & as.data.frame(atlas_sp)$date_index == 2, ],
                 aes(x = x, y = y, col = "Presence"), alpha = 0.8, pch = 15, cex = 3) +
      geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 0 & as.data.frame(ebird_sp)$date_index == 2, ],
                 aes(x = x, y = y, fill = "Absence"), pch = 21, alpha = 0.5) +
      geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 1 & as.data.frame(ebird_sp)$date_index == 2, ],
                 aes(x = x, y = y, fill = "Presence"), pch = 23, cex = 2) +
      theme_void() +
      coord_equal() +
      scale_fill_discrete(name = 'eBird') +
      scale_color_discrete(name = 'Atlas') +
      ggtitle("2000-2020") +
      theme(plot.title = element_text(colour = "white", hjust = 0.5),
            legend.text=element_text(colour = "white", size=13),
            legend.title=element_text(colour = "white", size=13))

raw_plots_comb <- ggpubr::ggarrange(plot_80s, plot_2000s, nrow = 2, common.legend = TRUE, legend = "right")

ggsave(plot = raw_plots_comb, filename = paste0("presentation_files/", sub(" ", "_", species), "_raw_ccurrence.png"),
       width = 18, height = 18, units = 'cm')


# ------------------------------------------------------------------------------------------------------------------------
# Mesh resolutions
# ------------------------------------------------------------------------------------------------------------------------
region.bdry <- inla.sp2segment(TZ_outline)

#create triangle based mesh from polygon boundary and mesh parameters set in the function call
mesh_fine <- inla.mesh.2d(boundary = region.bdry, cutoff = 0.6/2, 
                     max.edge = c(0.6, 0.6*4), offset = c(0.6, 0.6*5))
mesh_coarse <- inla.mesh.2d(boundary = region.bdry, cutoff = 1.2/2, 
                          max.edge = c(1.2, 1.2*4), offset = c(1.2, 0.6*5))

png("presentation_files/mesh_fine.png",
    bg = "transparent", res = 300, width = 20, height = 20, units = "cm")
plot(mesh_fine, vertex.color = "white", edge.color = "white")
dev.off()

png("presentation_files/mesh_coarse.png",
    bg = "transparent", res = 300, width = 20, height = 20, units = "cm")
plot(mesh_coarse, vertex.color = "white", edge.color = "white")
dev.off()

# ------------------------------------------------------------------------------------------------------------------------
# Effects plots
# ------------------------------------------------------------------------------------------------------------------------
covs_chocies <- c("TZ_BG")

original_values <- data.frame(orig_values = unlist(all.seq)) %>% 
      mutate(covariate = sub("*.seq\\d+", "", rownames(.)),
             sequence = as.numeric(gsub("\\D", "", rownames(.))))

scale_params <- median(model$summary.lincomb.derived$`0.5quant`[grep("TZ_HFP|TZ_ann_rain|TZ_BG|TZ_dryspell|TZ_max_temp|TZ_HFP", 
                                                                     rownames(model$summary.lincomb.derived))])
effect_combs <- data.frame(covariate = gsub('[[:digit:]]+', '', sub("*_lc\\d+", "", rownames(model$summary.lincomb.derived))),
                           sequence = as.numeric(gsub("\\D", "", rownames(model$summary.lincomb.derived))),
                           quant_05 = inla.link.cloglog((model$summary.lincomb.derived$`0.5quant` - scale_params) - 0.36, inverse = TRUE),
                           quant_0025 = inla.link.cloglog((model$summary.lincomb.derived$`0.025quant` - scale_params) - 0.36, inverse = TRUE),
                           quant_0975 = inla.link.cloglog((model$summary.lincomb.derived$`0.975quant` - scale_params) - 0.36, inverse = TRUE))
effect_combs_main <- effect_combs %>% 
      filter(covariate %in% covs_chocies)
effect_combs_m <- merge(original_values, effect_combs_main) %>%
      mutate(orig_values = if_else(str_detect(covariate, "log"), exp(orig_values) - 1, orig_values))

facet_labels <- c(
      TZ_ann_rain = "Annual rainfall",
      TZ_BG = "Bareground cover",
      TZ_dryspell = "Longest dryspell duration",
      TZ_max_temp = "Hottest temperature",
      TZ_HFP = "Human footprint"
)

# Make effect plots, using linear combinations
effects_plot <- ggplot(effect_combs_m) +
      geom_line(aes(x = orig_values, y = quant_05), col = "blue") +
      geom_line(aes(x = orig_values, y = quant_0025), col = "white", lty = 2, alpha = .5) +
      geom_line(aes(x = orig_values, y = quant_0975), col = "white", lty = 2, alpha = .5) +
      # geom_rug(data = atlas_filtered, aes(x = )) +
      # facet_wrap(~ covariate, scale = 'free_x', labeller = as_labeller(facet_labels)) +
      theme_minimal() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle("Kori bustard") +
      xlab("Bare ground cover") +
      ylab("Probability of occurrence") +
      theme(plot.title = element_text(color = "white", hjust = 0.5 ),
            axis.title = element_text(color = "white"),
            axis.text = element_text(color = "white"),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white"),
            panel.background = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent"),
            strip.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = NA, color = "white"),
            panel.grid.major = element_line(color = "white", size = 0.01),
            panel.grid.minor = element_line(color = "white", size = 0.01))

png(paste0("presentation_files/", sub(" ", "_", species), "_", paste0(covs_chocies, collapse = "_"), ".png"),
    bg = "transparent", res = 300, width = 15, height = 15, units = "cm")
effects_plot
dev.off()

# ------------------------------------------------------------------------------------------------------------------------
# Simulated non-linear plot
# ------------------------------------------------------------------------------------------------------------------------
# Set seed for reproducibility
set.seed(42)

# Generate simulated data
n <- 100
env_covariate <- seq(-10, 10, length.out = n)
gaussian <- function(x, mean, sd, exponent) {
      exp(-(abs(x - mean)^exponent) / (2 * sd^2))
}
species_occurrence <- gaussian(env_covariate, mean = 0, sd = 5, exponent = 3) + rnorm(n, sd = 0.2)
species_occurrence <- pmin(pmax(species_occurrence, 0), 1)
simulated_data <- data.frame(env_covariate, species_occurrence)

# Create a ggplot with points and a non-linear smoothed line
sim_plot <- ggplot(simulated_data, aes(x = env_covariate, y = species_occurrence)) +
      geom_point(col = "white") +
      geom_smooth(method = "loess", formula = "y ~ x", span = 0.6, se = FALSE, color = "blue") +
      labs(x = "Environmental Covariate Values",
           y = "Probability of Species Occurrence") +
      ylim(0, 1) +
      xlim(-8, 8) +
      theme_minimal() +
      theme(plot.title = element_text(color = "white", hjust = 0.5 ),
            axis.title = element_text(color = "white"),
            axis.text = element_text(color = "white"),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white"),
            panel.background = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent"),
            strip.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = NA, color = "white"),
            panel.grid.major = element_line(color = "white", size = 0.01),
            panel.grid.minor = element_line(color = "white", size = 0.01))

png("presentation_files/siim_plot_nonLinear.png",
    bg = "transparent", res = 300, width = 7, height = 10, units = "cm")
sim_plot
dev.off()




 # ------------------------------------------------------------------------------------------------------------------------
# Prediction plots
# ------------------------------------------------------------------------------------------------------------------------
TZ_no_lakes_vect <- vect(TZ_no_lakes)

pred.index <- inla.stack.index(stack = integrated_stack, tag = "pred")$data

IP_df <- data.frame(SP_Points_data) %>% dplyr::select(date_index, x, y)

IP_df$pred_median <- model$summary.fitted.values[pred.index, "0.5quant"]
IP_df$pred_sd <- model$summary.fitted.values[pred.index, "sd"]

pred_data <- data.frame(all_pred = inla.link.cloglog(IP_df$pred_median, inv = T),
                        sd = IP_df$pred_sd,
                        ind = rep(c(1,2), each = length(IP_df$pred_median)/2))
predcoordsGroup <- do.call(rbind, list(predcoords, predcoords))
pred_data_coords <- cbind(pred_data, predcoordsGroup) %>% 
      dplyr::select(x, y, everything())

preds_1 <- pred_data_coords %>% 
      filter(ind == 1) %>% 
      rast(., type = "xyz") %>% 
      mask(., TZ_no_lakes_vect)

preds_2 <- pred_data_coords %>% 
      filter(ind == 2) %>% 
      rast(., type = "xyz") %>% 
      mask(., TZ_no_lakes_vect)

median_plot_1 <- ggplot() +
      geom_spatraster(data = preds_1$all_pred) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "white") +
      coord_sf() +
      ggtitle("1980-1999") +
      scale_fill_viridis(name = "P(occurrence)", na.value="transparent") +
      theme_void() +
      theme(plot.title = element_text(color = "white"),
            plot.subtitle = element_text(color = "white"),
            axis.title = element_text(color = "white"),
            axis.text = element_text(color = "white"),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white"))

median_plot_2 <- ggplot() +
      geom_spatraster(data = preds_2$all_pred) +
      tidyterra::geom_spatvector(data = TZ_no_lakes_vect, fill = "transparent", colour = "white") +
      coord_sf() +
      ggtitle("2000-2020") +
      scale_fill_viridis(name = "P(occurrence)", na.value="transparent") +
      theme_void() +
      theme(plot.title = element_text(color = "white"),
            plot.subtitle = element_text(color = "white"),
            axis.title = element_text(color = "white"),
            axis.text = element_text(color = "white"),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white"))

maps_combined2 <- (median_plot_1 / median_plot_2) 

png(paste0("figures/median_", sub(" ", "_", species), "_E", max.edge, "_r",
           prior_range[1], "_", prior_range[2], "_s",
           prior_sigma[1], "_", prior_sigma[2], ".png"),
    bg = "transparent", res = 300, width = 12, height = 12, units = "cm")
median_plot_2
dev.off()

# ------------------------------------------------------------------------------------------------------------------------
# ROC plots
# ------------------------------------------------------------------------------------------------------------------------
all_obs <- rbind(atlas_sp_simple, ebird_sp_simple)
all_obs_80s <- all_obs %>% filter(ind == '1980-1999') 
all_obs_20s <- all_obs %>% filter(ind == '2000-2020')

all_obs_data <- all_pred_df %>% 
      dplyr::select(x, y, ind, all_pred, sd) %>% 
      mutate(pred_median = inla.link.cloglog(all_pred, inv = TRUE),
             pred_sd = sd)

pred_data_80s <- all_obs_data %>% filter(ind == '1980-1999') 
pred_data_20s <- all_obs_data %>% filter(ind == '2000-2020')

pred_raster_80s <- rasterFromXYZ(pred_data_80s[, c("x", "y", "pred_median")])
pred_raster_20s <- rasterFromXYZ(pred_data_20s[, c("x", "y", "pred_median")])

all_obs_80s_spdf <- SpatialPointsDataFrame(coords = data.frame(all_obs_80s$x, all_obs_80s$y), data = data.frame(presence = all_obs_80s$presence))
all_obs_20s_spdf <- SpatialPointsDataFrame(coords = data.frame(all_obs_20s$x, all_obs_20s$y), data = data.frame(presence = all_obs_20s$presence))

ext_80s <- extract(pred_raster_80s, all_obs_80s_spdf)
ext_20s <- extract(pred_raster_20s, all_obs_20s_spdf)

all_obs_80s$pred <- ext_80s
all_obs_20s$pred <- ext_20s

all_obs_80s <- all_obs_80s %>% 
      filter(!is.na(pred))
all_obs_20s <- all_obs_20s %>% 
      filter(!is.na(pred))

auc_80s <- round(as.numeric(auc(all_obs_80s$presence, all_obs_80s$pred)), digits = 3)
auc_20s <- round(as.numeric(auc(all_obs_20s$presence, all_obs_20s$pred)), digits = 3)

roc_80s <- roc(all_obs_80s$presence, all_obs_80s$pred)
roc_data_80s <- data.frame(
      sensitivity = roc_80s$sensitivities,
      fpr = 1 - roc_80s$specificities
)

roc_20s <- roc(all_obs_20s$presence, all_obs_20s$pred)
roc_data_20s <- data.frame(
      sensitivity = roc_20s$sensitivities,
      fpr = 1 - roc_20s$specificities
)

roc_plot_80s <- ggplot(roc_data_80s, aes(x = fpr, y = sensitivity)) +
      geom_line(color = "blue", lwd = 1.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", lwd = 0.5) +
      theme_classic() +
      xlab("False Positive Rate (1 - Specificity)") +
      ylab("True Positive Rate (Sensitivity)") +
      ggtitle("ROC Curve - 1980-1999") +
      theme(plot.title = element_text(hjust = 0.5)) + 
      annotate(
            "text",
            x = 0.75,
            y = 0.25,
            label = paste("AUC =", round(auc_80s, 3)),
            size = 5,
            color = "black"
      )

roc_plot_20s <- ggplot(roc_data_20s, aes(x = fpr, y = sensitivity)) +
      geom_line(color = "blue", lwd = 1.5) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray", lwd = 0.5) +
      theme_classic() +
      xlab("False Positive Rate (1 - Specificity)") +
      ylab("True Positive Rate (Sensitivity)") +
      ggtitle("ROC Curve - 2000-2020") +
      theme(plot.title = element_text(hjust = 0.5)) + 
      annotate(
            "text",
            x = 0.75,
            y = 0.25,
            label = paste("AUC =", round(auc_20s, 3)),
            size = 5,
            color = "black"
      )

roc_combined <- roc_plot_80s | roc_plot_20s

png(paste0(sub(" ", "_", species), "_ROC.png"),
    bg = "transparent", res = 300, width = 20, height = 10, units = "cm")
roc_combined
dev.off()



# -----------------------------------------------------------------------------------------------------------------
# Overview of range changes
# -----------------------------------------------------------------------------------------------------------------
cells_change_log <- ggplot(model_data, aes(x = reorder(species, log_prop_change))) +
      geom_bar(aes(y = log_prop_change, col = "Total range change"), stat = 'identity', alpha = 0.3) +
      geom_point(aes(y = log_relative_colonisation, col = "Relative colonisation"), pch = '-', size = 7) +
      geom_point(aes(y = log_relative_extinction, col = "Relative extinction"), pch = '-', size = 7) +
      geom_segment(aes(xend = reorder(species,  log_relative_colonisation), y = log_relative_extinction, yend = log_relative_colonisation), alpha = .3, lty = 2, lwd = .3) +
      geom_hline(aes(yintercept = 0), alpha = .5) + 
      theme_minimal() +
      theme(plot.title = element_text(color = "white", hjust = 0.5 ),
            axis.title = element_text(color = "white"),
            axis.text = element_text(color = "white"),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white"),
            panel.background = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent"),
            strip.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = NA, color = "white"),
            panel.grid.major = element_line(color = "white", size = 0.01),
            panel.grid.minor = element_line(color = "white", size = 0.01)) +
      xlab("Species") +
      ylab("Relative change factor") +
      theme(legend.title = element_blank(),
            legend.position = c(.45, .8),
            legend.key=element_blank()) +
      scale_x_discrete(labels = element_blank()) +
      scale_colour_colorblind() +
      scale_y_continuous(breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40))

cor_2_rsq <- round(rsq(model_data$log_prop_change, model_data$log_relative_colonisation), digits = 2)
cor_3_rsq <- round(rsq(model_data$log_prop_change, model_data$log_relative_extinction), digits = 2)
coeff <- 2
cor_plot <- ggplot(model_data, aes(log_prop_change)) +
      geom_point(aes(y = log_relative_colonisation), col = "white") +
      geom_point(aes(y = log_relative_extinction*coeff), col = "#E69F00") +
      theme_classic() +
      xlab("Total range change") +
      ylab("Relative colonisation") +
      scale_colour_colorblind(guide = "none") +
      geom_text(aes(x = Inf, y = -Inf, hjust = 1.8, vjust = -2, label = paste0('R²= ', cor_2_rsq)), col = "white")  +
      geom_text(aes(x = Inf, y = -Inf, hjust = 3, vjust = -2, label = paste0('R²= ', cor_3_rsq)), col = "#E69F00") +
      scale_y_continuous(name = "Relative colonisation", breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40),
                         sec.axis = sec_axis(~./coeff, name="Relative extinction",
                                             breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)))  +
      scale_x_continuous(breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)) + 
      theme_minimal() +
      theme(plot.title = element_text(color = "white", hjust = 0.5 ),
            axis.title = element_text(color = "white"),
            axis.text = element_text(color = "white"),
            legend.text = element_text(color = "white"),
            legend.title = element_text(color = "white"),
            panel.background = element_rect(fill = "transparent"),
            plot.background = element_rect(fill = "transparent"),
            strip.background = element_rect(fill = "transparent"),
            panel.border = element_rect(fill = NA, color = "white"),
            panel.grid.major = element_line(color = "white", size = 0.01),
            panel.grid.minor = element_line(color = "white", size = 0.01),
            axis.line.y.right = element_line(color = "#E69F00"), 
            axis.ticks.y.right = element_line(color = "#E69F00"),
            axis.text.y.right = element_text(color = "#E69F00"),
            axis.title.y.right = element_text(color = "#E69F00"))

png("presentation_files/cells_change_log.png",
    bg = "transparent", res = 300, width = 10, height = 10, units = "cm")
cells_change_log
dev.off()

png("presentation_files/cor_plot.png",
    bg = "transparent", res = 300, width = 10, height = 10, units = "cm")
cor_plot
dev.off()


# -----------------------------------------------------------------------------------------------------------------
# Forest plots
# -----------------------------------------------------------------------------------------------------------------

make_forest_plot <- function(inla_model, label, col){
      model_fixed <- data.frame(inla_model$summary.fixed) %>% 
            mutate(ID = rownames(.))
      
      if(label == "Niche breadth"){
            model_fixed <- model_fixed %>% 
                  mutate(ID = replace(ID, ID == "(Intercept)", "Intercept"),
                         ID = replace(ID, ID == "Trophic.LevelOmnivore", "Trophic level: Omnivore"),
                         ID = replace(ID, ID == "Trophic.LevelHerbivore", "Trophic level: Herbivore"),
                         ID = replace(ID, ID == "Primary.LifestyleGeneralist", "Primary lifestyle: Generalist"),
                         ID = replace(ID, ID == "Primary.LifestyleAerial", "Primary lifestyle: Aerial"),
                         ID = replace(ID, ID == "Primary.LifestyleTerrestrial", "Primary lifestyle: Terrestrial"),
                         ID = replace(ID, ID == "Migratory_abilityhigh", "Migratory ability: High"),
                         ID = replace(ID, ID == "Migratory_abilitymoderate", "Migratory ability: Moderate"),
                         ID = replace(ID, ID == "BG_breadth", "Bare ground cover"),
                         ID = replace(ID, ID == "rain_breadth", "Annual rainfall"),
                         ID = replace(ID, ID == "temp_breadth", "Hottest temperature"),
                         ID = replace(ID, ID == "dry_breadth", "Dry spell duration"),
                         ID = replace(ID, ID == "HFP_breadth", "Human footprint"),
                         ID = replace(ID, ID == "HWI", "Hand-wing Index"),
                         ID = replace(ID, ID == "Mass", "Body mass"),
                         ID = replace(ID, ID == "avg.r", "Dorsal reflectance"))
      }
      
      if(label == "Sensitivity"){
            model_fixed <- model_fixed %>% 
                  mutate(ID = replace(ID, ID == "(Intercept)", "Intercept"),
                         ID = replace(ID, ID == "Trophic.LevelOmnivore", "Trophic level: Omnivore"),
                         ID = replace(ID, ID == "Trophic.LevelHerbivore", "Trophic level: Herbivore"),
                         ID = replace(ID, ID == "Primary.LifestyleGeneralist", "Primary lifestyle: Generalist"),
                         ID = replace(ID, ID == "Primary.LifestyleAerial", "Primary lifestyle: Aerial"),
                         ID = replace(ID, ID == "Primary.LifestyleTerrestrial", "Primary lifestyle: Terrestrial"),
                         ID = replace(ID, ID == "Migratory_abilityhigh", "Migratory ability: High"),
                         ID = replace(ID, ID == "Migratory_abilitymoderate", "Migratory ability: Moderate"),
                         ID = replace(ID, ID == "BG_imp", "Bare ground cover"),
                         ID = replace(ID, ID == "rain_imp", "Annual rainfall"),
                         ID = replace(ID, ID == "temp_imp", "Hottest temperature"),
                         ID = replace(ID, ID == "dry_imp", "Dry spell duration"),
                         ID = replace(ID, ID == "HFP_imp", "Human footprint"),
                         ID = replace(ID, ID == "HWI", "Hand-wing Index"),
                         ID = replace(ID, ID == "Mass", "Body mass"),
                         ID = replace(ID, ID == "avg.r", "Dorsal reflectance"))
      }
      
      model_fixed$ID <- factor(model_fixed$ID, levels = c("Trophic level: Omnivore", "Trophic level: Herbivore", 
                                                          "Primary lifestyle: Generalist", "Primary lifestyle: Aerial", "Primary lifestyle: Terrestrial",
                                                          "Migratory ability: High", "Migratory ability: Moderate", "Hand-wing Index", "Body mass", "Dorsal reflectance",
                                                          "Bare ground cover", "Annual rainfall", "Hottest temperature", "Dry spell duration", "Human footprint",
                                                          "Intercept"))
      model_fixed$significant <- ifelse((model_fixed$X0.025quant > 0 & model_fixed$X0.975quant > 0)|(model_fixed$X0.025quant < 0 & model_fixed$X0.975quant < 0), "yes", "no")
      
      forest_plot <- ggplot() + 
            annotate("text", y = 1.6, x = max(inla_model$summary.fixed$`0.975quant`) - 0.1*max(inla_model$summary.fixed$`0.975quant`), 
                     col = col, label = label,
                     hjust = 0,
                     angle = 90,
                     size = 3.5) +
            geom_rect(aes(ymin = 1.5, ymax = 6.5, xmin = min(inla_model$summary.fixed$`0.025quant`), 
                          xmax = max(inla_model$summary.fixed$`0.975quant`)),
                      fill = col, col = col, alpha = .3) +
            geom_point(data = model_fixed, aes(y = ID, x = X0.5quant, col = significant)) +
            geom_errorbar(data = model_fixed, aes(y = ID, xmin = X0.025quant, xmax = X0.975quant, col = significant), width = 0.1) +
            geom_vline(aes(xintercept = 0), lty = 2, alpha = .7, col = "white") +
            # geom_text(data = model_fixed %>% filter(significant == "yes"), aes(x = ID, y = X0.975quant), label = "*", nudge_x = 0.1, nudge_y = 0.01) +
            xlab("Posterior estimates") +
            ylab(element_blank()) +
            theme_minimal() +
            theme(plot.title = element_text(color = "white", hjust = 0.5 ),
                  axis.title = element_text(color = "white"),
                  axis.text = element_text(color = "white"),
                  legend.text = element_text(color = "white"),
                  legend.title = element_text(color = "white"),
                  panel.background = element_rect(fill = "transparent"),
                  plot.background = element_rect(fill = "transparent"),
                  strip.background = element_rect(fill = "transparent"),
                  panel.border = element_rect(fill = NA, color = "white"),
                  panel.grid.major = element_line(color = "white", size = 0.01),
                  panel.grid.minor = element_line(color = "white", size = 0.01)) +
            scale_y_discrete(limits=rev) +
            scale_colour_manual(values = c("white", "#D55E00"), guide = "none") 
      return(forest_plot)
}

prop_change_forest <- make_forest_plot(model_prop_change, "Sensitivity", "#0072B2") + 
      ggtitle("Total range change")

loss_forest <- make_forest_plot(model_loss_log, "Sensitivity", "#0072B2") + 
      ggtitle("Extinction") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank())

col_forest <- make_forest_plot(model_col_log, "Sensitivity", "#0072B2")  + 
      ggtitle("Colonisation") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank()) 

png("presentation_files/prop_change_forest.png",
    bg = "transparent", res = 300, width = 10, height = 10, units = "cm")
prop_change_forest
dev.off()

png("presentation_files/loss_forest.png",
    bg = "transparent", res = 300, width = 6, height = 10, units = "cm")
loss_forest
dev.off()

png("presentation_files/col_forest.png",
    bg = "transparent", res = 300, width = 6, height = 10, units = "cm")
col_forest
dev.off()
