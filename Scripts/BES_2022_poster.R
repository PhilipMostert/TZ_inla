library(INLA)
library(INLAutils)
library(brinla)
library(GGally)
library(raster)
library(readxl)
library(ggpubr)
library(ggthemes)
library(ztable)
library(magrittr)


library(tidyverse)
library(patchwork)
library(terra)
library(tidyterra)
library(INLA)
        
# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
source("source/misc_functions.R")
load('model_data/regression_data.RData')

# -----------------------------------------------------------------------------------------------------------------
# Range figure
# -----------------------------------------------------------------------------------------------------------------
load("model_output/egyptian_vulture_output.RData")
load("model_output/thrush_output.RData")
atlas_vult_short <- atlas_pres_vulture %>% 
      dplyr::select(x, y, presence, ind) %>% 
      mutate(source = "Atlas")
      
eBird_vult_short <- eBird_pres_vulture %>% 
      dplyr::select(x, y, presence, ind) %>% 
      mutate(source = "eBird")   
      
vulture_comb <- rbind(atlas_vult_short, eBird_vult_short) %>% 
      mutate(test = paste0(source, "-", presence))

atlas_thrush_short <- atlas_pres_thrush %>% 
      dplyr::select(x, y, presence, ind) %>% 
      mutate(source = "Atlas")

eBird_thrush_short <- eBird_pres_thrush %>% 
      dplyr::select(x, y, presence, ind) %>% 
      mutate(source = "eBird")   

thrush_comb <- rbind(atlas_thrush_short, eBird_thrush_short) %>% 
      mutate(test = paste0(source, "-", presence))

thrush_dist_20s <- cloglog_inv(all_obs_data_thrush$all_pred[all_obs_data_thrush$ind ==  "2000-2020"])
thrush_dist_80s <- cloglog_inv(all_obs_data_thrush$all_pred[all_obs_data_thrush$ind ==  "1980-1999"])
thrush_col <- ((1-thrush_dist_80s) * thrush_dist_20s)
thrush_ext <-  (thrush_dist_80s * (1-thrush_dist_20s))

thrush_transition <- data.frame(x = all_obs_data_thrush$x[all_obs_data_thrush$ind == "2000-2020"],
                                y = all_obs_data_thrush$y[all_obs_data_thrush$ind == "2000-2020"],
                                col = thrush_col,
                                ext = thrush_ext)


vulture_dist_20s <- cloglog_inv(all_obs_data_vulture$all_pred[all_obs_data_vulture$ind ==  "2000-2020"])
vulture_dist_80s <- cloglog_inv(all_obs_data_vulture$all_pred[all_obs_data_vulture$ind ==  "1980-1999"])
vulture_col <- ((1-vulture_dist_80s) * vulture_dist_20s)
vulture_ext <-  (vulture_dist_80s * (1-vulture_dist_20s))

vulture_transition <- data.frame(x = all_obs_data_vulture$x[all_obs_data_vulture$ind == "2000-2020"],
                                y = all_obs_data_vulture$y[all_obs_data_vulture$ind == "2000-2020"],
                                col = vulture_col,
                                ext = vulture_ext)


legend <- ggplot() +
      geom_tile(data = all_obs_data_vulture, aes(x = x, y = y, fill = cloglog_inv(all_pred))) +
      geom_point(data = vulture_comb,
                 aes(x = x, y = y, col = test, alpha = as.factor(presence), shape = test, cex = source)) +
      geom_point(data = vulture_comb %>% filter(presence == 1, source == "eBird"),
                 aes(x = x, y = y), alpha = 0.9, cex = 0.5) +
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("Probability", na.value="white", option = "D") +
      scale_color_manual(name = "Input data", values = c("orangered", "black", "orangered", "black"),
                         label = c("Atlas absence", "Atlas presence", "eBird absence", "eBird presence")) +
      scale_alpha_manual(values = c(0.3, 0.9), guide = "none") +
      scale_shape_manual(name = "Input data", values = c(15, 15, 19, 19),
                         label = c("Atlas absence", "Atlas presence", "eBird absence", "eBird presence")) +
      scale_size_manual(values = c(1.5, 0.5), guide = "none") +
      theme_void() +
      theme(legend.position = "bottom") +
      ggtitle("Egyptian vulture - Estimated occurrence")

median_plot_egyptian_vulture <- ggplot() + 
      geom_tile(data = all_obs_data_vulture, aes(x = x, y = y, fill = pred_median)) + 
      geom_point(data = eBird_pres_vulture %>% filter(presence == 0),
                 aes(x = x, y = y), col = "orangered", alpha = 0.1, cex = .2) +
      geom_point(data = atlas_pres_vulture,
                 aes(x = x, y = y, col = as.factor(presence), alpha = as.factor(presence)), pch = 15, cex = 2) +
      geom_point(data = eBird_pres_vulture %>% filter(presence == 1),
                 aes(x = x, y = y), col = "black", alpha = 1, cex = 2) +
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("Probability", na.value="white", option = "D") +
      scale_color_manual("Occurrence", values = c("orangered", "black")) +
      scale_alpha_manual(values = c(0.3, 1), guide = "none") +
      theme_void() +
      theme(legend.position = "none") +
      ggtitle("Egyptian vulture - Estimated occurrence")

median_plot_trush <- ggplot() + 
      geom_tile(data = all_obs_data_thrush, aes(x = x, y = y, fill = cloglog_inv(all_pred))) + 
      geom_point(data = thrush_comb,
                 aes(x = x, y = y, col = test, alpha = as.factor(presence), shape = test, cex = source)) +
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("Probability", na.value="white", option = "D") +
      scale_color_manual(name = "Occurrence", values = c("orangered", "black", "orangered", "black"), 
                         label = c("Atlas-abs", "Atlas-pres", "eBird-abs", "eBird-pres")) +
      scale_alpha_manual(values = c(0.3, 0.8), guide = "none") +
      scale_shape_manual(name = "Occurrence", values = c(15, 15, 19, 19), 
                         label = c("Atlas-abs", "Atlas-pres", "eBird-abs", "eBird-pres")) +
      scale_size_manual(values = c(1.5, 0.2), guide = "none") +
      theme_void() +
      theme(legend.position = "none") +
      ggtitle("Bare-eyed thrush - Estimated occurrence")

vulture_colonisations <- ggplot() + 
      geom_tile(data = vulture_transition, aes(x = x, y = y, fill = col)) +
      coord_equal() +
      viridis::scale_fill_viridis("Probability", na.value="white", option = "D")  +
      theme_void() +
      theme(legend.position = "none") +
      ggtitle("Colonisation")

vulture_extinctions <- ggplot() + 
      geom_tile(data = vulture_transition, aes(x = x, y = y, fill = ext)) +
      coord_equal() +
      viridis::scale_fill_viridis("Probability", na.value="white", option = "D")  +
      theme_void() +
      theme(legend.position = "none") +
      ggtitle("Extinction")

thrush_colonisations <- ggplot() + 
      geom_tile(data = thrush_transition, aes(x = x, y = y, fill = col)) +
      coord_equal() +
      viridis::scale_fill_viridis("Probability", na.value="white", option = "D")  +
      theme_void() +
      theme(legend.position = "none") +
      ggtitle("Colonisation")

thrush_extinctions <- ggplot() + 
      geom_tile(data = thrush_transition, aes(x = x, y = y, fill = ext)) +
      coord_equal() +
      viridis::scale_fill_viridis("Probability", na.value="white", option = "D")  +
      theme_void() +
      theme(legend.position = "none") +
      ggtitle("Extinction")

vulture_trans_plots <- (vulture_colonisations | vulture_extinctions)
thrush_trans_plots <- (thrush_colonisations | thrush_extinctions)


ggsave(filename = "/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Conferences/BES_annual_2022/median_plot_egyptian_vulture.png", 
       plot = median_plot_egyptian_vulture, width = 15, height = 20,
       units = "cm", dpi = 300)
ggsave(filename = "/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Conferences/BES_annual_2022/legend.png", 
       plot = legend, width = 25, height = 20,
       units = "cm", dpi = 300)
ggsave(filename = "/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Conferences/BES_annual_2022/median_plot_trush.png", 
       plot = median_plot_trush, width = 15, height = 20,
       units = "cm", dpi = 300)
ggsave(filename = "/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Conferences/BES_annual_2022/vulture_trans_plots.png", 
       plot = vulture_trans_plots, width = 15, height = 20,
       units = "cm", dpi = 300)
ggsave(filename = "/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Conferences/BES_annual_2022/thrush_trans_plots.png", 
       plot = thrush_trans_plots, width = 15, height = 20,
       units = "cm", dpi = 300)


form_prop_change <- log_prop_change ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + Wing.Length + avg.r

model_prop_change_robust <- inla(form_prop_change, data = model_data, 
                                 family = "T",  # t distributed errors regression model
                                 lincomb = all_lc_sens,
                                 control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                 control.predictor=list(compute=TRUE))

form_col_log <- relative_colonisation ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + Wing.Length + avg.r
model_col_robust_log <- inla(form_col_log, data = model_data, 
                             family = "T",  
                             lincomb = all_lc_sens,
                             control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                             control.predictor=list(compute=TRUE))

form_loss_log <- relative_extinction ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + Wing.Length + avg.r
model_loss_robust_log <- inla(form_loss_log, data = model_data, 
                              family = "T", 
                              lincomb = all_lc_sens,
                              control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                              control.predictor=list(compute=TRUE))

form_prop_change_breadth <- log_prop_change ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      rain_breadth + temp_breadth + dry_breadth + HFP_breadth + BG_breadth + Wing.Length + avg.r
model_prop_change_breadth_robust <- inla(form_prop_change_breadth, data = model_data, 
                                         family = "T", 
                                         lincomb = all_lc_breadth,
                                         control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                         control.predictor=list(compute=TRUE))

form_col_log_breadth <- relative_colonisation ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      rain_breadth + temp_breadth + dry_breadth + HFP_breadth + BG_breadth + Wing.Length + avg.r
model_col_log_breadth_robust <- inla(form_col_log_breadth, data = model_data, 
                                     family = "T",  
                                     lincomb = all_lc_breadth,
                                     control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                     control.predictor=list(compute=TRUE))

form_loss_log_breadth <- relative_extinction ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      rain_breadth + temp_breadth + dry_breadth + HFP_breadth + BG_breadth + Wing.Length + avg.r
model_loss_log_breadth_robust <- inla(form_loss_log_breadth, data = model_data, 
                                      family = "T",  
                                      lincomb = all_lc_breadth,
                                      control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                      control.predictor=list(compute=TRUE))

# -----------------------------------------------------------------------------------------------------------------
# Model output - forest plot
# -----------------------------------------------------------------------------------------------------------------
prop_change_forest <- make_forest_plot(model_prop_change_robust, "Sensitivity", "#0072B2") + 
      ggtitle("Total range change")
loss_forest <- make_forest_plot(model_loss_robust_log, "Sensitivity", "#0072B2") + 
      ggtitle("Extinction") +
      ylab(element_blank())
col_forest <- make_forest_plot(model_col_robust_log, "Sensitivity", "#0072B2")  + 
      ggtitle("Colonisation") +
      ylab(element_blank())  
prop_change_forest_niche <- make_forest_plot(model_prop_change_breadth_robust, "Niche breadth", "#E69F00") + 
      ggtitle("Total range change") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank()) 
loss_forest_niche <- make_forest_plot(model_loss_log_breadth_robust, "Niche breadth", "#E69F00") + 
      ggtitle("Extinction") +
      ylab(element_blank()) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank())  
col_forest_niche <- make_forest_plot(model_col_log_breadth_robust, "Niche breadth", "#E69F00") + 
      ggtitle("Colonisation") +
      ylab(element_blank()) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank())    
forest_patch <- prop_change_forest_niche + loss_forest_niche + col_forest_niche + 
      prop_change_forest + loss_forest + col_forest 

ggsave(filename = "/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Conferences/BES_annual_2022/forest_plot.png", plot = forest_patch, width = 45, height = 27,
       units = "cm", dpi = 300)
