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

# -----------------------------------------------------------------------------------------------------------------
# Data import
# -----------------------------------------------------------------------------------------------------------------
source("source/misc_functions.R")
load('model_data/regression_data.RData')
load("results/traits_model_list.RData")

# -----------------------------------------------------------------------------------------------------------------
# Overview of range changes
# -----------------------------------------------------------------------------------------------------------------
cells_change_log <- ggplot(model_data, aes(x = reorder(species, log_prop_change))) +
      geom_bar(aes(y = log_prop_change, col = "Total range change"), stat = 'identity', alpha = 0.3) +
      geom_point(aes(y = log_relative_colonisation, col = "Relative colonisation"), pch = '-', size = 7) +
      geom_point(aes(y = log_relative_extinction, col = "Relative extinction"), pch = '-', size = 7) +
      geom_segment(aes(xend = reorder(species,  log_relative_colonisation), y = log_relative_extinction, yend = log_relative_colonisation), alpha = .3, lty = 2, lwd = .3) +
      geom_hline(aes(yintercept = 0), alpha = .5) + 
      theme_classic() +
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
      geom_point(aes(y = log_relative_colonisation), col = "#000000") +
      geom_point(aes(y = log_relative_extinction*coeff), col = "#E69F00") +
      theme_classic() +
      xlab("Total range change") +
      ylab("Relative colonisation") +
      scale_colour_colorblind(guide = "none") +
      geom_text(aes(x = Inf, y = -Inf, hjust = 1.8, vjust = -2, label = paste0('R²= ', cor_2_rsq)), col = "#000000")  +
      geom_text(aes(x = Inf, y = -Inf, hjust = 3, vjust = -2, label = paste0('R²= ', cor_3_rsq)), col = "#E69F00") +
      scale_y_continuous(name = "Relative colonisation", breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40),
                         sec.axis = sec_axis(~./coeff, name="Relative extinction",
                                             breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)))  +
      scale_x_continuous(breaks = log(c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)), labels = c(0.2, 0.5, seq(0, 5, by = 1), 4, 5, 10, 20, 30, 40)) + 
      theme(axis.line.y.right = element_line(color = "#E69F00"), 
            axis.ticks.y.right = element_line(color = "#E69F00"),
            axis.text.y.right = element_text(color = "#E69F00"),
            axis.title.y.right = element_text(color = "#E69F00"))

# load("model_output/egyptian_vulture_output.RData")
# load("model_output/thrush_output.RData")
# 
# median_plot_egyptian_vulture <- ggplot() + 
#       geom_tile(data = all_obs_data_vulture, aes(x = x, y = y, fill = pred_median)) + 
#       geom_point(data = eBird_pres_vulture %>% filter(presence == 0),
#                  aes(x = x, y = y), col = "orangered", alpha = 0.1, cex = .2) +
#       geom_point(data = atlas_pres_vulture,
#                  aes(x = x, y = y, col = as.factor(presence), alpha = as.factor(presence)), pch = 15, cex = 2) +
#       geom_point(data = eBird_pres_vulture %>% filter(presence == 1),
#                  aes(x = x, y = y), col = "black", alpha = 1, cex = 2) +
#       facet_grid( ~ ind) +
#       coord_equal() +
#       viridis::scale_fill_viridis("Probability", na.value="white", option = "D") +
#       scale_color_manual("Occurrence", values = c("orangered", "black")) +
#       scale_alpha_manual(values = c(0.3, 1), guide = "none") +
#       theme_void() +
#       theme(legend.position = "none") +
#       ggtitle("Egyptian vulture range - Model estimate")
# 
# median_plot_trush <- ggplot() + 
#       geom_tile(data = all_obs_data_thrush, aes(x = x, y = y, fill = pred_median)) + 
#       geom_point(data = eBird_pres_thrush %>% filter(presence == 0),
#                  aes(x = x, y = y), col = "orangered", alpha = 0.1, cex = .2) +
#       geom_point(data = atlas_pres_thrush,
#                  aes(x = x, y = y, col = as.factor(presence), alpha = as.factor(presence)), pch = 15, cex = 2) +
#       geom_point(data = eBird_pres_thrush %>% filter(presence == 1),
#                  aes(x = x, y = y), col = "black", alpha = 1, cex = 2) +
#       facet_grid( ~ ind) +
#       coord_equal() +
#       viridis::scale_fill_viridis("Probability", na.value="white", option = "D") +
#       scale_color_manual("Occurrence", values = c("orangered", "black")) +
#       scale_alpha_manual(values = c(0.3, 1), guide = "none") +
#       theme_void() +
#       theme(legend.position = "none") +
#       ggtitle("Bare-eyed thrush range - Model estimate")

# combined_range_plots <- (cells_change_log + cor_plot) / median_plot_egyptian_vulture / median_plot_trush + plot_annotation(tag_levels = 'A')
# ggsave(filename = "figures/final_figures/combined_range_plots.png", plot = combined_range_plots, width = 15, height = 20,
#        units = "cm", dpi = 300)

combined_range_plots_2 <- (cells_change_log | cor_plot)
png("figures/final_figures/combined_range_plots.png",
    bg = "transparent", res = 300, width = 20, height = 10, units = "cm")
combined_range_plots_2
dev.off()

# -----------------------------------------------------------------------------------------------------------------
# Model output - forest plot
# -----------------------------------------------------------------------------------------------------------------
prop_change_forest <- make_forest_plot(model_prop_change, "Sensitivity", "#0072B2") + 
      ggtitle("Total range change")
# geom_text(aes(x = Inf, y = -Inf, hjust = 1, vjust = -0.8, 
#               label = paste0('R²= ', round(rsq(model_prop_change$summary.fitted.values[, "0.5quant"], model_data$log_prop_change), digits = 2))))

loss_forest <- make_forest_plot(model_loss_log, "Sensitivity", "#0072B2") + 
      ggtitle("Extinction") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank())
# geom_text(aes(x = Inf, y = -Inf, hjust = 1, vjust = -0.8, 
#               label = paste0('R²= ', round(rsq(model_loss_log$summary.fitted.values[, "0.5quant"], model_data$relative_extinction), digits = 2))))

col_forest <- make_forest_plot(model_col_log, "Sensitivity", "#0072B2")  + 
      ggtitle("Colonisation") +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank()) 
# geom_text(aes(x = Inf, y = -Inf, hjust = 1, vjust = -0.8, 
#               label = paste0('R²= ', round(rsq(model_col_log$summary.fitted.values[, "0.5quant"], model_data$relative_colonisation), digits = 2))))

prop_change_forest_niche <- make_forest_plot(model_prop_change_breadth, "Niche breadth", "#E69F00") + 
      ggtitle("Total range change") +
      xlab(element_blank()) 
# geom_text(aes(x = Inf, y = -Inf, hjust = 1, vjust = -0.8, 
#               label = paste0('R²= ', round(rsq(model_prop_change_breadth$summary.fitted.values[, "0.5quant"], model_data$log_prop_change), digits = 2))))

loss_forest_niche <- make_forest_plot(model_loss_log_breadth, "Niche breadth", "#E69F00") + 
      ggtitle("Extinction") +
      xlab(element_blank()) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank())
# geom_text(aes(x = Inf, y = -Inf, hjust = 1, vjust = -0.8, 
#               label = paste0('R²= ', round(rsq(model_loss_log_breadth$summary.fitted.values[, "0.5quant"], model_data$relative_extinction), digits = 2))))

col_forest_niche <- make_forest_plot(model_col_log_breadth, "Niche breadth", "#E69F00") + 
      ggtitle("Colonisation") +
      xlab(element_blank()) +
      theme(axis.title.y=element_blank(),
            axis.text.y=element_blank())
# geom_text(aes(x = Inf, y = -Inf, hjust = 1, vjust = -0.8, 
#               label = paste0('R²= ', round(rsq(model_col_log_breadth$summary.fitted.values[, "0.5quant"], model_data$relative_colonisation), digits = 2))))

forest_patch <- prop_change_forest_niche + loss_forest_niche + col_forest_niche + 
      prop_change_forest + loss_forest + col_forest 
forest_annotated <- forest_patch + plot_annotation(tag_levels = 'A')

ggsave(filename = "figures/final_figures/forest_plots_comb.png", plot = forest_annotated, width = 22, height = 20,
       units = "cm", dpi = 400)

# -----------------------------------------------------------------------------------------------------------------
# Selected plots showing raw data with model results
# -----------------------------------------------------------------------------------------------------------------
figure_data <- model_data
figure_data$Mass <- exp(unscale_fun(model_data$Mass))
col_list <- c('BG_imp', 'rain_imp', 'temp_imp', 'dry_imp', 'HFP_imp', 
              'rain_breadth', 'temp_breadth', 'dry_breadth', 'HFP_breadth', 'BG_breadth',
              'avg.r', 'HWI')
for (i in 1:length(col_list)){
      figure_data[, col_list[i]] <- unscale_fun(figure_data[, col_list[i]])
}

relative_scale <- scale_y_continuous(breaks = log(c(seq(0.5, 3, by = 0.5), 4, 5, 10, 20, 30, 40)), labels = c(seq(0.5, 3, by = 0.5), 4, 5, 10, 20, 30, 40))

# Need niche rain and dry total change, total change wing niche, hfp sens extinct
dry_niche_total <-  make_raw_data_plots(model_prop_change_breadth, figure_data, "dry_breadth", newscale = relative_scale)
rain_niche_total <-  make_raw_data_plots(model_prop_change_breadth, figure_data, "rain_breadth", newscale = relative_scale)
mass_niche_total <-  make_raw_data_plots(model_prop_change_breadth, figure_data, "Mass", newscale = relative_scale)
HFP_sens_total <-  make_raw_data_plots(model_loss_log, figure_data, "HFP_imp", newscale = relative_scale)
make_raw_data_plots(model_loss_log, figure_data, "avg.r", newscale = relative_scale)


# Ecological generalisation: Bare ground niche breadth - colonisation
bareground_niche_plot <- make_raw_data_plots(model_col_log_breadth, figure_data, "BG_breadth", newscale = relative_scale)

# Ecological generalisation: Annual rainfall niche breadth - colonisation
rain_niche_plot <- make_raw_data_plots(model_col_log_breadth, figure_data, "rain_breadth", newscale = relative_scale)


lincomb_data3 <- data.frame(model_col_log_breadth$summary.lincomb.derived) %>% 
      mutate(ID = rownames(.)) %>% 
      dplyr::select(-mode, -kld, -mean, -sd)
lincomb_data3 <- lincomb_data3[grepl("Migratory", rownames(lincomb_data3)) == TRUE, ] %>% 
      mutate(Migratory_ability = c("low", "moderate", "high"),
             Variable = "Colonisation")

lincomb_data4 <- data.frame(model_loss_log_breadth$summary.lincomb.derived) %>% 
      mutate(ID = rownames(.)) %>% 
      dplyr::select(-mode, -kld, -mean, -sd)
lincomb_data4 <- lincomb_data4[grepl("Migratory", rownames(lincomb_data4)) == TRUE, ] %>% 
      mutate(Migratory_ability = c("low", "moderate", "high"),
             Variable = "Extinction")
lincomb_migratory <- rbind(lincomb_data3, lincomb_data4) %>% 
      dplyr::select(Migratory_ability, Variable, X0.5quant)

migration_data <- model_data %>% 
      dplyr::select(Migratory_ability, relative_extinction, relative_colonisation) %>% 
      rename(Extinction = relative_extinction, Colonisation = relative_colonisation) %>% 
      pivot_longer(2:3, names_to = "Variable", values_to = "Score")
migration_data <- merge(migration_data, lincomb_migratory, all.x = TRUE)

migrationPlot <- ggplot(migration_data) +
      geom_boxplot(aes(Migratory_ability, Score, col = Variable)) +
      geom_point(aes(Migratory_ability, X0.5quant, group = Variable), 
                 col = "red", cex = 2, position = position_dodge(width = 0.75)) +
      theme_linedraw() +
      xlab("Migratory ability") +
      ylab("Relative score") +
      theme(legend.position = c(.5, .9),
            legend.title = element_blank()) +
      scale_colour_manual(values = c("#000000", "#E69F00")) +
      relative_scale

raw_data_plots <- dry_niche_total + rain_niche_total + mass_niche_total + HFP_sens_total + plot_annotation(tag_levels = 'A') 

ggsave(filename = "figures/final_figures/raw_data_plots.png", plot = raw_data_plots, width = 22, height = 15,
       units = "cm", dpi = 400)



# -----------------------------------------------------------------------------------------------------------------
# Supplementary table of species data
# -----------------------------------------------------------------------------------------------------------------
species_table <- figure_data %>% 
      dplyr::select(species, log_prop_change, relative_extinction, relative_colonisation, everything()) 

options(ztable.type="viewer")
z = ztable(species_table) %>% makeHeatmap(margin = 2)
z
# -----------------------------------------------------------------------------------------------------------------
# Supplementary figures: All raw data plots
# -----------------------------------------------------------------------------------------------------------------
list_change_sens <- make_raw_data_plots(model_prop_change, figure_data, newscale = relative_scale)
plot_change_sens <- wrap_plots(list_change_sens, ncol = 3) + plot_annotation(tag_levels = 'A', title = "Total range change - Sensitivity") 
png(filename = "figures/final_figures/plot_change_sens.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_change_sens
dev.off()

list_loss_sens <- make_raw_data_plots(model_loss_log, figure_data, newscale = relative_scale)
plot_loss_sens <- wrap_plots(list_loss_sens) + plot_annotation(tag_levels = 'A', title = "Meaningful extinction - Sensitivity") 
png(filename = "figures/final_figures/plot_loss_sens.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_loss_sens
dev.off()

list_col_sens <- make_raw_data_plots(model_col_log, figure_data, newscale = relative_scale)
plot_col_sens <- wrap_plots(list_col_sens) + plot_annotation(tag_levels = 'A', title = "Meaningful colonisation - Sensitivity") 
png(filename = "figures/final_figures/plot_col_sens.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_col_sens
dev.off()

list_change_spec <- make_raw_data_plots(model_prop_change_breadth, figure_data, newscale = relative_scale)
plot_change_spec <- wrap_plots(list_change_spec) + plot_annotation(tag_levels = 'A', title = "Total range change - Niche breadth") 
png(filename = "figures/final_figures/plot_change_spec.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_change_spec
dev.off()

list_loss_spec <- make_raw_data_plots(model_loss_log_breadth, figure_data, newscale = relative_scale)
plot_loss_spec <- wrap_plots(list_loss_spec) + plot_annotation(tag_levels = 'A', title = "Meaningful extinction - Niche breadth") 
png(filename = "figures/final_figures/plot_loss_spec.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_loss_spec
dev.off()

list_col_spec <- make_raw_data_plots(model_col_log_breadth, figure_data, newscale = relative_scale)
plot_col_spec <- wrap_plots(list_col_spec) + plot_annotation(tag_levels = 'A', title = "Meaningful colonisation - Niche breadth") 
png(filename = "figures/final_figures/plot_col_spec.png", width = 30, height = 24,
    units = "cm", res = 400)
plot_col_spec
dev.off()
