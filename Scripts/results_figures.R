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
# Data import & preparation
# -----------------------------------------------------------------------------------------------------------------
source("source/misc_functions.R")
species_results <- read.csv('model_output/species_results.csv') %>% dplyr::select(-X)
traits <- readxl::read_xlsx('model_data/ELEData/TraitData/AVONET2_eBird.xlsx', sheet = "AVONET2_eBird")
colours <-  read.csv('model_data/plumage_lightness.csv') %>%
      mutate(species = gsub(pattern = "_", replacement = " ", species)) %>% 
      filter(sex == "F")
# Higher reflectance = lighter colors = like hotter areas, lesser heat load 
# Darker birds are limited from occupying areas with high temperatures

species_results$species[species_results$species == "Estrilda erythronotos"] <- "Brunhilda erythronotos"
species_results$species[species_results$species == "Riparia cincta"] <- "Neophedina cincta"
species_results$species[species_results$species == "Sylvia boehmi"] <- "Curruca boehmi"
species_results$species[species_results$species == "Sylvia communis"] <- "Curruca communis"
species_results$species[species_results$species == "Turdoides aylmeri"] <- "Argya aylmeri"

species_results$temp_breadth <- species_results$temp_max - species_results$temp_min
species_results$rain_breadth <- species_results$rain_max - species_results$rain_min
species_results$dry_breadth <- species_results$dry_max - species_results$dry_min
species_results$HFP_breadth <- species_results$HFP_max - species_results$HFP_min
species_results$BG_breadth <- species_results$BG_max - species_results$BG_min

# Add trait data from Tobias et al. 2022
df_merged <- merge(species_results, traits, all.x = TRUE, all.y = FALSE, by.x = "species", by.y = "Species2") %>% 
      mutate(Trophic.Level = factor(Trophic.Level, levels = c("Omnivore", "Carnivore", "Herbivore")), 
             Primary.Lifestyle = as.factor(Primary.Lifestyle))
df_merged <- merge(df_merged, colours, all.x = TRUE,by = "species")

df_merged$Migratory_ability[df_merged$Migration == 1] <- "low"
df_merged$Migratory_ability[df_merged$Migration == 2] <- "moderate"
df_merged$Migratory_ability[df_merged$Migration == 3] <- "high"
df_merged <- df_merged %>% 
      mutate(Migratory_ability = as.factor(Migratory_ability))

species_long <- df_merged %>% 
      mutate_if(is.numeric, round, digits = 3) %>% 
      pivot_longer(cols = c("rain_imp", "temp_imp", "dry_imp", "BG_imp", "HFP_imp" ), names_to = 'covariate', values_to = 'perc_explained') %>% 
      drop_na() %>% 
      mutate(covariate = factor(covariate, levels = c("temp_imp", "rain_imp", "BG_imp", "dry_imp", "HFP_imp")))

model_data <- df_merged %>% 
      mutate(BG_imp = scale(BG_imp),
             rain_imp = scale(rain_imp),
             temp_imp = scale(temp_imp),
             dry_imp = scale(dry_imp),
             HFP_imp = scale(HFP_imp),
             Wing.Length = scale(Wing.Length),
             avg.r = scale(avg.r),
             rain_breadth = scale(rain_breadth),
             temp_breadth = scale(temp_breadth),
             dry_breadth = scale(dry_breadth),
             HFP_breadth = scale(HFP_breadth),
             BG_breadth = scale(BG_breadth),
             log_cells_lost_norm = log(cells_lost_norm),
             log_cells_colonised_norm = log(cells_colonised_norm)) %>% 
      filter(Trophic.Level != "Scavenger",
             species != "Neophron percnopterus")

# -----------------------------------------------------------------------------------------------------------------
# Overview of range changes
# -----------------------------------------------------------------------------------------------------------------
cells_change_log <- ggplot(model_data, aes(x = reorder(species, log_prop_change))) +
      geom_bar(aes(y = log_prop_change, col = "Total range change"), stat = 'identity', alpha = 0.3) +
      geom_point(aes(y = log(cells_colonised_norm), col = "Colonisation"), pch = '-', size = 7) +
      geom_point(aes(y = log(cells_lost_norm), col = "Extinction"), pch = '-', size = 7) +
      geom_segment(aes(xend = reorder(species,  log(cells_colonised_norm)), y = log(cells_lost_norm), yend =  log(cells_colonised_norm)), alpha = .3, lty = 2, lwd = .3) +
      geom_hline(aes(yintercept = 0), alpha = .5) + 
      theme_classic() +
      xlab("Species") +
      ylab("Score") +
      theme(legend.title = element_blank(),
            legend.position = c(.3, .9)) +
      scale_x_discrete(labels = element_blank()) +
      scale_colour_colorblind()

cor_2_rsq <- round(rsq(model_data$log_prop_change, model_data$log_cells_colonised_norm), digits = 2)
cor_3_rsq <- round(rsq(model_data$log_prop_change, model_data$log_cells_lost_norm), digits = 2)

cor_2 <- ggplot(model_data, aes(log_prop_change, log_cells_colonised_norm)) +
      geom_point() +
      theme_classic() +
      xlab("Total range change") +
      ylab("Colonisation") +
      geom_text(aes(x = Inf, y = -Inf, hjust = 1.5, vjust = -2, label = paste0('R²= ', cor_2_rsq))) 

cor_3 <- ggplot(model_data, aes(log_prop_change, log_cells_lost_norm)) +
      geom_point() +
      theme_classic() +
      xlab("Total range change") +
      ylab("Extinction") +
      geom_text(aes(x = Inf, y = -Inf, hjust = 1.5, vjust = -2, label = paste0('R²= ', cor_3_rsq))) 


combined_range_plots <- cells_change_log | (cor_2/cor_3)  
annotated <- combined_range_plots + plot_annotation(tag_levels = 'A')
ggsave(filename = "figures/final_figures/combined_range_plots.png", plot = annotated, width = 15, height = 10,
       units = "cm", dpi = 300)


# -----------------------------------------------------------------------------------------------------------------
# Check specialization vs. sensitivities
# -----------------------------------------------------------------------------------------------------------------
ggplot(df_merged, aes(x = temp_breadth, y = temp_imp)) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic()

ggplot(df_merged, aes(x = rain_breadth, y = rain_imp)) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic() 

ggplot(df_merged, aes(x = dry_breadth, y = dry_imp)) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic() 

ggplot(df_merged, aes(x = HFP_breadth, y = HFP_imp)) +
      geom_point() +
      geom_smooth(method = "lm") +
      theme_classic() 

# -----------------------------------------------------------------------------------------------------------------
# Statistical test
# -----------------------------------------------------------------------------------------------------------------

BG_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                            BG_imp = seq(min(model_data$BG_imp), max(model_data$BG_imp), len = 100))
names(BG_lc) <- paste0("BG_imp", 1:100)
rain_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              rain_imp = seq(min(model_data$rain_imp), max(model_data$rain_imp), len = 100))
names(rain_lc) <- paste0("rain_imp", 1:100)
temp_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              temp_imp = seq(min(model_data$temp_imp), max(model_data$temp_imp), len = 100))
names(temp_lc) <- paste0("temp_imp", 1:100)
dry_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             dry_imp = seq(min(model_data$dry_imp), max(model_data$dry_imp), len = 100))
names(dry_lc) <- paste0("dry_imp", 1:100)
HFP_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             HFP_imp = seq(min(model_data$HFP_imp), max(model_data$HFP_imp), len = 100))
names(HFP_lc) <- paste0("HFP_imp", 1:100)
wing_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              Wing.Length = seq(min(model_data$Wing.Length), max(model_data$Wing.Length), len = 100))
names(wing_lc) <- paste0("Wing.Length", 1:100)
reflect_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                 avg.r = seq(min(model_data$avg.r, na.rm = TRUE), max(model_data$avg.r, na.rm = TRUE), len = 100))
names(reflect_lc) <- paste0("avg.r", 1:100)

Migratory_ability_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
                                           Migratory_abilitylow = c(0, 1, 0),
                                           Migratory_abilitymoderate = c(0, 0, 1))

all_lc_sens <- c(BG_lc, rain_lc, temp_lc, dry_lc, HFP_lc, wing_lc, reflect_lc)

BG_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                            BG_breadth = seq(min(model_data$BG_breadth), max(model_data$BG_breadth), len = 100))
names(BG_lc_niche) <- paste0("BG_breadth", 1:100)
rain_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              rain_breadth = seq(min(model_data$rain_breadth), max(model_data$rain_breadth), len = 100))
names(rain_lc_niche) <- paste0("rain_breadth", 1:100)
temp_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              temp_breadth = seq(min(model_data$temp_breadth), max(model_data$temp_breadth), len = 100))
names(temp_lc_niche) <- paste0("temp_breadth", 1:100)
dry_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             dry_breadth = seq(min(model_data$dry_breadth), max(model_data$dry_breadth), len = 100))
names(dry_lc_niche) <- paste0("dry_breadth", 1:100)
HFP_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             HFP_breadth = seq(min(model_data$HFP_breadth), max(model_data$HFP_breadth), len = 100))
names(HFP_lc_niche) <- paste0("HFP_breadth", 1:100)
wing_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              Wing.Length = seq(min(model_data$Wing.Length), max(model_data$Wing.Length), len = 100))
names(wing_lc) <- paste0("Wing.Length", 1:100)
reflect_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                 avg.r = seq(min(model_data$avg.r, na.rm = TRUE), max(model_data$avg.r, na.rm = TRUE), len = 100))
names(reflect_lc) <- paste0("avg.r", 1:100)

Migratory_ability_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
                                           Migratory_abilitylow = c(0, 1, 0),
                                           Migratory_abilitymoderate = c(0, 0, 1))

all_lc_niche <- c(BG_lc_niche, rain_lc_niche, temp_lc_niche, dry_lc_niche, HFP_lc_niche, wing_lc, reflect_lc)


# Wang et al. 2018: Experience suggests, that unless we specify a sharp prior that it is in contradiction with the data, 
# this will not make much difference to the posterior distribution. One is usually satisfied with the default flat prior 
# for the intercept in such models.

form_prop_change <- log_prop_change ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + Wing.Length + avg.r

# model_prop_change_gauss <- inla(form_prop_change, data = model_data, 
#                family = "gaussian", 
#                lincomb = all_lc_sens,
#                control.compute = list(dic = TRUE, cpo = TRUE), 
#                control.predictor=list(compute=TRUE))
# brinla::bri.lmresid.plot(model_prop_change_gauss)  # Outliers in residuals
# hist(model_gauss$cpo$pit, main="", breaks = 10, xlab = "PIT")
# qqplot(qunif(ppoints(length(model_gauss$cpo$pit))), 
#        model_gauss$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
# qqline(model_gauss$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
# Outliers in residuals, and skewed Q-Q plot. Try robust model, reducing effect of outliers.

model_prop_change_robust <- inla(form_prop_change, data = model_data, 
                    family = "T",  # t distributed errors regression model
                    lincomb = all_lc_sens,
                    control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                    control.predictor=list(compute=TRUE))

# hist(model_prop_change_robust$cpo$pit, main="", breaks = 10, xlab = "PIT")
# qqplot(qunif(ppoints(length(model_prop_change_robust$cpo$pit))), model_prop_change_robust$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
# qqline(model_prop_change_robust$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
# Q-Q plot looks much better. Compare the two models using CPO and PIT

# sum(model_gauss$cpo$failure)
# sum(model_prop_change_robust$cpo$failure)
# c(-sum(log(model_gauss$cpo$cpo)), -sum(log(model_prop_change_robust$cpo$cpo)))
# c(-sum(log(model_gauss$cpo$pit)), -sum(log(model_prop_change_robust$cpo$cpo)))
# Robust model is preferred

form_loss_log <- log_cells_lost_norm ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + Wing.Length + avg.r
model_loss_robust_log <- inla(form_loss_log, data = model_data, 
                          family = "T", 
                          lincomb = all_lc_sens,
                          control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                          control.predictor=list(compute=TRUE))

# hist(model_loss_robust_log$cpo$pit, main="", breaks = 10, xlab = "PIT")
# title("Log normalized extinction probability")
# qqplot(qunif(ppoints(length(model_loss_robust_log$cpo$pit))), model_loss_robust_log$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
# qqline(model_loss_robust_log$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))

form_col_log <- log_cells_colonised_norm ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + Wing.Length + avg.r
model_col_robust_log <- inla(form_col_log, data = model_data, 
                         family = "T",  
                         lincomb = all_lc_sens,
                         control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                         control.predictor=list(compute=TRUE))

# hist(model_col_robust_log$cpo$pit, main="", breaks = 10, xlab = "PIT")
# qqplot(qunif(ppoints(length(model_col_robust_log$cpo$pit))), model_col_robust_log$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
# qqline(model_col_robust_log$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))


form_prop_change_breadth <- log_prop_change ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      rain_breadth + temp_breadth + dry_breadth + HFP_breadth + BG_breadth + Wing.Length + avg.r
model_prop_change_breadth_robust <- inla(form_prop_change_breadth, data = model_data, 
                                 family = "T", 
                                 lincomb = all_lc_niche,
                                 control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                 control.predictor=list(compute=TRUE))


form_col_log_breadth <- log_cells_colonised_norm ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      rain_breadth + temp_breadth + dry_breadth + HFP_breadth + BG_breadth + Wing.Length + avg.r
model_col_log_breadth_robust <- inla(form_col_log_breadth, data = model_data, 
                                 family = "T",  
                                 lincomb = all_lc_niche,
                                 control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                 control.predictor=list(compute=TRUE))

form_loss_log_breadth <- log_cells_lost_norm ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      rain_breadth + temp_breadth + dry_breadth + HFP_breadth + BG_breadth + Wing.Length + avg.r
model_loss_log_breadth_robust <- inla(form_loss_log_breadth, data = model_data, 
                                      family = "T",  
                                      lincomb = all_lc_niche,
                                      control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                                      control.predictor=list(compute=TRUE))

# -----------------------------------------------------------------------------------------------------------------
# Model output - forest plot
# -----------------------------------------------------------------------------------------------------------------
prop_change_forest <- make_forest_plot(model_prop_change_robust, "Sensitivity", "#0072B2") + 
      ggtitle("Total range change") +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank())
loss_forest <- make_forest_plot(model_loss_robust_log, "Sensitivity", "#0072B2") + 
      ggtitle("Extinction") +
      ylab(element_blank()) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank())
col_forest <- make_forest_plot(model_col_robust_log, "Sensitivity", "#0072B2")  + 
      ggtitle("Colonisation") +
      ylab(element_blank()) +
      theme(axis.title.x=element_blank(),
            axis.text.x=element_blank())  

prop_change_forest_niche <- make_forest_plot(model_prop_change_breadth_robust, "Specialisation", "#E69F00") + ggtitle("Total range change") 
loss_forest_niche <- make_forest_plot(model_loss_log_breadth_robust, "Specialisation", "#E69F00") + ggtitle("Extinction") +
      ylab(element_blank()) 
col_forest_niche <- make_forest_plot(model_col_log_breadth_robust, "Specialisation", "#E69F00") + ggtitle("Colonisation") +
      ylab(element_blank())  

forest_patch <- prop_change_forest + loss_forest + col_forest +
      prop_change_forest_niche + loss_forest_niche + col_forest_niche
forest_annotated <- forest_patch + plot_annotation(tag_levels = 'A')

ggsave(filename = "figures/final_figures/forest_plots_comb_3.png", plot = forest_annotated, width = 22, height = 20,
       units = "cm", dpi = 400)

# -----------------------------------------------------------------------------------------------------------------
# Supplementary table of species data
# -----------------------------------------------------------------------------------------------------------------
species_table <- model_data %>% 
      dplyr::select(species, log_prop_change, cells_lost_norm, cells_colonised_norm,
                    Trophic.Level, Primary.Lifestyle, Migratory_ability, 
                    Wing.Length, avg.r,
                    BG_imp, rain_imp, temp_imp, dry_imp, HFP_imp,
                    BG_breadth, rain_breadth, temp_breadth, dry_breadth, HFP_breadth) %>% 
      mutate(log_cells_lost_norm = log(cells_lost_norm),
             log_cells_colonised_norm = log(cells_colonised_norm)) %>% 
      dplyr::select(species, log_prop_change, log_cells_lost_norm, log_cells_colonised_norm, everything(), -cells_lost_norm, -cells_colonised_norm) 

options(ztable.type="html")
z = ztable(species_table) %>% makeHeatmap(margin = 2)

# -----------------------------------------------------------------------------------------------------------------
# Supplementary figures: Total range change
# -----------------------------------------------------------------------------------------------------------------
list_change_sens <- make_raw_data_plots(model_prop_change_robust, model_data, "Score")
plot_change_sens <- wrap_plots(list_change_sens, ncol = 3) + plot_annotation(tag_levels = 'A', title = "Total range change - Sensitivity") 
png(filename = "figures/final_figures/plot_change_sens.png", width = 24, height = 22,
       units = "cm", res = 400)
plot_change_sens
dev.off()

list_loss_sens <- make_raw_data_plots(model_loss_robust_log, model_data, "Score")
plot_loss_sens <- wrap_plots(list_loss_sens) + plot_annotation(tag_levels = 'A', title = "Extinction - Sensitivity") 
png(filename = "figures/final_figures/plot_loss_sens.png", width = 24, height = 22,
    units = "cm", res = 400)
plot_loss_sens
dev.off()

list_col_sens <- make_raw_data_plots(model_col_robust_log, model_data, "Score")
plot_col_sens <- wrap_plots(list_col_sens) + plot_annotation(tag_levels = 'A', title = "Colonisation - Sensitivity") 
png(filename = "figures/final_figures/plot_col_sens.png", width = 24, height = 22,
    units = "cm", res = 400)
plot_col_sens
dev.off()

list_change_spec <- make_raw_data_plots(model_prop_change_breadth_robust, model_data, "Score")
plot_change_spec <- wrap_plots(list_change_spec) + plot_annotation(tag_levels = 'A', title = "Total range change - Specialisation") 
png(filename = "figures/final_figures/plot_change_spec.png", width = 24, height = 22,
    units = "cm", res = 400)
plot_change_spec
dev.off()

list_loss_spec <- make_raw_data_plots(model_loss_log_breadth_robust, model_data, "Score")
plot_loss_spec <- wrap_plots(list_loss_spec) + plot_annotation(tag_levels = 'A', title = "Extinction - Specialisation") 
png(filename = "figures/final_figures/plot_loss_spec.png", width = 24, height = 22,
    units = "cm", res = 400)
plot_loss_spec
dev.off()

list_col_spec <- make_raw_data_plots(model_col_log_breadth_robust, model_data, "Score")
plot_col_spec <- wrap_plots(list_col_spec) + plot_annotation(tag_levels = 'A', title = "Colonisation - Specialisation") 
png(filename = "figures/final_figures/plot_col_spec.png", width = 24, height = 22,
    units = "cm", res = 400)
plot_col_spec
dev.off()





