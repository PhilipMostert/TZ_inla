library(INLA)
library(INLAutils)
library(brinla)
library(tidyverse)
library(GGally)
library(raster)
library("wesanderson")
library(readxl)


# -----------------------------------------------------------------------------------------------------------------
# Data import & preparation
# -----------------------------------------------------------------------------------------------------------------
species_results <- read.csv('model_output/species_results.csv') %>% dplyr::select(-X)
traits <- readxl::read_xlsx('model_data/ELEData/TraitData/AVONET2_eBird.xlsx', sheet = "AVONET2_eBird")

species_results$species[species_results$species == "Estrilda erythronotos"] <- "Brunhilda erythronotos"
species_results$species[species_results$species == "Riparia cincta"] <- "Neophedina cincta"
species_results$species[species_results$species == "Sylvia boehmi"] <- "Curruca boehmi"
species_results$species[species_results$species == "Sylvia communis"] <- "Curruca communis"
species_results$species[species_results$species == "Turdoides aylmeri"] <- "Argya aylmeri"

# Add trait data from Tobias et al. 2022
df_merged <- merge(species_results, traits, all.x = TRUE, all.y = FALSE, by.x = "species", by.y = "Species2") %>% 
      filter(Trophic.Level != "Scavenger") %>% 
      mutate(Trophic.Level = factor(Trophic.Level, levels = c("Omnivore", "Carnivore", "Herbivore")), 
             Primary.Lifestyle = as.factor(Primary.Lifestyle),
             intercept = 1)

wing_quantiles <- quantile(df_merged$Wing.Length, probs = c(0.25, 0.75))

df_merged$wing_length_quantile <- NA
df_merged$wing_length_quantile[df_merged$Wing.Length < wing_quantiles[[1]] ] <- "small"
df_merged$wing_length_quantile[df_merged$Wing.Length >= wing_quantiles[[1]] & df_merged$Wing.Length <= wing_quantiles[[2]]] <- "moderate"
df_merged$wing_length_quantile[df_merged$Wing.Length > wing_quantiles[[2]]] <- "large"
df_merged$wing_length_quantile <- factor(df_merged$wing_length_quantile, levels = c("small", "moderate", "large"))

median_range_change <- median(df_merged$abs_range_change, na.rm = TRUE)
df_merged$range_change_group <- NA
df_merged$range_change_group[df_merged$abs_range_change < 0 ] <- "negative"
df_merged$range_change_group[df_merged$abs_range_change >= 0 & df_merged$abs_range_change <= median_range_change] <- "moderate positive"
df_merged$range_change_group[df_merged$abs_range_change > median_range_change] <- "large positive"
df_merged$range_change_group <- factor(df_merged$range_change_group, levels = c("negative", "moderate positive", "large positive"))

df_merged <- df_merged %>% 
      mutate(Wing.Length.s = scale(Wing.Length))
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

# -----------------------------------------------------------------------------------------------------------------
# Overview of range changes
# -----------------------------------------------------------------------------------------------------------------
coeff = 500

range_change_overview <- ggplot(species_results, aes(x = reorder(species, log_prop_change))) +
      geom_bar(aes(y = abs_range_change/coeff), stat = 'identity', fill = "#005AB5", alpha = 0.5) +
      geom_bar(aes(y = log_prop_change), stat = 'identity', fill = "#DC3220", alpha = 0.5) +
      geom_hline(aes(yintercept = 0), alpha = .5) +
      theme_classic() +
      xlab(element_blank()) +
      scale_x_discrete(labels = element_blank()) +
      scale_y_continuous(
            name = "Log proportional range change",
            sec.axis = sec_axis(~.*coeff, name="Absolute range change")) +
      theme(axis.line.y.right = element_line(color = "#005AB5"),
            axis.text.y.right = element_text(color = "#005AB5"),
            axis.ticks.y.right = element_line(color = "#005AB5"),
            axis.title.y.right = element_text(color = "#005AB5"),
            axis.line.y.left = element_line(color = "#DC3220"),
            axis.text.y.left = element_text(color = "#DC3220"),
            axis.ticks.y.left = element_line(color = "#DC3220"),
            axis.title.y.left = element_text(color = "#DC3220")); range_change_overview


species_results$cells_change <- species_results$cells_colonised_norm - species_results$cells_lost_norm

cells_change <- ggplot(species_results, aes(x = reorder(species, log_prop_change))) +
      geom_point(aes(y = cells_change, col = "Difference"), size = 2, alpha = .5) +
      geom_point(aes(y = cells_colonised_norm, col = "Cells colonised"), pch = '-', size = 7) +
      geom_point(aes(y = cells_lost_norm, col = "Cells lost"), pch = '-', size = 7) +
      geom_segment(aes(xend = reorder(species, cells_colonised_norm), y = cells_lost_norm, yend = cells_colonised_norm), alpha = .3, lty = 2, lwd = .3) +
      geom_hline(aes(yintercept = 0), alpha = .5) +
      theme_classic() +
      xlab(element_blank()) +
      ylab("N(cells) normalized") +
      theme(legend.title = element_blank()) +
      scale_x_discrete(labels = element_blank()); cells_change

cells_change_log <- ggplot(species_results, aes(x = reorder(species, log_prop_change))) +
      geom_point(aes(y = log(cells_change+1.2), col = "Difference"), size = 2, alpha = .5) +
      geom_point(aes(y = log(cells_colonised_norm+1.2), col = "Cells colonised"), pch = '-', size = 7) +
      geom_point(aes(y = log(cells_lost_norm+1.2), col = "Cells lost"), pch = '-', size = 7) +
      geom_segment(aes(xend = reorder(species,  log(cells_colonised_norm+1.2)), y = log(cells_lost_norm+1.2), yend =  log(cells_colonised_norm+1.2)), alpha = .3, lty = 2, lwd = .3) +
      geom_hline(aes(yintercept = log(1.2)), alpha = .5) +
      theme_classic() +
      xlab("Species") +
      ylab("log(N(cells) normalized + 1.2)") +
      theme(legend.title = element_blank()) +
      scale_x_discrete(labels = element_blank()); cells_change_log

overview_combined <- ggpubr::ggarrange(range_change_overview, cells_change_log,
                                       nrow = 2, ncol = 1, common.legend = T, legend = "bottom", align = "v"); overview_combined
ggsave(filename = "figures/final_figures/overview_combined.png", plot = overview_combined, width = 18, height = 20,
       units = "cm", dpi = 300)
# -----------------------------------------------------------------------------------------------------------------
# Statistical test
# -----------------------------------------------------------------------------------------------------------------
form <- log_prop_change ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + Wing.Length.s

form_loss <- cells_lost_norm ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + Wing.Length.s

form_col <- cells_colonised_norm ~ Trophic.Level + Primary.Lifestyle + Migratory_ability + 
      BG_imp + rain_imp + temp_imp + dry_imp + HFP_imp + Wing.Length.s

# Trophic.Level_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
#                                        Trophic.LevelCarnivore = c(0, 1, 0),
#                                        Trophic.LevelHerbivore = c(0, 0, 1))

BG_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                            BG_imp = seq(min(df_merged$BG_imp), max(df_merged$BG_imp), len = 100))
names(BG_lc) <- paste0("BG", 1:100)
rain_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              rain_imp = seq(min(df_merged$rain_imp), max(df_merged$rain_imp), len = 100))
names(rain_lc) <- paste0("rain", 1:100)
temp_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              temp_imp = seq(min(df_merged$temp_imp), max(df_merged$temp_imp), len = 100))
names(temp_lc) <- paste0("temp", 1:100)
dry_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             dry_imp = seq(min(df_merged$dry_imp), max(df_merged$dry_imp), len = 100))
names(dry_lc) <- paste0("dry", 1:100)
HFP_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             HFP_imp = seq(min(df_merged$HFP_imp), max(df_merged$HFP_imp), len = 100))
names(HFP_lc) <- paste0("HFP", 1:100)
wing_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              Wing.Length.s = seq(min(df_merged$Wing.Length.s), max(df_merged$Wing.Length.s), len = 100))
names(wing_lc) <- paste0("wing", 1:100)
all_lc <- c(BG_lc, rain_lc, temp_lc, dry_lc, HFP_lc, wing_lc)

hist(df_merged$log_prop_change)

# Wang et al. 2018: Experience suggests, that unless we specify a sharp prior that it is in contradiction with the data, 
# this will not make much difference to the posterior distribution. One is usually satisfied with the default flat prior 
# for the intercept in such models.

model_gauss <- inla(form, data = df_merged, 
                    family = "gaussian",  # t distributed errors regression model
                    lincomb = all_lc,
                    control.compute = list(dic = TRUE, cpo = TRUE), 
                    control.predictor=list(compute=TRUE))
bri.lmresid.plot(model_gauss)  
hist(model_gauss$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_gauss$cpo$pit))), 
       model_gauss$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_gauss$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
# Outliers in residuals, and skewed Q-Q plot. Try robust model, reducing effect of outliers.

model_robust <- inla(form, data = df_merged, 
                     family = "T",  # t distributed errors regression model
                     lincomb = all_lc,
                     control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                     control.predictor=list(compute=TRUE))
ggplot_inla_residuals(model_robust, df_merged$log_prop_change, binwidth = 0.1)

hist(model_robust$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_robust$cpo$pit))), model_robust$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_robust$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
# Q-Q plot looks much better. Compare the two models using CPO and PIT

sum(model_gauss$cpo$failure)
sum(model_robust$cpo$failure)
c(-sum(log(model_gauss$cpo$cpo)), -sum(log(model_robust$cpo$cpo)))
c(-sum(log(model_gauss$cpo$pit)), -sum(log(model_robust$cpo$cpo)))
# Robust model is preferred

summary(model_robust)

lincomb_out <- data.frame(model_robust$summary.lincomb.derived) %>% 
      mutate(ID = gsub('[[:digit:]]+', '', rownames(.))) %>% 
      dplyr::select(-mode, -kld, -mean, -sd) %>% 
      pivot_wider(names_from = ID, values_from = c(2:4),
                  names_glue = "{ID}_{.value}") %>% 
      unnest() %>% 
      mutate(
            "BG" = seq(min(df_merged$BG_imp), max(df_merged$BG_imp), len = 100),
            "rain" = seq(min(df_merged$rain_imp), max(df_merged$rain_imp), len = 100), 
            "temp" = seq(min(df_merged$temp_imp), max(df_merged$temp_imp), len = 100), 
            "dry" = seq(min(df_merged$dry_imp), max(df_merged$dry_imp), len = 100), 
            "HFP" = seq(min(df_merged$HFP_imp), max(df_merged$HFP_imp), len = 100),
            "wing" = seq(min(df_merged$Wing.Length), max(df_merged$Wing.Length), len = 100)
      )


model_loss_robust <- inla(form_loss, data = df_merged, 
                          family = "T",  # t distributed errors regression model
                          lincomb = all_lc,
                          control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                          control.predictor=list(compute=TRUE))

hist(model_loss_robust$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_loss_robust$cpo$pit))), model_loss_robust$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_loss_robust$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
summary(model_loss_robust)

model_col_robust <- inla(form_col, data = df_merged, 
                         family = "T",  # t distributed errors regression model
                         lincomb = all_lc,
                         control.compute = list(dic = TRUE, cpo = TRUE, return.marginals.predictor=TRUE), 
                         control.predictor=list(compute=TRUE))

hist(model_col_robust$cpo$pit, main="", breaks = 10, xlab = "PIT")
qqplot(qunif(ppoints(length(model_col_robust$cpo$pit))), model_col_robust$cpo$pit, main = "Q-Q plot for Unif(0,1)", xlab = "Theoretical Quantiles", ylab = "Sample Quantiles")
qqline(model_col_robust$cpo$pit, distribution = function(p) qunif(p), prob = c(0.1, 0.9))
summary(model_col_robust)





# -----------------------------------------------------------------------------------------------------------------
# Model output - forest plot
# -----------------------------------------------------------------------------------------------------------------
fixed_out <- data.frame(model_col_robust$summary.fixed) %>% 
      mutate(ID = c("Intercept", "Trophic level: Carnivore", "Trophic level: Herbivore", 
                    "Primary lifestyle: Generalist", "Primary lifestyle: Insessorial", "Primary lifestyle: Terrestrial",
                    "Migratory ability: Moderate", "Migratory ability: High",
                    "Bare ground cover", "Annual rainfall", "Hottest temperature", "Dryspell duration", "Human footprint", "Wing length")
      )
fixed_out$ID <- factor(fixed_out$ID, levels = c(fixed_out$ID))

forest_plot <- ggplot() +
      geom_point(data = fixed_out %>% filter(ID != "Intercept"), aes(y = ID, x = X0.5quant)) +
      geom_errorbar(data = fixed_out %>% filter(ID != "Intercept"), aes(y = ID, xmin = X0.025quant, xmax = X0.975quant), width = 0.1) +
      geom_vline(aes(xintercept = 0), lty = 2) +
      geom_vline(data = fixed_out %>% filter(ID == "Intercept"), aes(xintercept =  X0.5quant), col = "orange") +
      xlab("Posterior estimates") +
      ylab(element_blank()) +
      theme_classic() +
      scale_y_discrete(limits=rev)


ggsave(filename = "figures/final_figures/forest_plot.png", plot = forest_plot, width = 18, height = 20,
       units = "cm", dpi = 400)

# -----------------------------------------------------------------------------------------------------------------
# Overall covariate importance
# -----------------------------------------------------------------------------------------------------------------
cov_importance <- ggplot(species_long, aes(x = covariate, y = perc_explained*100)) +
      stat_boxplot(geom ='errorbar', width = 0.1) + 
      geom_boxplot(width = 0.2) +
      xlab("") +
      ylab("Relative importance [%]") +
      theme_classic() +
      scale_x_discrete(labels = c("Temperature",  "Rainfall", "Bare ground", "Dryspell", "Human footprint")); cov_importance

ggsave(filename = "figures/covariate_importance.png", plot = cov_importance, width = 15, height = 10,
       units = "cm", dpi = 400)

# -----------------------------------------------------------------------------------------------------------------
# Covariate importance by trait
# -----------------------------------------------------------------------------------------------------------------
migration_imp <- ggplot(species_long) +
      geom_boxplot(aes(x = Migratory_ability, y = perc_explained*100, fill = covariate)) +
      ylab("Relative importance [%]") +
      xlab(element_blank()) +
      theme_classic() +
      theme(legend.title=element_blank()) +
      scale_fill_manual(labels = c("Temperature",  "Rainfall", "Bare ground", "Dryspell", "Human footprint"),
                        values = c(wes_palettes$Rushmore1)) +
      ggtitle("Migratory ability") +
      scale_x_discrete(labels = c("Low", "Moderate", "High"))

wing_imp <- ggplot(species_long) +
      geom_boxplot(aes(x = wing_length_quantile, y = perc_explained*100, fill = covariate)) +
      ylab(element_blank()) +
      xlab(element_blank()) +
      theme_classic() +
      theme(legend.title=element_blank()) +
      scale_fill_manual(labels = c("Temperature",  "Rainfall", "Bare ground", "Dryspell", "Human footprint"),
                        values = c(wes_palettes$Rushmore1)) +
      ggtitle("Wing length") 

density_imp <- ggplot(species_long) +
      geom_boxplot(aes(x = as.factor(Habitat.Density), y = perc_explained*100, fill = covariate)) +
      ylab("Relative importance [%]") +
      xlab(element_blank()) +
      theme_classic() +
      theme(legend.title=element_blank()) +
      scale_fill_manual(labels = c("Temperature",  "Rainfall", "Bare ground", "Dryspell", "Human footprint"),
                        values = c(wes_palettes$Rushmore1)) +
      ggtitle("Habitat density") +
      scale_x_discrete(labels = c("Dense", "Semi-open", "Open"))

trophic_imp <- ggplot(species_long) +
      geom_boxplot(aes(x = Trophic.Level, y = perc_explained*100, fill = covariate)) +
      ylab("Relative importance [%]") +
      xlab(element_blank()) +
      theme_classic() +
      theme(legend.title=element_blank()) +
      scale_fill_manual(labels = c("Temperature",  "Rainfall", "Bare ground", "Dryspell", "Human footprint"),
                        values = c(wes_palettes$Rushmore1)) +
      ggtitle("Trophic level") 

lifestyle_imp <- ggplot(species_long) +
      geom_boxplot(aes(x = Primary.Lifestyle, y = perc_explained*100, fill = covariate)) +
      ylab(element_blank()) +
      xlab(element_blank()) +
      theme_classic() +
      theme(legend.title=element_blank()) +
      scale_fill_manual(labels = c("Temperature",  "Rainfall", "Bare ground", "Dryspell", "Human footprint"),
                        values = c(wes_palettes$Rushmore1)) +
      ggtitle("Primary lifestyle") 

importance_plots <- ggpubr::ggarrange(migration_imp, wing_imp, trophic_imp, lifestyle_imp,
                                      nrow = 2, ncol = 2, common.legend = TRUE, legend="bottom")

ggsave(filename = "figures/final_figures/importance_plots.png", plot = importance_plots, width = 18, height = 20,
       units = "cm", dpi = 400)

# -----------------------------------------------------------------------------------------------------------------
# Absolute range change by trait
# -----------------------------------------------------------------------------------------------------------------
migration <- ggplot(df_merged) +
      geom_boxplot(aes(Migratory_ability, abs_range_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ylab("Absolute range change") +
      ggtitle("Migratory ability") +
      scale_x_discrete(labels = c("Low", "Moderate", "High"))

wing_length <- ggplot(df_merged, aes(wing_length_quantile, abs_range_change)) +
      geom_boxplot() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      # ylab("Absolute range change") +
      ggtitle("Wing length") +
      theme_classic()

habitat_density <- ggplot(df_merged) +
      geom_boxplot(aes(as.factor(Habitat.Density), abs_range_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ylab("Absolute range change") +
      ggtitle("Habitat density") +
      scale_x_discrete(labels = c("Dense", "Semi-open", "Open"))

primary_lifestyle <- ggplot(df_merged) +
      geom_boxplot(aes(Primary.Lifestyle, abs_range_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Primary lifestyle") 
# ylab("Absolute range change")

trophic_level <- ggplot(df_merged) +
      geom_boxplot(aes(Trophic.Level, abs_range_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ggtitle("Trophic level") +
      ylab("Absolute range change")

trophic_niche <- ggplot(df_merged) +
      geom_boxplot(aes(Trophic.Niche, abs_range_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Trophic niche") 
# ylab("Absolute range change")

absChange_trait_plots <- ggpubr::ggarrange(migration, wing_length, trophic_level, primary_lifestyle,
                                           nrow = 2, ncol = 2)
ggsave(filename = "figures/final_figures/absChange_trait_plots.png", plot = absChange_trait_plots, width = 18, height = 20,
       units = "cm", dpi = 400)

# -----------------------------------------------------------------------------------------------------------------
# Log proportional range change by trait. 
# -----------------------------------------------------------------------------------------------------------------
# Zero is no change, same positive value for a doubling of area as the negative value for a halving of area.
prop_migration <- ggplot(df_merged) +
      geom_boxplot(aes(Migratory_ability, log_prop_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ylab("Log proportional range change") +
      ggtitle("Migratory ability") +
      scale_x_discrete(labels = c("Low", "Moderate", "High"))

# prop_wing_length <- ggplot(df_merged, aes(wing_length_quantile, log_prop_change)) +
#       geom_boxplot() +
#       xlab(element_blank()) +
#       ylab(element_blank()) +
#       # ylab("Log proportional range change") +
#       ggtitle("Wing length") +
#       theme_classic()

prop_wing_length <- ggplot() +
      geom_point(data = df_merged, aes(x = Wing.Length, y = log_prop_change)) +
      geom_line(data = lincomb_out, aes(x = wing, y = wing_X0.5quant), lty = 1.3, col = "orange") +
      geom_ribbon(data = lincomb_out, aes(x = wing, ymin = wing_X0.025quant, ymax = wing_X0.975quant), alpha = 0.2) +
      ylab(element_blank()) +
      xlab("Wing length [mm]") +
      ggtitle("Wing length") +
      theme_classic()

prop_habitat_density <- ggplot(df_merged) +
      geom_boxplot(aes(as.factor(Habitat.Density), log_prop_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ylab("Log proportional range change") +
      ggtitle("Habitat density") +
      scale_x_discrete(labels = c("Dense", "Semi-open", "Open"))

prop_primary_lifestyle <- ggplot(df_merged) +
      geom_boxplot(aes(Primary.Lifestyle, log_prop_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Primary lifestyle") 
# ylab("Log proportional range change") 

prop_trophic_level <- ggplot(df_merged %>% filter(Trophic.Level != "Scavenger")) +
      geom_boxplot(aes(Trophic.Level, log_prop_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ggtitle("Trophic level") +
      ylab("Log proportional range change") 

prop_trophic_niche <- ggplot(df_merged) +
      geom_boxplot(aes(Trophic.Niche, log_prop_change)) +
      theme_classic() +
      xlab(element_blank()) +
      ylab(element_blank()) +
      ggtitle("Trophic niche") 
# ylab("Log proportional range change") 

propChange_trait_plots <-  ggpubr::ggarrange(prop_migration, prop_wing_length, prop_trophic_level, prop_primary_lifestyle,
                                             nrow = 2, ncol = 2)
ggsave(filename = "figures/final_figures/propChange_trait_plots.png", plot = propChange_trait_plots, width = 18, height = 20,
       units = "cm", dpi = 400)

# -----------------------------------------------------------------------------------------------------------------
# Covariate importance against range change
# -----------------------------------------------------------------------------------------------------------------
BG_trend <- ggplot() +
      geom_point(data = df_merged, aes(x = BG_imp, y = log_prop_change)) +
      geom_line(data = lincomb_out, aes(x = BG, y = BG_X0.5quant), lty = 2, alpha = 0.5) +
      geom_ribbon(data = lincomb_out, aes(x = BG, ymin = BG_X0.025quant, ymax = BG_X0.975quant), alpha = 0.2) +
      ylab("Log proportional range change") +
      xlab("Relative importance") +
      ggtitle("Bare ground sensitivity") +
      theme_classic()

rain_trend <- ggplot() +
      geom_point(data = df_merged, aes(x = rain_imp, y = log_prop_change)) +
      geom_line(data = lincomb_out, aes(x = rain, y = rain_X0.5quant), lty = 2, alpha = 0.5) +
      geom_ribbon(data = lincomb_out, aes(x = rain, ymin = rain_X0.025quant, ymax = rain_X0.975quant), alpha = 0.3) +
      ylab(element_blank()) +
      xlab("Relative importance") +
      ggtitle("Annual rainfall sensitivity") +
      theme_classic()

temp_trend <- ggplot() +
      geom_point(data = df_merged, aes(x = temp_imp, y = log_prop_change)) +
      geom_line(data = lincomb_out, aes(x = temp, y = temp_X0.5quant), lty = 2, alpha = 0.5) +
      geom_ribbon(data = lincomb_out, aes(x = temp, ymin = temp_X0.025quant, ymax = temp_X0.975quant), alpha = 0.2) +
      ylab("Log proportional range change") +
      xlab("Relative importance") +
      ggtitle("Hottest temperature sensitivity") +
      theme_classic()

dry_trend <- ggplot() +
      geom_point(data = df_merged, aes(x = dry_imp, y = log_prop_change)) +
      geom_line(data = lincomb_out, aes(x = dry, y = dry_X0.5quant), lty = 2, alpha = 0.5) +
      geom_ribbon(data = lincomb_out, aes(x = dry, ymin = dry_X0.025quant, ymax = dry_X0.975quant), alpha = 0.2) +
      ylab(element_blank()) +
      xlab("Relative importance") +
      ggtitle("Dryspell duration sensitivity") +
      theme_classic()

HFP_trend <- ggplot() +
      geom_point(data = df_merged, aes(x = HFP_imp, y = log_prop_change)) +
      geom_line(data = lincomb_out, aes(x = HFP, y = HFP_X0.5quant), lty = 2, alpha = 0.5) +
      geom_ribbon(data = lincomb_out, aes(x = HFP, ymin = HFP_X0.025quant, ymax = HFP_X0.975quant), alpha = 0.2) +
      ylab("Log proportional range change") +
      xlab("Relative importance") +
      ggtitle("Human footprint sensitivity") +
      theme_classic()

trend_plots <-  ggpubr::ggarrange(BG_trend, rain_trend, temp_trend, dry_trend, HFP_trend, 
                                  nrow = 3, ncol = 2)
ggsave(filename = "figures/final_figures/trend_plots.png", plot = trend_plots, width = 18, height = 20,
       units = "cm", dpi = 400)

coeff = 500

range_change_overview <- ggplot(species_results, aes(x = reorder(species, log_prop_change))) +
      geom_bar(aes(y = abs_range_change/coeff), stat = 'identity', fill = "#005AB5", alpha = 0.5) +
      geom_bar(aes(y = log_prop_change), stat = 'identity', fill = "#DC3220", alpha = 0.5) +
      geom_hline(aes(yintercept = 0), alpha = .5)
theme_classic() +
      xlab(element_blank()) +
      scale_x_discrete(labels = element_blank()) +
      scale_y_continuous(
            name = "Log proportional range change",
            sec.axis = sec_axis(~.*coeff, name="Absolute range change")) +
      theme(axis.line.y.right = element_line(color = "#005AB5"),
            axis.text.y.right = element_text(color = "#005AB5"),
            axis.ticks.y.right = element_line(color = "#005AB5"),
            axis.title.y.right = element_text(color = "#005AB5"),
            axis.line.y.left = element_line(color = "#DC3220"),
            axis.text.y.left = element_text(color = "#DC3220"),
            axis.ticks.y.left = element_line(color = "#DC3220"),
            axis.title.y.left = element_text(color = "#DC3220")); range_change_overview


species_results$cells_change <- species_results$cells_colonised_norm - species_results$cells_lost_norm

wing_quantiles <- quantile(df_merged$Wing.Length, probs = c(0.25, 0.75))

df_merged$wing_length_quantile <- NA
df_merged$wing_length_quantile[df_merged$Wing.Length < wing_quantiles[[1]] ] <- "small"
df_merged$wing_length_quantile[df_merged$Wing.Length >= wing_quantiles[[1]] & df_merged$Wing.Length <= wing_quantiles[[2]]] <- "moderate"
df_merged$wing_length_quantile[df_merged$Wing.Length > wing_quantiles[[2]]] <- "large"
df_merged$wing_length_quantile <- factor(df_merged$wing_length_quantile, levels = c("small", "moderate", "large"))

median_range_change <- median(df_merged$abs_range_change, na.rm = TRUE)
df_merged$range_change_group <- NA
df_merged$range_change_group[df_merged$abs_range_change < 0 ] <- "negative"
df_merged$range_change_group[df_merged$abs_range_change >= 0 & df_merged$abs_range_change <= median_range_change] <- "moderate positive"
df_merged$range_change_group[df_merged$abs_range_change > median_range_change] <- "large positive"
df_merged$range_change_group <- factor(df_merged$range_change_group, levels = c("negative", "moderate positive", "large positive"))


# forest_plots_comb_1 <- ggarrange(prop_change_forest, loss_forest, col_forest,
#                                ncol = 3, nrow = 1, align = "v", 
#                                labels = c("A", "B", "C"))
# forest_plots_comb_2 <- ggarrange(prop_change_forest_niche, loss_forest_niche, col_forest_niche,
#                                  ncol = 3, nrow = 1, align = "v", 
#                                  labels = c("D", "E", "F"))
#                                  
#                                  
# ggsave(filename = "figures/final_figures/forest_plots_comb_1.png", plot = forest_plots_comb_1, width = 22, height = 12.5,
#        units = "cm", dpi = 400)
# ggsave(filename = "figures/final_figures/forest_plots_comb_2.png", plot = forest_plots_comb_2, width = 22, height = 12.5,
#        units = "cm", dpi = 400)


species_long <- df_merged %>% 
      mutate_if(is.numeric, round, digits = 3) %>% 
      pivot_longer(cols = c("rain_imp", "temp_imp", "dry_imp", "BG_imp", "HFP_imp" ), names_to = 'covariate', values_to = 'perc_explained') %>% 
      drop_na() %>% 
      mutate(covariate = factor(covariate, levels = c("temp_imp", "rain_imp", "BG_imp", "dry_imp", "HFP_imp")))


Migratory_ability_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
                                           Migratory_abilitylow = c(0, 1, 0),
                                           Migratory_abilitymoderate = c(0, 0, 1))



LC_no_lifestyle_sens <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                           Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                           Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                           Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                           Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                           BG_imp = model_data$BG_imp[!is.na(model_data$avg.r)],
                                           rain_imp = model_data$rain_imp[!is.na(model_data$avg.r)],
                                           temp_imp = model_data$temp_imp[!is.na(model_data$avg.r)],
                                           dry_imp = model_data$dry_imp[!is.na(model_data$avg.r)],
                                           HFP_imp = model_data$HFP_imp[!is.na(model_data$avg.r)],
                                           Wing.Length = model_data$Wing.Length[!is.na(model_data$avg.r)],
                                           avg.r = model_data$avg.r[!is.na(model_data$avg.r)])
names(LC_no_lifestyle_sens) <- paste0("LC_no_lifestyle_sens", 1:NROW(model_data))

LC_no_trophic_sens <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                         Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                         Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                         Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                         Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                         Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                         BG_imp = model_data$BG_imp,
                                         rain_imp = model_data$rain_imp,
                                         temp_imp = model_data$temp_imp,
                                         dry_imp = model_data$dry_imp,
                                         HFP_imp = model_data$HFP_imp,
                                         Wing.Length = model_data$Wing.Length)
names(LC_no_trophic_sens) <- paste0("LC_no_trophic_sens", 1:NROW(model_data))

LC_no_migratory_sens <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                           Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                           Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                           Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                           Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                           Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                           BG_imp = model_data$BG_imp,
                                           rain_imp = model_data$rain_imp,
                                           temp_imp = model_data$temp_imp,
                                           dry_imp = model_data$dry_imp,
                                           HFP_imp = model_data$HFP_imp,
                                           Wing.Length = model_data$Wing.Length)
names(LC_no_migratory_sens) <- paste0("LC_no_migratory_sens", 1:NROW(model_data))

LC_no_BG_imp <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                   Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                   Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                   Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                   Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                   Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                   Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                   Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                   rain_imp = model_data$rain_imp,
                                   temp_imp = model_data$temp_imp,
                                   dry_imp = model_data$dry_imp,
                                   HFP_imp = model_data$HFP_imp,
                                   Wing.Length = model_data$Wing.Length)
names(LC_no_BG_imp) <- paste0("LC_no_BG_imp", 1:NROW(model_data))

LC_no_rain_imp <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                     Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                     Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                     Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                     Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                     Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                     Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                     Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                     BG_imp = model_data$BG_imp,
                                     temp_imp = model_data$temp_imp,
                                     dry_imp = model_data$dry_imp,
                                     HFP_imp = model_data$HFP_imp,
                                     Wing.Length = model_data$Wing.Length)
names(LC_no_rain_imp) <- paste0("LC_no_rain_imp", 1:NROW(model_data))

LC_no_temp_imp <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                     Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                     Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                     Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                     Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                     Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                     Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                     Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                     BG_imp = model_data$BG_imp,
                                     rain_imp = model_data$rain_imp,
                                     dry_imp = model_data$dry_imp,
                                     HFP_imp = model_data$HFP_imp,
                                     Wing.Length = model_data$Wing.Length)
names(LC_no_temp_imp) <- paste0("LC_no_temp_imp", 1:NROW(model_data))

LC_no_dry_imp <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                    Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                    Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                    Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                    Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                    Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                    Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                    Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                    BG_imp = model_data$BG_imp,
                                    rain_imp = model_data$rain_imp,
                                    temp_imp = model_data$temp_imp,
                                    HFP_imp = model_data$HFP_imp,
                                    Wing.Length = model_data$Wing.Length)
names(LC_no_dry_imp) <- paste0("LC_no_dry_imp", 1:NROW(model_data))

LC_no_HFP_imp <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                    Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                    Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                    Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                    Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                    Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                    Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                    Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                    BG_imp = model_data$BG_imp,
                                    rain_imp = model_data$rain_imp,
                                    temp_imp = model_data$temp_imp,
                                    dry_imp = model_data$dry_imp,
                                    Wing.Length = model_data$Wing.Length)
names(LC_no_HFP_imp) <- paste0("LC_no_HFP_imp", 1:NROW(model_data))

LC_no_Wing_imp <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                     Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                     Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                     Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                     Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                     Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                     Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                     Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                     BG_imp = model_data$BG_imp,
                                     rain_imp = model_data$rain_imp,
                                     temp_imp = model_data$temp_imp,
                                     dry_imp = model_data$dry_imp,
                                     HFP_imp = model_data$HFP_imp)
names(LC_no_Wing_imp) <- paste0("LC_no_Wing_imp", 1:NROW(model_data))

# LC_no_avg.r_imp <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
#                                    Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
#                                    Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
#                                    Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
#                                    Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
#                                    Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
#                                    Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
#                                    Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
#                                    BG_imp = model_data$BG_imp,
#                                    rain_imp = model_data$rain_imp,
#                                    temp_imp = model_data$temp_imp,
#                                    dry_imp = model_data$dry_imp,
#                                    HFP_imp = model_data$HFP_imp,
#                                    Wing.Length = model_data$Wing.Length)
# names(LC_no_avg.r_imp) <- paste0("LC_no_avg.r_imp", 1:NROW(model_data))


LC_no_lifestyle_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                              Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                              Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                              Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                              Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                              BG_breadth = model_data$BG_breadth,
                                              rain_breadth = model_data$rain_breadth,
                                              temp_breadth = model_data$temp_breadth,
                                              dry_breadth = model_data$dry_breadth,
                                              HFP_breadth = model_data$HFP_breadth,
                                              Wing.Length = model_data$Wing.Length)
names(LC_no_lifestyle_breadth) <- paste0("LC_no_lifestyle_breadth", 1:NROW(model_data))

LC_no_trophic_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                            Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                            Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                            Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                            Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                            Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                            BG_breadth = model_data$BG_breadth,
                                            rain_breadth = model_data$rain_breadth,
                                            temp_breadth = model_data$temp_breadth,
                                            dry_breadth = model_data$dry_breadth,
                                            HFP_breadth = model_data$HFP_breadth,
                                            Wing.Length = model_data$Wing.Length)
names(LC_no_trophic_breadth) <- paste0("LC_no_trophic_breadth", 1:NROW(model_data))

LC_no_migratory_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                              Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                              Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                              Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                              Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                              Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                              BG_breadth = model_data$BG_breadth,
                                              rain_breadth = model_data$rain_breadth,
                                              temp_breadth = model_data$temp_breadth,
                                              dry_breadth = model_data$dry_breadth,
                                              HFP_breadth = model_data$HFP_breadth,
                                              Wing.Length = model_data$Wing.Length)
names(LC_no_migratory_breadth) <- paste0("LC_no_migratory_breadth", 1:NROW(model_data))

LC_no_BG_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                       Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                       Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                       Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                       Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                       Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                       Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                       Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                       rain_breadth = model_data$rain_breadth,
                                       temp_breadth = model_data$temp_breadth,
                                       dry_breadth = model_data$dry_breadth,
                                       HFP_breadth = model_data$HFP_breadth,
                                       Wing.Length = model_data$Wing.Length)
names(LC_no_BG_breadth) <- paste0("LC_no_BG_breadth", 1:NROW(model_data))

LC_no_rain_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                         Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                         Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                         Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                         Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                         Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                         Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                         Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                         BG_breadth = model_data$BG_breadth,
                                         temp_breadth = model_data$temp_breadth,
                                         dry_breadth = model_data$dry_breadth,
                                         HFP_breadth = model_data$HFP_breadth,
                                         Wing.Length = model_data$Wing.Length)
names(LC_no_rain_breadth) <- paste0("LC_no_rain_breadth", 1:NROW(model_data))

LC_no_temp_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                         Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                         Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                         Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                         Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                         Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                         Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                         Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                         BG_breadth = model_data$BG_breadth,
                                         rain_breadth = model_data$rain_breadth,
                                         dry_breadth = model_data$dry_breadth,
                                         HFP_breadth = model_data$HFP_breadth,
                                         Wing.Length = model_data$Wing.Length)
names(LC_no_temp_breadth) <- paste0("LC_no_temp_breadth", 1:NROW(model_data))

LC_no_dry_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                        Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                        Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                        Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                        Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                        Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                        Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                        Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                        BG_breadth = model_data$BG_breadth,
                                        rain_breadth = model_data$rain_breadth,
                                        temp_breadth = model_data$temp_breadth,
                                        HFP_breadth = model_data$HFP_breadth,
                                        Wing.Length = model_data$Wing.Length)
names(LC_no_dry_breadth) <- paste0("LC_no_dry_breadth", 1:NROW(model_data))

LC_no_HFP_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                        Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                        Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                        Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                        Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                        Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                        Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                        Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                        BG_breadth = model_data$BG_breadth,
                                        rain_breadth = model_data$rain_breadth,
                                        temp_breadth = model_data$temp_breadth,
                                        dry_breadth = model_data$dry_breadth,
                                        Wing.Length = model_data$Wing.Length)
names(LC_no_HFP_breadth) <- paste0("LC_no_HFP_breadth", 1:NROW(model_data))

LC_no_Wing_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                         Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                         Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                         Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                         Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                         Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                         Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                         Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                         BG_breadth = model_data$BG_breadth,
                                         rain_breadth = model_data$rain_breadth,
                                         temp_breadth = model_data$temp_breadth,
                                         dry_breadth = model_data$dry_breadth,
                                         HFP_breadth = model_data$HFP_breadth)
names(LC_no_Wing_breadth) <- paste0("LC_no_Wing_breadth", 1:NROW(model_data))

# LC_no_avg.r_breadth <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
#                                       Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
#                                       Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
#                                       Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
#                                       Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
#                                       Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
#                                       Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
#                                       Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
#                                       BG_breadth = model_data$BG_breadth,
#                                       rain_breadth = model_data$rain_breadth,
#                                       temp_breadth = model_data$temp_breadth,
#                                       dry_breadth = model_data$dry_breadth,
#                                       HFP_breadth = model_data$HFP_breadth,
#                                       Wing.Length = model_data$Wing.Length)
# names(LC_no_avg.r_breadth) <- paste0("LC_no_avg.r_breadth", 1:NROW(model_data))
# 
# 
# lincomb_data1 <- data.frame(model_col_robust_log$summary.lincomb.derived) %>% 
#       mutate(ID = rownames(.)) %>% 
#       dplyr::select(-mode, -kld, -mean, -sd)
# lincomb_data1 <- lincomb_data1[grepl("Wing", rownames(lincomb_data1)) == TRUE, ]
# 
# lincomb_data2 <- data.frame(model_prop_change_robust$summary.lincomb.derived) %>% 
#       mutate(ID = rownames(.)) %>% 
#       dplyr::select(-mode, -kld, -mean, -sd)
# lincomb_data2 <- lincomb_data2[grepl("Wing", rownames(lincomb_data2)) == TRUE, ]
# 
# wingPlot <- ggplot(model_data) +
#       geom_point(aes(Wing.Length, relative_colonisation, col = "Relative colonisation"), alpha = .2) +
#       geom_point(aes(Wing.Length, log_prop_change, col = "Total range change"), alpha = .2) +
#       geom_ribbon(data = lincomb_data1, 
#                  aes(x = seq(min(model_data$Wing.Length), max(model_data$Wing.Length), len = 100), 
#                      ymin = X0.025quant, 
#                      ymax = X0.975quant, 
#                      group = "Relative colonisation"), alpha = .1) +
#       geom_ribbon(data = lincomb_data2, 
#                  aes(x = seq(min(model_data$Wing.Length), max(model_data$Wing.Length), len = 100),
#                      ymin = X0.025quant, 
#                      ymax = X0.975quant,
#                      group = "Total range change"), alpha = .1) +
#       geom_line(data = lincomb_data1, 
#                 aes(x = seq(min(model_data$Wing.Length), max(model_data$Wing.Length), len = 100), 
#                     y = X0.5quant, col = "Relative colonisation"), lty = 2) +
#       geom_line(data = lincomb_data2, 
#                 aes(x = seq(min(model_data$Wing.Length), max(model_data$Wing.Length), len = 100),
#                     y = X0.5quant, col = "Total range change"), lty = 2) +
#       ggthemes::theme_calc() +
#       xlab("Wing length") +
#       ylab("Score") +
#       theme(legend.position = c(.5, .9),
#             legend.title = element_blank()) +
#       scale_colour_manual(values = c("#000000", "#56B4E9")) +
#       relative_scale

# migrationPlot <- make_raw_data_plots(model_loss_log_breadth_robust, model_data, "Extinction score")[[3]]
# migrationPlot <- migrationPlot + 
#       xlab("Migratory ability")

