library(INLA)
library(tidyverse)
library(raster)
library("wesanderson")

modelData_list <- list.files(path = "/Users/joriswiethase/Downloads/files_to_move", pattern = ".RData",
                      full.names = TRUE, recursive = TRUE)

species_results <- data.frame(matrix(NA, nrow = length(modelData_list), ncol = 7))
colnames(species_results) <- c('species', 
                               'rain_imp', 'temp_imp', 'dry_imp', 'BG_imp', 'HFP_imp',
                               'range_change')

add_quantiles <- function(df, column, new_var){
      med_value <- median(df[, column], na.rm = TRUE)
      df[, new_var] <- NA
      df[, new_var][df[, column] <= 0] <- "lower"
      df[, new_var][df[, column] > 0 & df[, column] <= med_value] <- "middle"
      df[, new_var][df[, column] > med_value] <- "upper"
      return(df)
}

for (j in 1:length(modelData_list)) {
      load(modelData_list[j])
      species_results[j, 'species'] <- species
      
      pred.index <- inla.stack.index(stack = integrated_stack, tag = "pred")$data
      lincomb.index.temp <- grep('TZ_lc_noTempMax', rownames(model$summary.lincomb.derived))
      lincomb.index.rain <- grep('TZ_lc_noRain', rownames(model$summary.lincomb.derived))
      lincomb.index.dry <- grep('TZ_lc_noDryspell', rownames(model$summary.lincomb.derived))
      lincomb.index.BG <- grep('TZ_lc_noBG', rownames(model$summary.lincomb.derived))
      lincomb.index.HFP <- grep('TZ_lc_noHFP', rownames(model$summary.lincomb.derived))
      lincomb.index.allFixed <- grep('TZ_lc_allFixed', rownames(model$summary.lincomb.derived))
      
      IP_df <- data.frame(IP_sp) %>% dplyr::select(date_index, x, y)
      
      IP_df$pred_median <- model$summary.fitted.values[pred.index, "0.5quant"]
      IP_df$pred_median_P <- inla.link.cloglog(IP_df$pred_median, inverse = TRUE)
      
      IP_df$pred_sd <- model$summary.fitted.values[pred.index, "sd"]
      
      IP_df$no_temp_median <- model$summary.lincomb.derived[lincomb.index.temp, "0.5quant"]
      IP_df$no_rain_median <- model$summary.lincomb.derived[lincomb.index.rain, "0.5quant"]
      IP_df$no_dry_median <- model$summary.lincomb.derived[lincomb.index.dry, "0.5quant"]
      IP_df$no_BG_median <- model$summary.lincomb.derived[lincomb.index.BG, "0.5quant"]
      IP_df$no_HFP_median <- model$summary.lincomb.derived[lincomb.index.HFP, "0.5quant"]
      IP_df$allFixed_median <- model$summary.lincomb.derived[TZ_lc_allFixed, "0.5quant"]
      
      pred_data <- data.frame(median = IP_df$pred_median,
                              median_P = IP_df$pred_median_P,
                              no_temp_median = IP_df$no_temp_median,
                              no_rain_median = IP_df$no_rain_median,
                              no_dry_median = IP_df$no_dry_median,
                              no_BG_median = IP_df$no_BG_median,
                              no_HFP_median = IP_df$no_HFP_median,
                              allFixed_median = IP_df$allFixed_median,
                              sd = IP_df$pred_sd,
                              ind = rep(c(1,2), each = length(IP_df$pred_median)/2))
      
      predcoordsGroup <- do.call(rbind, list(predcoords, predcoords))
      pred_data_spdf <- sp::SpatialPixelsDataFrame(points = predcoordsGroup, 
                                                   data = pred_data, 
                                                   proj4string = proj)
      
      pred_data_spdf <- raster::crop(pred_data_spdf, TZ_no_lakes)
      pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 1] <- "1980-1999"
      pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 2] <- "2000-2020"
      
      rsq <- function(x, y) summary(lm(y~x))$r.squared
      
      species_results[j, 'rain_imp'] <- 1- rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_rain_median)
      species_results[j, 'temp_imp'] <- 1 - rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_temp_total)
      species_results[j, 'dry_imp'] <- 1 - rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_dry_total)
      species_results[j, 'BG_imp'] <- 1 - rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_BG_total)
      species_results[j, 'HFP_imp'] <- 1 - rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_HFP_total)
      
      range_diff <- pred_data_spdf
      dist_20s <- range_diff@data[["median_P"]][range_diff@data[["ind"]] ==  "2000-2020"]
      dist_80s <- range_diff@data[["median_P"]][range_diff@data[["ind"]] ==  "1980-1999"]

      # net range expansion:
      range_exp <- sum(dist_20s - dist_80s, na.rm = TRUE) # expected number of pixels that have changed
      species_results[j, 'range_change'] <- range_exp
      
      print(paste0("Finished ", j, " of ", length(modelData_list), ": ", species))
      
}

results_unique <- species_results[!duplicated(species_results$species),] 
results_unique <- add_quantiles(results_unique, "range_change", "r_change_quant")

species_long <- results_unique %>% 
      mutate_if(is.numeric, round, digits = 3) %>% 
      pivot_longer(cols = 2:6, names_to = 'covariate', values_to = 'perc_explained') 

write.csv(species_long, 'model_output/species_long.csv')
species_long <- read.csv('model_output/species_long.csv') 

# 
# species_wide <- species_long %>% 
#       pivot_wider(names_from = covariate, values_from = perc_explained) 
# 
# species_wide$rain_sum <- paste0("S:", species_wide$rain_imp_quant, "- E:", species_wide$rain_exp_quant)
# species_wide$temp_sum <- paste0("S:", species_wide$temp_imp_quant, "- E:", species_wide$temp_exp_quant)
# species_wide$dry_sum <- paste0("S:", species_wide$dry_imp_quant, "- E:", species_wide$dry_exp_quant)
# species_wide$BG_sum <- paste0("S:", species_wide$BG_imp_quant, "- E:", species_wide$BG_exp_quant)
# 
# species_wide$vulnerability <- NA
# species_wide$vulnerability[species_wide$BG_imp_quant == "upper" & species_wide$BG_exp_quant == "upper"] <- "high"

species_long$covariate <- factor(species_long$covariate, levels = c("temp_imp", "rain_imp", "BG_imp", "dry_imp", "HFP_imp"))

cov_importance <- ggplot(species_long, aes(x = covariate, y = perc_explained*100)) +
      stat_boxplot(geom ='errorbar', width = 0.1) + 
      geom_boxplot(width = 0.2) +
      xlab("") +
      ylab("Relative importance [%]") +
      theme_classic() +
      scale_x_discrete(labels = c("Temperature",  "Rainfall", "Bare ground", "Dryspell", "Human footprint"))
ggsave(filename = "figures/covariate_importance.png", plot = cov_importance, width = 15, height = 10,
       units = "cm", dpi = 400)



range_change <- ggplot(species_long) +
      geom_boxplot(aes(x = r_change_quant, y = perc_explained*100, fill = covariate)) +
      ylab("Relative importance [%]") +
      xlab("Range change group") +
      theme_classic() +
      theme(legend.title=element_blank()) +
      scale_fill_manual(labels = c("Bare ground", "Dryspell", "Rainfall", "Temperature", "Human footprint"),
                        values = c(wes_palettes$Rushmore1))
ggsave(filename = "figures/range_change.png", plot = range_change, width = 15, height = 10,
       units = "cm", dpi = 400)

ggplot() +
      geom_point(data = species_long[species_long$covariate == 'rain_imp', ], 
                 aes(x = perc_explained, y = abs(range_change)), col = "blue") +
      geom_smooth(data = species_long[species_long$covariate == 'rain_imp', ], 
                 aes(x = perc_explained, y = abs(range_change)), col = "blue", method = "lm") +
      geom_point(data = species_long[species_long$covariate == 'temp_imp', ], 
                 aes(x = perc_explained, y = abs(range_change)), col = "red") +
      geom_smooth(data = species_long[species_long$covariate == 'temp_imp', ], 
                 aes(x = perc_explained, y = abs(range_change)), col = "red", method = "lm") +
      geom_point(data = species_long[species_long$covariate == 'dry_imp', ], 
                 aes(x = perc_explained, y = abs(range_change)), col = "green") +
      geom_smooth(data = species_long[species_long$covariate == 'dry_imp', ], 
                 aes(x = perc_explained, y = abs(range_change)), col = "green", method = "lm") +
      geom_point(data = species_long[species_long$covariate == 'BG_imp', ], 
                 aes(x = perc_explained, y = abs(range_change)), col = "yellow") +
      geom_smooth(data = species_long[species_long$covariate == 'BG_imp', ], 
                 aes(x = perc_explained, y = abs(range_change)), col = "yellow", method = "lm") +
      theme_classic() 

raindf <- species_wide %>% filter(rain_imp_quant != "middle", rain_exp_quant != "middle") %>% 
      mutate(var = "rain",
             sum = rain_sum)
tempdf <- species_wide %>% filter(temp_imp_quant  != "middle", temp_exp_quant != "middle") %>% 
      mutate(var = "temp",
             sum = temp_sum)
drydf <- species_wide %>% filter(dry_imp_quant  != "middle", dry_exp_quant != "middle") %>% 
      mutate(var = "dry",
             sum = dry_sum)
BGdf <- species_wide %>% filter(BG_imp_quant  != "middle", BG_exp_quant != "middle") %>% 
      mutate(var = "BG",
             sum = BG_sum)

df_total <- rbind(raindf, tempdf, drydf, BGdf)

ggplot(df_total, aes(x = sum, y = range_change, fill = var)) +
      geom_boxplot() +
      geom_point(aes(group = sum), pch = 21, 
                 position = position_dodge(1)) +
      theme_classic()



