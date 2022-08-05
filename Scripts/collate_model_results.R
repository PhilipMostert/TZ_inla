library(INLA)
library(tidyverse)
library(raster)

modelData_list <- list.files(path = "/Users/joriswiethase/Downloads/files_to_move/", pattern = ".RData",
                      full.names = TRUE, recursive = TRUE)

species_results <- data.frame(matrix(NA, nrow = length(modelData_list), ncol = 10))
colnames(species_results) <- c('species',
                               'rain_imp', 'temp_imp', 'dry_imp', 'BG_imp', 'HFP_imp',
                               'abs_range_change', 'log_prop_change', "cells_colonised", "cells_lost",
                               'rain_min', 'rain_max', 'temp_min', 'temp_max', 'BG_min', 'BG_max',
                               'dry_min', 'dry_max', 'HFP_min', 'HFP_max')

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
      IP_df$allFixed_median <- model$summary.lincomb.derived[lincomb.index.allFixed, "0.5quant"]

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

      species_results[j, 'rain_imp'] <- 1 - rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_rain_median)
      species_results[j, 'temp_imp'] <- 1 - rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_temp_median)
      species_results[j, 'dry_imp'] <- 1 - rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_dry_median)
      species_results[j, 'BG_imp'] <- 1 - rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_BG_median)
      species_results[j, 'HFP_imp'] <- 1 - rsq(pred_data_spdf$allFixed_median, pred_data_spdf$no_HFP_median)
      
      original_values <- data.frame(orig_values = unlist(all.seq)) %>% 
            mutate(covariate = sub("*.seq\\d+", "", rownames(.)),
                   sequence = as.numeric(gsub("\\D", "", rownames(.))))
      
      scale_params <- median(model$summary.lincomb.derived$`0.5quant`[grep("TZ_ann_rain|TZ_BG|TZ_dryspell|TZ_max_temp|TZ_HFP", 
                                                                           rownames(model$summary.lincomb.derived))])
      effect_combs <- data.frame(covariate = gsub('[[:digit:]]+', '', sub("*_lc\\d+", "", rownames(model$summary.lincomb.derived))),
                                 sequence = as.numeric(gsub("\\D", "", rownames(model$summary.lincomb.derived))),
                                 quant_05 = inla.link.cloglog((model$summary.lincomb.derived$`0.5quant` - scale_params) - 0.36, inverse = TRUE),
                                 quant_0025 = inla.link.cloglog((model$summary.lincomb.derived$`0.025quant` - scale_params) - 0.36, inverse = TRUE),
                                 quant_0975 = inla.link.cloglog((model$summary.lincomb.derived$`0.975quant` - scale_params) - 0.36, inverse = TRUE))
      effect_combs_main <- effect_combs %>% filter(covariate %in% c("TZ_ann_rain", "TZ_BG", "TZ_dryspell", "TZ_max_temp", "TZ_HFP"))
      effect_combs_m <- merge(original_values, effect_combs_main)
      
      df <- data.frame(covariate = unique(effect_combs_m$covariate))
      
      min_max <- effect_combs_m %>% 
            group_by(covariate) %>% 
            filter(quant_05 > 0.5) %>% 
            summarize(min = min(orig_values),
                      max = max(orig_values))
      df_merged <- merge(df, min_max, all.x = TRUE)
      
      species_results[j, 'rain_min'] <- df_merged$min[df_merged$covariate == "TZ_ann_rain"]
      species_results[j, 'rain_max'] <- df_merged$max[df_merged$covariate == "TZ_ann_rain"]
      species_results[j, 'temp_min'] <- df_merged$min[df_merged$covariate == "TZ_max_temp"]
      species_results[j, 'temp_max'] <- df_merged$max[df_merged$covariate == "TZ_max_temp"]
      species_results[j, 'BG_min'] <- df_merged$min[df_merged$covariate == "TZ_BG"]
      species_results[j, 'BG_min'] <- df_merged$max[df_merged$covariate == "TZ_BG"]
      species_results[j, 'dry_min'] <- df_merged$min[df_merged$covariate == "TZ_dryspell"]
      species_results[j, 'dry_max'] <- df_merged$max[df_merged$covariate == "TZ_dryspell"]
      species_results[j, 'HFP_min'] <- df_merged$min[df_merged$covariate == "TZ_HFP"]
      species_results[j, 'HFP_max'] <- df_merged$max[df_merged$covariate == "TZ_HFP"]
      

      range_diff <- pred_data_spdf
      dist_20s <- range_diff@data[["median_P"]][range_diff@data[["ind"]] ==  "2000-2020"]
      dist_80s <- range_diff@data[["median_P"]][range_diff@data[["ind"]] ==  "1980-1999"]

      # Absolute range change
      range_exp <- sum(dist_20s - dist_80s, na.rm = TRUE) # expected number of pixels that have changed
      species_results[j, 'abs_range_change'] <- range_exp
      
      # Log proportional range change
      log_prop_change <- log(sum(dist_20s)/sum(dist_80s))
      species_results[j, 'log_prop_change'] <- log_prop_change
      
      # Number of cells colonised and lost
      range_diff@data[["colonisation"]] <- (1-dist_80s) * dist_20s
      range_diff@data[["extinction"]] <- dist_80s * (1-dist_20s)
      
      species_results[j, 'cells_colonised'] <- sum(range_diff@data[["colonisation"]])
      species_results[j, 'cells_lost'] <- sum(range_diff@data[["extinction"]])
                                      

      print(paste0("Finished ", j, " of ", length(modelData_list), ": ", species))

}

write.csv(species_results, 'model_output/species_results.csv')