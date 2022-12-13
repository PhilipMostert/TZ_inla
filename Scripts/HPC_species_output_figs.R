library(tidyverse)
library(INLA)
library(patchwork)
library(ggthemes)
library(raster)

setwd("/users/jhw538/scratch/TZ_spatio_temporal_INLA/model_output/species_figures")

modelData_list <- list.files(path = "/users/jhw538/scratch/TZ_spatio_temporal_INLA/model_output/good/", pattern = ".RData",
                             full.names = TRUE, recursive = TRUE)

for (j in 1:length(modelData_list)) {
      load(modelData_list[j])
      
      dir.create(gsub(" ", "_", species))
      setwd(gsub(" ", "_", species))
      
      pred.index <- inla.stack.index(stack = integrated_stack, tag = "pred")$data
      
      IP_df <- data.frame(SP_Points_data) %>% dplyr::select(date_index, x, y)
      
      IP_df$pred_median <- model$summary.fitted.values[pred.index, "0.5quant"]
      IP_df$pred_sd <- model$summary.fitted.values[pred.index, "sd"]
      
      pred_data <- data.frame(all_pred = IP_df$pred_median,
                              sd = IP_df$pred_sd,
                              ind = rep(c(1,2), each = length(IP_df$pred_median)/2))
      
      predcoordsGroup <- do.call(rbind, list(predcoords, predcoords))
      pred_data_spdf <- sp::SpatialPixelsDataFrame(points = predcoordsGroup, 
                                                   data = pred_data, 
                                                   proj4string = proj)
      
      pred_data_spdf <- crop(pred_data_spdf, TZ_no_lakes)
      pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 1] <- "1980-1999"
      pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 2] <- "2000-2020"
      
      
      all_obs_data <- data.frame(pred_data_spdf) %>% 
            dplyr::select(x, y, ind, all_pred, sd) %>% 
            mutate(pred_median = inla.link.cloglog(all_pred, inv = TRUE),
                   pred_sd = inla.link.cloglog(sd, inv = TRUE))
      
      atlas_pres <- as.data.frame(atlas_sp) 
      atlas_pres$ind <- ifelse(atlas_pres$date_index == 1, "1980-1999", "2000-2020")
      
      eBird_pres <- as.data.frame(ebird_sp) 
      eBird_pres$ind <- ifelse(eBird_pres$date_index == 1, "1980-1999", "2000-2020")
      
      median_plot <- ggplot() + 
            geom_tile(data = all_obs_data, aes(x = x, y = y, fill = pred_median)) + 
            geom_point(data = atlas_pres,
                       aes(x = x, y = y, col = as.factor(presence), alpha = as.factor(presence)), pch = 15, cex = 3) +
            geom_point(data = eBird_pres %>% filter(presence == 0),
                       aes(x = x, y = y), col =  "orangered", alpha = 0.3) +
            geom_point(data = eBird_pres %>% filter(presence == 1),
                       aes(x = x, y = y), col =  "black", alpha = 1) +
            facet_grid( ~ ind) +
            theme(legend.position = "bottom") +
            coord_equal() +
            viridis::scale_fill_viridis("Median", na.value="white") +
            scale_color_manual("Occurrence", values = c("orangered", "black")) +
            scale_alpha_manual(values = c(0.3, 0.8), guide = "none") +
            theme_void() +
            theme(legend.position = "bottom") 
      
      sd_plot <- ggplot() + 
            geom_tile(data = all_obs_data, aes(x = x, y = y, fill = pred_sd)) + 
            facet_grid( ~ ind) +
            coord_equal() +
            theme(legend.position = "none") +
            viridis::scale_fill_viridis("Median", na.value="white") +
            scale_color_distiller(palette = "Reds", direction = 1) +
            theme_void() +
            theme(legend.position = "none") 
      
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
      
      facet_labels <- c(
            TZ_ann_rain = "Annual rainfall",
            TZ_BG = "Bareground cover",
            TZ_dryspell = "Longest dryspell duration",
            TZ_max_temp = "Hottest temperature",
            TZ_HFP = "Human footprint"
      )
      
      effects_plot <- ggplot(effect_combs_m) +
            geom_line(aes(x = orig_values, y = quant_05)) +
            geom_line(aes(x = orig_values, y = quant_0025), lty = 2, alpha = .5) +
            geom_line(aes(x = orig_values, y = quant_0975), lty = 2, alpha = .5) +
            facet_wrap(~ covariate, scale = 'free_x', labeller = as_labeller(facet_labels), ncol = 3) +
            theme_few() +
            xlab("Covariate value") +
            ylab("Probability of occurrence")
      
      species_plot <- median_plot + sd_plot + effects_plot + 
            plot_layout(nrow = 3, heights = c(1,1,1.2), widths = 1) +
            plot_annotation(title = paste0(species, ": Presence data and model predictions"),
                            theme = theme(plot.title = element_text(size = 20)))
      
      pdf(paste0(gsub(" ", "_", species), "_results_figures.pdf"), width = 10, height = 15)
      species_plot
      dev.off()
      
      setwd('..')
}


