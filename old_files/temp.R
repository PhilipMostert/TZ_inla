suppressPackageStartupMessages(library(INLA, quietly=TRUE))

library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(rlang)
library(inlabru)
library(ggthemes)
library(gridExtra)
library(ggpubr)
library(patchwork)



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

pred_data_spdf <- crop(pred_data_spdf, TZ_no_lakes)
pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 1] <- "1980-1999"
pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 2] <- "2000-2020"


range_diff <- pred_data_spdf
dist_20s <- range_diff@data[["median_P"]][range_diff@data[["ind"]] ==  "2000-2020"]
dist_80s <- range_diff@data[["median_P"]][range_diff@data[["ind"]] ==  "1980-1999"]
range_diff@data[["median_P"]] <- dist_20s - dist_80s
range_diff@data[["colonisation"]][range_diff@data[["ind"]] ==  "2000-2020"] <- (1-dist_80s) * dist_20s # probability of colonisation
range_diff@data[["extinction"]][range_diff@data[["ind"]] ==  "2000-2020"] <- (dist_80s) * (1-dist_20s) # probability of extinction
range_diff@data[["con_absence"]][range_diff@data[["ind"]] ==  "2000-2020"] <- (1-dist_80s) * (1-dist_20s) # probability of continued absence
range_diff@data[["con_presence"]][range_diff@data[["ind"]] ==  "2000-2020"] <- (dist_80s) * (dist_20s) # probability of continued presence

colonisation <- (1-dist_80s) * dist_20s
extinction <- (dist_80s) * (1-dist_20s)
sum_dist_80s_exp <-  sum((dist_80s) * (1-dist_80s))




median_plot <- ggplot() + 
      gg(pred_data_spdf, aes(x = x, y = y, fill = median_P)) + 
      geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), fill = NA, col = "black") +
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("Median", na.value="white") +
      theme_void(); median_plot

col_plot <- ggplot() +
      gg(data = range_diff, aes(x = x, y = y, fill = colonisation)) +
      geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), fill = NA, col = "black") +
      coord_equal() +
      viridis::scale_fill_viridis("P(Colonisation)", na.value="white") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, vjust = 5)); col_plot

ext_plot <- ggplot() +
      gg(range_diff, aes(x = x, y = y, fill = extinction)) +
      geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), fill = NA, col = "black") +
      coord_equal() +
      viridis::scale_fill_viridis("P(Extinction)", na.value="white") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, vjust = 5)); ext_plot

dist_maps <- median_plot + 
      ggtitle(paste0("Range 80s: ", round(sum(dist_80s)), ", Range 2000s: ",  round(sum(dist_20s)), ", Area of high uncertainty: ", round(sum_dist_80s_exp)))
col_plot2 <- col_plot + 
      ggtitle(paste0("Range colonised: ", round(sum(colonisation), digits = 0),
                     "\nMeaningful range colonised: ", round((sum(colonisation) - sum_dist_80s_exp) / sum(dist_80s), digits = 3)))
ext_plot2 <- ext_plot + 
      ggtitle(paste0("Range lost: ", round(sum(extinction), digits = 0),
                     "\nMeaningful range lost: ", round((sum(extinction) - sum_dist_80s_exp) / sum(dist_80s), digits = 3)))
out_map <- dist_maps / (col_plot2 | ext_plot2) + plot_annotation(
      title = species,
      subtitle = paste0("Proportional change: ", round(sum(dist_20s) / sum(dist_80s), digits = 2))
)

ggsave(plot = out_map, filename = paste0('/Users/joriswiethase/Downloads/Colin_species_examples/', gsub(' ', '_', species), ".png"),
       height = 20, width = 20, units = 'cm')

