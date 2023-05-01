rm(list = ls())
setwd("~/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_inla_spatial_temporal")
args = commandArgs(trailingOnly = TRUE)
options(scipen=999)

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

# Random field should not be too smooth (range > maximum study area extent (~10 degrees for TZ outline without buffer)), or very different between time periods.
# Filter 1: Visually implausible
# Filter 2: DIC or CPO

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Scripts')

#setwd('/home/ahomec/p/philism/Joris_work/scripts')
sapply(list.files(path = '/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_inla_spatial_temporal/source',
                  pattern="*.R", full.names = T), source, .GlobalEnv)

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')

#setwd('/home/ahomec/p/philism/Joris_work/Philip_data')
load("model_data/TZ_INLA_model_file_temporal.RData")

temporal_variables <- temporal_variables[!is.na(rowSums(temporal_variables@data)),]

# ------------------------------------------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------------------------------------------
# Species choices
# Based on list of savannah bird species given in Beale et al. 2013. Species were excluded if there were less than 20 
# records in all date_index - data source combinations

species_list = c('Oenanthe pileata')  
species <- species_list[1]

estimated_range = 3
max.edge = estimated_range/8

max.edge = 0.3

# Approach outlined in Simpson et al. (2017), weakly informative.
# Range: Distance at which spatial autocorrelation is small, or practical range.
# If prior.range = c(10, 0.01), this means that the probability that the range 
# of the spatial effect is less than 10 degrees is presumed to be very small.
# Range defines how far you need to go for two points to be independent, the κ parameter.
# A smaller value gives a rougher texture to the spde, a larger one a smoother.
# These priors shrink the spatial model towards the base model (without a 
# spatial term).
# Sigma: Variation in the spatial effect (how variable is the field from a point to the next, the δ parameter)

# sd(integrated_stack$effects$data$presence, na.rm = T)/10
# Need some knowledge of the scale of response and predictor, we just need to get the general 
# scale correct. 
# Set range so it represents a relatively small distance on the predictor scale.
# Set sigma quite large, we would be surprised if SD exceeded this.
# Higher range values means stricter prior, and stronger smoothing.

# Range guidelines:
# No larger than diff(range(Mesh$mesh$loc[, 1]))/2: 8
# No smaller than max.edge: 0.375
# Good starting point: c(10 * max.edge, 0.5)
# Range and sigma could be chosen so that random field explains up to 25% of variations is
# model predictions

# Sigma guidelines:
# After running an initial model, look at Stdev, compare with chosen sigma prior.
# If far apart, adjust. Keep a bit vague, with P 0.5

# Penalized complexity priors for spde
prior_range = c(1.05, 0.5)  
prior_sigma = c(1.05, 0.1)   

# Gaussian priors for fixed effects
# Default is "non-informative" (or "vague"): 0 mean and precision of 0.001.
# For models with log link, this is no longer non-informative
fixed_mean = 0
fixed_precision = 1

# Model output string
model_name <- paste0(sub(" ", "_", species), "_form2_r", 
prior_range[1], "_", prior_range[2], "_s", prior_sigma[1], "_", prior_sigma[2], 
"_mean", fixed_mean, "_prec", fixed_precision, ".RDS")

# ------------------------------------------------------------------------------------------------------------------------
# Data preparation
# ------------------------------------------------------------------------------------------------------------------------
if (!file.exists(paste0("model_data/", gsub(" ", "_", species), "_model_data.RData"))) {
      
      ebird_full <- ebird_full %>%
            mutate(date_index = ifelse(date > '2000-01-01',2,1))
      
      ebird_filtered <- ebird_full %>% 
            filter(APPROVED == 1,  # Only keep reviewed and approved records
                   `ALL SPECIES REPORTED` == 1,           # Only keep complete checklists
                   `EFFORT DISTANCE KM` < 15,
                   `DURATION MINUTES` >= 5,
                   `DURATION MINUTES` <= 5*60,
                   `NUMBER OBSERVERS` <= 10)   

      ebird_filtered <- ebird_filtered %>% 
            group_by(LATITUDE, LONGITUDE, `SAMPLING EVENT IDENTIFIER`, `DURATION MINUTES`, 
                     `EFFORT DISTANCE KM`, `NUMBER OBSERVERS`, `OBSERVATION DATE`, `LOCALITY`, date_index) %>% 
            summarise(occurrence = ifelse(species %in% `SCIENTIFIC NAME`, TRUE, FALSE)) %>% 
            ungroup()  %>%
            group_by(LATITUDE, LONGITUDE, `DURATION MINUTES`, 
                     `EFFORT DISTANCE KM`, `NUMBER OBSERVERS`, `OBSERVATION DATE`, `LOCALITY`) %>% 
            slice_head() %>%    # Where duplicate checklists occurred, keep only the first
            ungroup() %>% 
            rename(duration_minutes = `DURATION MINUTES`,
                   effort_distance_km = `EFFORT DISTANCE KM`,
                   number_observers = `NUMBER OBSERVERS`,
                   x = LONGITUDE,
                   y = LATITUDE) %>% 
            mutate(ebird_effort = number_observers*duration_minutes)   # Following Humphreys et al. 2018, Diversity and Distributions
      
      ebird_sp <- SpatialPointsDataFrame(
            coords = ebird_filtered[, c("x", "y")],
            data = data.frame(presence = ebird_filtered$occurrence, 
                              ebird_effort = ebird_filtered$ebird_effort,
                              duration_minutes = ebird_filtered$duration_minutes, 
                              number_observers = ebird_filtered$number_observers, 
                              date_index = ebird_filtered$date_index),
            proj4string = crs(proj))
      
      # Only include eBird data points for Tanzania
      # Get intersecting points
      in_sp <- rgeos::gIntersection(ebird_sp, TZ_outline)
      
      # Only keep intersecting points in original spdf
      ebird_sp <- ebird_sp[in_sp, ]
      
      atlas_full <- atlas_full %>%
            filter(time_period != 'x') %>%
            mutate(date_index = ifelse(time_period == '20s',2,1))
      
      atlas_filtered <- atlas_full %>% 
            mutate(Scientific = trimws(Scientific, which = 'both')) %>% 
            filter(Scientific == species) %>% 
            rename(x = Long,
                   y = Lat) %>% 
            mutate(presence = ifelse(occurrence == 1, TRUE, FALSE)) %>% 
            dplyr::select(-V1); if(is_empty(atlas_filtered$presence)){print("ERROR: No Atlas data available")}
      
      atlas_sp <- SpatialPointsDataFrame(
            coords = data.frame(atlas_filtered[, c("x", "y")]),
            data = data.frame(presence = atlas_filtered$presence, effort = atlas_filtered$effort,
                              date_index = atlas_filtered$date_index),
            proj4string = crs(proj))
  
      # Log transform effort, then range 0-1
      range01 <- function(x){(x - min(x))/(max(x) - min(x))}
      ebird_sp$duration_minutes <- range01(log(ebird_sp$duration_minutes+1))
      atlas_sp$effort <- range01(log(atlas_sp$effort+1))

      filtered_covs <- temporal_variables[, c('TZ_ann_rain_log_1980s_1.s', 'TZ_ann_rain_log_2000s_1.s', 
                                              'TZ_ann_rain_log_1980s_2.s', 'TZ_ann_rain_log_2000s_2.s', 
                                              'TZ_max_temp_1980s_1.s', 'TZ_max_temp_2000s_1.s',
                                              'TZ_max_temp_1980s_2.s', 'TZ_max_temp_2000s_2.s',
                                              'TZ_dryspell_1980s_1.s', 'TZ_dryspell_2000s_1.s',
                                              'TZ_dryspell_1980s_2.s', 'TZ_dryspell_2000s_2.s',
                                              'TZ_BG_1980s_1.s', 'TZ_BG_2000s_1.s',
                                              'TZ_BG_1980s_2.s', 'TZ_BG_2000s_2.s',
                                              'indicator_90s', 'indicator_2010s',
                                              'TZ_HFP_1980s_1.s', 'TZ_HFP_2000s_1.s',
                                              'TZ_HFP_1980s_2.s', 'TZ_HFP_2000s_2.s')]
      
            
      # Samples the covariates spatially closest to the occurrence data points    
      Nearest_covs_ebird <- GetNearestCovariate(ebird_sp, filtered_covs)  
      Nearest_covs_atlas <- GetNearestCovariate(atlas_sp, filtered_covs)
      
      # Add sampled covariates to occurrence data
      ebird_sp@data[, names(Nearest_covs_ebird@data)] <- Nearest_covs_ebird@data
      ebird_sp <- as(ebird_sp, 'data.frame')
      
      # Combine covariates from different time periods into single variable, different times identified by separate date_index  
      ebird_sp <- ebird_sp %>% mutate(annual_rain_log_1 = ifelse(date_index == 1, TZ_ann_rain_log_1980s_1.s, TZ_ann_rain_log_2000s_1.s),
                                      annual_rain_log_2 = ifelse(date_index == 1, TZ_ann_rain_log_1980s_2.s, TZ_ann_rain_log_2000s_2.s),
                                      hottest_temp_1 = ifelse(date_index == 1, TZ_max_temp_1980s_1.s, TZ_max_temp_2000s_1.s),
                                      hottest_temp_2 = ifelse(date_index == 1, TZ_max_temp_1980s_2.s, TZ_max_temp_2000s_2.s),
                                      max_dryspell_1 = ifelse(date_index == 1, TZ_dryspell_1980s_1.s, TZ_dryspell_2000s_1.s),
                                      max_dryspell_2 = ifelse(date_index == 1, TZ_dryspell_1980s_2.s, TZ_dryspell_2000s_2.s),
                                      BG_1 = ifelse(date_index == 1, TZ_BG_1980s_1.s, TZ_BG_2000s_1.s),
                                      BG_2 = ifelse(date_index == 1, TZ_BG_1980s_2.s, TZ_BG_2000s_2.s),
                                      indicator = ifelse(date_index == 1, indicator_90s, indicator_2010s),
                                      HFP_1 = ifelse(date_index == 1, TZ_HFP_1980s_1.s, TZ_HFP_2000s_1.s),
                                      HFP_2 = ifelse(date_index == 1, TZ_HFP_1980s_2.s, TZ_HFP_2000s_2.s))
      
      # Make spdf, add intercept, so model can distinguish eBird presence from Atlas presence
      ebird_sp <- SpatialPointsDataFrame(coords = ebird_sp[, c("x", "y")],
                                         data = ebird_sp,
                                         proj4string = crs(proj))
      ebird_sp@data[, 'ebird_intercept'] <- 1
      ebird_sp$presence <- as.numeric(ebird_sp$presence)
      
      atlas_sp@data[, names(Nearest_covs_atlas@data)] <- Nearest_covs_atlas@data
      atlas_sp <- as(atlas_sp, 'data.frame')
      
      atlas_sp <- atlas_sp %>% mutate(annual_rain_log_1 = ifelse(date_index == 1, TZ_ann_rain_log_1980s_1.s, TZ_ann_rain_log_2000s_1.s),
                                      annual_rain_log_2 = ifelse(date_index == 1, TZ_ann_rain_log_1980s_2.s, TZ_ann_rain_log_2000s_2.s),
                                      hottest_temp_1 = ifelse(date_index == 1, TZ_max_temp_1980s_1.s, TZ_max_temp_2000s_1.s),
                                      hottest_temp_2 = ifelse(date_index == 1, TZ_max_temp_1980s_2.s, TZ_max_temp_2000s_2.s),
                                      max_dryspell_1 = ifelse(date_index == 1, TZ_dryspell_1980s_1.s, TZ_dryspell_2000s_1.s),
                                      max_dryspell_2 = ifelse(date_index == 1, TZ_dryspell_1980s_2.s, TZ_dryspell_2000s_2.s),
                                      BG_1 = ifelse(date_index == 1, TZ_BG_1980s_1.s, TZ_BG_2000s_1.s),
                                      BG_2 = ifelse(date_index == 1, TZ_BG_1980s_2.s, TZ_BG_2000s_2.s),
                                      indicator = ifelse(date_index == 1, indicator_90s, indicator_2010s),
                                      HFP_1 = ifelse(date_index == 1, TZ_HFP_1980s_1.s, TZ_HFP_2000s_1.s),
                                      HFP_2 = ifelse(date_index == 1, TZ_HFP_1980s_2.s, TZ_HFP_2000s_2.s))
      
      atlas_sp <- SpatialPointsDataFrame(coords = atlas_sp[, c('x', 'y')],
                                         data = atlas_sp,
                                         proj4string = crs(proj))
      atlas_sp@data[,'atlas_intercept'] <- 1
      atlas_sp$presence <- as.numeric(atlas_sp$presence)

      save(ebird_sp, atlas_sp, file = paste0("model_data/", gsub(" ", "_", species), "_model_data.RData"))

} else {
      
      # setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Covariate_data')
      load(paste0("model_data/", gsub(" ", "_", species), "_model_data.RData"))
}

print(species)
print(paste0('eBird presence records 1980s: ', length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 1, ]))))
print(paste0('eBird presence records 2000s: ', length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 2, ]))))
print(paste0('Atlas presence records 1980s: ', length(rownames(atlas_sp@data[atlas_sp@data$presence == TRUE & as.data.frame(atlas_sp)$date_index == 1, ]))))
print(paste0('Atlas presence records 2000s: ', length(rownames(atlas_sp@data[atlas_sp@data$presence == TRUE & as.data.frame(atlas_sp)$date_index == 2, ]))))


if (!file.exists(paste0("figures/", sub(" ", "_", species), "_raw_ccurrence.png"))){
      plot_80s <- ggplot() +
            geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), 
                         colour = "black", fill = NA) +
            geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 0 & as.data.frame(atlas_sp)$date_index == 1, ], 
                       aes(x = x, y = y, col = "Absence"), alpha = 0.8, pch = 15, cex = 2) +
            geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 1 & as.data.frame(atlas_sp)$date_index == 1, ], 
                       aes(x = x, y = y, col = "Presence"), alpha = 0.8, pch = 15, cex = 3) +
            # geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 0 & as.data.frame(ebird_sp)$date_index == 1, ],
            #            aes(x = x, y = y, fill = "Absence"), pch = 21, alpha = 0.5) +
            # geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 1 & as.data.frame(ebird_sp)$date_index == 1, ],
            #            aes(x = x, y = y, fill = "Presence"), pch = 23, cex = 2) +
            theme_void() +
            coord_equal() +
            # scale_fill_discrete(name = 'eBird') +
            # scale_color_discrete(name = 'Atlas') +
            # ggtitle(paste0("1980-1999, eBird presence: ", length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 1, ]))))  +
            ggtitle("Time period 1") +
            theme(plot.title = element_text(hjust = 0.5))
      
      plot_2000s <- ggplot() +
            geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), 
                         colour = "black", fill = NA) +
            geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 0 & as.data.frame(atlas_sp)$date_index == 2, ], 
                       aes(x = x, y = y, col = "Absence"), alpha = 0.8, pch = 15, cex = 2) +
            geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 1 & as.data.frame(atlas_sp)$date_index == 2, ], 
                       aes(x = x, y = y, col = "Presence"), alpha = 0.8, pch = 15, cex = 3) +
            # geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 0 & as.data.frame(ebird_sp)$date_index == 2, ], 
            #            aes(x = x, y = y, fill = "Absence"), pch = 21, alpha = 0.5) +
            # geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 1 & as.data.frame(ebird_sp)$date_index == 2, ], 
            #            aes(x = x, y = y, fill = "Presence"), pch = 23, cex = 2) +
            theme_void() +
            coord_equal() +
            # scale_fill_discrete(name = 'eBird') +
            # scale_color_discrete(name = 'Atlas') +
            # ggtitle(paste0("2000-2020, eBird presence: ", length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 2, ])))) +
            ggtitle("Time period 2") +
            theme(plot.title = element_text(hjust = 0.5))
      
      raw_plots_comb <- ggpubr::ggarrange(plot_80s, plot_2000s, nrow = 2, common.legend = TRUE, legend = "right")
      
      ggsave(plot = raw_plots_comb, filename = paste0("figures/", sub(" ", "_", species), "_raw_ccurrence.png"),
             width = 18, height = 18, units = 'cm')
}

# ------------------------------------------------------------------------------------------------------------------------
# Mesh
# ------------------------------------------------------------------------------------------------------------------------
if (!file.exists(paste0("model_data/mesh_E", max.edge,".RData"))) {
      Meshpars <- list(max.edge = c(max.edge, max.edge*4), 
                       offset = c(max.edge, max.edge*5), 
                       cutoff = max.edge/2)

      Mesh <- MakeSpatialRegion(
            data = NULL,
            bdry = ROI,    
            meshpars = Meshpars,
            proj = proj
      )
      Mesh$mesh$crs <- proj
      plot(Mesh$mesh)
      save(Mesh, file = paste0("model_data/mesh_E", max.edge,".RData"))
} else {
      load(paste0("model_data/mesh_E", max.edge,".RData"))
}

# ------------------------------------------------------------------------------------------------------------------------
# SPDE model on the mesh
# ------------------------------------------------------------------------------------------------------------------------
# Use Penalized Complexity prior, this determines the smoothness of the spatial effect.
pcspde <- inla.spde2.pcmatern(
      mesh = Mesh$mesh,
      alpha = 2,
      prior.range = prior_range,   
      prior.sigma = prior_sigma)

eBird_spde <- inla.spde2.pcmatern(
      mesh = Mesh$mesh,
      alpha = 2,
      prior.range = prior_range,
      prior.sigma = prior_sigma)

# To use default priors: inla.spde2.matern(mesh = Mesh$mesh, alpha = 2)

# ------------------------------------------------------------------------------------------------------------------------
# Space-time index set
# ------------------------------------------------------------------------------------------------------------------------
# Set the number of time periods
n_time_layers = 2

# Set index for the latent field, taking into account the number of mesh points in the SPDE model, 
# and the number of time groups.
# This index does not depend on any data set locations, only SPDE model size and time dimensions.
# A double index, identifying both the spatial location and associated time point. 1 to n.spde index points 
# for n.group times.
index_set <- inla.spde.make.index(name ='shared.field',
                                  n.spde = pcspde$n.spde,
                                  n.group = n_time_layers)
lengths(index_set)

index_set_eBird <- inla.spde.make.index(name ='eBird.field',
                                  n.spde = eBird_spde$n.spde,
                                  n.group = n_time_layers)

# ------------------------------------------------------------------------------------------------------------------------
# Projection matrices
# ------------------------------------------------------------------------------------------------------------------------
# Make the projector matrix, which projects the spatio-temporal continuous Gaussian
# random field from observations to mesh nodes.
# Uses coordinates of the observed data

projmat_eBird <- inla.spde.make.A(mesh = Mesh$mesh,
                                  loc = as.matrix(ebird_sp@coords),
                                  group = ebird_sp$date_index)

projmat_atlas <- inla.spde.make.A(mesh = Mesh$mesh,
                                  loc = as.matrix(atlas_sp@coords),
                                  group = atlas_sp$date_index)

dim(projmat_eBird)  # Matrix equal to number of observations by number of indices
dim(projmat_atlas)

# ------------------------------------------------------------------------------------------------------------------------
# Estimation stacks
# ------------------------------------------------------------------------------------------------------------------------
stk.eBird <- inla.stack(data = list(resp = ebird_sp@data[, 'presence']),
                        A = list(1, projmat_eBird, projmat_eBird),
                        tag = 'eBird',
                        effects = list(ebird_sp@data, eBird.field = index_set_eBird, shared.field = index_set))

stk.atlas <- inla.stack(data = list(resp = atlas_sp@data[, 'presence']),
                        A = list(1, projmat_atlas),
                        tag = 'atlas',
                        effects = list(atlas_sp@data, shared.field = index_set))


# ------------------------------------------------------------------------------------------------------------------------
# Prediction data
# ------------------------------------------------------------------------------------------------------------------------
# Contains the locations and times where we want to make predictions. Code for this is adapted from the 
# 'MakeProjectionGrid' function.

if(!file.exists(paste0("model_data/pred_files_E", max.edge,".RData"))){
      # Set grid locations for predictions
      Nxy.scale <- 0.1  # about 10km resolution
      
      Boundary <- Mesh$mesh$loc[Mesh$mesh$segm$int$idx[, 2], ]
      Nxy.size <- c(diff(range(Boundary[, 1])), diff(range(Boundary[, 2])))
      Nxy <- round(Nxy.size / Nxy.scale)
      
      projgrid <- inla.mesh.projector(Mesh$mesh, 
                                      xlim = range(Boundary[, 1]), 
                                      ylim = range(Boundary[, 2]), 
                                      dims = Nxy)
      
      # Get the index of points on the grid within the boundary
      xy.in <- splancs::inout(projgrid$lattice$loc, Boundary)
      
      # Select only points on the grid that fall within the boundary
      predcoords <- projgrid$lattice$loc[which(xy.in), ]
      colnames(predcoords) <- c('x','y')
      Apred <- projgrid$proj$A[which(xy.in), ]
      
      
      # Construct the prediction data, which is spatial points by temporal layers
      spatial_points <- cbind(c(Mesh$mesh$loc[,1]), c(Mesh$mesh$loc[,2]))
      spatial_points <- projgrid$lattice$loc[which(xy.in), ]
      colnames(spatial_points) <- c("x", "y")
      
      NearestPredCovs <- GetNearestCovariate(points = spatial_points, covs = temporal_variables)
      
      spatiotemporal_points <- rbind(spatial_points, spatial_points)
      SP_Points_data <- data.frame(date_index = rep(c(1,2), each = nrow(spatiotemporal_points)/2))
      SP_Points_data[, 'annual_rain_log_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_ann_rain_log_1980s_1.s, NearestPredCovs@data$TZ_ann_rain_log_2000s_1.s)
      SP_Points_data[, 'annual_rain_log_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_ann_rain_log_1980s_2.s, NearestPredCovs@data$TZ_ann_rain_log_2000s_2.s)
      
      SP_Points_data[, 'hottest_temp_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_max_temp_1980s_1.s, NearestPredCovs@data$TZ_max_temp_2000s_1.s)
      SP_Points_data[, 'hottest_temp_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_max_temp_1980s_2.s, NearestPredCovs@data$TZ_max_temp_2000s_2.s)
      
      SP_Points_data[, 'max_dryspell_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_dryspell_1980s_1.s, NearestPredCovs@data$TZ_dryspell_2000s_1.s)
      SP_Points_data[, 'max_dryspell_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_dryspell_1980s_2.s, NearestPredCovs@data$TZ_dryspell_2000s_2.s)
      
      SP_Points_data[, 'BG_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_BG_1980s_1.s, NearestPredCovs@data$TZ_BG_2000s_1.s)
      SP_Points_data[, 'BG_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_BG_1980s_2.s, NearestPredCovs@data$TZ_BG_2000s_2.s)
      
      SP_Points_data[, 'indicator'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$indicator_90s, NearestPredCovs@data$indicator_2010s)
      
      SP_Points_data[, 'HFP_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_HFP_1980s_1.s, NearestPredCovs@data$TZ_HFP_2000s_1.s)
      SP_Points_data[, 'HFP_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_HFP_1980s_2.s, NearestPredCovs@data$TZ_HFP_2000s_2.s)
      
      # Add effort variables as constantly high effort across the projection grid
      SP_Points_data[, 'duration_minutes'] <- 0.9
      SP_Points_data[, 'effort'] <- 0.9
      
      SP_Points_data <- cbind(SP_Points_data, spatiotemporal_points)
      
      IP_sp <- sp::SpatialPointsDataFrame(coords = spatiotemporal_points, data = SP_Points_data, proj4string = proj)
      IP_sp@data$Intercept <- 1
      IP_sp@data[, c("x", "y")] <- IP_sp@coords
      
      projmat.pred <- inla.spde.make.A(mesh = Mesh$mesh, 
                                       loc = IP_sp@coords, 
                                       group = IP_sp@data$date_index)

# ------------------------------------------------------------------------------------------------------------------------
# Prediction stack
# ------------------------------------------------------------------------------------------------------------------------
stk.pred <- inla.stack(tag='pred',
                       data = list(resp = NA),
                       A = list(1, projmat.pred),
                       effects = list(IP_sp@data,
                                      shared.field = index_set))

      
      save(stk.pred, SP_Points_data, projmat.pred, IP_sp, predcoords, projgrid, xy.in, file = paste0("model_data/CHELSA_pred_files_E", max.edge,"_1km.RData"))
} else {
      load(paste0("model_data/pred_files_E", max.edge,".RData"))
}

# ------------------------------------------------------------------------------------------------------------------------
# Integration stack
# ------------------------------------------------------------------------------------------------------------------------
# LGCP, spatially varying intensity. 'MakeIntegrationStack' adds area of mesh 
# polygon as 'e'. Probability of observing a certain number of points in that 
# area follows a Poisson distribution with defined intensity.

if (!file.exists(paste0("model_data/stk_ip_E", max.edge,".RData"))) {
      #get coordindate and covariate column names
      coordnames <- colnames(temporal_variables@coords)
      Names <- colnames(temporal_variables@data)
      
      Names  <- c(coordnames, Names)
      
      #get mesh triangle centroids
      Points <- cbind(c(Mesh$mesh$loc[,1]), c(Mesh$mesh$loc[,2]))
      
      #set column names of centroids to cordinate names from covariate data
      colnames(Points) <- coordnames
      
      #get value of nearest covariate to each mesh centroids
      NearestCovs <- GetNearestCovariate(points=Points, covs=temporal_variables)
      
      #set intercept to 1 for each mesh element
      NearestCovs$Intercept <- rep(1,nrow(NearestCovs))
      
      #add coordinate to data part of spatialpoints object if they are given
      NearestCovs@data[, colnames(NearestCovs@coords)] <- NearestCovs@coords
      covs_duplicated <- rbind(NearestCovs@data, NearestCovs@data)
      covs_duplicated$index <- rep(c(1,2), each = nrow(covs_duplicated)/2)
            
      # Projector matrix for integration points.
      projmat.ip <- Matrix::Diagonal(n = n_time_layers*pcspde$n.spde)  # from mesh to integration points
      
      # repeat the i-points for the number of time periods under consideration. 
      # So if you have 1000 i-points and 2 time periods, you'll end up with a 
      # data.frame with 2000 points, with either the index 1 or 2.
      
      #create inla stack
      stk.ip <- inla.stack(tag = "ip",
                           data = list(resp = NA, e = Mesh$w),  # 0 for count, NA for incidental, e is area of mesh polygon
                           A = list(1, projmat.ip, projmat.ip),
                           effects = list(covs_duplicated,
                                          shared.field = index_set,
                                          eBird.field = index_set_eBird))
      
      save(stk.ip, file = paste0("model_data/CHELSA_stk_ip_E", max.edge,"_1km.RData"))
} else {
      load(paste0("model_data/stk_ip_E", max.edge,".RData"))
}



integrated_stack <- inla.stack(stk.eBird, stk.atlas, stk.pred, stk.ip)


# ------------------------------------------------------------------------------------------------------------------------
# Model formula
# ------------------------------------------------------------------------------------------------------------------------
index_list <- paste0("date_index", 1:6)

form_1_stack <- integrated_stack
for (i in 1:length(index_list)){
      new_var <- index_list[i]
      form_1_stack[["effects"]][["data"]][[new_var]] <- form_1_stack[["effects"]][["data"]][["date_index"]]
      form_1_stack[["effects"]][["ncol"]][[new_var]] <- form_1_stack[["effects"]][["ncol"]][["date_index"]]
      form_1_stack[["effects"]][["names"]][[new_var]] <- form_1_stack[["effects"]][["names"]][["date_index"]]
}

# Define a PC-prior for the temporal autoregressive parameter. This is the probability of the 
# standard deviation being higher than the given number. If 'param = c(4, 0.01)', 
# P(sd > 4) = 0.01
sdres <- sd(integrated_stack$effects$data$presence, na.rm = T)
h.spec <- list(rho = list(prior="pc.prec", param = c(0.5, 0.01)))

# At each time point, spatial locations are linked through the spde.
# Across time, the process evolves according to an AR(1) process.

# MODEL 1 -----------------------------------------------------------------------------------------------------------------
form_1 <- resp ~ 0 +
      ebird_intercept +
      atlas_intercept +
      ebird_effort + effort + 
      f(date_index, annual_rain_log_1, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +  
      f(date_index1, annual_rain_log_2, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +   
      f(date_index2, hottest_temp_1, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +   
      f(date_index3, hottest_temp_2, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +
      f(date_index4, max_dryspell_1, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +
      f(date_index5, max_dryspell_2, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +
      f(i, model = pcspde, group = i.group, control.group = list(model = 'ar1'))

# MODEL 2 -----------------------------------------------------------------------------------------------------------------
# f(annual_rain_log_1, model = "linear", mean.linear = 0, prec.linear = 0.001)
form_2 <- resp ~ 0 + ebird_intercept + atlas_intercept + 
      annual_rain_log_1 + annual_rain_log_2 + hottest_temp_1 + hottest_temp_2 + 
      max_dryspell_1 + max_dryspell_2 + 
      BG_1 + BG_2 + indicator + HFP_1 + HFP_2 +
      x + y +
      duration_minutes + effort + 
      f(shared.field, model = pcspde, group = shared.field.group,  # In each year, spatial locations are linked by spde
        control.group = list(model = 'ar1', hyper = h.spec)) +  # Across time, process evolves according to AR(1) process
      f(eBird.field, model = eBird_spde, group = eBird.field.group,
        control.group = list(model = 'ar1', hyper = h.spec))

# ------------------------------------------------------------------------------------------------------------------------
# Linear combinations
# ------------------------------------------------------------------------------------------------------------------------

# Add linear combinations to calculate pseudo-R-Squared for covariates.
# This computes linear combinations on some effects without altering model fitting, estimating their
# posterior marginals (derived using model$summary.lincomb.derived). Linear combinations effectively
# allow us to isolate the effect of different covariates (or combinations of covariates).

TZ_lc_allFixed <-  inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                      atlas_intercept = rep(1, NROW(SP_Points_data)),
                                      annual_rain_log_1 = SP_Points_data$annual_rain_log_1,
                                      annual_rain_log_2 = SP_Points_data$annual_rain_log_2,
                                      max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                      max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                      hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                      hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                      HFP_1 = SP_Points_data$HFP_1,
                                      HFP_2 = SP_Points_data$HFP_2,
                                      BG_1 = SP_Points_data$BG_1,
                                      BG_2 = SP_Points_data$BG_2,
                                      duration_minutes = SP_Points_data$duration_minutes,
                                      effort = SP_Points_data$effort,
                                      x = SP_Points_data$x,
                                      y = SP_Points_data$y,
                                      indicator = SP_Points_data$indicator)
names(TZ_lc_allFixed) <- paste0("TZ_lc_allFixed", 1:NROW(SP_Points_data))
TZ_lc_noTempMax <-  inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                       atlas_intercept = rep(1, NROW(SP_Points_data)),
                                       annual_rain_log_1 = SP_Points_data$annual_rain_log_1,
                                       annual_rain_log_2 = SP_Points_data$annual_rain_log_2,
                                       max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                       max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                       HFP_1 = SP_Points_data$HFP_1,
                                       HFP_2 = SP_Points_data$HFP_2,
                                       BG_1 = SP_Points_data$BG_1,
                                       BG_2 = SP_Points_data$BG_2,
                                       duration_minutes = SP_Points_data$duration_minutes,
                                       effort = SP_Points_data$effort,
                                       x = SP_Points_data$x,
                                       y = SP_Points_data$y,
                                       indicator = SP_Points_data$indicator)
names(TZ_lc_noTempMax) <- paste0("TZ_lc_noTempMax", 1:NROW(SP_Points_data))
TZ_lc_noRain <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                   atlas_intercept = rep(1, NROW(SP_Points_data)),
                                   max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                   max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                   hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                   hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                   HFP_1 = SP_Points_data$HFP_1,
                                   HFP_2 = SP_Points_data$HFP_2,
                                   BG_1 = SP_Points_data$BG_1,
                                   BG_2 = SP_Points_data$BG_2,
                                   duration_minutes = SP_Points_data$duration_minutes,
                                   effort = SP_Points_data$effort,
                                   x = SP_Points_data$x,
                                   y = SP_Points_data$y,
                                   indicator = SP_Points_data$indicator)
names(TZ_lc_noRain) <- paste0("TZ_lc_noRain", 1:NROW(SP_Points_data))
TZ_lc_noDryspell <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                       atlas_intercept = rep(1, NROW(SP_Points_data)),
                                       annual_rain_log_1 = SP_Points_data$annual_rain_log_1,
                                       annual_rain_log_2 = SP_Points_data$annual_rain_log_2,
                                       hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                       hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                       HFP_1 = SP_Points_data$HFP_1,
                                       HFP_2 = SP_Points_data$HFP_2,
                                       BG_1 = SP_Points_data$BG_1,
                                       BG_2 = SP_Points_data$BG_2,
                                       duration_minutes = SP_Points_data$duration_minutes,
                                       effort = SP_Points_data$effort,
                                       x = SP_Points_data$x,
                                       y = SP_Points_data$y,
                                       indicator = SP_Points_data$indicator)
names(TZ_lc_noDryspell) <- paste0("TZ_lc_noDryspell", 1:NROW(SP_Points_data))
TZ_lc_noBG <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                 atlas_intercept = rep(1, NROW(SP_Points_data)),
                                 annual_rain_log_1 = SP_Points_data$annual_rain_log_1,
                                 annual_rain_log_2 = SP_Points_data$annual_rain_log_2,
                                 max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                 max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                 hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                 hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                 HFP_1 = SP_Points_data$HFP_1,
                                 HFP_2 = SP_Points_data$HFP_2,
                                 duration_minutes = SP_Points_data$duration_minutes,
                                 effort = SP_Points_data$effort,
                                 x = SP_Points_data$x,
                                 y = SP_Points_data$y)
names(TZ_lc_noBG) <- paste0("TZ_lc_noBG", 1:NROW(SP_Points_data))
TZ_lc_noHFP <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                  atlas_intercept = rep(1, NROW(SP_Points_data)),
                                  annual_rain_log_1 = SP_Points_data$annual_rain_log_1,
                                  annual_rain_log_2 = SP_Points_data$annual_rain_log_2,
                                  max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                  max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                  hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                  hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                  BG_1 = SP_Points_data$BG_1,
                                  BG_2 = SP_Points_data$BG_2,
                                  duration_minutes = SP_Points_data$duration_minutes,
                                  effort = SP_Points_data$effort,
                                  x = SP_Points_data$x,
                                  y = SP_Points_data$y,
                                  indicator = SP_Points_data$indicator)
names(TZ_lc_noHFP) <- paste0("TZ_lc_noHFP", 1:NROW(SP_Points_data))

lc_combined <- c(all_lc, TZ_lc_noRain, TZ_lc_noTempMax, TZ_lc_noDryspell, TZ_lc_noBG, TZ_lc_noHFP, TZ_lc_allFixed)
save(lc_combined, file = "model_data/CHELSA_lc_combined.RData")
load("model_data/lc_combined.RData")
      
# ------------------------------------------------------------------------------------------------------------------------
# INLA model
# ------------------------------------------------------------------------------------------------------------------------

# Gaussian priors
C.F. <- list(
      mean = fixed_mean,
      prec = fixed_precision   # Precision for all fixed effects except intercept
)

if(!file.exists(paste0("model_output/", model_name))){
      # lc_no_BG <- all_lc[grepl("BG", all_lc) == FALSE]
      model <- inla(form_2, family = "binomial", control.family = list(link = "cloglog"), 
                    # lincomb = allFixed_lc,
                    data = inla.stack.data(integrated_stack), 
                    verbose = FALSE,
                    control.predictor = list(A = inla.stack.A(integrated_stack), 
                                             link = NULL, compute = TRUE), 
                    control.fixed = C.F.,
                    E = inla.stack.data(integrated_stack)$e, 
                    control.compute = list(waic = FALSE, dic = TRUE, cpo = TRUE))
      saveRDS(model, file = '/Users/joriswiethase/Downloads/model_output/good/Iduna_pallida/test_model_with_effort.RDS')
      saveRDS(model, file = paste0("model_output/", model_name))
} else {
      model <- readRDS(paste0("model_output/", model_name))
}

# Look at hyperparameters. sd should be somewhat smaller than mean,
# Stdev should not be very small. Range should not be larger than study extent.
model[["summary.hyperpar"]]

# ------------------------------------------------------------------------------------------------------------------------
# Effect plots 
# --------------------------------------------------------------------- ---------------------------------------------------
original_values <- data.frame(orig_values = unlist(all.seq)) %>% 
      mutate(covariate = sub("*.seq\\d+", "", rownames(.)),
             sequence = as.numeric(gsub("\\D", "", rownames(.))))

# Make effects data frame. Add constant equally to all derived values, to better visualize effects
# Fitted intercepts are related to average intensity of sampling, don't mean much ecologically.
# Value of intercept, which one is sensible? None of them, just pick one that gives good visualisation.
# Accept it doesn't mean much ecologically.

# Scale effect plots so that flat trends lie around the median
scale_params <- median(model$summary.lincomb.derived$`0.5quant`[grep("TZ_ann_rain_log|TZ_BG|TZ_dryspell|TZ_max_temp|TZ_HFP", 
                                                                     rownames(model$summary.lincomb.derived))])
effect_combs <- data.frame(covariate = gsub('[[:digit:]]+', '', sub("*_lc\\d+", "", rownames(model$summary.lincomb.derived))),
                           sequence = as.numeric(gsub("\\D", "", rownames(model$summary.lincomb.derived))),
                           quant_05 = inla.link.cloglog((model$summary.lincomb.derived$`0.5quant` - scale_params) - 0.36, inverse = TRUE),
                           quant_0025 = inla.link.cloglog((model$summary.lincomb.derived$`0.025quant` - scale_params) - 0.36, inverse = TRUE),
                           quant_0975 = inla.link.cloglog((model$summary.lincomb.derived$`0.975quant` - scale_params) - 0.36, inverse = TRUE))
effect_combs_main <- effect_combs %>% filter(covariate %in% c("TZ_ann_rain_log", "TZ_BG", "TZ_dryspell", "TZ_max_temp", "TZ_HFP"))
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
      geom_line(aes(x = orig_values, y = quant_05)) +
      geom_line(aes(x = orig_values, y = quant_0025), lty = 2, alpha = .5) +
      geom_line(aes(x = orig_values, y = quant_0975), lty = 2, alpha = .5) +
      # geom_rug(data = atlas_filtered, aes(x = )) +
      facet_wrap(~ covariate, scale = 'free_x', labeller = as_labeller(facet_labels)) +
      ggthemes::theme_few() +
      # theme(plot.title = element_text(hjust = 0.5)) +
      # ggtitle(paste0(species, " - Effect plots")) +
      xlab("Covariate value") +
      ylab("Probability of occurrence"); effects_plot

ggsave(plot = effects_plot, filename = paste0("figures/effects_", sub(" ", "_", species), "_r", 
                                              prior_range[1], "_", prior_range[2], "_s", 
                                              prior_sigma[1], "_", prior_sigma[2], 
                                              "_mean", fixed_mean, "_prec", fixed_precision,  
                                              ".png"))

# ------------------------------------------------------------------------------------------------------------------------
# Prediction plots
# ------------------------------------------------------------------------------------------------------------------------
pred.index <- inla.stack.index(stack = integrated_stack, tag = "pred")$data
lincomb.index.temp <- grep('TZ_lc_noTempMax', rownames(model$summary.lincomb.derived))
lincomb.index.rain <- grep('TZ_lc_noRain', rownames(model$summary.lincomb.derived))
lincomb.index.dry <- grep('TZ_lc_noDryspell', rownames(model$summary.lincomb.derived))
lincomb.index.BG <- grep('TZ_lc_noBG', rownames(model$summary.lincomb.derived))
lincomb.index.HFP <- grep('TZ_lc_noHFP', rownames(model$summary.lincomb.derived))
lincomb.index.allFixed <- grep('TZ_lc_allFixed', rownames(model$summary.lincomb.derived))

IP_df <- data.frame(SP_Points_data) %>% dplyr::select(date_index, x, y)

IP_df$pred_median <- model$summary.fitted.values[pred.index, "0.5quant"]
IP_df$pred_sd <- model$summary.fitted.values[pred.index, "sd"]
IP_df$no_temp_median <- model$summary.lincomb.derived[lincomb.index.temp, "0.5quant"]
IP_df$no_rain_median <- model$summary.lincomb.derived[lincomb.index.rain, "0.5quant"]
IP_df$no_dry_median <- model$summary.lincomb.derived[lincomb.index.dry, "0.5quant"]
IP_df$no_BG_median <- model$summary.lincomb.derived[lincomb.index.BG, "0.5quant"]
IP_df$no_HFP_median <- model$summary.lincomb.derived[lincomb.index.HFP, "0.5quant"]
IP_df$allFixed_median <- model$summary.lincomb.derived[lincomb.index.allFixed, "0.5quant"]

pred_data <- data.frame(all_pred = IP_df$pred_median,
                        # no_temp_median = IP_df$no_temp_median,
                        # no_rain_median = IP_df$no_rain_median,
                        # no_dry_median = IP_df$no_dry_median,
                        # no_BG_median = IP_df$no_BG_median,
                        # no_HFP_median = IP_df$no_HFP_median,
                        # allFixed_median = IP_df$allFixed_median,
                        sd = IP_df$pred_sd,
                        ind = rep(c(1,2), each = length(IP_df$pred_median)/2))

predcoordsGroup <- do.call(rbind, list(predcoords, predcoords))
pred_data_spdf <- sp::SpatialPixelsDataFrame(points = predcoordsGroup, 
                                             data = pred_data, 
                                             proj4string = proj)
test2 <- stack(pred_data_spdf)
median <- test2$all_pred
plot(median)
writeRaster(median, "testrast.tif")

pred_data_spdf <- crop(pred_data_spdf, TZ_no_lakes)
pred_data_spdf <- as(pred_data_spdf, "SpatialPixelsDataFrame")
test <- mask(pred_data_spdf, TZ_no_lakes)
pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 1] <- "1980-1999"
pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 2] <- "2000-2020"


median_plot <- ggplot() + 
      inlabru::gg(pred_data_spdf, aes(x = x, y = y, fill = inla.link.cloglog(all_pred, inv = TRUE))) + 
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("Median", na.value="white") +
      theme_void(); median_plot

sd_plot <- ggplot() + 
      gg(pred_data_spdf, aes(x = x, y = y, fill = sd)) + 
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("SD") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, vjust = 5)) 

# ------------------------------------------------------------------------------------------------------------------------
# Posterior median of the space-time random field 
# ------------------------------------------------------------------------------------------------------------------------
# Plot the posterior median of the space-time random field = the latent field (not directly observed)
# This shows us the variation in the spatial effect, as well as spatial dependence.

xmedian <- list()
for (j in 1:2) {
      xmedian[[j]] <- inla.mesh.project(
            projgrid,  model$summary.random$shared.field$`0.5quant`[index_set$shared.field.group == j])
}

xy.inGroup <- c(xy.in,xy.in)
xmedian <- unlist(xmedian)
xmedian <- xmedian[xy.inGroup]

dataObj <- data.frame(random_shared = xmedian,
                      ind = rep(c(1,2), each = length(xmedian)/2))

spatObj <- sp::SpatialPixelsDataFrame(points = predcoordsGroup, 
                                      data = dataObj, 
                                      proj4string = proj)
spatObj <- crop(spatObj, TZ_no_lakes)
spatObj@data[["ind"]][spatObj@data[["ind"]] == 1] <- "1980-1999"
spatObj@data[["ind"]][spatObj@data[["ind"]] == 2] <- "2000-2020"

# random_plot <- ggplot() + 
#       gg(spatObj, aes(x = x, y = y, fill = inla.link.cloglog(random_shared, inverse = T))) + 
#       facet_grid( ~ ind) +
#       coord_equal() +
#       viridis::scale_fill_viridis("Median") +
#       theme_void() 

xmedian2 <- list()
for (j in 1:2) {
      xmedian2[[j]] <- inla.mesh.project(
            projgrid,  model$summary.random$eBird.field$`0.5quant`[index_set_eBird$eBird.field.group == j])
}

xmedian2 <- unlist(xmedian2)
xmedian2 <- xmedian2[xy.inGroup]

dataObj2 <- data.frame(random_eBird = xmedian2,
                       ind = rep(c(1,2), each = length(xmedian2)/2))

spatObj2 <- sp::SpatialPixelsDataFrame(points = predcoordsGroup, 
                                       data = dataObj2, 
                                       proj4string = proj)
spatObj2 <- crop(spatObj2, TZ_no_lakes)
spatObj2@data[["ind"]][spatObj2@data[["ind"]] == 1] <- "1980-1999"
spatObj2@data[["ind"]][spatObj2@data[["ind"]] == 2] <- "2000-2020"

# random_plot_ebird <- ggplot() + 
#       gg(spatObj2, aes(x = x, y = y, fill = inla.link.cloglog(random_eBird, inverse = T))) + 
#       facet_grid( ~ ind) +
#       coord_equal() +
#       viridis::scale_fill_viridis("Median") +
#       theme_void() 
# 
# maps_combined <- gridExtra::grid.arrange(median_plot, random_plot)
# 
# ggsave(plot = maps_combined, filename =paste0("figures/maps_", sub(" ", "_", species), "_r",
#                                               prior_range[1], "_", prior_range[2], "_s",
#                                               prior_sigma[1], "_", prior_sigma[2],
#                                               "_mean", fixed_mean, "_prec", fixed_precision,
#                                               ".png"),
#        width = 20, height = 20, units = 'cm')

# ------------------------------------------------------------------------------------------------------------------------
# Make model output data  
# ------------------------------------------------------------------------------------------------------------------------
all_pred_df <- as.data.frame(pred_data_spdf)
random_eBird_df <- as.data.frame(spatObj2)
random_shared_df <- as.data.frame(spatObj)
all_merged_1 <- merge(all_pred_df, random_eBird_df, by = c("ind", "x", "y"))
all_merged <- merge(all_merged_1, random_shared_df, by = c("ind", "x", "y"))
all_merged$all_preds_est <- all_merged$random_eBird + all_merged$random_shared + all_merged$allFixed_median
all_merged$all_preds_no_rain <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_rain_median
all_merged$all_preds_no_temp <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_temp_median
all_merged$all_preds_no_dry <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_dry_median
all_merged$all_preds_no_BG <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_BG_median
all_merged$all_preds_no_HFP <- all_merged$random_eBird + all_merged$random_shared + all_merged$no_HFP_median
all_merged$all_random <- all_merged$random_eBird + all_merged$random_shared

all_merged[, 4:20] <- lapply(all_merged[, 4:20], inla.link.cloglog, inverse = TRUE) 


# ------------------------------------------------------------------------------------------------------------------------
# Diagnostic plot
# ------------------------------------------------------------------------------------------------------------------------
# Find out how much of the variation in linear predictions is explained by the random effects alone
rsq <- function(x, y) summary(lm(y~x))$r.squared
rsquared <- rsq(all_merged$all_preds_est, all_merged$all_random)
r2_random_linear <- ggplot(all_merged, aes(x = all_random, y = all_preds_est, col = ind)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, color = "darkred") +
      annotate("text", -Inf, Inf, hjust = -0.2, vjust = 3, label = paste0("R-squared: ", round(rsquared, digits = 3))) +
      annotate("text", -Inf, Inf, hjust = -0.25, vjust = 5, label = paste0("Range: ", round(model[["summary.hyperpar"]]$mean[1], digits = 3))) +
      annotate("text", -Inf, Inf, hjust = -0.25, vjust = 7, label = paste0("Sigma: ", round(model[["summary.hyperpar"]]$mean[2], digits = 3))) +
      facet_grid(~ind) +
      theme_few() +
      xlab("Random field estimate") +
      ylab("Linear prediction"); r2_random_linear

ggsave(plot = r2_random_linear, filename =paste0("figures/r2_randomLinear_", sub(" ", "_", species), "_r",
                                              prior_range[1], "_", prior_range[2], "_s",
                                              prior_sigma[1], "_", prior_sigma[2],
                                              "_mean", fixed_mean, "_prec", fixed_precision,
                                              ".png"),
       width = 20, height = 20, units = 'cm')

# ------------------------------------------------------------------------------------------------------------------------
#  R-squared of the full model
# ------------------------------------------------------------------------------------------------------------------------ 
all_obs_data <- as.data.frame(pred_data_spdf) %>% 
      dplyr::select(x, y, ind, all_pred) %>% 
      mutate(all_pred = inla.link.cloglog(all_pred, inv = TRUE))
atlas_pres <- as.data.frame(atlas_sp) %>% 
      filter(presence == 1)
atlas_pres$ind <- ifelse(atlas_pres$date_index == 1, "1980-1999", "2000-2020")

eBird_pres <- as.data.frame(ebird_sp) %>% 
      filter(presence == 1)
eBird_pres$ind <- ifelse(eBird_pres$date_index == 1, "1980-1999", "2000-2020")

median_plot <- ggplot() + 
      gg(pred_data_spdf, aes(x = x, y = y, fill = inla.link.cloglog(all_pred, inv = TRUE))) + 
      geom_point(data = atlas_pres,
                 aes(x = x, y = y, col = effort), alpha = 0.8, pch = 15, cex = 2) +
      geom_point(data = eBird_pres,
                 aes(x = x, y = y, col = duration_minutes), alpha = 0.6) +
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("Median", na.value="white") +
      scale_color_distiller(palette = "Reds", direction = 1) +
      theme_void(); median_plot

atlas_sp@data$ind <- atlas_sp@data$date_index
atlas_sp@data[["ind"]][atlas_sp@data[["ind"]] == 1] <- "1980-1999"
atlas_sp@data[["ind"]][atlas_sp@data[["ind"]] == 2] <- "2000-2020"
atlas_sp_simple <- atlas_sp[, names(atlas_sp) %in% c("presence", "ind", "x", 'y')]

ebird_sp@data$ind <- ebird_sp@data$date_index
ebird_sp@data[["ind"]][ebird_sp@data[["ind"]] == 1] <- "1980-1999"
ebird_sp@data[["ind"]][ebird_sp@data[["ind"]] == 2] <- "2000-2020"
ebird_sp_simple <- ebird_sp[, names(ebird_sp) %in% c("presence", "ind", "x", 'y')]

all_obs <- rbind(atlas_sp_simple, ebird_sp_simple)
all_obs_80s <- all_obs[all_obs@data[["ind"]] == "1980-1999", ]
all_obs_20s <- all_obs[all_obs@data[["ind"]] == "2000-2020", ]
      
pred_data_80s <- as.data.frame(pred_data_spdf) %>% filter(ind == '1980-1999') %>% mutate(all_pred = inla.link.cloglog(all_pred, inverse = T))
pred_data_20s <- as.data.frame(pred_data_spdf) %>% filter(ind == '2000-2020') %>% mutate(all_pred = inla.link.cloglog(all_pred, inverse = T))

pred_raster_80s <- rasterFromXYZ(pred_data_80s[, c("x", "y", "all_pred")])
pred_raster_20s <- rasterFromXYZ(pred_data_20s[, c("x", "y", "all_pred")])

ext_80s <- extract(pred_raster_80s, all_obs_80s)
ext_20s <- extract(pred_raster_20s, all_obs_20s)

all_obs_80s$pred <- ext_80s
all_obs_20s$pred <- ext_20s

plot(inla.link.cloglog(all_obs_80s$pred, inverse = T), all_obs_80s$presence)
plot(inla.link.cloglog(all_obs_20s$pred, inverse = T), all_obs_20s$presence)


# ------------------------------------------------------------------------------------------------------------------------
# Relative importance of covariates
# ------------------------------------------------------------------------------------------------------------------------
# Proportion of explained fixed effect variance
# 1- (cor(pred_data_spdf$no_temp_new, pred_data_spdf$fixed_effect))^2
1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_rain)
1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_temp)
1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_dry)
1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_BG)
1 - rsq(all_merged$all_preds_est, all_merged$all_preds_no_HFP)

# ------------------------------------------------------------------------------------------------------------------------
# Range change
# ------------------------------------------------------------------------------------------------------------------------
range_diff <- pred_data_spdf
dist_20s <- range_diff@data[["all_pred"]][range_diff@data[["ind"]] ==  "2000-2020"]
dist_80s <- range_diff@data[["all_pred"]][range_diff@data[["ind"]] ==  "1980-1999"]

colonisation <- (1-dist_80s) * dist_20s
extinction <- (dist_80s) * (1-dist_20s)
area_chance_transitions <-  sum((dist_80s) * (1-dist_80s))
