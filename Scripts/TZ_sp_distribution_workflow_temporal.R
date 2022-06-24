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

# Random field should not be too smooth (range > 20), and very different between time periods.
# Filter 1: Visually implausible
# Filter 2: DIC or CPO




# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Scripts')

#setwd('/home/ahomec/p/philism/Joris_work/scripts')
sapply(list.files(path = '/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_inla_spatial_temporal/source',
                  pattern="*.R", full.names = T), source, .GlobalEnv)

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')

#setwd('/home/ahomec/p/philism/Joris_work/Philip_data')
load("model_data/TZ_INLA_model_file_temporal.RData")

# ------------------------------------------------------------------------------------------------------------------------
# Global variables
# ------------------------------------------------------------------------------------------------------------------------
# Species choices
# Based on list of savannah bird species given in Beale et al. 2013. Species were excluded if there were less than 20 
# records in all date_index - data source combinations
species_list = c('Gyps africanus')
species <- species_list[1]

#species_list = c('Passer domesticus', 'Cisticola juncidis', 'Estrilda astrild', 'Histurgops ruficauda', 'Ploceus nigricollis', 
#                 'Cisticola brunnescens', 'Chrysococcyx cupreus', 'Tauraco hartlaubi', 'Ploceus castaneiceps', 'Nigrita canicapilla', 
#                 'Nectarinia kilimensis', 'Lanius collaris', 'Terpsiphone viridis', 'Oriolus auratus', 'Bubo capensis', 'Bubo africanus', 'Eremopterix leucopareia')

estimated_range = 3
max.edge = estimated_range/8

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
prior_range = c(1, 0.1)  
prior_sigma = c(6, 0.5)   

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
# 1. Data preparation
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
            data = data.frame(presence = ebird_filtered$occurrence, ebird_effort = ebird_filtered$ebird_effort,
                              duration_minutes = ebird_filtered$duration_minutes, number_observers = ebird_filtered$number_observers, date_index = ebird_filtered$date_index),
            proj4string = crs(proj))
      
      # Only include eBird data points for Tanzania
      # Get intersecting points
      in_sp <- rgeos::gIntersection(ebird_sp, TZ_outline)
      
      # Only keep intersecting points in original spdf
      ebird_sp <- ebird_sp[in_sp, ]
      
      atlas_full <- atlas_full %>%
            filter(time_period != 'x') %>%
            mutate(date_index = ifelse(time_period == '20s',2,1))
      
      atlas_full$Scientific[atlas_full$Scientific == 'Nectarinia olivacea'] <- 'Cyanomitra olivacea'
      atlas_full$Scientific[atlas_full$Scientific == 'Andropadus virens'] <- 'Eurillas virens'
      atlas_full$Scientific[atlas_full$Scientific == 'Gyps rueppellii'] <- 'Gyps africanus'
      atlas_full$Scientific[atlas_full$Scientific == 'Otis kori'] <- 'Ardeotis kori'
      atlas_full$Scientific[atlas_full$Scientific == 'Eupodotis ruficristata'] <- 'Eupodotis gindiana'
      atlas_full$Scientific[atlas_full$Scientific == 'Eupodotis senegalensis'] <- 'Lissotis melanogaster'
      atlas_full$Scientific[atlas_full$Scientific == 'Eupodotis hartlaubii'] <- 'Lissotis hartlaubii'
      atlas_full$Scientific[atlas_full$Scientific == 'Rhinoptilus africanus'] <- 'Smutsornis africanus'
      atlas_full$Scientific[atlas_full$Scientific == 'Coracias naevia'] <- 'Coracias naevius'
      atlas_full$Scientific[atlas_full$Scientific == 'Phoeniculus minor'] <- 'Rhinopomastus minor'
      atlas_full$Scientific[atlas_full$Scientific == 'Tockus erythrorhyncus'] <- 'Tockus erythrorhynchus'
      atlas_full$Scientific[atlas_full$Scientific == 'Lybius diadematus'] <- 'Tricholaema diademata'
      atlas_full$Scientific[atlas_full$Scientific == 'Lybius melanocephalus'] <- 'Tricholaema melanocephala'
      atlas_full$Scientific[atlas_full$Scientific == 'Mirafra poecilosterna'] <- 'Calendulauda poecilosterna'
      
      
     
      
      
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
      # run vulture again
      # 
      filtered_covs <- temporal_variables[, c('TZ_ann_rain_1980s_1.s', 'TZ_ann_rain_2000s_1.s', 
                                             'TZ_ann_rain_1980s_2.s', 'TZ_ann_rain_2000s_2.s', 
                                             'TZ_max_temp_1980s_1.s', 'TZ_max_temp_2000s_1.s',
                                             'TZ_max_temp_1980s_2.s', 'TZ_max_temp_2000s_2.s',
                                             'TZ_dryspell_1980s_1.s', 'TZ_dryspell_2000s_1.s',
                                             'TZ_dryspell_1980s_2.s', 'TZ_dryspell_2000s_2.s',
                                             'TZ_BG_1980s_1.s', 'TZ_BG_2000s_1.s',
                                             'TZ_BG_1980s_2.s', 'TZ_BG_2000s_2.s',
                                             'indicator_90s', 'indicator_2010s')]
            
      # Samples the covariates spatially closest to the occurrence data points    
      Nearest_covs_ebird <- GetNearestCovariate(ebird_sp, filtered_covs)  
      Nearest_covs_atlas <- GetNearestCovariate(atlas_sp, filtered_covs)
      # Use convolution layer for Atlas here? Would make more sense that just 
      # extracting centroid..
      
      # Add sampled covariates to occurrence data
      ebird_sp@data[, names(Nearest_covs_ebird@data)] <- Nearest_covs_ebird@data
      ebird_sp <- as(ebird_sp, 'data.frame')
      
      # Combine covariates from different time periods into single variable, different times identified by separate date_index  
      ebird_sp <- ebird_sp %>% mutate(annual_rain_1 = ifelse(date_index == 1, TZ_ann_rain_1980s_1.s, TZ_ann_rain_2000s_1.s),
                                      annual_rain_2 = ifelse(date_index == 1, TZ_ann_rain_1980s_2.s, TZ_ann_rain_2000s_2.s),
                                      hottest_temp_1 = ifelse(date_index == 1, TZ_max_temp_1980s_1.s, TZ_max_temp_2000s_1.s),
                                      hottest_temp_2 = ifelse(date_index == 1, TZ_max_temp_1980s_2.s, TZ_max_temp_2000s_2.s),
                                      max_dryspell_1 = ifelse(date_index == 1, TZ_dryspell_1980s_1.s, TZ_dryspell_2000s_1.s),
                                      max_dryspell_2 = ifelse(date_index == 1, TZ_dryspell_1980s_2.s, TZ_dryspell_2000s_2.s),
                                      BG_1 = ifelse(date_index == 1, TZ_BG_1980s_1.s, TZ_BG_2000s_1.s),
                                      BG_2 = ifelse(date_index == 1, TZ_BG_1980s_2.s, TZ_BG_2000s_2.s),
                                      indicator = ifelse(date_index == 1, indicator_90s, indicator_2010s))
      
      # Make spdf, add intercept, so model can distinguish eBird presence from Atlas presence
      ebird_sp <- SpatialPointsDataFrame(coords = ebird_sp[, c("x", "y")],
                                         data = ebird_sp,
                                         proj4string = crs(proj))
      ebird_sp@data[, 'ebird_intercept'] <- 1
      ebird_sp$presence <- as.numeric(ebird_sp$presence)
      
      atlas_sp@data[, names(Nearest_covs_atlas@data)] <- Nearest_covs_atlas@data
      atlas_sp <- as(atlas_sp, 'data.frame')
      
      atlas_sp <- atlas_sp %>% mutate(annual_rain_1 = ifelse(date_index == 1, TZ_ann_rain_1980s_1.s, TZ_ann_rain_2000s_1.s),
                                      annual_rain_2 = ifelse(date_index == 1, TZ_ann_rain_1980s_2.s, TZ_ann_rain_2000s_2.s),
                                      hottest_temp_1 = ifelse(date_index == 1, TZ_max_temp_1980s_1.s, TZ_max_temp_2000s_1.s),
                                      hottest_temp_2 = ifelse(date_index == 1, TZ_max_temp_1980s_2.s, TZ_max_temp_2000s_2.s),
                                      max_dryspell_1 = ifelse(date_index == 1, TZ_dryspell_1980s_1.s, TZ_dryspell_2000s_1.s),
                                      max_dryspell_2 = ifelse(date_index == 1, TZ_dryspell_1980s_2.s, TZ_dryspell_2000s_2.s),
                                      BG_1 = ifelse(date_index == 1, TZ_BG_1980s_1.s, TZ_BG_2000s_1.s),
                                      BG_2 = ifelse(date_index == 1, TZ_BG_1980s_2.s, TZ_BG_2000s_2.s),
                                      indicator = ifelse(date_index == 1, indicator_90s, indicator_2010s))
      
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
                       aes(x = x, y = y, col = "Absence"), alpha = 0.4, pch = 15, cex = 2) +
            geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 1 & as.data.frame(atlas_sp)$date_index == 1, ], 
                       aes(x = x, y = y, col = "Presence"), alpha = 0.4, pch = 15, cex = 3) +
            geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 0 & as.data.frame(ebird_sp)$date_index == 1, ],
                       aes(x = x, y = y, fill = "Absence"), pch = 21, alpha = 0.5) +
            geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 1 & as.data.frame(ebird_sp)$date_index == 1, ],
                       aes(x = x, y = y, fill = "Presence"), pch = 23, cex = 2) +
            theme_void() +
            coord_equal() +
            scale_fill_discrete(name = 'eBird') +
            scale_color_discrete(name = 'Atlas') +
            # ggtitle(paste0("1980-1999, eBird presence: ", length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 1, ]))))  +
            ggtitle("1980-1999") +
            theme(plot.title = element_text(hjust = 0.5))
      
      plot_2000s <- ggplot() +
            geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), 
                         colour = "black", fill = NA) +
            geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 0 & as.data.frame(atlas_sp)$date_index == 2, ], 
                       aes(x = x, y = y, col = "Absence"), alpha = 0.4, pch = 15, cex = 2) +
            geom_point(data = as.data.frame(atlas_sp)[as.data.frame(atlas_sp)$presence == 1 & as.data.frame(atlas_sp)$date_index == 2, ], 
                       aes(x = x, y = y, col = "Presence"), alpha = 0.4, pch = 15, cex = 3) +
            geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 0 & as.data.frame(ebird_sp)$date_index == 2, ], 
                       aes(x = x, y = y, fill = "Absence"), pch = 21, alpha = 0.5) +
            geom_point(data = as.data.frame(ebird_sp)[as.data.frame(ebird_sp)$presence == 1 & as.data.frame(ebird_sp)$date_index == 2, ], 
                       aes(x = x, y = y, fill = "Presence"), pch = 23, cex = 2) +
            theme_void() +
            coord_equal() +
            scale_fill_discrete(name = 'eBird') +
            scale_color_discrete(name = 'Atlas') +
            # ggtitle(paste0("2000-2020, eBird presence: ", length(rownames(ebird_sp@data[ebird_sp@data$presence == TRUE & as.data.frame(ebird_sp)$date_index == 2, ])))) +
            ggtitle("2000-2020") +
            theme(plot.title = element_text(hjust = 0.5))
      
      # 
      # raw_plots_comb <- gridExtra::grid.arrange(plot_80s, plot_2000s, nrow = 1, 
      #                                           top = text_grob(paste0(species, " - Occurrence Data"), size = 16))
      raw_plots_comb <- ggpubr::ggarrange(plot_80s, plot_2000s, nrow = 1, common.legend = TRUE, legend = "right")
      
      ggsave(plot = raw_plots_comb, filename = paste0("figures/", sub(" ", "_", species), "_raw_ccurrence.png"),
             width = 18, height = 18, units = 'cm')
}

# ------------------------------------------------------------------------------------------------------------------------
# 2. Mesh
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
      
      plot(Mesh$mesh)
      save(Mesh, file = paste0("model_data/mesh_E", max.edge,".RData"))
} else {
      load(paste0("model_data/mesh_E", max.edge,".RData"))
}

# ------------------------------------------------------------------------------------------------------------------------
# 3. SPDE model on the mesh
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
# 4. Space-time index set
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
# 5. Projection matrices
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
# 6. Estimation stacks
# ------------------------------------------------------------------------------------------------------------------------
stk.eBird <- inla.stack(data = list(resp = ebird_sp@data[, 'presence']),
                        A = list(1, projmat_eBird, projmat_eBird),
                        tag = 'eBird',
                        effects = list(ebird_sp@data, eBird.field = index_set_eBird, shared.field = index_set))
# stk.eBird <- inla.stack(data = list(resp = ebird_sp@data[, 'presence']),
#                         A = list(1, projmat_eBird),
#                         tag = 'eBird',
#                         effects = list(ebird_sp@data, shared.field = index_set))

stk.atlas <- inla.stack(data = list(resp = atlas_sp@data[, 'presence']),
                        A = list(1, projmat_atlas),
                        tag = 'atlas',
                        effects = list(atlas_sp@data, shared.field = index_set))


# ------------------------------------------------------------------------------------------------------------------------
# 7. Prediction data
# ------------------------------------------------------------------------------------------------------------------------
# Contains the locations and times where we want to make predictions. Code for this is adapted from the 
# 'MakeProjectionGrid' function.

if(!file.exists(paste0("model_data/", gsub(" ", "_", species), "_pred_files.RData"))){
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
      SP_Points_data[, 'annual_rain_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_ann_rain_1980s_1.s, NearestPredCovs@data$TZ_ann_rain_2000s_1.s)
      SP_Points_data[, 'annual_rain_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_ann_rain_1980s_2.s, NearestPredCovs@data$TZ_ann_rain_2000s_2.s)
      
      SP_Points_data[, 'hottest_temp_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_max_temp_1980s_1.s, NearestPredCovs@data$TZ_max_temp_2000s_1.s)
      SP_Points_data[, 'hottest_temp_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_max_temp_1980s_2.s, NearestPredCovs@data$TZ_max_temp_2000s_2.s)
      
      SP_Points_data[, 'max_dryspell_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_dryspell_1980s_1.s, NearestPredCovs@data$TZ_dryspell_2000s_1.s)
      SP_Points_data[, 'max_dryspell_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_dryspell_1980s_2.s, NearestPredCovs@data$TZ_dryspell_2000s_2.s)
      
      SP_Points_data[, 'BG_1'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_BG_1980s_1.s, NearestPredCovs@data$TZ_BG_2000s_1.s)
      SP_Points_data[, 'BG_2'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$TZ_BG_1980s_2.s, NearestPredCovs@data$TZ_BG_2000s_2.s)
      
      SP_Points_data[, 'indicator'] <- ifelse(SP_Points_data$date_index == 1, NearestPredCovs@data$indicator_90s, NearestPredCovs@data$indicator_2010s)
      
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
      
      save(SP_Points_data, projmat.pred, IP_sp, predcoords, projgrid, xy.in, file = paste0("model_data/", gsub(" ", "_", species), "_pred_files.RData"))
} else {
      load(paste0("model_data/", gsub(" ", "_", species), "_pred_files.RData"))
}

# ------------------------------------------------------------------------------------------------------------------------
# 8. Prediction stack
# ------------------------------------------------------------------------------------------------------------------------
stk.pred <- inla.stack(tag='pred',
                       data = list(resp = NA),
                       A = list(1, projmat.pred, projmat.pred),
                       effects = list(IP_sp@data,
                                      shared.field = index_set,
                                      eBird.field = index_set_eBird))

# 9. Integration stack
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
      NearestCovs@data[,colnames(NearestCovs@coords)] <- NearestCovs@coords
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
      
      save(stk.ip, file = paste0("model_data/stk_ip_E", max.edge,".RData"))
} else {
      load(paste0("model_data/stk_ip_E", max.edge,".RData"))
}



integrated_stack <- inla.stack(stk.eBird, stk.atlas, stk.pred, stk.ip)


# ------------------------------------------------------------------------------------------------------------------------
# 9. Model formula
# ------------------------------------------------------------------------------------------------------------------------
# Not sure if this is correct, but need to somehow add a copy of 'date_index', with different names, for form_1
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
h.spec <- list(rho = list(prior="pc.prec", param = c(3*sdres, 0.01)))

# At each time point, spatial locations are linked through the spde.
# Across time, the process evolves according to an AR(1) process.

# MODEL 1 -----------------------------------------------------------------------------------------------------------------
form_1 <- resp ~ 0 +
      ebird_intercept +
      atlas_intercept +
      ebird_effort + effort + 
      f(date_index, annual_rain_1, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +  
      f(date_index1, annual_rain_2, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +   
      f(date_index2, hottest_temp_1, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +   
      f(date_index3, hottest_temp_2, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +
      f(date_index4, max_dryspell_1, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +
      f(date_index5, max_dryspell_2, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +
      f(i, model = pcspde, group = i.group, control.group = list(model = 'ar1'))

# MODEL 2 -----------------------------------------------------------------------------------------------------------------
# f(annual_rain_1, model = "linear", mean.linear = 0, prec.linear = 0.001)
form_2 <- resp ~ 0 +
      ebird_intercept +
      atlas_intercept +
      annual_rain_1 + annual_rain_2 + hottest_temp_1 + hottest_temp_2 + 
      max_dryspell_1 + max_dryspell_2 + 
      BG_1 + BG_2 + indicator +
      x + y +
      duration_minutes + effort + 
      f(shared.field, model = pcspde, group = shared.field.group,  # In each year, spatial locations are linked by spde
        control.group = list(model = 'ar1', hyper = h.spec)) +  # Across time, process evolves according to AR(1) process
      f(eBird.field, model = eBird_spde, group = eBird.field.group,
        control.group = list(model = 'ar1', hyper = h.spec))
# ------------------------------------------------------------------------------------------------------------------------
# 10. INLA model
# ------------------------------------------------------------------------------------------------------------------------
# Set parameters on the default priors of fixed effects.
# If default priors produce too much shrinkage of coefficients
# towards zero, consider a larger precision, e.g. prec = 0.001^2.
# The default: INLA::inla.set.control.fixed.default()

# Add linear combinations to identify covariate importance
TZ_lc_noTempMax <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                      atlas_intercept = rep(1, NROW(SP_Points_data)),
                                      annual_rain_1 = SP_Points_data$annual_rain_1,
                                      annual_rain_2 = SP_Points_data$annual_rain_2,
                                      max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                      max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                      BG_1 = SP_Points_data$BG_1,
                                      BG_2 = SP_Points_data$BG_2,
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
                                   BG_1 = SP_Points_data$BG_1,
                                   BG_2 = SP_Points_data$BG_2,
                                   x = SP_Points_data$x,
                                   y = SP_Points_data$y,
                                   indicator = SP_Points_data$indicator)
names(TZ_lc_noRain) <- paste0("TZ_lc_noRain", 1:NROW(SP_Points_data))
                              
TZ_lc_noDryspell <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                   atlas_intercept = rep(1, NROW(SP_Points_data)),
                                   annual_rain_1 = SP_Points_data$annual_rain_1,
                                   annual_rain_2 = SP_Points_data$annual_rain_2,
                                   hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                   hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                   BG_1 = SP_Points_data$BG_1,
                                   BG_2 = SP_Points_data$BG_2,
                                   x = SP_Points_data$x,
                                   y = SP_Points_data$y,
                                   indicator = SP_Points_data$indicator)
names(TZ_lc_noDryspell) <- paste0("TZ_lc_noDryspell", 1:NROW(SP_Points_data))

TZ_lc_noBG <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
                                       atlas_intercept = rep(1, NROW(SP_Points_data)),
                                       annual_rain_1 = SP_Points_data$annual_rain_1,
                                       annual_rain_2 = SP_Points_data$annual_rain_2,
                                       max_dryspell_1 =SP_Points_data$max_dryspell_1,
                                       max_dryspell_2 = SP_Points_data$max_dryspell_2,
                                       hottest_temp_1 = SP_Points_data$hottest_temp_1,
                                       hottest_temp_2 = SP_Points_data$hottest_temp_2,
                                 x = SP_Points_data$x,
                                 y = SP_Points_data$y)
names(TZ_lc_noBG) <- paste0("TZ_lc_noBG", 1:NROW(SP_Points_data))
# TZ_lc_noEffort <- inla.make.lincombs(ebird_intercept = rep(1, NROW(SP_Points_data)),
#                                        atlas_intercept = rep(1, NROW(SP_Points_data)),
#                                        annual_rain_1 = SP_Points_data$annual_rain_1,
#                                        annual_rain_2 = SP_Points_data$annual_rain_2,
#                                        max_dryspell_1 =SP_Points_data$max_dryspell_1,
#                                        max_dryspell_2 = SP_Points_data$max_dryspell_2,
#                                        hottest_temp_1 = SP_Points_data$hottest_temp_1,
#                                        hottest_temp_2 = SP_Points_data$hottest_temp_2,
#                                        BG_1 = SP_Points_data$BG_1,
#                                        BG_2 = SP_Points_data$BG_2,
#                                        x = SP_Points_data$x,
#                                        y = SP_Points_data$y,
#                                        indicator = SP_Points_data$indicator)
# names(TZ_lc_noEffort) <- paste0("TZ_lc_noEffort", 1:NROW(SP_Points_data))


lc_combined <- c(all_lc, TZ_lc_noDryspell, TZ_lc_noTempMax, TZ_lc_noBG, TZ_lc_noRain)

# 1: Get into same spatial points shape as random effect and full prediction
# 2: Add random effect to those spatial points

# Gaussian priors
C.F. <- list(
      mean = fixed_mean,
      prec = fixed_precision   # Precision for all fixed effects except intercept
)

if(!file.exists(paste0("model_output/", model_name))){
      # lc_no_BG <- all_lc[grepl("BG", all_lc) == FALSE]
      model <- inla(form_2, family = "binomial", control.family = list(link = "cloglog"), 
                    lincomb = lc_combined,
                    data = inla.stack.data(integrated_stack), 
                    verbose = FALSE,
                    control.predictor = list(A = inla.stack.A(integrated_stack), 
                                             link = NULL, compute = TRUE), 
                    control.fixed = C.F.,
                    E = inla.stack.data(integrated_stack)$e, 
                    control.compute = list(waic = FALSE, dic = TRUE, cpo = TRUE))
      saveRDS(model, file = paste0("model_output/", model_name))
} else {
      model <- readRDS(paste0("model_output/", model_name))
}

# Look at hyperparameters. sd should be somewhat smaller than mean,
# Stdev should not be very small.
model[["summary.hyperpar"]]

# ------------------------------------------------------------------------------------------------------------------------
# 11. Effect plots 
# ------------------------------------------------------------------------------------------------------------------------
original_values <- data.frame(orig_values = unlist(all.seq)) %>% 
      mutate(covariate = sub("*.seq\\d+", "", rownames(.)),
             sequence = as.numeric(gsub("\\D", "", rownames(.))))

# Make effects data frame. Add constant equally to all derived values, to better visualize effects
# Fitted intercepts are related to average intensity of sampling, don't mean much ecologically.
# Value of intercept, which one is sensible? None of them, just pick one that gives good visualisation.
# Accept it doesn't mean much ecologically, pick a number that makes the plots look nice.


scale_params <- median(model$summary.lincomb.derived$`0.5quant`[grep("TZ_ann_rain|TZ_BG|TZ_dryspell|TZ_max_temp", rownames(model$summary.lincomb.derived))])
effect_combs <- data.frame(covariate = gsub('[[:digit:]]+', '', sub("*_lc\\d+", "", rownames(model$summary.lincomb.derived))),
                           sequence = as.numeric(gsub("\\D", "", rownames(model$summary.lincomb.derived))),
                           quant_05 = inla.link.cloglog((model$summary.lincomb.derived$`0.5quant` - scale_params) - 0.36, inverse = TRUE),
                           quant_0025 = inla.link.cloglog((model$summary.lincomb.derived$`0.025quant` - scale_params) - 0.36, inverse = TRUE),
                           quant_0975 = inla.link.cloglog((model$summary.lincomb.derived$`0.975quant` - scale_params) - 0.36, inverse = TRUE))
effect_combs_main <- effect_combs %>% filter(covariate %in% c("TZ_ann_rain", "TZ_BG", "TZ_dryspell", "TZ_max_temp"))
effect_combs_m <- merge(original_values, effect_combs_main)

facet_labels <- c(
      TZ_ann_rain = "Annual rainfall",
      TZ_BG = "Bareground cover",
      TZ_dryspell = "Longest dryspell duration",
      TZ_max_temp = "Hottest temperature"
)

# Make effect plots, using linear combinations
effects_plot <- ggplot(effect_combs_m) +
      geom_line(aes(x = orig_values, y = quant_05)) +
      geom_line(aes(x = orig_values, y = quant_0025), lty = 2, alpha = .5) +
      geom_line(aes(x = orig_values, y = quant_0975), lty = 2, alpha = .5) +
      # geom_rug(data = atlas_filtered, aes(x = )) +
      facet_wrap(~ covariate, scale = 'free_x', labeller = as_labeller(facet_labels)) +
      theme_few() +
      theme(plot.title = element_text(hjust = 0.5)) +
      ggtitle(paste0(species, " - Effect plots")) +
      xlab("Covariate value") +
      ylab("Probability of occurrence"); effects_plot

ggsave(plot = effects_plot, filename = paste0("figures/effects_", sub(" ", "_", species), "_r", 
                                              prior_range[1], "_", prior_range[2], "_s", 
                                              prior_sigma[1], "_", prior_sigma[2], 
                                              "_mean", fixed_mean, "_prec", fixed_precision,  
                                              ".png"))

# ------------------------------------------------------------------------------------------------------------------------
# 12. Prediction plots
# ------------------------------------------------------------------------------------------------------------------------
pred.index <- inla.stack.index(stack = integrated_stack, tag = "pred")$data
lincomb.index.temp <- grep('TZ_lc_noTempMax', rownames(model$summary.lincomb.derived))
lincomb.index.rain <- grep('TZ_lc_noRain', rownames(model$summary.lincomb.derived))
lincomb.index.dry <- grep('TZ_lc_noDryspell', rownames(model$summary.lincomb.derived))
lincomb.index.BG <- grep('TZ_lc_noBG', rownames(model$summary.lincomb.derived))

IP_df <- data.frame(IP_sp) %>% dplyr::select(date_index, x, y)

IP_df$pred_median <- model$summary.fitted.values[pred.index, "0.5quant"]
IP_df$pred_median_P <- inla.link.cloglog(IP_df$pred_median, inverse = TRUE)

IP_df$pred_sd <- model$summary.fitted.values[pred.index, "sd"]

IP_df$no_temp_median <- model$summary.lincomb.derived[lincomb.index.temp, "0.5quant"]
IP_df$no_rain_median <- model$summary.lincomb.derived[lincomb.index.rain, "0.5quant"]
IP_df$no_dry_median <- model$summary.lincomb.derived[lincomb.index.dry, "0.5quant"]
IP_df$no_BG_median <- model$summary.lincomb.derived[lincomb.index.BG, "0.5quant"]

# IP_df$pred_ll <- model$summary.fitted.values[pred.index, "0.025quant"]
# IP_df$pred_ul <- model$summary.fitted.values[pred.index, "0.975quant"]

pred_data <- data.frame(median = IP_df$pred_median,
                        median_P = IP_df$pred_median_P,
                        no_temp_median = IP_df$no_temp_median,
                        no_rain_median = IP_df$no_rain_median,
                        no_dry_median = IP_df$no_dry_median,
                        no_BG_median = IP_df$no_BG_median,
                        sd = IP_df$pred_sd,
                        ind = rep(c(1,2), each = length(IP_df$pred_median)/2))

predcoordsGroup <- do.call(rbind, list(predcoords, predcoords))
pred_data_spdf <- sp::SpatialPixelsDataFrame(points = predcoordsGroup, 
                                             data = pred_data, 
                                             proj4string = proj)

pred_data_spdf <- crop(pred_data_spdf, TZ_no_lakes)
pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 1] <- "1980-1999"
pred_data_spdf@data[["ind"]][pred_data_spdf@data[["ind"]] == 2] <- "2000-2020"

median_plot <- ggplot() + 
      gg(pred_data_spdf, aes(x = x, y = y, fill = median_P)) + 
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("Median") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, vjust = 5)) 

sd_plot <- ggplot() + 
      gg(pred_data_spdf, aes(x = x, y = y, fill = sd)) + 
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("SD") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, vjust = 5)) 

# ------------------------------------------------------------------------------------------------------------------------
# 13. Posterior median of the space-time random field 
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

dataObj <- data.frame(median = xmedian,
                      ind = rep(c(1,2), each = length(xmedian)/2))

spatObj <- sp::SpatialPixelsDataFrame(points = predcoordsGroup, 
                                      data = dataObj, 
                                      proj4string = proj)
spatObj <- crop(spatObj, TZ_no_lakes)
spatObj@data[["ind"]][spatObj@data[["ind"]] == 1] <- "1980-1999"
spatObj@data[["ind"]][spatObj@data[["ind"]] == 2] <- "2000-2020"

random_plot <- ggplot() + 
      gg(spatObj, aes(x = x, y = y, fill = median)) + 
      facet_grid( ~ ind) +
      coord_equal() +
      viridis::scale_fill_viridis("Median") +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, vjust = 5)) +
      ggtitle(paste0(species, " - Random field"))

maps_combined <- gridExtra::grid.arrange(median_plot, sd_plot, random_plot, 
                                         top = ggpubr::text_grob(species))

ggsave(plot = maps_combined, filename =paste0("figures/maps_", sub(" ", "_", species), "_r",
                                              prior_range[1], "_", prior_range[2], "_s",
                                              prior_sigma[1], "_", prior_sigma[2],
                                              "_mean", fixed_mean, "_prec", fixed_precision,
                                              ".png"),
       width = 20, height = 20, units = 'cm')

# ------------------------------------------------------------------------------------------------------------------------
# 14. Diagnostic plot
# ------------------------------------------------------------------------------------------------------------------------
# Find out how much of the variation in linear predictions is explained by the random effects alone
pred_data_spdf$random <- spatObj$median
pred_data_spdf$fixed_effect <- pred_data_spdf$median - pred_data_spdf$random

rsq <- function(x, y) summary(lm(y~x))$r.squared
rsquared <- rsq(pred_data_spdf$median, pred_data_spdf$random)
r2_random_linear <- ggplot(as.data.frame(pred_data_spdf), aes(x = random, y = median)) +
      geom_point() +
      geom_smooth(method = "lm", se = FALSE, color = "darkred") +
      annotate("text", -Inf, Inf, hjust = -0.2, vjust = 3, label = paste0("R-squared: ", round(rsquared, digits = 3))) +
      annotate("text", -Inf, Inf, hjust = -0.25, vjust = 5, label = paste0("Range: ", round(model[["summary.hyperpar"]]$mean[1], digits = 3))) +
      annotate("text", -Inf, Inf, hjust = -0.25, vjust = 7, label = paste0("Sigma: ", round(model[["summary.hyperpar"]]$mean[2], digits = 3))) +
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
# 15. Relative importance of covariates
# ------------------------------------------------------------------------------------------------------------------------
# Get the total predictions for the model without a given covariate
pred_data_spdf$no_temp_total <- pred_data_spdf$no_temp_median + pred_data_spdf$random
pred_data_spdf$no_rain_total <- pred_data_spdf$no_rain_median + pred_data_spdf$random
pred_data_spdf$no_dry_total <- pred_data_spdf$no_dry_median + pred_data_spdf$random
pred_data_spdf$no_BG_total <- pred_data_spdf$no_BG_median + pred_data_spdf$random
spplot(pred_data_spdf)

# Overall covariate importance
# 1- (cor(pred_data_spdf$no_rain_total, pred_data_spdf$median))^2
1 - rsq(pred_data_spdf$median, pred_data_spdf$no_rain_total)
1 - rsq(pred_data_spdf$median, pred_data_spdf$no_temp_total)
1 - rsq(pred_data_spdf$median, pred_data_spdf$no_dry_total)
1 - rsq(pred_data_spdf$median, pred_data_spdf$no_BG_total)

# Proportion of explained fixed effect variance
# 1- (cor(pred_data_spdf$no_temp_new, pred_data_spdf$fixed_effect))^2
1- rsq(pred_data_spdf$fixed_effect, pred_data_spdf$no_rain_median)
1- rsq(pred_data_spdf$fixed_effect, pred_data_spdf$no_temp_median)
1- rsq(pred_data_spdf$fixed_effect, pred_data_spdf$no_dry_median)
1- rsq(pred_data_spdf$fixed_effect, pred_data_spdf$no_BG_median)

# What happens with effort?
# Try linear combination without effort, check if that equals fixed effects computed by
# subtracting random effects from fitted (not transformed to probability)

# ------------------------------------------------------------------------------------------------------------------------
# 16. Range change
# ------------------------------------------------------------------------------------------------------------------------
range_diff <- pred_data_spdf
#
dist_20s <- range_diff@data[["median_P"]][range_diff@data[["ind"]] ==  "2000-2020"]
dist_80s <- range_diff@data[["median_P"]][range_diff@data[["ind"]] ==  "1980-1999"]
range_diff@data[["median_P"]] <- dist_20s - dist_80s
range_diff@data[["colonisation"]] <- (1-dist_80s) * dist_20s # probability of colonisation
range_diff@data[["extinction"]] <- (dist_80s) * (1-dist_20s) # probability of extinction
range_diff@data[["con_absence"]] <- (1-dist_80s) * (1-dist_20s) # probability of continued absence
range_diff@data[["con_presence"]] <- (dist_80s) * (dist_20s) # probability of continued presence

# net range expansion:
range_exp <- sum(dist_20s - dist_80s, na.rm = TRUE) # expected number of pixels that have changed

col_plot <- ggplot() +
    gg(data = range_diff, aes(x = x, y = y, fill = colonisation)) +
    geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), fill = NA, col = "black") +
    coord_equal() +
    viridis::scale_fill_viridis("P(Colonisation)", limits = c(0.5, 1), na.value="white") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 5)); col_plot

ext_plot <- ggplot() +
    gg(range_diff, aes(x = x, y = y, fill = extinction)) +
    geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), fill = NA, col = "black") +
    coord_equal() +
    viridis::scale_fill_viridis("P(Extinction)", na.value="white") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 5))

con_absence_plot <- ggplot() +
    gg(range_diff, aes(x = x, y = y, fill = con_absence)) +
    geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), fill = NA, col = "black") +
    coord_equal() +
    viridis::scale_fill_viridis("P(Continued absence)", limits = c(0.5, 1), na.value="white") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 5))

con_presence_plot <- ggplot() +
    gg(range_diff, aes(x = x, y = y, fill = con_presence)) +
    geom_polygon(data = TZ_no_lakes, aes(x = long, y = lat, group = group), fill = NA, col = "black") +
    coord_equal() +
    viridis::scale_fill_viridis("P(Continued presence)", limits = c(0.5, 1), na.value="white") +
    theme_void() +
    theme(plot.title = element_text(hjust = 0.5, vjust = 5)); con_presence_plot

range_plots_comb <- gridExtra::grid.arrange(col_plot, con_absence_plot, con_presence_plot,
                                          nrow = 1)

ggsave(plot = range_plots_comb, filename =paste0("figures/difference_", sub(" ", "_", species), "_r",
                                               prior_range[1], "_", prior_range[2], "_s",
                                               prior_sigma[1], "_", prior_sigma[2],
                                               "_mean", fixed_mean, "_prec", fixed_precision,
                                               ".png"),
     width = 18, height = 18, units = 'cm')

