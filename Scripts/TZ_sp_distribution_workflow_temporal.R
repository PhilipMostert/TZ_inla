args = commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages(library(INLA, quietly=TRUE))

library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(rlang)
library(inlabru)
library(sf)
library(reshape2)

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Scripts')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_inla_spatial_temporal/source/')

#setwd('/home/ahomec/p/philism/Joris_work/scripts')
sapply(list.files(pattern="*.R"),source,.GlobalEnv)

species_list = c('Cisticola juncidis', 'Eremopterix leucopareia', 'Estrilda astrild', 'Histurgops ruficauda')
species <- species_list[1]
#species_list = c('Passer domesticus', 'Cisticola juncidis', 'Estrilda astrild', 'Histurgops ruficauda', 'Ploceus nigricollis', 
#                 'Cisticola brunnescens', 'Chrysococcyx cupreus', 'Tauraco hartlaubi', 'Ploceus castaneiceps', 'Nigrita canicapilla', 
#                 'Nectarinia kilimensis', 'Lanius collaris', 'Terpsiphone viridis', 'Oriolus auratus', 'Bubo capensis', 'Bubo africanus', 'Eremopterix leucopareia')

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/data_processed')

#setwd('/home/ahomec/p/philism/Joris_work/Philip_data')
estimated_range = 2
max.edge = estimated_range/8
load(paste0("TZ_INLA_model_file_temporal_E", round(max.edge, digits = 3), ".RData"))

##Set time 1960s -- 1990s as time 1;
##Set time 2000s as time 2.
n_time_layers = 2


# ------------------------------------------------------------------------------------------------------------------------
# 1. Data preparation
# ------------------------------------------------------------------------------------------------------------------------

ebird_full <- ebird_full %>%
      mutate(date_index = ifelse(date > '2000-01-01',2,1))

ebird_filtered <- ebird_full %>% 
      filter(APPROVED == 1,  # Only keep reviewed and approved records
             `ALL SPECIES REPORTED` == 1,           # Only keep complete checklists
             `EFFORT DISTANCE KM` < 15,
             `DURATION MINUTES` >= 5,
             `DURATION MINUTES` <= 240)   

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
             number_observers = `NUMBER OBSERVERS`)

ebird_sp <- SpatialPointsDataFrame(
      coords = ebird_filtered[, c("LONGITUDE", "LATITUDE")],
      data = data.frame(presence = ebird_filtered$occurrence, duration_minutes = ebird_filtered$duration_minutes,
                        effort_distance_km = ebird_filtered$effort_distance_km, number_observers = ebird_filtered$number_observers, date_index = ebird_filtered$date_index),
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
      mutate(presence = ifelse(occurrence == 1, TRUE, FALSE)) %>% 
      dplyr::select(-V1); if(is_empty(atlas_filtered$presence)){print("ERROR: No Atlas data available")}

atlas_sp <- SpatialPointsDataFrame(
      coords = data.frame(atlas_filtered[, c("Long", "Lat")]),
      data = data.frame(presence = atlas_filtered$presence, effort = atlas_filtered$effort,
                        date_index = atlas_filtered$date_index),
      proj4string = crs(proj))

range01 <- function(x){(x - min(x))/(max(x) - min(x))}
ebird_sp$duration_minutes <- range01(ebird_sp$duration_minutes)
atlas_sp$effort <- range01(atlas_sp$effort)

filtered_covs <- temporal_variables[,c('TZ_ann_rain_1980s_1.s', 'TZ_ann_rain_2000s_1.s', 
                                       'TZ_ann_rain_1980s_2.s', 'TZ_ann_rain_2000s_2.s', 
                                       'TZ_max_temp_1980s_1.s', 'TZ_max_temp_2000s_1.s',
                                       'TZ_max_temp_1980s_2.s', 'TZ_max_temp_2000s_2.s',
                                       'TZ_dryspell_1980s_1.s', 'TZ_dryspell_2000s_1.s',
                                       'TZ_dryspell_1980s_2.s', 'TZ_dryspell_2000s_2.s',
                                       'TZ_BG_1980s_1.s', 'TZ_BG_2000s_1.s',
                                       'TZ_BG_1980s_2.s', 'TZ_BG_2000s_2.s',
                                       'indicator_90s', 'indicator_2010s')]
calc_covs <- TRUE
if (calc_covs) {
      
      # Samples the covariates spatially closest to the occurrence data points    
      Nearest_covs_ebird <- GetNearestCovariate(ebird_sp,filtered_covs)  
      Nearest_covs_atlas <- GetNearestCovariate(atlas_sp, filtered_covs)
      
      
} else {
      
      #Nearest_covs_ebird <- readRDS('/Users/philism/Downloads/Nearest_covs_ebird.RDS')
      #Nearest_covs_atlas <-  readRDS('/Users/philism/Downloads/Nearest_covs_atlas.RDS')
      setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Covariate_data')
      Nearest_covs_ebird <- readRDS("Nearest_covs_ebird.RDS")
      Nearest_covs_atlas <- readRDS("Nearest_covs_atlas.RDS")
      # Add covariates to the bird data 
      ebird_sp@data[, names(Nearest_covs_ebird@data)] <- Nearest_covs_ebird@data
      atlas_sp@data[, names(Nearest_covs_atlas@data)] <- Nearest_covs_atlas@data
      
}

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
ebird_sp <- SpatialPointsDataFrame(coords = ebird_sp[, c("LONGITUDE", "LATITUDE")],
                                   data = ebird_sp[, !names(ebird_sp)%in%c('LONGITUDE', 'LATITUDE')],
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

atlas_sp <- SpatialPointsDataFrame(coords = atlas_sp[, c('Long', 'Lat')],
                                   data = atlas_sp[, !names(atlas_sp)%in%c('Long','Lat')],
                                   proj4string = crs(proj))
atlas_sp@data[,'atlas_intercept'] <- 1
atlas_sp$presence <- as.numeric(atlas_sp$presence)

# ------------------------------------------------------------------------------------------------------------------------
# 2. Mesh
# ------------------------------------------------------------------------------------------------------------------------

# Max.edge based on an estimated range
estimated_range = 2
max.edge = estimated_range/8

# Create the mesh 
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

# ------------------------------------------------------------------------------------------------------------------------
# 3. SPDE model on the mesh
# ------------------------------------------------------------------------------------------------------------------------

# spde <- inla.spde2.matern(mesh = Mesh$mesh)

pcspde <- inla.spde2.pcmatern(
      mesh = Mesh$mesh,
      alpha = 2,
      prior.range = c(5, 0.01),   
      prior.sigma = c(2, 0.01))


# ------------------------------------------------------------------------------------------------------------------------
# 4. Space-time index set
# ------------------------------------------------------------------------------------------------------------------------

# Set index for the latent field, taking into account the number of mesh points in the SPDE model, 
# and the number of time groups.
# This index does not depend on any data set locations, only SPDE model size and time dimensions.
# A double index, identifying both the spatial location and associated time point. 1 to n.spde index points 
# for n.group times.

index_set <- inla.spde.make.index(name ='i',
                            n.spde = pcspde$n.spde,
                            n.group = n_time_layers)
lengths(index_set)


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
                        A = list(1, projmat_eBird), 
                        tag = 'eBird',
                        effects = list(ebird_sp@data, i = index_set))

stk.atlas <- inla.stack(data = list(resp = atlas_sp@data[, 'presence']),
                        A = list(1, projmat_atlas),
                        tag = 'atlas',
                        effects = list(atlas_sp@data, i = index_set))


# ------------------------------------------------------------------------------------------------------------------------
# 7. Prediction data
# ------------------------------------------------------------------------------------------------------------------------

# Contains the locations and times where we want to make predictions. Code for this is adapted from the 
# 'MakeProjectionGrid' function.

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
colnames(predcoords) <- c('Long','Lat')
Apred <- projgrid$proj$A[which(xy.in), ]


# Construct the prediction data, which is spatial points by temporal layers
spatial_points <- cbind(c(Mesh$mesh$loc[,1]), c(Mesh$mesh$loc[,2]))
colnames(spatial_points) <- c("LONGITUDE", "LATITUDE")

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

IP_sp <- sp::SpatialPointsDataFrame(coords = spatiotemporal_points, data = SP_Points_data, proj4string = proj)
IP_sp@data$Intercept <- 1
IP_sp@data[, c("LONGITUDE", "LATITUDE")] <- IP_sp@coords

projmat.pred <- inla.spde.make.A(mesh = Mesh$mesh, loc = IP_sp@coords, group = IP_sp@data$date_index)

# ------------------------------------------------------------------------------------------------------------------------
# 8. Prediction stack
# ------------------------------------------------------------------------------------------------------------------------

stk.pred <- inla.stack(tag='pred',
                       data = list(resp = NA), 
                       A = list(1, projmat.pred), 
                       effects = list(IP_sp@data, 
                                      i = index_set))
# stk.pred[["effects"]][["data"]][["i.group"]] <- stk.pred[["effects"]][["data"]][["date_index"]]

integrated_stack <- inla.stack(stk.eBird, stk.atlas, stk.pred)


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

# Define a PC-prior for the temporal autoregressive parameter.
h.spec <- list(rho = list(prior = 'pc.prec', param = c(4,0.01)))

# At each time point, spatial locations are linked through the spde.
# Across time, the process evolves according to an AR(1) process.

# MODEL 1 -----------------------------------------------------------------------------------------------------------------
form_1 <- resp ~ 0 +
      ebird_intercept +
      atlas_intercept +
      duration_minutes + effort + 
      f(date_index, annual_rain_1, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +  
      f(date_index1, annual_rain_2, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +   
      f(date_index2, hottest_temp_1, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +   
      f(date_index3, hottest_temp_2, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +
      f(date_index4, max_dryspell_1, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +
      f(date_index5, max_dryspell_2, model="rw1", scale.model=TRUE, constr=FALSE, hyper = h.spec) +
      f(i, model = pcspde, group = i.group, control.group = list(model = 'ar1'))

# MODEL 2 -----------------------------------------------------------------------------------------------------------------
form_2 <- resp ~ 0 +
      ebird_intercept +
      atlas_intercept +
      duration_minutes + effort + 
      annual_rain_1 + annual_rain_2 + hottest_temp_1 + hottest_temp_2 + 
      max_dryspell_1 + max_dryspell_2 + BG_1 + BG_2 +
      date_index + indicator +
      f(i, model = pcspde, group = i.group, control.group = list(model = 'ar1', hyper = h.spec))


# ------------------------------------------------------------------------------------------------------------------------
# 10. INLA model
# ------------------------------------------------------------------------------------------------------------------------

model <- inla(form_2, family = "binomial", control.family = list(link = "cloglog"), # Backtransform for probability scale
              lincomb = all_lc,
              data = inla.stack.data(integrated_stack), 
              verbose = FALSE,
              control.predictor = list(A = inla.stack.A(integrated_stack), 
                                       link = NULL, compute = TRUE), 
              E = inla.stack.data(integrated_stack)$e, 
              control.compute = list(waic = FALSE, dic = FALSE, cpo = FALSE))

setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_inla_spatial_temporal/model_output')
saveRDS(model, file = "model_all_lc_GAM_form2.RDS")
# model <- readRDS('model_all_lc_GAM_form2.RDS')

# ------------------------------------------------------------------------------------------------------------------------
# 11. Effect plots 
# ------------------------------------------------------------------------------------------------------------------------

# for (i in 1:length(res.bits$marginals.fixed)) {
#       tmp = inla.tmarginal(function(x) x, res.bits$marginals.fixed[[i]]) ## not sure how to fix?
#       plot(tmp, type = "l", xlab = paste("Fixed effect marginal", i, ":", names(res.bits[["marginals.fixed"]])[i]), ylab = "Density")
#       abline(v = 0, lty = 2)
# }

# Collect model results
res.bits <- list("summary.lincomb" = model$summary.lincomb, 
                 "summary.lincomb.derived" = model$summary.lincomb.derived,
                 "marginals.lincomb.derived" = model$marginals.lincomb.derived,
                 "summary.fixed"  = model$summary.fixed,
                 "summary.hyperpar" = model$summary.hyperpar,
                 "marginals.fixed" = model$marginals.fixed)


# Make effect plots, using linear combinations
par(mfrow = c(2,2))
plot(all.seq$TZ_ann_rain.seq, cloglog_inv(res.bits$summary.lincomb.derived$`0.5quant`[grep("rain", rownames(res.bits$summary.lincomb.derived))]), 
     lwd = 2 , type = 'l', main = 'Annual rainfall', ylab = '')
plot(all.seq$TZ_max_temp.seq, cloglog_inv(res.bits$summary.lincomb.derived$`0.5quant`[grep("temp", rownames(res.bits$summary.lincomb.derived))]), 
     lwd = 2 , type = 'l', main = 'Maximum temperature', ylab = '')
plot(all.seq$TZ_dryspell.seq, cloglog_inv(res.bits$summary.lincomb.derived$`0.5quant`[grep("dry", rownames(res.bits$summary.lincomb.derived))]), 
     lwd = 2 , type = 'l', main = 'Dryspell duration', ylab = '')
plot(all.seq$TZ_dryspell.seq, cloglog_inv(res.bits$summary.lincomb.derived$`0.5quant`[grep("BG", rownames(res.bits$summary.lincomb.derived))]), 
     lwd = 2 , type = 'l', main = 'Bareground cover', ylab = '')
par(mfrow = c(1,1))

# ------------------------------------------------------------------------------------------------------------------------
# 12. Prediction plots
# ------------------------------------------------------------------------------------------------------------------------
pred.index <- inla.stack.index(stack = integrated_stack, tag = "pred")$data

IP_df <- data.frame(IP_sp) %>% dplyr::select(date_index, LONGITUDE, LATITUDE)

IP_df$pred_mean <- model$summary.fitted.values[pred.index, "mean"]
IP_df$pred_ll <- model$summary.fitted.values[pred.index, "0.025quant"]
IP_df$pred_ul <- model$summary.fitted.values[pred.index, "0.975quant"]

IP_df_m <- melt(IP_df,
            id.vars = c("LONGITUDE", "LATITUDE", "date_index"),
            measure.vars = c("pred_mean", "pred_ll", "pred_ul")
)
IP_df_sp <- sp::SpatialPointsDataFrame(coords = data.frame(long = IP_df$LONGITUDE, lat = IP_df$LATITUDE), data = IP_df[, -c(2:3)], proj4string = proj)
ggplot() +
      geom_sf(data = IP_df_sp, aes(fill = pred_mean))

ggplot() + 
      geom_point(data = IP_df, aes(x = LONGITUDE, y = LATITUDE, col = pred_mean)) +
      labs(x = "", y = "") +
      facet_wrap(. ~ date_index) +
      theme_bw() +
      coord_equal()


# ------------------------------------------------------------------------------------------------------------------------
# 13. Posterior mean of the space-time random field 
# ------------------------------------------------------------------------------------------------------------------------

# Plot the posterior mean of the space-time random field = the latent field (not directly observed)
# This shows us the variation in the spatial effect, as well as spatial dependence.

xmean <- list()
for (j in 1:2) {
      xmean[[j]] <- inla.mesh.project(
            projgrid,  model$summary.random$i$mean[index_set$i.group == j])
}

xy.inGroup <- c(xy.in,xy.in)
xmean <- unlist(xmean)
xmean <- xmean[xy.inGroup]

dataObj <- data.frame(mean = xmean,
                      ind = rep(c(1,2), each = length(xmean)/2))

predcoordsGroup <- do.call(rbind, list(predcoords, predcoords))
spatObj <- sp::SpatialPixelsDataFrame(points  = predcoordsGroup, data = dataObj, proj4string = proj)
spatObj <- crop(spatObj, TZ_outline)
spatObj@data[["ind"]][spatObj@data[["ind"]] == 1] <- "1980-1999"
spatObj@data[["ind"]][spatObj@data[["ind"]] == 2] <- "2000-2020"

ggplot() + 
      gg(spatObj) + 
      facet_grid(~ind) +
      coord_equal() +
      viridis::scale_fill_viridis() +
      theme_void() +
      theme(plot.title = element_text(hjust = 0.5, vjust = 5)) +
      ggtitle(paste0("Posterior random field - ", species))



