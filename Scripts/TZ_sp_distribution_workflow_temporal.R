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
      summarise(occurrence = ifelse(species_list[1] %in% `SCIENTIFIC NAME`, TRUE, FALSE)) %>% 
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
      filter(Scientific == species_list[1]) %>% 
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

spde <- inla.spde2.matern(mesh = Mesh$mesh)

pcspde <- inla.spde2.pcmatern(
      mesh = Mesh$mesh,
      alpha = 2,
      prior.range = c(5, 0.01),   
      prior.sigma = c(2, 0.01))

ind <- inla.spde.make.index(name ='i',
                            n.spde = spde$n.spde,
                            n.group = 2)

#projmat <- inla.spde.make.A(Mesh$mesh, as.matrix(ebird_sp@coords)) 


##Change this somehow??
projmat_eBird <- inla.spde.make.A(mesh = Mesh$mesh,
                                  loc = as.matrix(ebird_sp@coords),
                                  group = ebird_sp$date_index)

# projmat <- inla.spde.make.A(mesh = Mesh$mesh, 
#                             loc = as.matrix(ebird_sp@coords),
#                             group = ebird_sp$date_index) 

stk.eBird <- inla.stack(data=list(resp=ebird_sp@data[,'presence']),
                        A=list(1,projmat_eBird), 
                        tag='eBird',
                        effects=list(ebird_sp@data, i = ind))

projmat_atlas <- inla.spde.make.A(mesh = Mesh$mesh,
                                  loc = as.matrix(atlas_sp@coords),
                                  group = atlas_sp$date_index)

stk.atlas <- inla.stack(data = list(resp = atlas_sp@data[,'presence']),
                        A = list(1,projmat_atlas),
                        tag = 'atlas',
                        effects = list(atlas_sp@data, i = ind))
##Make predictions grid
if(!exists("Nxy.scale")) Nxy.scale <- 0.1  # about 10km resolution

Boundary <- Mesh$mesh$loc[Mesh$mesh$segm$int$idx[, 2], ]
Nxy.size <- c(diff(range(Boundary[, 1])), diff(range(Boundary[, 2])))
Nxy <- round(Nxy.size / Nxy.scale)

# Code adapted from 'MakeProjectionGrid' function.
if(class(Boundary)=="SpatialPolygons") {
      #create grid based on inla mesh and number of cells specified by the nxy parameter
      projgrid <- inla.mesh.projector(Mesh$mesh, xlim=Boundary@bbox["x",], ylim=Boundary@bbox["y",], dims=Nxy)
      #get the index of points on the grid within the boundary
      xy.in <- !is.na(over(SpatialPoints(projgrid$lattice$loc, proj4string=Boundary@proj4string), Boundary))
} else {
      if(ncol(Boundary)<2) stop("Boundary should have at least 2 columns")
      #create grid based on inla mesh
      projgrid <- inla.mesh.projector(Mesh$mesh, xlim=range(Boundary[,1]), ylim=range(Boundary[,2]), dims=Nxy)
      #get the index of points on the grid within the boundary
      xy.in <- splancs::inout(projgrid$lattice$loc, Boundary)
}
#select only points on the grid that fall within the boundary
predcoords <- projgrid$lattice$loc[which(xy.in),]
colnames(predcoords) <- c('Long','Lat')
Apred <- projgrid$proj$A[which(xy.in), ]

# Extract covariates for points, add intercept and coordinates
NearestCovs=GetNearestCovariate(points=predcoords, covs=temporal_variables)
NearestCovs$Intercept=1
NearestCovs@data[,colnames(NearestCovs@coords)] <- NearestCovs@coords

#NewCoords <- do.call(rbind, list(NearestCovs@coords, NearestCovs@coords))
#NewData <- do.call(rbind, list(NearestCovs@data, NearestCovs@data))

#NearestCovsGroup <- sp::SpatialPointsDataFrame(coords = NewCoords,
#                                               data = NewData,
#                                               proj4string = proj)

#colnames(NearestCovsGroup@coords) <- c('Long', 'Lat')
# stack the predicted data
# Need to add an index?


## Why is effects two lists?
## Should I also just have effects = list(ind, Neearestcovs@data) ?

##Need new Apred
ApredGroup <- inla.spde.make.A(mesh = Mesh$mesh, loc = cbind(predcoords[,1], predcoords[,2]),
                               n.group = 2)

#Things to do here:: replicate NearestCovs@data + coords twice for the 2 time periods
#                 :: add date_index = 1,2 to data
#                 :: Do we need a joint intercept?

stk.predGroup <- inla.stack(list(resp = rep(NA, nrow(NearestCovs@data))),
                            A=list(1,ApredGroup), tag= 'pred.group', effects=list(NearestCovs@data, list(i.group = ind$i.group)))

integrated_stack <- inla.stack(stk.eBird, stk.atlas, stk.predGroup, stk.ip)

# Not sure if this is correct, but need to somehow add a copy of 'date_index', with different namez
# index_list <- paste0("date_index", 1:6)
# 
# for (i in 1:length(index_list)){
#       new_var <- index_list[i]
#       integrated_stack[["effects"]][["data"]][[new_var]] <- integrated_stack[["effects"]][["data"]][["date_index"]]
#       integrated_stack[["effects"]][["ncol"]][[new_var]] <- integrated_stack[["effects"]][["ncol"]][["date_index"]]
#       integrated_stack[["effects"]][["names"]][[new_var]] <- integrated_stack[["effects"]][["names"]][["date_index"]]
# }


#Add other covs here
#Do I add covs like effort to predstack?

# MODEL 1 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
form_1 <- resp ~ 0 +
      ebird_intercept +
      atlas_intercept +
      duration_minutes + effort + 
      f(date_index, annual_rain_1, model="rw1", scale.model=TRUE, constr=FALSE,
        hyper = list(theta = list(prior="pc.prec", param=c(4,0.01)))) +  
      f(date_index1, annual_rain_2, model="rw1", scale.model=TRUE, constr=FALSE,
        hyper = list(theta = list(prior="pc.prec", param=c(4,0.01)))) +   
      f(date_index2, hottest_temp_1, model="rw1", scale.model=TRUE, constr=FALSE,
        hyper = list(theta = list(prior="pc.prec", param=c(4,0.01)))) +   
      f(date_index3, hottest_temp_2, model="rw1", scale.model=TRUE, constr=FALSE,
        hyper = list(theta = list(prior="pc.prec", param=c(4,0.01)))) +
      f(date_index4, max_dryspell_1, model="rw1", scale.model=TRUE, constr=FALSE,
        hyper = list(theta = list(prior="pc.prec", param=c(4,0.01)))) +
      f(date_index5, max_dryspell_2, model="rw1", scale.model=TRUE, constr=FALSE,
        hyper = list(theta = list(prior="pc.prec", param=c(4,0.01)))) +
      f(i, model = spde, group = i.group, control.group = list(model = 'ar1'))

# MODEL 2 ----------------------------------------------------------------------------------------------------------------------------------------------------------------------
form_2 <- resp ~ 0 +
      ebird_intercept +
      atlas_intercept +
      duration_minutes + effort + 
      annual_rain_1 + annual_rain_2 + hottest_temp_1 + hottest_temp_2 + 
      max_dryspell_1 + max_dryspell_2 + BG_1 + BG_2 +
      date_index + indicator +
      f(i, model = spde, group = i.group, control.group = list(model = 'ar1'))


#form <- resp ~ 0 +
#  intercept + 
#  annual_rain +
#  f(time_index, annual_rain, model="rw1", scale.model=TRUE, constr=FALSE) +   # Accounts for temporal structure of the covariate
#  f(i, model = spde, group = i.group, control.group = list(model = 'ar1'))  
# At each time point, spatial locations are linked through the spde.
# Across time, the process evolves according to an AR(1) process.

##Things to do::
#Add annual rain to the stk.pred. Can't do predections without it
model <- inla(form_2, family = "binomial", control.family = list(link = "cloglog"), # Backtransform for probability scale
              lincomb = all_lc,
              data = inla.stack.data(integrated_stack), 
              verbose = FALSE,
              control.predictor = list(A = inla.stack.A(integrated_stack), 
                                       link = NULL, compute = TRUE), 
              E = inla.stack.data(integrated_stack)$e, 
              control.compute = list(waic = FALSE, dic = FALSE, cpo = FALSE))
summary(model)
model$summary.random

setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_inla_spatial_temporal/model_output')
saveRDS(model, file = "model_all_lc_GAM_form2.RDS")
model <- readRDS('model_all_lc_GAM_form2.RDS')

# Collect model results
res.bits <- list("summary.lincomb" = model$summary.lincomb, 
                 "summary.lincomb.derived" = model$summary.lincomb.derived,
                 "marginals.lincomb.derived" = model$marginals.lincomb.derived,
                 "summary.fixed"  = model$summary.fixed,
                 "summary.hyperpar" = model$summary.hyperpar,
                 "marginals.fixed" = model$marginals.fixed)

# for (i in 1:length(res.bits$marginals.fixed)) {
#       tmp = inla.tmarginal(function(x) x, res.bits$marginals.fixed[[i]]) ## not sure how to fix?
#       plot(tmp, type = "l", xlab = paste("Fixed effect marginal", i, ":", names(res.bits[["marginals.fixed"]])[i]), ylab = "Density")
#       abline(v = 0, lty = 2)
# }

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



# From bee example:
# ## Plot here is real data, we use empty box 0-1. Use rug() to still display real data.
# plot(cov.seq, rep(1, NROW(cov.seq), type = "n")
# plot(variables$et, variables$TD, axes = TRUE, bty = "l",
#      xlab = "ET",
#      ylab = "Bee species richness", pch = 20, col = rgb(0,0,0,0.05),
#      log = "y")
# box(bty = "o")
# sc.p <- attributes(scale(variables$et))
# xs <- et.seq
# polygon(c(xs, rev(xs)), exp(c(res.bits$summary.lincomb.derived[grep("et", rownames(res.bits$summary.lincomb.derived)), "0.025quant"],
#                               rev(res.bits$summary.lincomb.derived[grep("et", rownames(res.bits$summary.lincomb.derived)), "0.975quant"]))),
#         border = NA, col = rgb(0.7, 0, 0.1, 0.5))
# lines(xs, exp(res.bits$summary.lincomb.derived$`0.5quant`[grep("et", rownames(res.bits$summary.lincomb.derived))]), lwd = 2, col = rgb(0.7, 0, 0.1))
# rug()
index <- inla.stack.index(stack = integrated_stack, tag = "pred.group")$data

dp <- rbind(cbind(predcoords, 1), cbind(predcoords, 2))
dp <- data.frame(dp)
names(dp) <- c("x", "y", "time")

dp$pred_mean <- model$summary.fitted.values[index, "mean"]
dp$pred_ll <- model$summary.fitted.values[index, "0.025quant"]
dp$pred_ul <- model$summary.fitted.values[index, "0.975quant"]

dpm <- melt(dp,
            id.vars = c("x", "y", "time"),
            measure.vars = c("pred_mean", "pred_ll", "pred_ul")
)

ggplot() + 
      geom_tile(data = dpm, aes(x = x, y = y, fill = value)) +
      labs(x = "", y = "") +
      facet_wrap(variable ~ time) +
      theme_bw() +
      coord_equal()





model$summary.random$i$mean_inv <- inla.link.cloglog(model$summary.random$i$mean, inverse = TRUE)
xmean <- list()
for (j in 1:2) {
      xmean[[j]] <- inla.mesh.project(
            projgrid,  model$summary.random$i$mean_inv[ind$i.group == j])
}

xy.inGroup <- c(xy.in,xy.in)
xmean <- unlist(xmean)
xmean <- xmean[xy.inGroup]

dataObj <- data.frame(mean = xmean,
                      ind = rep(c(1,2),each = length(xmean)/2))

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
      ggtitle("Distribution")


#saveRDS(model, 'model.RDS')

