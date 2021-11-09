args = commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages(library(INLA, quietly=TRUE))

library(data.table)
library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)
library(dplyr)
library(rlang)
library(inlabru)

setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Scripts')
#setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/source')

#setwd('/home/ahomec/p/philism/Joris_work/scripts')
sapply(list.files(pattern="*.R"),source,.GlobalEnv)

species_list = c('Cisticola juncidis', 'Eremopterix leucopareia', 'Estrilda astrild', 'Histurgops ruficauda')
#species_list = c('Passer domesticus', 'Cisticola juncidis', 'Estrilda astrild', 'Histurgops ruficauda', 'Ploceus nigricollis', 
#                 'Cisticola brunnescens', 'Chrysococcyx cupreus', 'Tauraco hartlaubi', 'Ploceus castaneiceps', 'Nigrita canicapilla', 
#                 'Nectarinia kilimensis', 'Lanius collaris', 'Terpsiphone viridis', 'Oriolus auratus', 'Bubo capensis', 'Bubo africanus', 'Eremopterix leucopareia')

setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
#setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/data_processed')

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
print(summary(ebird_filtered))

ebird_filtered <- ebird_filtered %>% 
  group_by(LATITUDE, LONGITUDE, `SAMPLING EVENT IDENTIFIER`, `DURATION MINUTES`, 
           `EFFORT DISTANCE KM`, `NUMBER OBSERVERS`, `OBSERVATION DATE`, `LOCALITY`, date_index) %>% 
  summarise(occurrence = ifelse(species_list[2] %in% `SCIENTIFIC NAME`, TRUE, FALSE)) %>% 
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

colnames(ebird_sp@coords) <- c('Long','Lat') ##Standardizing colnames to match atlas names

# Only include eBird data points for the region of interest
# Get intersecting points

in_sp <- rgeos::gIntersection(ebird_sp, TZ_outline)
#in_sp <- rgeos::gIntersection(df_sp, ROI)
in_sp <- rgeos::gIntersection(ebird_sp, ROI)

# Only keep intersecting points in original spdf
ebird_sp <- ebird_sp[in_sp, ]

atlas_full <- atlas_full %>%
  filter(time_period != 'x') %>%
  mutate(date_index = ifelse(time_period == '20s',2,1))

atlas_filtered <- atlas_full %>% 
  mutate(Scientific = trimws(Scientific, which = 'both')) %>% 
  filter(Scientific == species_list[2]) %>% 
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

# Take only non-GAM data for now
filtered_covs <- temporal_variables_no_BG[,1:6]

calc_covs <- FALSE

if (calc_covs) {
  
  Nearest_covs_ebird <- GetNearestCovariate(ebird_sp,filtered_covs)  
  Nearest_covs_atlas <- GetNearestCovariate(atlas_sp, filtered_covs)
  
  
} else {
  
  Nearest_covs_ebird <- readRDS('/Users/philism/Downloads/Nearest_covs_ebird.RDS')
  Nearest_covs_atlas <-  readRDS('/Users/philism/Downloads/Nearest_covs_atlas.RDS')
  # Add covariates to the bird data 
  ebird_sp@data[, names(Nearest_covs_ebird@data)] <- Nearest_covs_ebird@data
  atlas_sp@data[, names(Nearest_covs_atlas@data)] <- Nearest_covs_atlas@data
  
}

ebird_sp@data[, names(Nearest_covs_ebird@data)] <- Nearest_covs_ebird@data
ebird_sp <- as(ebird_sp, 'data.frame')

ebird_sp <- ebird_sp %>% mutate(annual_rain = ifelse(date_index == 1, TZ_ann_rain_1980s, TZ_ann_rain_2000s))

ebird_sp <- SpatialPointsDataFrame(coords = ebird_sp[, c("Long", "Lat")],
                                   data = ebird_sp[, !names(ebird_sp)%in%c('Long', 'Lat')],
                                   proj4string = crs(proj))
ebird_sp$presence <- as.numeric(ebird_sp$presence)

atlas_sp@data[, names(Nearest_covs_atlas@data)] <- Nearest_covs_atlas@data
atlas_sp <- as(atlas_sp, 'data.frame')

atlas_sp <- atlas_sp %>% mutate(annual_rain = ifelse(date_index == 1, TZ_ann_rain_1980s, TZ_ann_rain_2000s))
atlas_sp <- SpatialPointsDataFrame(coords = atlas_sp[, c('Long', 'Lat')],
                                   data = atlas_sp[, !names(atlas_sp)%in%c('Long','Lat')],
                                   proj4string = crs(proj))
atlas_sp$presence <- as.numeric(atlas_sp$presence)


ebird_likelihood <- inlabru::like(formula = presence ~ ebird_intercept + TZ_ann_rain_1960s + spatial,
                                  mesh = Mesh$mesh,
                                  family = 'binomial',
                                  data = ebird_sp,
                                  )

atlas_likelihood <- inlabru::like(formula = presence ~ atlas_intercept + TZ_ann_rain_1960s + spatial,
                                  mesh = Mesh$mesh,
                                  family = 'binomial',
                                  data = atlas_sp)


temporal_variables_no_BG <- as(temporal_variables_no_BG, 'SpatialPixelsDataFrame')

spde <- inla.spde2.matern(mesh = Mesh$mesh)
pcspde <- inla.spde2.pcmatern(
  mesh = Mesh$mesh,
  alpha = 2,
  prior.range = c(5, 0.01),   
  prior.sigma = c(2, 0.01))

components <- precence + date_index ~ atlas_intercept(1) +
                                      ebird_intercept(1) +
                                      TZ_ann_rain_1960s(main = temporal_variables_no_BG, model = 'linear') +
                                      spatial(main = coordinates, model = pcspde, group = date_index, control.group = list(model = 'ar1')) - 1

joint_model <- bru(components, ebird_likelihood,
                               atlas_likelihood)

