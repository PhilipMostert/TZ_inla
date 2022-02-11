suppressPackageStartupMessages(library(INLA, quietly=TRUE))

library(data.table)
library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)

setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Scripts')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/source')
sapply(list.files(pattern="*.R"),source,.GlobalEnv)

# Import and prepare data-----------------------------------------------------------------------------------

setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/data')

# Set Region of Interest
TZ_outline <- readOGR('TZ_simpler.shp')   # Full extent
TZ_buffered <- readOGR('TZ_simpler_buffered.shp')
NTRI <- readOGR("NTRI_outline.shp")  # A small subset 
ROI <- TZ_buffered

# Import bird data
ebird_full <- fread("ebd_TZ_relMay-2021.txt") %>% 
  mutate(date = ymd(get("OBSERVATION DATE"))) %>%
  filter(!is.na(`DURATION MINUTES`), !is.na(`EFFORT DISTANCE KM`))

atlas_full <- fread("TZ_bird_atlas_data.csv") %>%
  filter(!is.na(effort))

proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Import and prepare the temporally varying covariates
#setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Temporal_variables')
setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Joris_tif_files/Philip_additional_files')
TZ_annual_median_rain_80_00 <- raster('TZ_annual_median_rain_80_00.tif') %>% mask(., ROI) 
TZ_annual_median_rain_00_20 <- raster('TZ_annual_median_rain_00_20.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_coldest_80_00 <- raster('TZ_ERA5_coldest_temperature_1980_2000.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_coldest_00_20 <- raster('TZ_ERA5_coldest_temperature_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_hottest_80_00 <- raster('TZ_ERA5_hottest_temperature_1980_2000.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_hottest_00_20 <- raster('TZ_ERA5_hottest_temperature_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)

setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Temporal_variables')

TZ_annual_median_rain_80_00 <- raster('TZbuff_annual_median_rain_1981_1999.tif') %>% mask(., ROI) 
TZ_annual_median_rain_00_20 <- raster('TZbuff_annual_median_rain_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_coldest_80_00 <- raster('TZbuff_ERA5_coldest_temperature_1981_1999.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_coldest_00_20 <- raster('TZbuff_ERA5_coldest_temperature_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_hottest_80_00 <- raster('TZbuff_ERA5_hottest_temperature_1981_1999.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_hottest_00_20 <- raster('TZbuff_ERA5_hottest_temperature_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_dryspell_80_00 <- raster('TZbuff_median_annual_dryspell_length_1981_1999.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_dryspell_00_20 <- raster('TZbuff_median_annual_dryspell_length_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)

temporal_variables <- stack(TZ_annual_median_rain_80_00, TZ_annual_median_rain_00_20, 
                            TZ_ERA5_coldest_80_00, TZ_ERA5_coldest_00_20, 
                            TZ_ERA5_hottest_80_00, TZ_ERA5_hottest_00_20,
                            TZ_dryspell_80_00, TZ_dryspell_00_20)

names(temporal_variables) <- c('TZ_ann_rain_1980s', 'TZ_ann_rain_2000s', 
                               'TZ_min_temp_1980s', 'TZ_min_temp_2000s', 
                               'TZ_max_temp_1980s', 'TZ_max_temp_2000s',
                               'TZ_dryspell_1980s', 'TZ_dryspell_2000s')

temporal_variables <- as(temporal_variables, 'SpatialPointsDataFrame')

temporal_variables <- prepare_GAM(temporal_variables, 'TZ_ann_rain_1980s')
temporal_variables <- prepare_GAM(temporal_variables, 'TZ_ann_rain_2000s')
temporal_variables <- prepare_GAM(temporal_variables, 'TZ_min_temp_1980s')
temporal_variables <- prepare_GAM(temporal_variables, 'TZ_min_temp_2000s')
temporal_variables <- prepare_GAM(temporal_variables, 'TZ_max_temp_1980s')
temporal_variables <- prepare_GAM(temporal_variables, 'TZ_max_temp_2000s')

# Prepare model parameters-----------------------------------------------------------------------------
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

# Make projection stack stack for background mesh, the prediction stack with NA as response.
# Projection stack = Integration stack
stk.ip <- MakeIntegrationStack(
  mesh = Mesh$mesh,
  data = temporal_variables,
  area = Mesh$w,
  tag = "ip",
  InclCoords = TRUE
)

# Make data for projections
if(!exists("Nxy.scale")) Nxy.scale <- 0.1  # about 10km resolution

Boundary <- Mesh$mesh$loc[Mesh$mesh$segm$int$idx[, 2], ]
Nxy.size <- c(diff(range(Boundary[, 1])), diff(range(Boundary[, 2])))
Nxy <- round(Nxy.size / Nxy.scale)

# Make stack for projections
#Rewrite ProjectionGrid to incorporate multiple time periods::
# Look at MakeProjectionGrid script + read https://becarioprecario.bitbucket.io/spde-gitbook/ch-spacetime.html again
# Will need to load in some other objects ie the ind = inla.make.index from TZ_sp...
stk.pred <- MakeProjectionGrid(
  nxy = Nxy,
  mesh = Mesh$mesh,
  data = temporal_variables,
  tag = "pred",
  boundary = Boundary
)


setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/data_processed')

save(proj, ROI, TZ_outline, ebird_full, atlas_full, Mesh, stk.ip, stk.pred, temporal_variables_no_BG, file = paste0("TZ_INLA_model_file_temporal_E", round(max.edge, digits = 3), ".RData"))
save(proj, ROI, ebird_full, atlas_full, Mesh, stk.ip, stk.pred, temporal_variables, TZ_outline, file = paste0("TZ_INLA_model_file_temporal_E", round(max.edge, digits = 3), ".RData"))
