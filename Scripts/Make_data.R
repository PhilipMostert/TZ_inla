suppressPackageStartupMessages(library(INLA, quietly=TRUE))

library(data.table)
library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Scripts')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_inla_spatial_temporal/source')
sapply(list.files(pattern="*.R"),source,.GlobalEnv)

# Import and prepare data-----------------------------------------------------------------------------------

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/data')

# Set Region of Interest
# setwd('~/Documents/TZ_inla_spatial_temporal/Data/New_files_INLA')
TZ_outline <- readOGR('TZ_simpler.shp')
TZ_buffered <- readOGR('TZ_simpler_buffered.shp')
NTRI <- readOGR("NTRI_outline.shp")  # A small subset 
ROI <- TZ_buffered

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
# Import bird data, initial filtering
ebird_full <- fread("ebd_TZ_relMay-2021.txt") %>% 
  mutate(date = ymd(get("OBSERVATION DATE"))) %>%
  filter(!is.na(`DURATION MINUTES`), !is.na(`EFFORT DISTANCE KM`))

atlas_full <- fread("TZ_bird_atlas_data.csv") %>%
  filter(!is.na(effort))

proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# Import and prepare the temporally varying covariates
#setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Temporal_variables')
# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Joris_tif_files/Philip_additional_files')
# setwd('~/Documents/TZ_inla_spatial_temporal/Data/New_files_INLA')
# TZ_annual_median_rain_80_00 <- raster('TZ_annual_median_rain_80_00.tif') %>% mask(., ROI) 
# TZ_annual_median_rain_00_20 <- raster('TZ_annual_median_rain_00_20.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
# TZ_ERA5_coldest_80_00 <- raster('TZ_ERA5_coldest_temperature_1980_2000.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
# TZ_ERA5_coldest_00_20 <- raster('TZ_ERA5_coldest_temperature_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
# TZ_ERA5_hottest_80_00 <- raster('TZ_ERA5_hottest_temperature_1980_2000.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
# TZ_ERA5_hottest_00_20 <- raster('TZ_ERA5_hottest_temperature_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
# 
# setwd('~/Documents/TZ_inla_spatial_temporal/Data/New_files_INLA')
# 
# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Temporal_variables')

# Import and process covariate data
TZ_annual_median_rain_80_00 <- raster('TZbuff_annual_median_rain_1981_1999.tif') %>% mask(., ROI)
TZ_annual_median_rain_80_00[is.nan(TZ_annual_median_rain_80_00)] <- NA
TZ_annual_median_rain_00_20 <- raster('TZbuff_annual_median_rain_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_annual_median_rain_00_20[is.nan(TZ_annual_median_rain_00_20)] <- NA
TZ_ERA5_hottest_80_00 <- raster('TZbuff_ERA5_hottest_temperature_1981_1999.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_hottest_80_00[is.nan(TZ_ERA5_hottest_80_00)] <- NA
TZ_ERA5_hottest_00_20 <- raster('TZbuff_ERA5_hottest_temperature_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_ERA5_hottest_00_20[is.nan(TZ_ERA5_hottest_00_20)] <- NA
TZ_dryspell_80_00 <- raster('TZbuff_median_annual_dryspell_length_1981_1999.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_dryspell_80_00[is.nan(TZ_dryspell_80_00)] <- NA
TZ_dryspell_00_20 <- raster('TZbuff_median_annual_dryspell_length_2000_2020.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_dryspell_00_20[is.nan(TZ_dryspell_00_20)] <- NA
TZ_BG_90_99 <- raster('BG_1990_1999_1000m.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_BG_90_99[is.nan(TZ_BG_90_99)] <- NA
TZ_BG_10_19 <- raster('BG_2010_2019_1000m.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_BG_10_19[is.nan(TZ_BG_10_19)] <- NA

# BG layer has large gaps in data, this needs to be accounted for. Manually create indicator layer and BG interaction layer,
# to make sure model ignores areas with NA BG
# Make indicator layer where 0 is NA in BG, and 1 is value in BG
indicator_90s <- TZ_BG_90_99 
values(indicator_90s)[!is.na(values(indicator_90s))] <- 1
values(indicator_90s)[is.na(values(indicator_90s))] <- 0

indicator_2010s <- TZ_BG_10_19 
values(indicator_2010s)[!is.na(values(indicator_2010s))] <- 1
values(indicator_2010s)[is.na(values(indicator_2010s))] <- 0

# Stack all covariate layers
variables <- stack(TZ_annual_median_rain_80_00, TZ_annual_median_rain_00_20, 
                            TZ_ERA5_hottest_80_00, TZ_ERA5_hottest_00_20,
                            TZ_dryspell_80_00, TZ_dryspell_00_20,
                            TZ_BG_90_99, TZ_BG_10_19, indicator_90s, indicator_2010s)

names(variables) <- c('TZ_ann_rain_1980s', 'TZ_ann_rain_2000s', 
                               'TZ_max_temp_1980s', 'TZ_max_temp_2000s',
                               'TZ_dryspell_1980s', 'TZ_dryspell_2000s',
                               'TZ_BG_90_99', 'TZ_BG_10_19', 'indicator_90s', 'indicator_2010s')

variables <- as(variables, 'SpatialPointsDataFrame')

# Make spdf without NAs, do transformation and GAM prep
variables_no_NA <- variables[!is.na(rowSums(variables@data)),]
variables_no_NA <- prepare_GAM(variables_no_NA, 'TZ_BG_90_99')
variables_no_NA <- prepare_GAM(variables_no_NA, 'TZ_BG_10_19')

# Put NAs back, replace with any other value, here 0 (will be ignored)
variables$z.TZ_BG_90_991.s <- variables$TZ_BG_90_99
variables$z.TZ_BG_90_991.s[!is.na(variables$z.TZ_BG_90_991.s)] <- variables_no_NA$z.TZ_BG_90_991.s
variables$z.TZ_BG_90_991.s[is.na(variables$z.TZ_BG_90_991.s)] <- 0

variables$z.TZ_BG_90_992.s <- variables$TZ_BG_90_99
variables$z.TZ_BG_90_992.s[!is.na(variables$z.TZ_BG_90_992.s)] <- variables_no_NA$z.TZ_BG_90_992.s
variables$z.TZ_BG_90_992.s[is.na(variables$z.TZ_BG_90_992.s)] <- 0

variables$z.TZ_BG_10_191.s <- variables$TZ_BG_10_19
variables$z.TZ_BG_10_191.s[!is.na(variables$z.TZ_BG_10_191.s)] <- variables_no_NA$z.TZ_BG_10_191.s
variables$z.TZ_BG_10_191.s[is.na(variables$z.TZ_BG_10_191.s)] <- 0

variables$z.TZ_BG_10_192.s <- variables$TZ_BG_10_19
variables$z.TZ_BG_10_192.s[!is.na(variables$z.TZ_BG_10_192.s)] <- variables_no_NA$z.TZ_BG_10_192.s
variables$z.TZ_BG_10_192.s[is.na(variables$z.TZ_BG_10_192.s)] <- 0

# Remove BG layer, no longer needed
temporal_variables <- variables[ , !names(variables) %in% c("TZ_BG_90_99", "TZ_BG_10_19")]
temporal_variables <- temporal_variables[!is.na(rowSums(temporal_variables@data)), ]

# Due to the model structure when using form_2, we don't have to separate the different time periods for the linear combinations.
# Consider them together:
TZ_ann_rain <- c(temporal_variables@data[["TZ_ann_rain_1980s"]], temporal_variables@data[["TZ_ann_rain_2000s"]])
TZ_max_temp <- c(temporal_variables@data[["TZ_max_temp_1980s"]], temporal_variables@data[["TZ_max_temp_2000s"]])
TZ_dryspell <- c(temporal_variables@data[["TZ_dryspell_1980s"]], temporal_variables@data[["TZ_dryspell_2000s"]])
TZ_BG <- c(variables@data[["TZ_BG_90_99"]], variables@data[["TZ_BG_10_19"]])

# Make sequences for linear combinations, based on combined covariates data for both time periods
prepare_lincombs(TZ_ann_rain)
prepare_lincombs(TZ_max_temp)
prepare_lincombs(TZ_dryspell)
prepare_lincombs(TZ_BG)

# Prepare remaining data for a GAM model
remaining_vars <- c('TZ_ann_rain_1980s', 'TZ_ann_rain_2000s', 
  'TZ_max_temp_1980s', 'TZ_max_temp_2000s',
  'TZ_dryspell_1980s', 'TZ_dryspell_2000s')

for (i in 1:length(remaining_vars)){
      temporal_variables <- prepare_GAM(temporal_variables, remaining_vars[i])
      assign("temporal_variables", temporal_variables, envir = .GlobalEnv)
}

# Combine the unscaled prediction vectors
all.seq <- mget(ls(pattern = "TZ_.*.seq"))

# Make linear combinations, needed to make effect plots. Variable names have to match model term names
# Effect plot scaled according to eBird observations, if only eBird intercept included
# Covars for two time periods have to be transformed together. 
# Make seq min to max for all covars values, two time periods combined, add time 1 and time 2 data after that (var.s step)

TZ_max_temp_lc <- inla.make.lincombs(ebird_intercept = rep(1, 100),  
                                     hottest_temp_1 = z.TZ_max_temp1.s,   
                                     hottest_temp_2 = z.TZ_max_temp2.s); names(TZ_max_temp_lc) <- paste0("TZ_max_temp_lc", 1:100)
TZ_ann_rain_lc <- inla.make.lincombs(ebird_intercept = rep(1, 100),  
                                     annual_rain_1 = z.TZ_ann_rain1.s,   
                                     annual_rain_2 = z.TZ_ann_rain2.s); names(TZ_ann_rain_lc) <- paste0("TZ_ann_rain_lc", 1:100)
TZ_dryspell_lc <- inla.make.lincombs(ebird_intercept = rep(1, 100),
                                     max_dryspell_1 = z.TZ_dryspell1.s,   
                                     max_dryspell_2 = z.TZ_dryspell2.s); names(TZ_dryspell_lc) <- paste0("TZ_dryspell_lc", 1:100)
TZ_BG_lc       <- inla.make.lincombs(ebird_intercept = rep(1, 100), 
                                     BG_1 = z.TZ_BG1.s,   
                                     BG_2 = z.TZ_BG2.s); names(TZ_BG_lc) <- paste0("TZ_BG_lc", 1:100)

# Combine all linear combinations, to include in the final model.
all_lc <- c(TZ_max_temp_lc, TZ_ann_rain_lc, TZ_dryspell_lc, TZ_BG_lc)


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
# Get mesh triangle centroids
Points <- cbind(c(Mesh$mesh$loc[,1]), c(Mesh$mesh$loc[,2]))
colnames(Points) <- c("LONGITUDE", "LATITUDE")
# Get value of nearest covariate to each mesh centroids
NearestCovs <- GetNearestCovariate(points=Points, covs=temporal_variables)

Points <- rbind(Points, Points)
Points_data <- data.frame(date_index = rep(c(1,2), each = nrow(Points)/2))
Points_data[, 'annual_rain_1'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_ann_rain_1980s1.s, NearestCovs@data$z.TZ_ann_rain_2000s1.s)
Points_data[, 'annual_rain_2'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_ann_rain_1980s2.s, NearestCovs@data$z.TZ_ann_rain_2000s2.s)

Points_data[, 'hottest_temp_1'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_max_temp_1980s1.s, NearestCovs@data$z.TZ_max_temp_2000s1.s)
Points_data[, 'hottest_temp_2'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_max_temp_1980s2.s, NearestCovs@data$z.TZ_max_temp_2000s2.s)

Points_data[, 'max_dryspell_1'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_dryspell_1980s1.s, NearestCovs@data$z.TZ_dryspell_2000s1.s)
Points_data[, 'max_dryspell_2'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_dryspell_1980s2.s, NearestCovs@data$z.TZ_dryspell_2000s2.s)

Points_data[, 'BG_1'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_BG_90_991.s, NearestCovs@data$z.TZ_BG_10_191.s)
Points_data[, 'BG_2'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_BG_90_992.s, NearestCovs@data$z.TZ_BG_10_192.s)

sequence <- seq(min(Points_data[, 'annual_rain_1'], na.rm = TRUE), max(Points_data[, 'annual_rain_1'], na.rm = TRUE), length = 100)
var.s <- c(scale(c(sequence, Points_data[, 'annual_rain_1'])))  
comb_vector <- var.s[1:100]
      
      
IP_sp <- sp::SpatialPointsDataFrame(coords = Points, data = Points_data, proj4string = proj)
IP_sp@data$Intercept <- 1
IP_sp@data[, c("LONGITUDE", "LATITUDE")] <- IP_sp@coords
projmat.ip <- Matrix::Diagonal(2 * Mesh$mesh$n, rep(1, Mesh$mesh$n * 2)) 

ind <- inla.spde.make.index(name ='i',
                            n.spde = Mesh$mesh$n,
                            n.group = 2)

stk.ip <- inla.stack(data=list(resp= NA, e=c(Mesh$w, Mesh$w)), ##Removed the rep(0, ...)?
                     A=list(1,projmat.ip), tag='ip',
                     effects=list(IP_sp@data, i = ind))

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/data_processed')

save(proj, ROI, ebird_full, atlas_full, Mesh, stk.ip,temporal_variables, TZ_outline, all_lc, all.seq, file = paste0("TZ_INLA_model_file_temporal_E", round(max.edge, digits = 3), ".RData"))
