suppressPackageStartupMessages(library(INLA, quietly=TRUE))

library(data.table)
library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Scripts')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/source')
sapply(list.files(pattern="*.R"),source,.GlobalEnv)

# Import and prepare data-----------------------------------------------------------------------------------

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/data')

# Set Region of Interest
TZ_outline <- readOGR('TZ_simpler.shp')   # Full extent
# setwd('~/Documents/TZ_inla_spatial_temporal/Data/New_files_INLA')
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
TZ_BG_10_19 <- raster('BG_2010_19_1000m.tif') %>% mask(., ROI) %>% projectRaster(., TZ_annual_median_rain_80_00)
TZ_BG_10_19[is.nan(TZ_BG_10_19)] <- NA

temporal_variables <- stack(TZ_annual_median_rain_80_00, TZ_annual_median_rain_00_20, 
                            TZ_ERA5_hottest_80_00, TZ_ERA5_hottest_00_20,
                            TZ_dryspell_80_00, TZ_dryspell_00_20,
                            TZ_BG_90_99, TZ_BG_10_19)

names(temporal_variables) <- c('TZ_ann_rain_1980s', 'TZ_ann_rain_2000s', 
                               'TZ_max_temp_1980s', 'TZ_max_temp_2000s',
                               'TZ_dryspell_1980s', 'TZ_dryspell_2000s')

temporal_variables <- as(temporal_variables, 'SpatialPointsDataFrame')
temporal_variables <- temporal_variables[!is.na(rowSums(temporal_variables@data)),]

# Prepare data for a GAM model, keep simple to avoid overfitting (two spline bases, limited flexibility GAM). 
# Separates out and scales a prediction vector
for (i in 1:length(names(temporal_variables))){
      temporal_variables <- prepare_GAM(temporal_variables, names(temporal_variables)[i])
      assign("temporal_variables", temporal_variables, envir = .GlobalEnv)
}

prepare_GAM <- function(df, colname){
      # Create prediction vector, as sequence of 100 values, spanning full range of values. 
      seq <- seq(min(df@data[, colname], na.rm = TRUE), max(df@data[, colname], na.rm = TRUE), length = 100)
      # Add prediction vector to beginning of original values, scale all to mean 0 and sd 1
      var.s <- c(scale(c(seq, df@data[, colname])))  
      # Split the full vector into two spline bases, and scale the new vectors to mean 0 and sd 1.
      # Explain:
      z.var.s <- scale(sortBase(var.s, n.knots = 2))  
      
      # Add the two scaled spline bases back to the original dataset. Omit the first 100 values,
      # this is the prediction vector, not original data
      df@data[, paste0("z.", colname, "1.s")] <- c(z.var.s[101:NROW(z.var.s), 1])
      df@data[, paste0("z.", colname, "2.s")] <- c(z.var.s[101:NROW(z.var.s), 2])
      
      # Add original, but scaled values back to the dataset
      df@data[, paste0(colname, ".s")]  <- var.s[-c(1:100)]
      
      ## Separate out prediction vectors for the two spline bases, and assign to environment as
      ## separate objects
      assign(paste0("z.", colname, "1.s"), c(z.var.s[1:100, 1]), envir = .GlobalEnv)
      assign(paste0("z.", colname, "2.s"), c(z.var.s[1:100, 2]), envir = .GlobalEnv)
      
      return(df)
}

# Make linear combinations, needed to make effect plots
TZ_max_temp_1980s_lc <- inla.make.lincombs(ebird_intercept = rep(1, 100),  
                                           atlas_intercept = rep(1, 100),
                                           date_index = rep(1, 100),
                                           hottest_temp_1 = z.TZ_max_temp_1980s1.s,   
                                           hottest_temp_2 = z.TZ_max_temp_1980s2.s); names(TZ_max_temp_1980s_lc) <- paste0("TZ_max_temp_1980s_lc", 1:100)
TZ_max_temp_2000s_lc <- inla.make.lincombs(ebird_intercept = rep(1, 100),  
                                           atlas_intercept = rep(1, 100),
                                           date_index = rep(2, 100),
                                           hottest_temp_1 = z.TZ_max_temp_2000s1.s,   
                                           hottest_temp_2 = z.TZ_max_temp_2000s2.s); names(TZ_max_temp_2000s_lc) <- paste0("TZ_max_temp_2000s_lc", 1:100)
TZ_ann_rain_1980s_lc <- inla.make.lincombs(ebird_intercept = rep(1, 100),  
                                           atlas_intercept = rep(1, 100),
                                           date_index = rep(1, 100),
                                           annual_rain_1 = z.TZ_ann_rain_1980s1.s,   
                                           annual_rain_2 = z.TZ_ann_rain_1980s2.s); names(TZ_ann_rain_1980s_lc) <- paste0("TZ_ann_rain_1980s_lc", 1:100)
TZ_ann_rain_2000s_lc <- inla.make.lincombs(ebird_intercept = rep(1, 100),  
                                           atlas_intercept = rep(1, 100),
                                           date_index = rep(2, 100),
                                           annual_rain_1 = z.TZ_ann_rain_2000s1.s,   
                                           annual_rain_2 = z.TZ_max_temp_2000s2.s); names(TZ_ann_rain_2000s_lc) <- paste0("TZ_ann_rain_2000s_lc", 1:100)
TZ_dryspell_1980s_lc <- inla.make.lincombs(ebird_intercept = rep(1, 100),  
                                           atlas_intercept = rep(1, 100),
                                           date_index = rep(1, 100),
                                           max_dryspell_1 = z.TZ_dryspell_1980s1.s,   
                                           max_dryspell_2 = z.TZ_dryspell_1980s2.s); names(TZ_dryspell_1980s_lc) <- paste0("TZ_dryspell_1980s_lc", 1:100)
TZ_dryspell_2000s_lc <- inla.make.lincombs(ebird_intercept = rep(1, 100),  
                                           atlas_intercept = rep(1, 100),
                                           date_index = rep(2, 100),
                                           max_dryspell_1 = z.TZ_dryspell_2000s1.s,   
                                           max_dryspell_2 = z.TZ_dryspell_2000s2.s); names(TZ_dryspell_2000s_lc) <- paste0("TZ_dryspell_2000s_lc", 1:100)
all_lc <- c(TZ_max_temp_1980s_lc, TZ_max_temp_2000s_lc, TZ_ann_rain_1980s_lc, TZ_ann_rain_2000s_lc, TZ_dryspell_1980s_lc, TZ_dryspell_2000s_lc)


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
#stk.ip <- MakeIntegrationStack(
#mesh = Mesh$mesh,
#data = temporal_variables,
#area = Mesh$w,
#tag = "ip",
#InclCoords = TRUE
#)

Points <- cbind(c(Mesh$mesh$loc[,1]), c(Mesh$mesh$loc[,2]))
colnames(Points) <- c("LONGITUDE", "LATITUDE")
NearestCovs <- GetNearestCovariate(points=Points, covs=temporal_variables)

Points <- rbind(Points, Points)
Points_data <- data.frame(date_index = rep(c(1,2), each = nrow(Points)/2))
Points_data[, 'annual_rain_s1'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_ann_rain_1980s1.s, NearestCovs@data$z.TZ_ann_rain_2000s1.s)
Points_data[, 'annual_rain_s2'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_ann_rain_1980s2.s, NearestCovs@data$z.TZ_ann_rain_2000s2.s)

Points_data[, 'hottest_temp_s1'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_max_temp_1980s1.s, NearestCovs@data$z.TZ_max_temp_2000s1.s)
Points_data[, 'hottest_temp_s2'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_max_temp_1980s2.s, NearestCovs@data$z.TZ_max_temp_2000s2.s)

Points_data[, 'max_dryspell_s1'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_dryspell_1980s1.s, NearestCovs@data$z.TZ_dryspell_2000s1.s)
Points_data[, 'max_dryspell_s2'] <- ifelse(Points_data$date_index == 1, NearestCovs@data$z.TZ_dryspell_1980s2.s, NearestCovs@data$z.TZ_dryspell_2000s2.s)

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

# Make data for projections
if(!exists("Nxy.scale")) Nxy.scale <- 0.1  # about 10km resolution

Boundary <- Mesh$mesh$loc[Mesh$mesh$segm$int$idx[, 2], ]
Nxy.size <- c(diff(range(Boundary[, 1])), diff(range(Boundary[, 2])))
Nxy <- round(Nxy.size / Nxy.scale)

# Make stack for projections
stk.pred <- MakeProjectionGrid(
  nxy = Nxy,
  mesh = Mesh$mesh,
  data = temporal_variables,
  tag = "pred",
  boundary = Boundary
)

# setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
setwd('/Users/joriswiethase/Google Drive (jhw538@york.ac.uk)/Work/PhD_York/Chapter3/TZ_INLA/data_processed')

save(proj, ROI, ebird_full, atlas_full, Mesh, stk.ip, stk.pred, temporal_variables, TZ_outline, all_lc, file = paste0("TZ_INLA_model_file_temporal_E", round(max.edge, digits = 3), ".RData"))
