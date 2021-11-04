args = commandArgs(trailingOnly = TRUE)

suppressPackageStartupMessages(library(INLA, quietly=TRUE))

library(data.table)
library(tidyverse)
library(lubridate)
library(raster)
library(rgdal)

setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Scripts')
#setwd('/home/ahomec/p/philism/Joris_work/scripts')
sapply(list.files(pattern="*.R"),source,.GlobalEnv)

species_list = c('Cisticola juncidis', 'Eremopterix leucopareia', 'Estrilda astrild', 'Histurgops ruficauda')
#species_list = c('Passer domesticus', 'Cisticola juncidis', 'Estrilda astrild', 'Histurgops ruficauda', 'Ploceus nigricollis', 
#                 'Cisticola brunnescens', 'Chrysococcyx cupreus', 'Tauraco hartlaubi', 'Ploceus castaneiceps', 'Nigrita canicapilla', 
#                 'Nectarinia kilimensis', 'Lanius collaris', 'Terpsiphone viridis', 'Oriolus auratus', 'Bubo capensis', 'Bubo africanus', 'Eremopterix leucopareia')
i <- as.numeric(args[1])

time_list = c('20s')

# Max.edge based on estimated range from pilot model
estimated_range = 2
max.edge = estimated_range/8

setwd('/Users/philism/OneDrive - NTNU/PhD/Joris_work/Philip_data')
#setwd('/home/ahomec/p/philism/Joris_work/Philip_data')
load(paste0("TZ_INLA_model_file_temporal_E", round(max.edge, digits = 3), ".RData"))

# The prior ranges to use for the model
# Prior range could be 10*max.edge. Should not be smaller than mesh resolution.
range0_ls = c(5) 
Prange_ls = c(0.5)
range_combs <- crossing(range0_ls, Prange_ls)

sigma0_ls = c(3)
Psigma_ls =  c(0.01)
sigma_combs <- crossing(sigma0_ls, Psigma_ls)

all_combs <- crossing(range_combs, sigma_combs)


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
#ebird_sp <- make_ebird_sp(species_list[1], TZ_outline)

df <- ebird_filtered %>% 
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

df_sp <- SpatialPointsDataFrame(
  coords = df[, c("LONGITUDE", "LATITUDE")],
  data = data.frame(presence = df$occurrence, duration_minutes = df$duration_minutes,
                    effort_distance_km = df$effort_distance_km, number_observers = df$number_observers, time_index = df$date_index),
  proj4string = crs(proj)
)

# Only include eBird data points for the region of interest
# Get intersecting points
in_sp <- rgeos::gIntersection(df_sp, TZ_outline)

# Only keep intersecting points in original spdf
ebird_sp <- df_sp[in_sp, ]



range01 <- function(x){(x - min(x))/(max(x) - min(x))}
ebird_sp$duration_minutes <- range01(ebird_sp$duration_minutes)
#atlas_sp$effort <- range01(atlas_sp$effort)

# Create data stacks for each of the different data types
#  stk.eBird <- MakeBinomStack(
#    observs = ebird_sp,
#    data = temporal_variables_no_BG,
#    mesh = Mesh$mesh,
#    presname = "presence",
#    tag = "eBird",
#    InclCoords = TRUE)#,
#add_effort = TRUE, 
#effort_names = names(ebird_sp@data[ , -which(names(ebird_sp) %in% c("presence"))]))

filtered_covs <- temporal_variables_no_BG[,1:2]

calc_covs <- FALSE

if (calc_covs) {

Nearest_covs <- GetNearestCovariate(ebird_sp,filtered_covs)  

} else Nearest_covs <- readRDS('/Users/philism/Downloads/Nearest_covs.RDS')

ebird_sp@data[, names(Nearest_covs@data)] <- Nearest_covs@data

ebird_sp <- as(ebird_sp, 'data.frame')

## Form new covariate called annual rain which gives value based on time

ebird_sp <- ebird_sp %>% mutate(annual_rain = ifelse(time_index == 1, TZ_ann_rain_1960s, TZ_ann_rain_2000s))

ebird_sp <- SpatialPointsDataFrame(coords = ebird_sp[, c("LONGITUDE", "LATITUDE")],
                                   data = ebird_sp[, !names(ebird_sp)%in%c('LONGITUDE', 'LATITUDE')],
                                   proj4string = crs(proj))
ebird_sp@data[, 'intercept'] <- 1
ebird_sp$presence <- as.numeric(ebird_sp$presence)

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

projmat <- inla.spde.make.A(mesh = Mesh$mesh, 
                            loc = as.matrix(ebird_sp@coords),
                            group = ebird_sp$time_index) 


#stk.eBird <- inla.stack(data=list(resp=ebird_sp@data[,'presence']), A=list(1,projmat), tag='eBird',
#                        effects=list(ebird_sp@data, list(i=1:Mesh$mesh$n)))

stk.eBird <- inla.stack(data=list(resp=ebird_sp@data[,'presence']),
                        A=list(1,projmat), 
                        tag='eBird',
                        effects=list(ebird_sp@data, i = ind))


form <- resp ~ 0 +
  intercept + 
  annual_field +
  f(time_index, model = pcspde, covariates = annual_rain) +   # Accounts for temporal structure of the covariate
  f(i, model = spde, group = i.group, control.group = list(model = 'ar1'))

model <- inla(form, family = "binomial", control.family = list(link = "cloglog"), data = inla.stack.data(stk.eBird), 
            verbose = FALSE,
            control.predictor = list(A = inla.stack.A(stk.eBird), 
            link = NULL, compute = TRUE), 
            E = inla.stack.data(stk.eBird)$e, 
            control.compute = list(waic = FALSE, dic = FALSE, cpo = FALSE))
summary(model)


#saveRDS(model, 'model.RDS')

