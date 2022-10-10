library(tidyverse)
library(readxl)
library(INLA)





# -----------------------------------------------------------------------------------------------------------------
# Data import & preparation
# -----------------------------------------------------------------------------------------------------------------
source("source/misc_functions.R")
species_results <- read.csv('model_output/species_results.csv') %>% dplyr::select(-X)
traits <- readxl::read_xlsx('model_data/ELEData/TraitData/AVONET2_eBird.xlsx', sheet = "AVONET2_eBird")
colours <-  read.csv('model_data/plumage_lightness.csv') %>%
      mutate(species = gsub(pattern = "_", replacement = " ", species)) %>% 
      filter(sex == "F")
# Higher reflectance = lighter colors = like hotter areas, lesser heat load 
# Darker birds are limited from occupying areas with high temperatures

species_results$species[species_results$species == "Estrilda erythronotos"] <- "Brunhilda erythronotos"
species_results$species[species_results$species == "Riparia cincta"] <- "Neophedina cincta"
species_results$species[species_results$species == "Sylvia boehmi"] <- "Curruca boehmi"
species_results$species[species_results$species == "Sylvia communis"] <- "Curruca communis"
species_results$species[species_results$species == "Turdoides aylmeri"] <- "Argya aylmeri"

species_results$temp_breadth <- species_results$temp_max - species_results$temp_min
species_results$rain_breadth <- species_results$rain_max - species_results$rain_min
species_results$dry_breadth <- species_results$dry_max - species_results$dry_min
species_results$HFP_breadth <- species_results$HFP_max - species_results$HFP_min
species_results$BG_breadth <- species_results$BG_max - species_results$BG_min

# Add trait data from Tobias et al. 2022
df_merged <- merge(species_results, traits, all.x = TRUE, all.y = FALSE, by.x = "species", by.y = "Species2") %>% 
      mutate(Trophic.Level = factor(Trophic.Level, levels = c("Omnivore", "Carnivore", "Herbivore")), 
             Primary.Lifestyle = as.factor(Primary.Lifestyle))
df_merged <- merge(df_merged, colours, all.x = TRUE,by = "species")

df_merged$Migratory_ability[df_merged$Migration == 1] <- "low"
df_merged$Migratory_ability[df_merged$Migration == 2] <- "moderate"
df_merged$Migratory_ability[df_merged$Migration == 3] <- "high"
df_merged <- df_merged %>% 
      mutate(Migratory_ability = as.factor(Migratory_ability))

model_data <- df_merged %>% 
      mutate(BG_imp = scale(BG_imp),
             rain_imp = scale(rain_imp),
             temp_imp = scale(temp_imp),
             dry_imp = scale(dry_imp),
             HFP_imp = scale(HFP_imp),
             Wing.Length = scale(Wing.Length),
             avg.r = scale(avg.r),
             rain_breadth = scale(rain_breadth),
             temp_breadth = scale(temp_breadth),
             dry_breadth = scale(dry_breadth),
             HFP_breadth = scale(HFP_breadth),
             BG_breadth = scale(BG_breadth),
             log_cells_lost_norm = log(cells_lost_norm),
             log_cells_colonised_norm = log(cells_colonised_norm)) %>% 
      filter(Trophic.Level != "Scavenger",
             species != "Neophron percnopterus") %>% 
      mutate(Migratory_ability = factor(Migratory_ability, levels = c("low", "moderate", "high")),
             Primary.Lifestyle = factor(Primary.Lifestyle, levels = c("Insessorial", "Terrestrial", "Aerial", "Generalist")),
             Trophic.Level = factor(Trophic.Level, levels = c("Carnivore", "Herbivore", "Omnivore"))) %>% 
      dplyr::select(Trophic.Level, Primary.Lifestyle, Migratory_ability,
                    BG_imp, rain_imp, temp_imp,
                    dry_imp, HFP_imp, Wing.Length, avg.r,
                    rain_breadth, temp_breadth, dry_breadth,
                    HFP_breadth, BG_breadth, 
                    log_cells_lost_norm, log_cells_colonised_norm, log_prop_change)

# -----------------------------------------------------------------------------------------------------------------
# Linear combinations for INLA model
# -----------------------------------------------------------------------------------------------------------------

BG_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                            BG_imp = seq(min(model_data$BG_imp), max(model_data$BG_imp), len = 100))
names(BG_lc) <- paste0("BG_imp", 1:100)
rain_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              rain_imp = seq(min(model_data$rain_imp), max(model_data$rain_imp), len = 100))
names(rain_lc) <- paste0("rain_imp", 1:100)
temp_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              temp_imp = seq(min(model_data$temp_imp), max(model_data$temp_imp), len = 100))
names(temp_lc) <- paste0("temp_imp", 1:100)
dry_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             dry_imp = seq(min(model_data$dry_imp), max(model_data$dry_imp), len = 100))
names(dry_lc) <- paste0("dry_imp", 1:100)
HFP_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                             HFP_imp = seq(min(model_data$HFP_imp), max(model_data$HFP_imp), len = 100))
names(HFP_lc) <- paste0("HFP_imp", 1:100)
wing_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              Wing.Length = seq(min(model_data$Wing.Length), max(model_data$Wing.Length), len = 100))
names(wing_lc) <- paste0("Wing.Length", 1:100)
reflect_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                 avg.r = seq(min(model_data$avg.r, na.rm = TRUE), max(model_data$avg.r, na.rm = TRUE), len = 100))
names(reflect_lc) <- paste0("avg.r", 1:100)
Trophic_level_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
                                       Trophic.LevelHerbivore = c(0, 1, 0),
                                       Trophic.LevelOmnivore = c(0, 0, 1))
names(Trophic_level_lc) <- paste0("Trophic.Level", 1:3)
Migratory_ability_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
                                           Migratory_abilitymoderate = c(0, 1, 0),
                                           Migratory_abilityhigh = c(0, 0, 1))
names(Migratory_ability_lc) <- paste0("Migratory_ability", 1:3)
Primary.Lifestyle_lc <- inla.make.lincombs("(Intercept)" = rep(1, 4),
                                           Primary.LifestyleTerrestrial = c(0, 1, 0, 0),
                                           Primary.LifestyleAerial = c(0, 0, 1, 0),
                                           Primary.LifestyleGeneralist = c(0, 0, 0, 1))
names(Primary.Lifestyle_lc) <- paste0("Primary.Lifestyle", 1:4)


LC_no_lifestyle_sens <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                           Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                           Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                           Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                           Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                           BG_imp = model_data$BG_imp,
                                           rain_imp = model_data$rain_imp,
                                           temp_imp = model_data$temp_imp,
                                           dry_imp = model_data$dry_imp,
                                           HFP_imp = model_data$HFP_imp,
                                           Wing.Length = model_data$Wing.Length)
names(LC_no_lifestyle_sens) <- paste0("LC_no_lifestyle_sens", 1:NROW(model_data))

LC_no_trophic_sens <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                         Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                         Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                         Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                         Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                         Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                         BG_imp = model_data$BG_imp,
                                         rain_imp = model_data$rain_imp,
                                         temp_imp = model_data$temp_imp,
                                         dry_imp = model_data$dry_imp,
                                         HFP_imp = model_data$HFP_imp,
                                         Wing.Length = model_data$Wing.Length)
names(LC_no_trophic_sens) <- paste0("LC_no_trophic_sens", 1:NROW(model_data))

LC_no_migratory_sens <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                           Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                           Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                           Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                           Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                           Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                           BG_imp = model_data$BG_imp,
                                           rain_imp = model_data$rain_imp,
                                           temp_imp = model_data$temp_imp,
                                           dry_imp = model_data$dry_imp,
                                           HFP_imp = model_data$HFP_imp,
                                           Wing.Length = model_data$Wing.Length)
names(LC_no_migratory_sens) <- paste0("LC_no_migratory_sens", 1:NROW(model_data))

LC_no_BG_imp <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                   Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                   Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                   Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                   Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                   Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                   Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                   Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                   rain_imp = model_data$rain_imp,
                                   temp_imp = model_data$temp_imp,
                                   dry_imp = model_data$dry_imp,
                                   HFP_imp = model_data$HFP_imp,
                                   Wing.Length = model_data$Wing.Length)
names(LC_no_BG_imp) <- paste0("LC_no_BG_imp", 1:NROW(model_data))

all_lc_sens <- c(BG_lc, rain_lc, temp_lc, dry_lc, HFP_lc, wing_lc, reflect_lc, 
                 Trophic_level_lc, Migratory_ability_lc, Primary.Lifestyle_lc,
                 LC_no_lifestyle_sens, LC_no_trophic_sens, LC_no_migratory_sens,
                 LC_no_BG_imp)


BG_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                  BG_breadth = seq(min(model_data$BG_breadth), max(model_data$BG_breadth), len = 100))
names(BG_lc_niche) <- paste0("BG_breadth", 1:100)
rain_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                    rain_breadth = seq(min(model_data$rain_breadth), max(model_data$rain_breadth), len = 100))
names(rain_lc_niche) <- paste0("rain_breadth", 1:100)
temp_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                    temp_breadth = seq(min(model_data$temp_breadth), max(model_data$temp_breadth), len = 100))
names(temp_lc_niche) <- paste0("temp_breadth", 1:100)
dry_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                   dry_breadth = seq(min(model_data$dry_breadth), max(model_data$dry_breadth), len = 100))
names(dry_lc_niche) <- paste0("dry_breadth", 1:100)
HFP_lc_niche <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                   HFP_breadth = seq(min(model_data$HFP_breadth), max(model_data$HFP_breadth), len = 100))
names(HFP_lc_niche) <- paste0("HFP_breadth", 1:100)
wing_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                              Wing.Length = seq(min(model_data$Wing.Length), max(model_data$Wing.Length), len = 100))
names(wing_lc) <- paste0("Wing.Length", 1:100)
reflect_lc <- inla.make.lincombs("(Intercept)" = rep(1, 100),
                                 avg.r = seq(min(model_data$avg.r, na.rm = TRUE), max(model_data$avg.r, na.rm = TRUE), len = 100))
names(reflect_lc) <- paste0("avg.r", 1:100)
Trophic_level_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
                                       Trophic.LevelHerbivore = c(0, 1, 0),
                                       Trophic.LevelOmnivore = c(0, 0, 1))
names(Trophic_level_lc) <- paste0("Trophic.Level", 1:3)
Migratory_ability_lc <- inla.make.lincombs("(Intercept)" = rep(1, 3),
                                           Migratory_abilitymoderate = c(0, 1, 0),
                                           Migratory_abilityhigh = c(0, 0, 1))
names(Migratory_ability_lc) <- paste0("Migratory_ability", 1:3)
Primary.Lifestyle_lc <- inla.make.lincombs("(Intercept)" = rep(1, 4),
                                           Primary.LifestyleTerrestrial = c(0, 1, 0, 0),
                                           Primary.LifestyleAerial = c(0, 0, 1, 0),
                                           Primary.LifestyleGeneralist = c(0, 0, 0, 1))
names(Primary.Lifestyle_lc) <- paste0("Primary.Lifestyle", 1:4)

LC_no_lifestyle_niche <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                            Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                            Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                            Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                            Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                            BG_breadth = model_data$BG_breadth,
                                            rain_breadth = model_data$rain_breadth,
                                            temp_breadth = model_data$temp_breadth,
                                            dry_breadth = model_data$dry_breadth,
                                            HFP_breadth = model_data$HFP_breadth,
                                            Wing.Length = model_data$Wing.Length)
names(LC_no_lifestyle_niche) <- paste0("LC_no_lifestyle_niche", 1:NROW(model_data))

LC_no_trophic_niche <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                          Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                          Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                          Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                          Migratory_abilitymoderate = (model_data$Migratory_ability == "moderate")*1,
                                          Migratory_abilityhigh = (model_data$Migratory_ability == "high")*1,
                                          BG_breadth = model_data$BG_breadth,
                                          rain_breadth = model_data$rain_breadth,
                                          temp_breadth = model_data$temp_breadth,
                                          dry_breadth = model_data$dry_breadth,
                                          HFP_breadth = model_data$HFP_breadth,
                                          Wing.Length = model_data$Wing.Length)
names(LC_no_trophic_niche) <- paste0("LC_no_trophic_niche", 1:NROW(model_data))

LC_no_migratory_niche <- inla.make.lincombs("(Intercept)" = rep(1, NROW(model_data)),
                                            Trophic.LevelHerbivore = (model_data$Trophic.Level == "Herbivore")*1,
                                            Trophic.LevelOmnivore = (model_data$Trophic.Level == "Omnivore")*1,
                                            Primary.LifestyleTerrestrial = (model_data$Primary.Lifestyle == "Terrestrial")*1,
                                            Primary.LifestyleAerial = (model_data$Primary.Lifestyle == "Aerial")*1,
                                            Primary.LifestyleGeneralist = (model_data$Primary.Lifestyle == "Generalist")*1,
                                            BG_breadth = model_data$BG_breadth,
                                            rain_breadth = model_data$rain_breadth,
                                            temp_breadth = model_data$temp_breadth,
                                            dry_breadth = model_data$dry_breadth,
                                            HFP_breadth = model_data$HFP_breadth,
                                            Wing.Length = model_data$Wing.Length)
names(LC_no_migratory_niche) <- paste0("LC_no_migratory_niche", 1:NROW(model_data))

all_lc_niche <- c(BG_lc_niche, rain_lc_niche, temp_lc_niche, dry_lc_niche, HFP_lc_niche, wing_lc, reflect_lc, 
                  Trophic_level_lc, Migratory_ability_lc, Primary.Lifestyle_lc,
                  LC_no_lifestyle_niche, LC_no_trophic_niche, LC_no_migratory_niche)

save(model_data, all_lc_sens, all_lc_niche, file = 'model_data/regression_data.RData')