sortBase <- function(vec, n.knots = 2) {
   ## Function to calculate bases for regression splines. Modified from code
   ## provided in Crainiceanu, C., Ruppert, D. & Wand, M.P. Bayesian analysis for
   ## penalized spline regression using WinBUGS. J. Stat. Soft. 14, 1 24(2005).
   ## Parameter vec is a vector defining the raw data vector, n.knots defines the
   ## number of knots in the GAM.
   ## Creates transformation of covariates with weights, full length vectors
   # Define sample size
   N         <- length(vec)   
   x.time    <- c(vec)  
   # Define the X matrix of fixed effects for the thin-plate spline
   zFE       <- cbind(rep(1,N), x.time)   
   # Define knot values as the sample quantiles of the covariate, based on the number of knots.
   # With two knots, defines knot values as covariate values at the 33.3% and 66.6% quantile.
   x.knots   <- quantile(unique(x.time), seq(0, 1, length = (n.knots+2))[-c(1, (n.knots+2))], na.rm = TRUE)
   
   z_K       <- (abs(outer(x.time,x.knots,"-")))^3
   OMEGA.all <- (abs(outer(x.knots,x.knots,"-")))^3
   svd.OMEGA.all  <- svd(OMEGA.all)
   sqrt.OMEGA.all <- t(svd.OMEGA.all$v %*% (t(svd.OMEGA.all$u) *
                                               sqrt(svd.OMEGA.all$d)))
   z.out     <- t(solve(sqrt.OMEGA.all, t(z_K)))
   return(z.out)
}


make_ebird_sp <- function(scientific_name, ROI){
   require(spatialEco)
      df <- ebird_filtered %>% 
            group_by(LATITUDE, LONGITUDE, `SAMPLING EVENT IDENTIFIER`, `DURATION MINUTES`, 
                     `EFFORT DISTANCE KM`, `NUMBER OBSERVERS`, `OBSERVATION DATE`, `LOCALITY`) %>% 
            summarise(occurrence = ifelse(scientific_name %in% `SCIENTIFIC NAME`, TRUE, FALSE)) %>% 
            ungroup() %>% 
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
                              effort_distance_km = df$effort_distance_km, number_observers = df$number_observers),
            proj4string = crs(proj)
      )
      
      # Only include eBird data points for the region of interest
      # Get intersecting points
      in_sp <- rgeos::gIntersection(df_sp, ROI)
      
      # Only keep intersecting points in original spdf
      df_sp_ROI <- df_sp[in_sp, ]
      cat('\nNumber of eBird presence records: ', length(rownames(df_sp_ROI@data[df_sp_ROI@data$presence == TRUE, ])), '\n')
      
      
      return(df_sp_ROI)
}

get_hyperpars <- function(path){
   out = NULL
   list <- list.files(path = path)
   
   for(i in 1:NROW(list)){
      file = list[i]
      
      load(paste0(path, file))
      
      # Check hyperparameters. sd ideally should be somewhat smaller than mean. Stdev should not be very small
      outdf <- as.data.frame(model[["model"]][["summary.hyperpar"]])
      outdf$par <- c("Range", "Stdev")
      outdf <- cbind(outdf, rbind(data, data))
      out <-  rbind(out, outdf)
   }
   return(out)
}

plot_INLA_predictions <- function(model_out){
   load(model_out)
   cloglog_inv <- function(x){
      1-exp(-exp(x))
   }
   
   proj <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
   
   Pred <- SpatialPixelsDataFrame( 
      points = stk.pred$predcoords,
      data = model$predictions,
      proj4string = crs(proj)
   )
   
   Pred@data$inv_mean <- cloglog_inv(Pred@data$mean)
   # test <- raster(Pred)
   # writeRaster(test, "TZ_fsl_predictions.tif")
   # Plot of predictions
   plot_title <- str_split(model_out, pattern = "/")[[1]][length(str_split(model_out, pattern = "/")[[1]])]

   p1 <- ggplot() +
      gg(Pred, aes(fill = inv_mean)) +
      coord_equal() +
      theme_void() +
      ggtitle(plot_title) + 
      labs(fill = "Probability") 
   return(p1)
}



plot_INLA_random <- function(model_out, ROI){
   load(model_out)
   # Visualize the random field, the unexplained spatial error
   projector <- inla.mesh.projector(Mesh$mesh)
   projection <- inla.mesh.project(projector,
                                   model[["model"]][["summary.random"]]$i$mean)
   plot_title <- str_split(model_out, pattern = "/")[[1]][length(str_split(model_out, pattern = "/")[[1]])]
   p1 <- INLAutils::ggplot_projection_shapefile(projection, projector) +
      geom_polygon(data = ROI, aes(long, lat, group = group), fill = NA, col = 'black') +
      coord_equal() +
      ggtitle(plot_title)
   
   return(p1)
}

plot_INLA_fit <- function(model_out){
   require(cowplot)
   require(gridGraphics)
   load(model_out)
   
   failure <- sum(model$model$cpo$failure) 
   print(paste0("Number of CPO failures: ", failure))
   data <- data.frame(pit = model$model$cpo$pit)
   p1 <- ggplot(data) +
      geom_histogram(aes(x = pit))
   # hist(model$model$cpo$pit, main="", breaks = 10, xlab = "PIT")
   # p1 <- recordPlot()  
   # 
   # plot.ecdf(model$model$cpo$pit); abline(0, 1)
   # p2 <- recordPlot()  
   # 
   # p3 <- plot_grid(p1, p2, ncol = 2)
   return(p1)
}

prepare_GAM <- function(df, vector){
      # Create prediction vector, as sequence of 100 values, spanning full range of values. 
      sequence <- seq(min(vector, na.rm = TRUE), max(vector, na.rm = TRUE), length = 100)
      # Save this vector, we use it later to label the x axis in the model effect plots
      assign(paste0(deparse(substitute(vector)), ".seq"), sequence, envir = .GlobalEnv)
      
      # Add prediction vector to beginning of original values, scale all to mean 0 and sd 1
      var.s <- c(scale(c(sequence, vector)))  
      # Split the full vector into two spline bases, and scale the new vectors to mean 0 and sd 1.
      # Scale again, covariates should be scaled for the INLA model
      z.var.s <- scale(sortBase(var.s, n.knots = 2))  
      
      # Separate out prediction vectors for the two spline bases, and assign to environment as
      # separate objects. These will be used for linear combinations
      assign(paste0(deparse(substitute(vector)), "_1.s"), c(z.var.s[1:100, 1]), envir = .GlobalEnv)
      assign(paste0(deparse(substitute(vector)), "_2.s"), c(z.var.s[1:100, 2]), envir = .GlobalEnv)
      
      # Need to separate out the two time periods, and add to the main data, for the model.
      # Get the length of each time period vector
      n_data <- NROW(vector)/2
      
      df@data[, paste0(deparse(substitute(vector)), "_1980s_1.s")] <- c(z.var.s[101:(100+n_data), 1])
      df@data[, paste0(deparse(substitute(vector)), "_1980s_2.s")] <- c(z.var.s[101:(100+n_data), 2])
      
      df@data[, paste0(deparse(substitute(vector)), "_2000s_1.s")] <- c(z.var.s[(n_data+101):NROW(z.var.s), 1])
      df@data[, paste0(deparse(substitute(vector)), "_2000s_2.s")] <- c(z.var.s[(n_data+101):NROW(z.var.s), 2])
      return(df)
}

unscale <- function(x, scale.params = sc.p) {
      return((x * scale.params$'scaled:scale') + scale.params$'scaled:center')
}

cloglog_inv <- function(x){
      1-exp(-exp(x))
}

# make_lincomb_contrasts <- function(data, colname){
#       k <- nlevels(data[, colname])
#       if(k < 3 | k > 4){
#             stop("Currently only works for factors with 3 or 4 levels.")
#       }
#       if (k == 4) {
#             ncombs = 6
#             LC <- inla.make.lincombs(lc1 = c( 1,  1,  1,  0,  0, 0),
#                                      lc2 = c(-1,  0,  0,  1,  1, 0),
#                                      lc3 = c( 0, -1,  0, -1,  0, 1),
#                                      lc4 = c( 0,  0, -1,  0, -1,-1))
#       }
#       
#       if (k == 3) {
#             ncombs = 3
#             LC <- inla.make.lincombs(lc1 = c( 1,  1,  0),
#                                      lc2 = c(-1,  0, -1),
#                                      lc3 = c( 0, -1,  1))
#       }
# 
#       for(i in 1:ncombs){
#             for (j in 1:k){
#                   names(LC[[paste0("lc", i)]][[j]]) <- paste0(colname, levels(data[, colname])[j])
#             }
#       }
#       
#       # Assign names to combinations
#       t <- 0
#       for(i in 1:(k-1)) {
#             for(j in (i+1):k) {
#                   tryCatch({
#                         t <- t + 1
#                         names(LC)[t] <- paste0(levels(data[, colname])[i], "-", levels(data[, colname])[j])
#                   }, error=function(e){})
#             }
#       }
#       return(LC)
# }
# 
# 
# 
# lincomb_output <- function(model, data){
#       lincomb_out_df <- data.frame(model$summary.lincomb.derived) %>% 
#             mutate(ID = gsub('[[:digit:]]+', '', rownames(.))) %>% 
#             dplyr::select(-mode, -kld, -mean, -sd) %>% 
#             pivot_wider(names_from = ID, values_from = c(2:4),
#                         names_glue = "{ID}_{.value}",
#                         values_fn = list) %>% 
#             unnest()
#       
#       lincomb_out_df[, strsplit(colnames(lincomb_out_df)[grep("BG", colnames(lincomb_out_df))], split = "_X")[[1]][1]] <- 
#             seq(min(data[, strsplit(colnames(lincomb_out_df)[grep("BG", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 max(data[, strsplit(colnames(lincomb_out_df)[grep("BG", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 len = 100)
#       lincomb_out_df[, strsplit(colnames(lincomb_out_df)[grep("rain", colnames(lincomb_out_df))], split = "_X")[[1]][1]] <- 
#             seq(min(data[, strsplit(colnames(lincomb_out_df)[grep("rain", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 max(data[, strsplit(colnames(lincomb_out_df)[grep("rain", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 len = 100)
#       lincomb_out_df[, strsplit(colnames(lincomb_out_df)[grep("temp", colnames(lincomb_out_df))], split = "_X")[[1]][1]] <- 
#             seq(min(data[, strsplit(colnames(lincomb_out_df)[grep("temp", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 max(data[, strsplit(colnames(lincomb_out_df)[grep("temp", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 len = 100)
#       lincomb_out_df[, strsplit(colnames(lincomb_out_df)[grep("dry", colnames(lincomb_out_df))], split = "_X")[[1]][1]] <- 
#             seq(min(data[, strsplit(colnames(lincomb_out_df)[grep("dry", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 max(data[, strsplit(colnames(lincomb_out_df)[grep("dry", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 len = 100)
#       lincomb_out_df[, strsplit(colnames(lincomb_out_df)[grep("HFP", colnames(lincomb_out_df))], split = "_X")[[1]][1]] <- 
#             seq(min(data[, strsplit(colnames(lincomb_out_df)[grep("HFP", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 max(data[, strsplit(colnames(lincomb_out_df)[grep("HFP", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 len = 100)
#       lincomb_out_df[, strsplit(colnames(lincomb_out_df)[grep("Wing.Length", colnames(lincomb_out_df))], split = "_X")[[1]][1]] <- 
#             seq(min(data[, strsplit(colnames(lincomb_out_df)[grep("Wing.Length", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 max(data[, strsplit(colnames(lincomb_out_df)[grep("Wing.Length", colnames(lincomb_out_df))], split = "_X")[[1]][1]]),
#                 len = 100)
#       lincomb_out_df[, strsplit(colnames(lincomb_out_df)[grep("avg.r", colnames(lincomb_out_df))], split = "_X")[[1]][1]] <- 
#             seq(min(data[, strsplit(colnames(lincomb_out_df)[grep("avg.r", colnames(lincomb_out_df))], split = "_X")[[1]][1]], na.rm = TRUE),
#                 max(data[, strsplit(colnames(lincomb_out_df)[grep("avg.r", colnames(lincomb_out_df))], split = "_X")[[1]][1]], na.rm = TRUE),
#                 len = 100)
#       lincomb_out_df[, strsplit(colnames(lincomb_out_df)[grep("Trophic", colnames(lincomb_out_df))], split = "_X")[[1]][1]] <- 
#             seq(min(data[, strsplit(colnames(lincomb_out_df)[grep("Trophic", colnames(lincomb_out_df))], split = "_X")[[1]][1]], na.rm = TRUE),
#                 max(data[, strsplit(colnames(lincomb_out_df)[grep("Trophic", colnames(lincomb_out_df))], split = "_X")[[1]][1]], na.rm = TRUE),
#                 len = 100)
#       
#       return(lincomb_out_df)
# }

make_forest_plot <- function(inla_model, label, col){
      model_fixed <- data.frame(inla_model$summary.fixed) %>% 
            mutate(ID = rownames(.))
      
      if(label == "Niche breadth"){
            model_fixed <- model_fixed %>% 
                  mutate(ID = replace(ID, ID == "(Intercept)", "Intercept"),
                         ID = replace(ID, ID == "Trophic.LevelOmnivore", "Trophic level: Omnivore"),
                         ID = replace(ID, ID == "Trophic.LevelHerbivore", "Trophic level: Herbivore"),
                         ID = replace(ID, ID == "Primary.LifestyleGeneralist", "Primary lifestyle: Generalist"),
                         ID = replace(ID, ID == "Primary.LifestyleAerial", "Primary lifestyle: Aerial"),
                         ID = replace(ID, ID == "Primary.LifestyleTerrestrial", "Primary lifestyle: Terrestrial"),
                         ID = replace(ID, ID == "Migratory_abilityhigh", "Migratory ability: High"),
                         ID = replace(ID, ID == "Migratory_abilitymoderate", "Migratory ability: Moderate"),
                         ID = replace(ID, ID == "BG_breadth", "Bare ground cover"),
                         ID = replace(ID, ID == "rain_breadth", "Annual rainfall"),
                         ID = replace(ID, ID == "temp_breadth", "Hottest temperature"),
                         ID = replace(ID, ID == "dry_breadth", "Dry spell duration"),
                         ID = replace(ID, ID == "HFP_breadth", "Human footprint"),
                         ID = replace(ID, ID == "Wing.Length", "Wing length"),
                         ID = replace(ID, ID == "avg.r", "Dorsal reflectance"))
      }
      
      if(label == "Sensitivity"){
            model_fixed <- model_fixed %>% 
                  mutate(ID = replace(ID, ID == "(Intercept)", "Intercept"),
                         ID = replace(ID, ID == "Trophic.LevelOmnivore", "Trophic level: Omnivore"),
                         ID = replace(ID, ID == "Trophic.LevelHerbivore", "Trophic level: Herbivore"),
                         ID = replace(ID, ID == "Primary.LifestyleGeneralist", "Primary lifestyle: Generalist"),
                         ID = replace(ID, ID == "Primary.LifestyleAerial", "Primary lifestyle: Aerial"),
                         ID = replace(ID, ID == "Primary.LifestyleTerrestrial", "Primary lifestyle: Terrestrial"),
                         ID = replace(ID, ID == "Migratory_abilityhigh", "Migratory ability: High"),
                         ID = replace(ID, ID == "Migratory_abilitymoderate", "Migratory ability: Moderate"),
                         ID = replace(ID, ID == "BG_imp", "Bare ground cover"),
                         ID = replace(ID, ID == "rain_imp", "Annual rainfall"),
                         ID = replace(ID, ID == "temp_imp", "Hottest temperature"),
                         ID = replace(ID, ID == "dry_imp", "Dry spell duration"),
                         ID = replace(ID, ID == "HFP_imp", "Human footprint"),
                         ID = replace(ID, ID == "Wing.Length", "Wing length"),
                         ID = replace(ID, ID == "avg.r", "Dorsal reflectance"))
      }
            
      model_fixed$ID <- factor(model_fixed$ID, levels = c("Trophic level: Omnivore", "Trophic level: Herbivore", 
                                                          "Primary lifestyle: Generalist", "Primary lifestyle: Aerial", "Primary lifestyle: Terrestrial",
                                                          "Migratory ability: High", "Migratory ability: Moderate", "Wing length", "Dorsal reflectance",
                                                          "Bare ground cover", "Annual rainfall", "Hottest temperature", "Dry spell duration", "Human footprint",
                                                          "Intercept"))
      model_fixed$significant <- ifelse((model_fixed$X0.025quant > 0 & model_fixed$X0.975quant > 0)|(model_fixed$X0.025quant < 0 & model_fixed$X0.975quant < 0), "yes", "no")
      
      forest_plot <- ggplot() + 
            annotate("text", x = 7, y = max(inla_model$summary.fixed$`0.975quant`) + 0.1*max(inla_model$summary.fixed$`0.975quant`), 
                     col = col, label = label,
                     hjust = 0) +
            geom_rect(aes(xmin = 1.5, xmax = 6.5, ymin = min(inla_model$summary.fixed$`0.025quant`), 
                          ymax = max(inla_model$summary.fixed$`0.975quant`)),
                      fill = col, col = col, alpha = .3) +
            geom_point(data = model_fixed, aes(x = ID, y = X0.5quant, col = significant)) +
            geom_errorbar(data = model_fixed, aes(x = ID, ymin = X0.025quant, ymax = X0.975quant, col = significant), width = 0.1) +
            geom_hline(aes(yintercept = 0), lty = 2, alpha = .3) +
            # geom_text(data = model_fixed %>% filter(significant == "yes"), aes(x = ID, y = X0.975quant), label = "*", nudge_x = 0.1, nudge_y = 0.01) +
            ylab("Posterior estimates") +
            xlab(element_blank()) +
            theme_minimal() +
            scale_x_discrete(limits=rev) +
            scale_colour_manual(values = c("black", "#D55E00"), guide = "none") +
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
      return(forest_plot)
}

rsq <- function(x, y) summary(lm(y~x))$r.squared

make_raw_data_plots <- function(inla_model, data, xlabel){
      covariates <- strsplit(as.character(inla_model$.args$formula[3]), split = " \\+ ")[[1]]
      response <- as.character(inla_model$.args$formula[2])
      plot_list <- vector('list', length(covariates))
      
      lincomb_data <- data.frame(inla_model$summary.lincomb.derived) %>% 
            mutate(ID = rownames(.)) %>% 
            dplyr::select(-mode, -kld, -mean, -sd)
      
      for (i in 1:length(covariates)){
            if(is.numeric(data[, covariates[i]])){
                  plot_list[[i]] <- local({
                        i <- i
                        p1 <- ggplot() +
                              geom_point(data = data, aes(x = get(covariates[i]), y = get(response))) +
                              geom_line(data = lincomb_data[grepl(covariates[i], rownames(lincomb_data)) == TRUE, ],
                                        aes(x = seq(min(data[, covariates[i]], na.rm = T), max(data[, covariates[i]], na.rm = T), len = 100),
                                            y = X0.5quant), lty = 2, alpha = 0.5, col = "red") +
                              geom_ribbon(data = lincomb_data[grepl(covariates[i], rownames(lincomb_data)) == TRUE, ],
                                          aes(x = seq(min(data[, covariates[i]], na.rm = T), max(data[, covariates[i]], na.rm = T), len = 100), 
                                              ymin = X0.025quant, 
                                              ymax = X0.975quant), 
                                          alpha = 0.1) +
                              theme_classic() +
                              xlab(covariates[i]) +
                              ylab(xlabel) + 
                              geom_rug(data = data, aes(x = get(covariates[i]), y = get(response)), col=rgb(.5,0,0, alpha=.4))
                        print(p1)})
            } else {
                  plot_list[[i]] <- local({
                        i <- i
                        p2 <- ggplot(data) +
                              geom_boxplot(aes(get(covariates[i]), get(response))) +
                              geom_point(data = lincomb_data[grepl(covariates[i], rownames(lincomb_data)) == TRUE, ], 
                                         aes(x = levels(data[, covariates[i]]), y = X0.5quant), col = "red", cex = 2) +
                              theme_classic() +
                              xlab(covariates[i]) +
                              ylab(xlabel)
                        print(p2)})
            }
      }
      
      return(plot_list)
}
