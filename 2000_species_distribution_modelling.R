options(echo=TRUE)

fit_maxent2 <- function(occ, bg, name, predictors, outdir, template, max_bg_size = 100000, 
                        bg_buffer_width = 1000000, shapefiles = TRUE, features, 
                        replicates, responsecurves = TRUE, rep_args, full_args, random_args,
                        evaluation, sensitivity = 0.9) {
  require(ff)
  require(rgeos)
  require(sp)
  require(raster)
  require(rJava)
  require(dismo)
  require(things)
  
  if (!file.exists(outdir)) stop('outdir does not exist :(', call. = FALSE)
  outdir_sp <- file.path(outdir, gsub(' ', '_', name))
  if (!file.exists(outdir_sp)) dir.create(outdir_sp)
  features <- unlist(strsplit(features, ''))
  if (length(setdiff(features, c('l', 'p', 'q', 'h', 't'))) > 1)
    stop("features must be a vector of one or more of ',
         'l', 'p', 'q', 'h', and 't'.")
  
  b <- aggregate(gBuffer(occ, width = bg_buffer_width, byid = TRUE))
  #  proj4string(b) <- proj4string(bg)
  
  if (length(occ) < 30) {
    warning('Fewer occurrence records than the number of cross-validation ',
            'replicates for species ', name, 
            '. Model not fit for this species.')
    return(NULL)
  }
  
  # subset bg to the buffer poly (bg = 200 km, bg_poly = IBRA poly)
  
  koppen_occ <- gIntersects(koppen,occ,byid = TRUE)
  n <- colSums(koppen_occ)
  koppen_subset <- koppen[n > 5,]
  koppen_subset_clean <- clgeo_Clean(koppen_subset)
  kop_buf <- gIntersection(koppen_subset_clean,b)
  bg <- spTransform(bg, crs(kop_buf))
  bg_kop_buf <- over(bg, kop_buf) # SLOW
  bg_kop_buf_nona <- data.frame(bg[!is.na(bg_kop_buf),])
  crp <- raster::crop(template, b)
  crp1 <- raster::mask(crp, b)
  random <- dismo::randomPoints(crp1, n=nrow(bg_kop_buf_nona))
  random_swd <- raster::extract(s, random)
  speciesName <- rep(gsub(' ', '_', name), nrow(random_swd))
  random_swd1 <- cbind.data.frame(speciesName, random, random_swd)
  random_swd_spdf <- SpatialPointsDataFrame(random_swd1[,2:3], random_swd1, proj4string=crs(s.crs))
  
  # Reduce background sample if it's larger than max_bg_size
  if (nrow(bg_kop_buf_nona) > max_bg_size) {
    message(nrow(bg_kop_buf_nona), ' target species background records, reduced to random ',
            max_bg_size, '.')
    bg_2 <- bg_kop_buf_nona[sample(nrow(bg_kop_buf_nona), max_bg_size), ]} else {
      message(nrow(bg_kop_buf_nona), ' target species background records.')
      bg_2 <-bg_kop_buf_nona
    }
  
  
  # Save objects for future reference
  if (shapefiles) {
    suppressWarnings({
      writeOGR(SpatialPolygonsDataFrame(b, data.frame(ID = seq_len(length(b)))), 
               outdir_sp, 'bg_buffer', 'ESRI Shapefile', overwrite_layer = TRUE)  
      writeOGR(occ, outdir_sp, 'occ', 'ESRI Shapefile', overwrite_layer = TRUE)
      writeOGR(koppen_subset_clean, outdir_sp, 'koppen', 'ESRI Shapefile', overwrite_layer = TRUE)
      writeOGR(SpatialPointsDataFrame(random_swd1[,2:3], random_swd1, proj4string=s.crs),
               outdir_sp, 'random_bkgr', 'ESRI Shapefile', overwrite_layer = TRUE)
    })
  }
  saveRDS(bg_2,  file.path(outdir_sp, 'bg.rds'))
  saveRDS(occ, file.path(outdir_sp, 'occ.rds'))
  
  # Sample predictors at occurrence and background points
  all.predictors <- names(occ@data[,-c(1:4,20)])
  swd_occ <- as.data.frame(occ[,all.predictors])
  swd_occ <- swd_occ[,-c(16:17)]
  swd_occ <-swd_occ[,predictors]
  swd_bg <- as.data.frame(bg_2[,all.predictors])
  swd_bg <- swd_bg[,predictors]
  saveRDS(swd_occ, file.path(outdir_sp, 'occ_swd.rds'))
  saveRDS(swd_bg, file.path(outdir_sp, 'bg_swd.rds'))
  swd_bg_spdf <- SpatialPointsDataFrame(bg_2[,2:3], bg_2, proj4string = crs(s.crs))
  
  if (shapefiles) {
    writeOGR(occ, outdir_sp, 'occ_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
    writeOGR(swd_bg_spdf, outdir_sp, 'bg_swd', 'ESRI Shapefile', overwrite_layer = TRUE)
  }
  
  # Combine occ and bg SWD data
  swd <- as.data.frame(rbind(swd_occ, swd_bg))
  saveRDS(swd, file.path(outdir_sp, 'swd.rds'))
  pa <- rep(1:0, c(nrow(swd_occ), nrow(swd_bg)))
  
  random_swd2 <- random_swd[,predictors]
  swd1 <- as.data.frame(rbind(swd_occ, random_swd2))
  saveRDS(swd1, file.path(outdir_sp, 'swd_random.rds'))
  pa1 <- rep(1:0, c(nrow(swd_occ), nrow(random_swd2)))
  
  # Fit model
  off <- setdiff(c('l', 'p', 'q', 't', 'h'), features)
  if (length(off) > 0) {
    off <- c(l = 'linear=FALSE', p = 'product=FALSE', q = 'quadratic=FALSE',
             t = 'threshold=FALSE', h = 'hinge=FALSE')[off]
  }
  off <- unname(off)
  if (replicates > 1) {
    if (missing(rep_args)) rep_args <- NULL
    me_xval <- maxent(swd, pa, path = file.path(outdir_sp, 'xval'), 
                      args = c(paste0('replicates=', replicates),
                               'responsecurves=TRUE',
                               off, paste(names(rep_args), rep_args, sep = '=')))
  }
  if (missing(full_args)) full_args <- NULL
  me_full <- maxent(swd, pa, path = file.path(outdir_sp, 'full'), 
                    args = c(off, paste(names(full_args), full_args, sep = '='),
                             'responsecurves=TRUE'))
  
  if (missing(random_args)) random_args <- NULL
  me_random <- maxent(swd1, pa1, path = file.path(outdir_sp, 'rndm_bkgr'), 
                      args = c(off, paste(names(random_args), random_args, sep = '='),
                               'responsecurves=TRUE'))
  
  if (missing(evaluation)) evaluation <- NULL
  eval_mod <- dismo::evaluate(occ, swd_bg_spdf, me_full, predictors = predictors, type = "logistic")
  eval_df <- data.frame(threshold = eval_mod@t, 
                        prevalence = eval_mod@prevalence, 
                        TPR = eval_mod@TPR, TNR = eval_mod@TNR, FPR = eval_mod@FPR, FNR = eval_mod@FNR, 
                        PPP = eval_mod@PPP, NPP = eval_mod@NPP, MCR = eval_mod@MCR, 
                        kappa = eval_mod@kappa, TSS = (eval_mod@TPR + eval_mod@TNR) - 1)
  th_table <- dismo::threshold(eval_mod, sensitivity = sensitivity)
  th_table$species_name <- name
  th_table$n.presencen <- eval_mod@np
  th_table$n.absencen <- eval_mod@na
  th_table$AUC <- eval_mod@auc
  th_table$TSS10th <- min(eval_df$TSS[which(abs(eval_df$TPR - 0.9)==min(abs(eval_df$TPR - 0.9)))])
  th_table$TSSmax <- max(eval_df$TSS)
  th_table$Kappa10th <- min(eval_mod@kappa[which(abs(eval_mod@TPR - 0.9)==min(abs(eval_mod@TPR - 0.9)))])
  th_table$KAPPAmax <- max(eval_df$kappa)
  colnames(th_table) <- c("th_max.KAPPA", "th_max.spec+sens", "th_no.omission","prevalence",
                          "th_spec=sens", "th_sensitivity90", "species_name","n.presencen",
                          "n.absencen","AUC", "TSS10th", "max.TSS", "Kappa10th", "max.KAPPA")
  th_table <- th_table[,c(7,8,9,4,10,11,12,13,14,1,2,3,5,6)]
  rownames(th_table) <- ""
  message("writing evaluation table")
  write.csv(th_table, file = paste0(outdir_sp, 
                                    "/evaluate_", name, ".csv"), row.names=F)
  
  eval_mod_random <- dismo::evaluate(occ, random_swd_spdf, me_random, predictors = predictors, type = "logistic")
  eval_df_random <- data.frame(threshold = eval_mod_random@t, 
                               prevalence = eval_mod_random@prevalence, 
                               TPR = eval_mod_random@TPR, TNR = eval_mod_random@TNR, FPR = eval_mod_random@FPR, FNR = eval_mod_random@FNR, 
                               PPP = eval_mod_random@PPP, NPP = eval_mod_random@NPP, MCR = eval_mod_random@MCR, 
                               kappa = eval_mod_random@kappa, TSS = (eval_mod_random@TPR + eval_mod_random@TNR) - 1)
  th_table_random <- dismo::threshold(eval_mod_random, sensitivity = sensitivity)
  th_table_random$species_name <- name
  th_table_random$n.presencen <- eval_mod_random@np
  th_table_random$n.absencen <- eval_mod_random@na
  th_table_random$AUC <- eval_mod_random@auc
  th_table_random$TSS10th <- min(eval_df_random$TSS[which(abs(eval_df_random$TPR - 0.9)==min(abs(eval_df_random$TPR - 0.9)))])
  th_table_random$TSSmax <- max(eval_df_random$TSS)
  th_table_random$Kappa10th <- min(eval_mod_random@kappa[which(abs(eval_mod_random@TPR - 0.9)==min(abs(eval_mod_random@TPR - 0.9)))])
  th_table_random$KAPPAmax <- max(eval_df_random$kappa)
  colnames(th_table_random) <- c("th_max.KAPPA", "th_max.spec+sens", "th_no.omission","prevalence",
                                 "th_spec=sens", "th_sensitivity90", "species_name","n.presencen",
                                 "n.absencen","AUC", "TSS10th", "max.TSS", "Kappa10th", "max.KAPPA")
  th_table_random <- th_table_random[,c(7,8,9,4,10,11,12,13,14,1,2,3,5,6)]
  rownames(th_table_random) <- ""
  message("writing evaluation table")
  write.csv(th_table_random, file = paste0(outdir_sp, 
                                           "/evaluate_random_", name, ".csv"), row.names=F)
  
  
  # Save fitted model object, and the model-fitting data.
  saveRDS(list(me_xval = me_xval, me_full = me_full, me_random = me_random, swd = swd, pa = pa), 
          file.path(outdir_sp, 'maxent_fitted.rds'))
  
  return(invisible(NULL))
}



options(bitmapType='cairo')
## Load libraries individually for debugging purposes
library('dismo')
library('sp')
library('maptools')
library('raster')
library('ff')
library('rgdal')
library('rgeos')
library('gdalUtils')
library('rmaxent')
library('readr')
library('dplyr')
library('tidyr')
library('cleangeo')
library('readr')
library('rnaturalearth')
library('rasterVis')
library('RColorBrewer')
library('latticeExtra')
library('data.table')
library('rJava')

## Load libraries
p <- c('ff', 'things', 'raster', 'dismo', 'sp', 'latticeExtra', 
       'rgdal', 'rgeos', 'gdalUtils', 'rmaxent', 'readr', 'dplyr', 'tidyr','cleangeo',
       'readr', 'rnaturalearth', 'rasterVis', 'RColorBrewer', 'latticeExtra', 'data.table')
sapply(p, require, character.only=TRUE)


## Load environmental data and Koppen zones
template <- raster('CHELSA_RASTERS_10km/Current/biol_01_mol.tif')
koppen   <- readOGR("Given/Koppen_single_mol.shp")
files_current  <- list.files("CHELSA_RASTERS_10km/Current/",full.names = TRUE, pattern = "tif$")
s <- stack(files_current)
s <- dropLayer(s, c(7,8,18,19))
s.crs <- crs(s[[1]])


#only if we want to do the projection for AUS
files_aus  <- list.files("CHELSA_RASTERS_1km/Current_cropped_AUS",full.names = TRUE, pattern = "tif$")
s_aus <- stack(files_aus)
s_aus <- dropLayer(s_aus, c(7,8,18,19))



## Load occs, bg, and names of target species
occ <- readRDS("Given/occ_swd_SPDF_single_species_1km_mol_chelsa.rds")
bg1 <- readRDS("Given/bg_swd_df_1rec_per_spp_Chelsa.rds")

spdf.bg <- SpatialPoints(bg1[,c(2:3)], s.crs)
bg <- SpatialPointsDataFrame(spdf.bg, data.frame(bg1))
colnames(bg@data)[1] = "speciesName"
bg <- spTransform(bg, crs(occ[[2]]))

occ_red <- occ[__RANGE__]
## Set parameters
bg_buffer_width <- 1000000
max_bg_size <- 50000
group <- 'plants'

# Load background 
names(s)    <- c("biol_01",
                 "biol_02",
                 "biol_03",
                 "biol_04",
                 "biol_05",
                 "biol_06",
                 "biol_09",
                 "biol_10",
                 "biol_11",
                 "biol_12",
                 "biol_13",
                 "biol_14",
                 "biol_15",
                 "biol_16",
                 "biol_17")

#only if we want to do for AUS
names(s_aus)    <- c("biol_01",
                     "biol_02",
                     "biol_03",
                     "biol_04",
                     "biol_05",
                     "biol_06",
                     "biol_09",
                     "biol_10",
                     "biol_11",
                     "biol_12",
                     "biol_13",
                     "biol_14",
                     "biol_15",
                     "biol_16",
                     "biol_17")



## Read in predictors
gc()


## Load the function to fit Maxent models
## Fit models
lapply(c('clim1'), function(x) {
  outdir            <- file.path('__OUTDIR__/', 
                                 paste(group, 
                                       x, sep = '_'))
  if (!dir.exists(outdir)) dir.create(outdir, recursive=TRUE)
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" ) 
  clim2 <- c("biol_01","biol_04","biol_06","biol_12","biol_13","biol_14","biol_15")
  preds <- clim1
  # switch(x,
  #               clim1 = clim1,
  #               clim2 = clim2)
  mapply(function(occ, name) {
    message('Doing ', name)
    name <- gsub('\\*', '', name) # remove asterisks from names
    if (!dir.exists(file.path(outdir, gsub(' ', '_', name)))) {
      fit_maxent2(occ = occ, bg = bg, predictors = preds, name = name, 
                  outdir = outdir, 
                  template = template, shapefiles = TRUE, features = 'lpq',
                  replicates = 2) ## Change replicates if necessary
    }
  }, occ_red, names(occ_red))
})




## project model onto current scenario
###    with ff_table

ff_models_FN <- grep(
  'full', 
  list.files('__OUTDIR__/', '\\.lambdas', recursive=TRUE, full=TRUE), 
  value=TRUE)

## prepare ff table for current



#locs <- which(!is.na(sum(s_aus)[]))
locs <- which(!is.na(sum(s)[]))
locs <- cbind(cell=locs, xyFromCell(s, locs))
cells <- locs[, 'cell']

#If we want AUS
locs_aus <- which(!is.na(sum(s_aus)[]))
locs_aus <- cbind(cell_aus=locs_aus, xyFromCell(s_aus, locs_aus))
cells_aus <- locs_aus[, 'cell_aus']


library(ff)



s_ff <- ff::ff(vmode="double", dim=c(length(cells), nlayers(s)),
               filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(s)) {
  s_ff[, i] <- s[[i]][][cells]
}
colnames(s_ff) <- names(s)

s_ff_aus <- ff::ff(vmode="double", dim=c(length(cells_aus), nlayers(s_aus)),
                   filename=ff_swd_aus <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(s_aus)) {
  s_ff_aus[, i] <- s_aus[[i]][][cells_aus]
}
colnames(s_ff_aus) <- names(s_aus)

gc()



### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_current.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- s[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, s_ff[, seq_len(ncol(s_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


### apply funciton to write the prediciton map and binary map for Australia

lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_current_aus.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- s_aus[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells_aus), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, s_ff_aus[, seq_len(ncol(s_ff_aus))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells_aus] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

# projection based on randomly-selected background points
ff_models_random <- grep(
  'rndm_bkgr', 
  list.files('__OUTDIR__//', '\\.lambdas', recursive=TRUE, full=TRUE), 
  value=TRUE)

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_current_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- s[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, s_ff[, seq_len(ncol(s_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

# projection to Australia for model fitted based on randomly-selected background points 
lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_current_aus_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- s_aus[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells_aus), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, s_ff_aus[, seq_len(ncol(s_ff_aus))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells_aus] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


#project onto future scenario
## 2030 RCP 8.5
## prepare ff table for RCP 8.5 ACCESS10 2030 (GCM1)


access_8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/ACCESS10/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_ACCESS_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_ACCESS_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})



## 2050 RCP 8.5
## prepare ff table for RCP 8.5 ACCESS10 2050 (GCM1)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/ACCESS10/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_ACCESS_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_ACCESS_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})



########
## prepare ff table for RCP 8.5 ACCESS10 2070 (GCM1)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/ACCESS10/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_ACCESS_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_ACCESS_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


#project onto future scenario
## 2030 RCP 8.5
## prepare ff table for RCP 8.5 CESM1BGC 2030 (GCM2)


CESM1BGC_30.8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/CESM1BGC/",full.names = TRUE)
f1              <- stack(CESM1BGC_30.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1BGC_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1BGC_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})




## prepare ff table for RCP 8.5 CESM1BGC 2050 (GCM2)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/CESM1BGC/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1BGC_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1BGC_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


########
## prepare ff table for RCP 8.5 CESM1BGC 2070 (GCM2)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/CESM1BGC/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1BGC_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1BGC_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


# GCM 3_CESM1CAM5
## prepare ff table for RCP 8.5 CESM1CAM5 2030 (GCM3)


CESM1CAM5_30.8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/CESM1CAM5/",full.names = TRUE)
f1              <- stack(CESM1CAM5_30.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1CAM5_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1CAM5_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})



## prepare ff table for RCP 8.5 CESM1CAM5 2050 (GCM3)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/CESM1CAM5/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1CAM5_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1CAM5_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


########
## prepare ff table for RCP 8.5 CESM1CAM5 2070 (GCM3)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/CESM1CAM5/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1CAM5_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CESM1CAM5_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

# GCM 4_CMCCCM
## prepare ff table for RCP 8.5 CMCCCM 2030 (GCM4)


CMCCCM_30.8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/CMCCCM/",full.names = TRUE)
f1              <- stack(CMCCCM_30.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CMCCCM_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CMCCCM_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


# GCM 4
## prepare ff table for RCP 8.5 CMCCCM 2050 (GCM4)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/CMCCCM/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CMCCCM_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CMCCCM_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


########
## prepare ff table for RCP 8.5 CMCCCM 2070 (GCM4)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/CMCCCM/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CMCCCM_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_CMCCCM_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

# GCM 5_FIOESM
## prepare ff table for RCP 8.5 FIOESM 2030 (GCM5)


FIOESM_30.8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/FIOESM/",full.names = TRUE)
f1              <- stack(FIOESM_30.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_FIOESM_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_FIOESM_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})








## prepare ff table for RCP 8.5 FIOESM 2050 (GCM5)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/FIOESM/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_FIOESM_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_FIOESM_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


########
## prepare ff table for RCP 8.5 FIOESM 2070 (GCM5)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/FIOESM/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_FIOESM_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_FIOESM_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


# GCM 6_GISSE2H
## prepare ff table for RCP 8.5 GISSE2H 2030 (GCM6)


GISSE2H_30.8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/GISSE2H/",full.names = TRUE)
f1              <- stack(GISSE2H_30.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_GISSE2H_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_GISSE2H_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})





## prepare ff table for RCP 8.5 GISSE2H 2050 (GCM6)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/GISSE2H/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_GISSE2H_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_GISSE2H_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


########
## prepare ff table for RCP 8.5 GISSE2H 2070 (GCM6)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/GISSE2H/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_GISSE2H_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_GISSE2H_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

# GCM 7_Inmcm4
## prepare ff table for RCP 8.5 Inmcm4 2030 (GCM7)


Inmcm4_30.8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/Inmcm4/",full.names = TRUE)
f1              <- stack(Inmcm4_30.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_Inmcm4_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_Inmcm4_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})



## prepare ff table for RCP 8.5 Inmcm4 2050 (GCM7)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/Inmcm4/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_Inmcm4_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_Inmcm4_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


########
## prepare ff table for RCP 8.5 Inmcm4 2070 (GCM7)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/Inmcm4/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_Inmcm4_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_Inmcm4_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})



# GCM 8_IPSLCM5AMR
## prepare ff table for RCP 8.5 IPSLCM5AMR 2030 (GCM8)


IPSLCM5AMR_30.8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/IPSLCM5AMR/",full.names = TRUE)
f1              <- stack(IPSLCM5AMR_30.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_IPSLCM5AMR_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_IPSLCM5AMR_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})





## prepare ff table for RCP 8.5 IPSLCM5AMR 2050 (GCM8)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/IPSLCM5AMR/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_IPSLCM5AMR_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_IPSLCM5AMR_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


########
## prepare ff table for RCP 8.5 IPSLCM5AMR 2070 (GCM8)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/IPSLCM5AMR/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_IPSLCM5AMR_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_IPSLCM5AMR_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


# GCM 9_MIROC5
## prepare ff table for RCP 8.5 MIROC5 2030 (GCM9)


MIROC5_30.8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/MIROC5/",full.names = TRUE)
f1              <- stack(MIROC5_30.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MIROC5_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MIROC5_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})





## prepare ff table for RCP 8.5 MIROC5 2050 (GCM9)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/MIROC5/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MIROC5_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MIROC5_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


########
## prepare ff table for RCP 8.5 MIROC5 2070 (GCM9)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/MIROC5/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MIROC5_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MIROC5_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

# GCM 10_MPIESMMR
## prepare ff table for RCP 8.5 MPIESMMR 2030 (GCM10)


MPIESMMR_30.8.5  <- list.files("CHELSA_RASTERS_1km/2030_correct/RCP85/MPIESMMR/",full.names = TRUE)
f1              <- stack(MPIESMMR_30.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MPIESMMR_30_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MPIESMMR_30_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})



## prepare ff table for RCP 8.5 MPIESMMR 2050 (GCM10)

access_8.5  <- list.files("CHELSA_RASTERS_1km/2050/RCP85/MPIESMMR/",full.names = TRUE)
f1              <- stack(access_8.5)
f1 <- dropLayer(f1, c(7,8,18,19))

names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
# fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MPIESMMR_50_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})

lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MPIESMMR_50_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


########
## prepare ff table for RCP 8.5 MPIESMMR 2070 (GCM10)

access_70.8.5  <- list.files("CHELSA_RASTERS_1km/2070/RCP85/MPIESMMR/",full.names = TRUE)
f1              <- stack(access_70.8.5)
f1 <- dropLayer(f1, c(7,8,18,19))
names(f1) = names(s)
locs1 <- which(!is.na(sum(f1)[]))
locs1 <- cbind(cell=locs1, xyFromCell(f1, locs1))
cells1 <- locs1[, 'cell']


library(ff)
f1_ff <- ff::ff(vmode="double", dim=c(length(cells1), nlayers(f1)),
                filename=ff_swd <- tempfile(fileext='.ff'))
#fill ff_matrix with data
for(i in 1:nlayers(f1)) {
  f1_ff[, i] <- f1[[i]][][cells1]
}
colnames(f1_ff) <- names(f1)
gc()


### apply funciton to write the prediciton map and binary map according to selected threshold
lapply(ff_models_FN, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MPIESMMR_70_8.5.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('full/species\\.lambdas', 'full/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})


lapply(ff_models_random, function(m) {
  species <- basename(dirname(dirname(m)))
  outfile <- sprintf('%s/%s_MPIESMMR_70_8.5_random.ff', 
                     (dirname(dirname(m))),
                     gsub(' ', '_', species))
  if(!file.exists(extension(outfile, '.tif'))){
    r_pred <- f1[[1]]
    message('Doing species: ', gsub(' ', '_', species))
    
    preds_ff <- ff(vmode="double", dim=c(length(cells1), 1),
                   filename=outfile)
    finalizer(preds_ff) <- 'close'
    
    preds_ff[, 1] <- 
      round(rmaxent::project(m, f1_ff[, seq_len(ncol(f1_ff))], 
                             quiet=TRUE)$prediction_logistic * 1000)
    r_pred[cells1] <- preds_ff[, 1]
    writeRaster(r_pred, extension(outfile, '.tif'), datatype='INT2S', NAflag=-9999, overwrite=TRUE)
    
    results <- suppressMessages(
      read_csv(sub('rndm_bkgr/species\\.lambdas', 'rndm_bkgr/maxentResults.csv', m))
    )
    
    thr <- results$`10 percentile training presence logistic threshold`[
      nrow(results)]
    r_pred_binary <- r_pred >= thr*1000
    outfile_binary <- sub('\\.ff$', '_binary', outfile)
    writeRaster(r_pred_binary, extension(outfile_binary, '.tif'), datatype='INT1U', overwrite=TRUE)
    close(preds_ff)
    delete(preds_ff)
    rm(preds_ff, r_pred, r_pred_binary)
    
  }
})






####################-------------------------------------------------####################
####################--->  Count number of point per koppen zone  <---####################
####################-------------------------------------------------####################

spp_occ1 <-  list.files('__OUTDIR__/plants_clim1/', 'occ\\.shp', recursive=TRUE, full=TRUE)

spp_occ <- lapply(spp_occ1, function(f) {
  readOGR(f)
})

spp_koppen <- list.files('__OUTDIR__/plants_clim1/', 'koppen\\.shp', recursive=TRUE, full=TRUE)

spp_koppen <- lapply(spp_koppen, function(f) {
  readOGR(f)
})

sp <- gsub('__OUTDIR__/plants_clim1//', '', spp_occ1)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)


library(ggplot2)
for (i in 1:length(spp_occ)) {
  spp_occ[[i]] <- spTransform(spp_occ[[i]], crs(spp_koppen[[i]]))
  over=over(spp_occ[[i]], spp_koppen[[i]])
  table=as.data.frame(table(over$newKop))
  table$percent <- round(table$Freq/sum(table$Freq), 3)
  colnames(table) <- c("koppen_class", "frequency", "percent")
  write.csv(table, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp.nm[[i]], '_point_koppen.csv')), row.names=F)
  
  table1 <- table %>%
    arrange(desc(koppen_class)) %>%
    mutate(lab.ypos = cumsum(percent) - 0.5*percent)
  
  pie <- ggplot(table1, aes(x = "", y = percent, fill = koppen_class)) +
    geom_bar(width = 1, stat = "identity", color = "white") +
    coord_polar("y", start = 0)+
    geom_text(aes(y = lab.ypos, label = percent), color = "black")+
    theme_void()
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp.nm[[i]], '_point_koppen.png'))
  png(f, 5, 5, units='in', res=300)
  print(pie)
  dev.off()
  
}



#PLOT THE MAPS

library(magrittr)
library(raster)
library(sp)
library(rasterVis)
library(rgdal)
library(dplyr)
library(rnaturalearth)
library(RColorBrewer)
library(gdalUtils)
library(gridExtra)
library(ggplot2)


#plot for OCC
spp_occ= list.files('__OUTDIR__/plants_clim1/', 'occ\\.shp', recursive=TRUE, full=TRUE)

sp <- gsub('__OUTDIR__/plants_clim1//', '', spp_occ)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)

#spp_occ<- lapply(spp_occ, function(f) {
#  readOGR(f)
#})


data("wrld_simpl", package = 'maptools')
shape <- wrld_simpl
shape_df <- fortify(shape, region = 'ISO3')


for (i in 1:length(spp_occ)) {
  sp_oc <- readOGR(spp_occ[[i]])
  sp_oc <- spTransform(sp_oc, CRS(proj4string(shape)))
  p <- ggplot() +
    geom_polygon(data = shape_df , aes(x = long, y = lat, group = group), fill = 'white', color = 'black')+
    geom_point(data= data.frame(sp_oc@coords), aes(x=coords.x1, y=coords.x2), color="#0A8064", size = 0.5)+
    ggtitle(paste0(sp.nm[[i]], ' occurrence points'))
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_occ.png'))
  png(f, 7, 4.5, units='in', res=200)
  print(p)
  dev.off()
  
}


##### thereshold 0.5

spp_occ= list.files('__OUTDIR__/plants_clim1/', 'occ\\.shp', recursive=TRUE, full=TRUE)

sp <- gsub('__OUTDIR__/plants_clim1//', '', spp_occ)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)

ras.list <- list.files('__OUTDIR__/plants_clim1', 'current_random.tif$', full=TRUE, recursive=TRUE)

for (i in 1:length(spp_occ)) {
  sp_oc <- readOGR(spp_occ[[i]])
  ras <- lapply(ras.list[i], raster)
  max.val <- raster::extract(ras[[1]], sp_oc)
  th.50 <- min(max.val)
  bin.50 <- reclassify(ras[[1]], c(-Inf,th.50,0, th.50,Inf,1))
  writeRaster(bin.50, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_current_random_binary_100%.tif')), overwrite=T)
}

aus <- readOGR('additional_data_for_modeling/Australia/Aust_states.shp')
aus.moll <- spTransform(aus, s.crs)

ras.list1 <- list.files('__OUTDIR__/plants_clim1', '_current_random_binary_100%.tif$', full=TRUE, recursive=TRUE)

for (i in 1:length(spp_occ)) {
  sp_oc <- readOGR(spp_occ[[i]])
  ras <- lapply(ras.list1[i], raster)
  crp.aus <- crop(ras[[1]], aus.moll)
  crp.aus <- mask(crp.aus, aus.moll)
  writeRaster(crp.aus, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_current_random_binary_100%_aus.tif')), overwrite=T)
}


####### run MESS

spp_occ= list.files('__OUTDIR__/plants_clim1/', 'occ\\.shp', recursive=TRUE, full=TRUE)
maxent_fitted= list.files('__OUTDIR__/plants_clim1/', 'maxent_fitted\\.rds', recursive=TRUE, full=TRUE)

aus.shp <- shapefile('additional_data_for_modeling/Australia/Aust_states.shp')
aus.moll <- spTransform(aus.shp, s.crs)

sp <- gsub('__OUTDIR__/plants_clim1//', '', spp_occ)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)

for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS= rmaxent::similarity(s[[(clim1)]], mx_fit$swd)
  writeRaster(MESS$similarity_min, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS.tiff')), overwrite=TRUE)
  crop <- crop(MESS$similarity_min, aus.moll)
  crop <- mask(crop, aus.moll)
  
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS$similarity_min, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS)'))
  dev.off()
  
  f1 <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_Aus.png'))
  png(f1, 9, 6, units='in', res=300)
  plot(crop, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS)'))
  dev.off()
  
  lim <- rmaxent::limiting(s[[(clim1)]], mx_fit$me_random)
  
  f2 <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_limiting_variables.png'))
  png(f2, 9, 6, units='in', res=300)
  plot(lim, legend = FALSE, main=paste0(sp.nm[[i]], ' - limiting variables'))
  legend("bottomleft", legend = c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15"), fill=rev(terrain.colors(6)))
  dev.off()
  
  crop1 <- crop(lim, aus.moll)
  crop1 <- mask(crop1, aus.moll)
  
  f3 <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_limiting_variables_aus.png'))
  png(f3, 9, 6, units='in', res=300)
  plot(crop1, legend = FALSE, main=paste0(sp.nm[[i]], ' - limiting variables'))
  legend("bottomleft", legend = c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15"), fill=rev(terrain.colors(6)))
  dev.off()
}

mes <- stack(list.files('__OUTDIR__/plants_clim1', 'MESS.tif$', full=TRUE, recursive=TRUE))


for (i in 1:nlayers(mes)) {
  crop <- crop(mes[[i]], aus.moll)
  crop <- mask(crop, aus.moll)
  writeRaster(crop, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_Aus.tiff')), overwrite=TRUE)
  
}


###################----------------------------------###################
###################--->   Generating Ensemble    <---###################
###################----------------------------------###################
spp_occ=  list.files('__OUTDIR__/plants_clim1/', 'occ\\.shp', recursive=TRUE, full=TRUE)
sp <- gsub('__OUTDIR__/plants_clim1//', '', spp_occ)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)
path <- '__OUTDIR__/plants_clim1/'


for (i in 1:length(spp_occ)) {
  GCM1_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_ACCESS_30_8.5_random.tif'))))
  GCM2_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_CESM1BGC_30_8.5_random.tif'))))
  GCM3_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_CESM1CAM5_30_8.5_random.tif'))))
  GCM4_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_CMCCCM_30_8.5_random.tif'))))
  GCM5_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FIOESM_30_8.5_random.tif'))))
  GCM6_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_GISSE2H_30_8.5_random.tif'))))
  GCM7_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_Inmcm4_30_8.5_random.tif'))))
  GCM8_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_IPSLCM5AMR_30_8.5_random.tif'))))
  GCM9_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_MIROC5_30_8.5_random.tif'))))
  GCM10_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_MPIESMMR_30_8.5_random.tif'))))
  
  Ens_10GCMs_30 <- (GCM1_30 + GCM2_30 + GCM3_30 + GCM4_30 + GCM5_30 + GCM6_30 + GCM7_30 + GCM8_30 + GCM9_30 +  GCM10_30)/10
  writeRaster(Ens_10GCMs_30, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_Ensem_30.tiff')), overwrite=TRUE)
  
}


for (i in 1:length(spp_occ)) {
  GCM1_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_ACCESS_50_8.5_random.tif'))))
  GCM2_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_CESM1BGC_50_8.5_random.tif'))))
  GCM3_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_CESM1CAM5_50_8.5_random.tif'))))
  GCM4_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_CMCCCM_50_8.5_random.tif'))))
  GCM5_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FIOESM_50_8.5_random.tif'))))
  GCM6_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_GISSE2H_50_8.5_random.tif'))))
  GCM7_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_Inmcm4_50_8.5_random.tif'))))
  GCM8_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_IPSLCM5AMR_50_8.5_random.tif'))))
  GCM9_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_MIROC5_50_8.5_random.tif'))))
  GCM10_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_MPIESMMR_50_8.5_random.tif'))))
  
  Ens_10GCMs_50 <- (GCM1_50 + GCM2_50 + GCM3_50 + GCM4_50 + GCM5_50 + GCM6_50 + GCM7_50 + GCM8_50 + GCM9_50 +  GCM10_50)/10
  writeRaster(Ens_10GCMs_50, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_Ensem_50.tiff')), overwrite=TRUE)
  
}


for (i in 1:length(spp_occ)) {
  GCM1_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_ACCESS_70_8.5_random.tif'))))
  GCM2_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_CESM1BGC_70_8.5_random.tif'))))
  GCM3_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_CESM1CAM5_70_8.5_random.tif'))))
  GCM4_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_CMCCCM_70_8.5_random.tif'))))
  GCM5_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FIOESM_70_8.5_random.tif'))))
  GCM6_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_GISSE2H_70_8.5_random.tif'))))
  GCM7_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_Inmcm4_70_8.5_random.tif'))))
  GCM8_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_IPSLCM5AMR_70_8.5_random.tif'))))
  GCM9_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_MIROC5_70_8.5_random.tif'))))
  GCM10_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_MPIESMMR_70_8.5_random.tif'))))
  
  Ens_10GCMs_70 <- (GCM1_70 + GCM2_70 + GCM3_70 + GCM4_70 + GCM5_70 + GCM6_70 + GCM7_70 + GCM8_70 + GCM9_70 +  GCM10_70)/10
  writeRaster(Ens_10GCMs_70, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_Ensem_70.tiff')), overwrite=TRUE)
  
}



###################----------------------------------###################
###################--->   plot all maps in one   <---###################
###################----------------------------------###################

spp_occ=  list.files('__OUTDIR__/plants_clim1/', 'occ\\.shp', recursive=TRUE, full=TRUE)
sp <- gsub('__OUTDIR__/plants_clim1//', '', spp_occ)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)

data("wrld_simpl", package = 'maptools')
shape <- wrld_simpl
shape_df <- fortify(shape, region = 'ISO3')

for (i in 1:length(spp_occ)) {
  ras <- raster(paste0('__OUTDIR__/plants_clim1/', sp[i], paste0("/", sp[i], '_current_random.tif')))
  ras.bin <- raster(paste0('__OUTDIR__/plants_clim1/', sp[i], paste0("/", sp[i], '_current_random_binary.tif')))
  mes <- raster(paste0('__OUTDIR__/plants_clim1/', sp[i], paste0("/", sp[i], '_MESS.tif')))
  rec.mes <- reclassify(mes, c(-Inf,0,0, 0,Inf,1))
  writeRaster(rec.mes, paste0('__OUTDIR__/plants_clim1/', sp[i], paste0("/", sp[i], '_current_random_binary_MESS.tif')),overwrite=TRUE)                                
  rr <- ras * rec.mes
  rr2 <- ras.bin * rec.mes
  writeRaster(rr, paste0('__OUTDIR__/plants_clim1/', sp[i], paste0("/", sp[i], '_current_random_Extracted_MESS.tif')),overwrite=TRUE)
  writeRaster(rr2, paste0('__OUTDIR__/plants_clim1/', sp[i], paste0("/", sp[i], '_current_random_binary_Extracted_MESS.tif')),overwrite=TRUE)
  crop3 <- crop(rec.mes, aus.moll)
  crop3 <- mask(crop3, aus.moll)
  writeRaster(crop3, file.path(paste0('__OUTDIR__/plants_clim1//', sp[i]), paste0(sp[[i]], '_MESS_binary_Aus.tiff')), overwrite=TRUE)
  crop4 <- crop(rr2, aus.moll)
  crop4 <- mask(crop4, aus.moll)
  writeRaster(crop4, file.path(paste0('__OUTDIR__/plants_clim1//', sp[i]), paste0(sp[[i]], '_MESS_binary_extract_Aus.tiff')), overwrite=TRUE)
  
}      

#cur <- stack(list.files('__OUTDIR__/plants_clim1', 'current.tif$', full=TRUE, recursive=TRUE))
#bin <- stack(list.files('__OUTDIR__/plants_clim1', 'current_binary.tif$', full=TRUE, recursive=TRUE))
#aus <- stack(list.files('__OUTDIR__/plants_clim1', 'current_aus_binary.tif$', full=TRUE, recursive=TRUE))

cur_random <- stack(list.files('__OUTDIR__/plants_clim1', 'current_random.tif$', full=TRUE, recursive=TRUE))
bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'current_random_binary.tif$', full=TRUE, recursive=TRUE))
aus_random <- stack(list.files('__OUTDIR__/plants_clim1', 'current_aus_random_binary.tif$', full=TRUE, recursive=TRUE))


mes <- stack(list.files('__OUTDIR__/plants_clim1/', 'MESS_Aus.tif$', full=TRUE, recursive=TRUE))
mes_bin <- stack(list.files('__OUTDIR__/plants_clim1/', '_current_random_binary_MESS.tif$', full=TRUE, recursive=TRUE))
mes_ext <- stack(list.files('__OUTDIR__/plants_clim1/', '_current_random_Extracted_MESS.tif$', full=TRUE, recursive=TRUE))
mes_ext.bin <- stack(list.files('__OUTDIR__/plants_clim1/', '_current_random_binary_Extracted_MESS.tif$', full=TRUE, recursive=TRUE))
mes_bin_aus <- stack(list.files('__OUTDIR__/plants_clim1/', '_MESS_binary_Aus.tif$', full=TRUE, recursive=TRUE))
mes_ext.bin_aus <- stack(list.files('__OUTDIR__/plants_clim1/', '_MESS_binary_extract_Aus.tif$', full=TRUE, recursive=TRUE))


bin_random50 <- stack(list.files('__OUTDIR__/plants_clim1', 'current_random_binary_100%.tif$', full=TRUE, recursive=TRUE))
bin_random50_aus <- stack(list.files('__OUTDIR__/plants_clim1', 'current_random_binary_100%_aus.tif$', full=TRUE, recursive=TRUE))


ac3085 <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_30_8.5.tif$', full=TRUE, recursive=TRUE))
ac3085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_30_8.5_binary.tif$', full=TRUE, recursive=TRUE))

ac3085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_30_8.5_random.tif$', full=TRUE, recursive=TRUE))
ac3085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_30_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


ac5085 <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_50_8.5.tif$', full=TRUE, recursive=TRUE))
ac5085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

ac5085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
ac5085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


ac7085 <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_70_8.5.tif$', full=TRUE, recursive=TRUE))
ac7085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

ac7085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
ac7085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'ACCESS_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))

#GCM2#CESM1BGC
CESM1BGC5085 <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1BGC_50_8.5.tif$', full=TRUE, recursive=TRUE))
CESM1BGC5085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1BGC_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

CESM1BGC5085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1BGC_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
CESM1BGC5085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1BGC_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


CESM1BGC7085 <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1BGC_70_8.5.tif$', full=TRUE, recursive=TRUE))
CESM1BGC7085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1BGC_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

CESM1BGC7085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1BGC_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
CESM1BGC7085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1BGC_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


#GCM3#CESM1CAM
CESM1CAM55085 <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1CAM5_50_8.5.tif$', full=TRUE, recursive=TRUE))
CESM1CAM55085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1CAM5_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

CESM1CAM55085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1CAM5_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
CESM1CAM55085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1CAM5_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


CESM1CAM57085 <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1CAM5_70_8.5.tif$', full=TRUE, recursive=TRUE))
CESM1CAM57085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1CAM5_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

CESM1CAM57085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1CAM5_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
CESM1CAM57085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CESM1CAM5_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))

#GCM4#CESM1CAM
CMCCCM5085 <- stack(list.files('__OUTDIR__/plants_clim1', 'CMCCCM_50_8.5.tif$', full=TRUE, recursive=TRUE))
CMCCCM5085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'CMCCCM_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

CMCCCM5085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CMCCCM_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
CMCCCM5085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CMCCCM_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


CMCCCM7085 <- stack(list.files('__OUTDIR__/plants_clim1', 'CMCCCM_70_8.5.tif$', full=TRUE, recursive=TRUE))
CMCCCM7085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'CMCCCM_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

CMCCCM7085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CMCCCM_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
CMCCCM7085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'CMCCCM_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))

#GCM5#FIOESM
FIOESM5085 <- stack(list.files('__OUTDIR__/plants_clim1', 'FIOESM_50_8.5.tif$', full=TRUE, recursive=TRUE))
FIOESM5085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'FIOESM_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

FIOESM5085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'FIOESM_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
FIOESM5085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'FIOESM_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


FIOESM7085 <- stack(list.files('__OUTDIR__/plants_clim1', 'FIOESM_70_8.5.tif$', full=TRUE, recursive=TRUE))
FIOESM7085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'FIOESM_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

FIOESM7085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'FIOESM_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
FIOESM7085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'FIOESM_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


#GCM6#GISSE2H
GISSE2H5085 <- stack(list.files('__OUTDIR__/plants_clim1', 'GISSE2H_50_8.5.tif$', full=TRUE, recursive=TRUE))
GISSE2H5085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'GISSE2H_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

GISSE2H5085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'GISSE2H_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
GISSE2H5085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'GISSE2H_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


GISSE2H7085 <- stack(list.files('__OUTDIR__/plants_clim1', 'GISSE2H_70_8.5.tif$', full=TRUE, recursive=TRUE))
GISSE2H7085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'GISSE2H_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

GISSE2H7085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'GISSE2H_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
GISSE2H7085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'GISSE2H_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))

#GCM7#Inmcm4
Inmcm45085 <- stack(list.files('__OUTDIR__/plants_clim1', 'Inmcm4_50_8.5.tif$', full=TRUE, recursive=TRUE))
Inmcm45085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'Inmcm4_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

Inmcm45085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'Inmcm4_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
Inmcm45085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'Inmcm4_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


Inmcm47085 <- stack(list.files('__OUTDIR__/plants_clim1', 'Inmcm4_70_8.5.tif$', full=TRUE, recursive=TRUE))
Inmcm47085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'Inmcm4_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

Inmcm47085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'Inmcm4_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
Inmcm47085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'Inmcm4_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))

#GCM8#IPSLCM5AMR
IPSLCM5AMR5085 <- stack(list.files('__OUTDIR__/plants_clim1', 'IPSLCM5AMR_50_8.5.tif$', full=TRUE, recursive=TRUE))
IPSLCM5AMR5085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'IPSLCM5AMR_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

IPSLCM5AMR5085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'IPSLCM5AMR_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
IPSLCM5AMR5085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'IPSLCM5AMR_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


IPSLCM5AMR7085 <- stack(list.files('__OUTDIR__/plants_clim1', 'IPSLCM5AMR_70_8.5.tif$', full=TRUE, recursive=TRUE))
IPSLCM5AMR7085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'IPSLCM5AMR_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

IPSLCM5AMR7085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'IPSLCM5AMR_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
IPSLCM5AMR7085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'IPSLCM5AMR_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))

#GCM9#MIROC5
MIROC55085 <- stack(list.files('__OUTDIR__/plants_clim1', 'MIROC5_50_8.5.tif$', full=TRUE, recursive=TRUE))
MIROC55085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'MIROC5_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

MIROC55085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'MIROC5_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
MIROC55085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'MIROC5_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


MIROC57085 <- stack(list.files('__OUTDIR__/plants_clim1', 'MIROC5_70_8.5.tif$', full=TRUE, recursive=TRUE))
MIROC57085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'MIROC5_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

MIROC57085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'MIROC5_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
MIROC57085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'MIROC5_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))

#GCM10#MPIESMMR
MPIESMMR5085 <- stack(list.files('__OUTDIR__/plants_clim1', 'MPIESMMR_50_8.5.tif$', full=TRUE, recursive=TRUE))
MPIESMMR5085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'MPIESMMR_50_8.5_binary.tif$', full=TRUE, recursive=TRUE))

MPIESMMR5085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'MPIESMMR_50_8.5_random.tif$', full=TRUE, recursive=TRUE))
MPIESMMR5085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'MPIESMMR_50_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


MPIESMMR7085 <- stack(list.files('__OUTDIR__/plants_clim1', 'MPIESMMR_70_8.5.tif$', full=TRUE, recursive=TRUE))
MPIESMMR7085.bin <- stack(list.files('__OUTDIR__/plants_clim1', 'MPIESMMR_70_8.5_binary.tif$', full=TRUE, recursive=TRUE))

MPIESMMR7085_random <- stack(list.files('__OUTDIR__/plants_clim1', 'MPIESMMR_70_8.5_random.tif$', full=TRUE, recursive=TRUE))
MPIESMMR7085.bin_random <- stack(list.files('__OUTDIR__/plants_clim1', 'MPIESMMR_70_8.5_random_binary.tif$', full=TRUE, recursive=TRUE))


diverge0 <- function(p, ramp) {
  # p: a trellis object resulting from rasterVis::levelplot
  # ramp: the name of an RColorBrewer palette (as character), a character 
  #       vector of colour names to interpolate, or a colorRampPalette.
  require(RColorBrewer)
  require(rasterVis)
  if(length(ramp)==1 && is.character(ramp) && ramp %in% 
     row.names(brewer.pal.info)) {
    ramp <- suppressWarnings(colorRampPalette(brewer.pal(11, ramp)))
  } else if(length(ramp) > 1 && is.character(ramp) && all(ramp %in% colors())) {
    ramp <- colorRampPalette(ramp)
  } else if(!is.function(ramp)) 
    stop('ramp should be either the name of a RColorBrewer palette, ', 
         'a vector of colours to be interpolated, or a colorRampPalette.')
  rng <- range(p$legend[[1]]$args$key$at)
  s <- seq(-max(abs(rng)), max(abs(rng)), len=1001)
  i <- findInterval(rng[which.min(abs(rng))], s)
  zlim <- switch(which.min(abs(rng)), `1`=i:(1000+1), `2`=1:(i+1))
  p$legend[[1]]$args$key$at <- s[zlim]
  p[[grep('^legend', names(p))]][[1]]$args$key$col <- ramp(1000)[zlim[-length(zlim)]]
  p$panel.args.common$col.regions <- ramp(1000)[zlim[-length(zlim)]]
  p
}


for (i in 1:length(spp_occ)) {
  sp_oc <- readOGR(spp_occ[[i]])
  
  p5 <- levelplot(cur_random[[i]], lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' current suitability random'), font=4))
  r2 <- bin_random[[i]]
  r2 <- as.factor(r2)
  ID <- c(0, 1)
  Prs_Abs2 <- c("Absence","Presence")
  rat2 <- data.frame(ID, Prs_Abs2)
  levels(r2) <- rat2
  p6 <- levelplot(r2, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' binary suitability random'), font=4))
  r3 <- aus_random[[i]]
  r3 <- as.factor(r3)
  ID <- c(0, 1)
  Prs_Abs3 <- c("Absence","Presence")
  rat3 <- data.frame(ID, Prs_Abs3)
  levels(r3) <- rat3
  p7 <- levelplot(r3, lwd=0,col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' Australia suitability random'), font=4))
  r4 <- bin_random50[[i]]
  r4 <- as.factor(r4)
  ID <- c(0, 1)
  Prs_Abs4 <- c("Absence","Presence")
  rat4 <- data.frame(ID, Prs_Abs4)
  levels(r4) <- rat4
  p8 <- levelplot(r4, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' 100% binary suitability random'), font=4))
  
  r5 <- bin_random50_aus[[i]]
  r5 <- as.factor(r5)
  ID <- c(0, 1)
  Prs_Abs5 <- c("Absence","Presence")
  rat5 <- data.frame(ID, Prs_Abs5)
  levels(r5) <- rat5
  p9 <- levelplot(r5, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' Australia 100% suitability random'), font=4))
  
  messs <- setValues(raster(mes[[i]]), mes[[i]][])
  p10 <- levelplot(messs, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                   margin=F, main=list(paste0(sp.nm[[i]], ' MESS'), font=4))
  
  r6 <- mes_bin_aus[[i]]
  r6 <- as.factor(r6)
  ID <- c(0, 1)
  Prs_Abs6 <- c("mess negative","mess possitive")
  rat6 <- data.frame(ID, Prs_Abs6)
  levels(r6) <- rat6
  p11 <- levelplot(r6, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' MESS binary Aus'), font=3))
  
  p12 <- levelplot(mes_ext[[i]], lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' current suitability with extracted MESS'), font=3))
  
  r7 <- mes_ext.bin[[i]]
  r7 <- as.factor(r7)
  ID <- c(0, 1)
  Prs_Abs7 <- c("Absence","Presence")
  rat7 <- data.frame(ID, Prs_Abs7)
  levels(r7) <- rat7
  p13 <- levelplot(r7, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                   margin=F, main=list(paste0(sp.nm[[i]], ' current binary with extracted MESS'), font=3))
  
  r8 <- mes_ext.bin_aus[[i]]
  r8 <- as.factor(r8)
  ID <- c(0, 1)
  Prs_Abs8 <- c("Absence","Presence")
  rat8 <- data.frame(ID, Prs_Abs8)
  levels(r8) <- rat8
  p14 <- levelplot(r8, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                   margin=F, main=list(paste0(sp.nm[[i]], ' current binary Australia with extracted MESS'), font=3))
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1//', sp[i]), paste0(sp[[i]], '_plots.png'))
  png(f, 21, 15, units='in', res=650)
  gridExtra::grid.arrange(diverge0(p10, 'RdBu'),diverge0(p10, 'RdBu'), p11, p5, p5, p12, p6, p8, p13, p7, p9, p14, ncol=3)
  dev.off()
}


for (i in 1:nlayers(ac5085)) {
  p1 <- levelplot(abs(ac5085[[i]]), lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' future suitability (2050 - ACCESS)'), font=4))
  r <-  ac5085.bin[[i]]
  r <- as.factor(r)
  ID <- c(0, 1)
  Prs_Abs <- c("Absence","Presence")
  rat <- data.frame(ID, Prs_Abs)
  levels(r) <- rat
  p2 <- levelplot(r, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' binary future suitability (2050 - ACCESS)'), font=4))
  p3 <- levelplot(abs(ac7085[[i]]), lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' future suitability (2070 - ACCESS)'), font=4))
  r1 <-  ac7085.bin[[i]]
  r1 <- as.factor(r1)
  ID <- c(0, 1)
  Prs_Abs1 <- c("Absence","Presence")
  rat1 <- data.frame(ID, Prs_Abs1)
  levels(r1) <- rat1
  p4 <- levelplot(r1, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' binary future suitability (2070 - ACCESS)'), font=4))
  
  p5 <- levelplot(abs(ac5085_random[[i]]), lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' future suitability (2050 - ACCESS) random'), font=4))
  r2 <-  ac5085.bin_random[[i]]
  r2 <- as.factor(r2)
  ID <- c(0, 1)
  Prs_Abs2 <- c("Absence","Presence")
  rat2 <- data.frame(ID, Prs_Abs2)
  levels(r2) <- rat2
  p6 <- levelplot(r2, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' binary future suitability (2050 - ACCESS) random'), font=4))
  p7 <- levelplot(abs(ac7085_random[[i]]), lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' future suitability (2070 - ACCESS) random'), font=4))
  r3 <-  ac7085.bin_random[[i]]
  r3 <- as.factor(r3)
  ID <- c(0, 1)
  Prs_Abs3 <- c("Absence","Presence")
  rat3 <- data.frame(ID, Prs_Abs3)
  levels(r3) <- rat3
  p8 <- levelplot(r3, lwd=0, col.regions=colorRampPalette(rev(terrain.colors(100))),
                  margin=F, main=list(paste0(sp.nm[[i]], ' binary future suitability (2070 - ACCESS) random'), font=4))
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_future projection plots.png'))
  png(f, 14, 15, units='in', res=650)
  gridExtra::grid.arrange(p1, p5, p2, p6, p3, p7, p4, p8, ncol=2)
  dev.off()
}




#paste current time plots
new.folder <- '__OUTDIR__/plants_clim1/1_Current_suitability_plots_occ_maps/'
dir.create(new.folder, recursive=TRUE)
for (i in 1:length(sp)) {
  plots <- list.files('__OUTDIR__/plants_clim1', '.png$', full=TRUE, recursive=TRUE) %>% 
    grep(paste0(sp[i], '_plot'), ., value=TRUE, invert=F)
  file.copy(plots, new.folder)
}

#paste future projection plots
plots1 <- list.files('__OUTDIR__/plants_clim1', '.png$', full=TRUE, recursive=TRUE) %>% 
  grep('future projection plots', ., value=TRUE, invert=F)
new.folder1 <- '__OUTDIR__/plants_clim1/2_future_suitability_plots/'
dir.create(new.folder1, recursive=TRUE)
file.copy(plots1, new.folder1)

#paste occ plots
for (i in 1:length(sp)) {
  plots3 <- list.files('__OUTDIR__/plants_clim1', '.png$', full=TRUE, recursive=TRUE) %>% 
    grep(paste0(sp[i], '_occ'), ., value=TRUE, invert=F)
  new.folder <- '__OUTDIR__/plants_clim1/1_Current_suitability_plots_occ_maps/'
  file.copy(plots3, new.folder)
}

#paste MESS_maps
new.folder2 <- '__OUTDIR__/plants_clim1/3_MESS_maps_Globe_and_AUS/'
dir.create(new.folder2, recursive=TRUE)
for (i in 1:length(sp)) {
  plots4 <- list.files('__OUTDIR__/plants_clim1', '.png$', full=TRUE, recursive=TRUE) %>% 
    grep(paste0(sp[i], '_MESS'), ., value=TRUE, invert=F)
  new.folder <- '__OUTDIR__/plants_clim1/3_MESS_maps_Globe_and_AUS/'
  file.copy(plots4, new.folder2)
}

#paste limiting_variables
new.folder3 <- '__OUTDIR__/plants_clim1/4_limiting_variables_Globe_and_AUS/'
dir.create(new.folder3, recursive=TRUE)
for (i in 1:length(sp)) {
  plots5 <- list.files('__OUTDIR__/plants_clim1', '.png$', full=TRUE, recursive=TRUE) %>% 
    grep(paste0(sp[i], '_limiting_variables'), ., value=TRUE, invert=F)
  new.folder <- '__OUTDIR__/plants_clim1/4_limiting_variables_Globe_and_AUS/'
  file.copy(plots5, new.folder3)
} 




###################----------------------------------###################
###################--->   Current mol to Alber   <---###################
###################----------------------------------###################

spp_occ=  list.files(path, 'occ\\.shp', recursive=TRUE, full=TRUE)
sp <- gsub('(.*/)?plants_clim1//', '', spp_occ)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)
crs.alber <- "+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs"

for (i in 1:length(sp)) {
  f_in <- file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_current_aus_random.tif')))
  f_out <- file.path(paste0(path, sp[i]), paste0(sp[[i]], '_current_aus_random_Alber.tif'))
  template_to_match <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_Ensem_50.tif')))) # crop to this extent
  
  # Either this...
  ras1 <- raster(f_in)
  ras.alber <- projectRaster(ras1, crs=crs.alber, res=1000, method="ngb")
  # the origins will not match exactly if crop is used, so we instead use resample
  #ras.alber <- crop(ras.alber, extent(template_to_crop))
  ras.alber <- raster::resample(ras.alber, template_to_match, 'ngb')
  writeRaster(ras.alber, f_out, overwrite=TRUE)
  
  # Or this... can replace projectRaster with gdalUtils::gdalwarp if you have GDAL system library installed
  # gdalUtils::gdalwarp(f_in, f_out, tr=c(1000, 1000), t_srs='EPSG:3577', r='near',
  #                     te=c(bbox(template_to_crop)), overwrite=TRUE)
  
}




###################----------------------------------###################
###################--->    MESS                  <---###################
###################----------------------------------###################

spp_occ= list.files('__OUTDIR__/plants_clim1/', 'occ\\.shp', recursive=TRUE, full=TRUE)
maxent_fitted= list.files('__OUTDIR__/plants_clim1/', 'maxent_fitted\\.rds', recursive=TRUE, full=TRUE)

aus.shp <- shapefile('additional_data_for_modeling/Australia/Aust_states.shp')
aus.moll <- spTransform(aus.shp, s.crs)

sp <- gsub('__OUTDIR__/plants_clim1//', '', spp_occ)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)


s_cur <- list.files('CHELSA_RASTERS_1km/Current_cropped_AUS_1km_alber/', '\\.tif$', full.names=TRUE)
s_current <- stack(s_cur)

names(s_current)    <- c("biol_01",
                         "biol_02",
                         "biol_03",
                         "biol_04",
                         "biol_05",
                         "biol_06",
                         "biol_07",
                         "biol_08",
                         "biol_09",
                         "biol_10",
                         "biol_11",
                         "biol_12",
                         "biol_13",
                         "biol_14",
                         "biol_15",
                         "biol_16",
                         "biol_17",
                         "biol_18",
                         "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(s_current[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_Aus.tiff')), overwrite=TRUE)
  suitability_current <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_current_aus_random_Alber.tif')))
  suitability_current[MESS < 0] <- NA
  suitability_current <- crop(suitability_current, extent(MESS))
  writeRaster(suitability_current, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_CURRENT.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_CURRENT.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) current'))
  dev.off()
}





GCM1_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/ACCESS10', '\\.tif$', full.names=TRUE)
GCM1_30  <- stack(GCM1_30)

names(GCM1_30)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM1_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM1_30.tiff')), overwrite=TRUE)
  suitability_current1 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_ACCESS_30_8.5_random.tif')))
  suitability_current1[MESS < 0] <- NA
  suitability_current1 <- crop(suitability_current1, extent(MESS))
  writeRaster(suitability_current1, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM1.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM1.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 ACCESS (GCM1)'))
  dev.off()
}

GCM2_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/CESM1BGC', '\\.tif$', full.names=TRUE)
GCM2_30  <- stack(GCM2_30)

names(GCM2_30)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM2_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM2_30.tiff')), overwrite=TRUE)
  suitability_current2 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_CESM1BGC_30_8.5_random.tif')))
  suitability_current2[MESS < 0] <- NA
  suitability_current2 <- crop(suitability_current2, extent(MESS))
  writeRaster(suitability_current2, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM2.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM2.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 CESM1BGC (GCM2)'))
  dev.off()
}


GCM3_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/CESM1CAM5', '\\.tif$', full.names=TRUE)
GCM3_30  <- stack(GCM3_30)

names(GCM3_30)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM3_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM3_30.tiff')), overwrite=TRUE)
  suitability_current3 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_CESM1CAM5_30_8.5_random.tif')))
  suitability_current3[MESS < 0] <- NA
  suitability_current3 <- crop(suitability_current3, extent(MESS))
  writeRaster(suitability_current3, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM3.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM3.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 CESM1CAM5'))
  dev.off()
}



GCM4_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/CMCCCM', '\\.tif$', full.names=TRUE)
GCM4_30  <- stack(GCM4_30)

names(GCM4_30)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM4_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM4_30.tiff')), overwrite=TRUE)
  suitability_current4 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_CMCCCM_30_8.5_random.tif')))
  suitability_current4[MESS < 0] <- NA
  suitability_current4 <- crop(suitability_current4, extent(MESS))
  writeRaster(suitability_current4, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM4.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM4.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 CMCCCM (GCM4)'))
  dev.off()
}

GCM5_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/FIOESM', '\\.tif$', full.names=TRUE)
GCM5_30  <- stack(GCM5_30)

names(GCM5_30)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM5_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM5_30.tiff')), overwrite=TRUE)
  suitability_current5 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FIOESM_30_8.5_random.tif')))
  suitability_current5[MESS < 0] <- NA
  suitability_current5 <- crop(suitability_current5, extent(MESS))
  writeRaster(suitability_current5, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM5.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM5.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 FIOESM (GCM5)'))
  dev.off()
}

GCM6_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/GISSE2H', '\\.tif$', full.names=TRUE)
GCM6_30  <- stack(GCM6_30)

names(GCM6_30)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM6_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM6_30.tiff')), overwrite=TRUE)
  suitability_current6 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_GISSE2H_30_8.5_random.tif')))
  suitability_current6[MESS < 0] <- NA
  suitability_current6 <- crop(suitability_current6, extent(MESS))
  writeRaster(suitability_current6, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM6.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM6.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 GISSE2H (GCM6)'))
  dev.off()
}


GCM7_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/Inmcm4', '\\.tif$', full.names=TRUE)
GCM7_30  <- stack(GCM7_30)

names(GCM7_30)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM7_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM7_30.tiff')), overwrite=TRUE)
  suitability_current7 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_Inmcm4_30_8.5_random.tif')))
  suitability_current7[MESS < 0] <- NA
  suitability_current7 <- crop(suitability_current7, extent(MESS))
  writeRaster(suitability_current7, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM7.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM7.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 Inmcm4 (GCM7)'))
  dev.off()
}

GCM8_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/IPSLCM5AMR', '\\.tif$', full.names=TRUE)
GCM8_30  <- stack(GCM8_30)

names(GCM8_30)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM8_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM8_30.tiff')), overwrite=TRUE)
  suitability_current8 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_IPSLCM5AMR_30_8.5_random.tif')))
  suitability_current8[MESS < 0] <- NA
  suitability_current8 <- crop(suitability_current8, extent(MESS))
  writeRaster(suitability_current8, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM8.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM8.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 IPSLCM5AMR (GCM8)'))
  dev.off()
}

GCM9_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/MIROC5', '\\.tif$', full.names=TRUE)
GCM9_30  <- stack(GCM9_30)

names(GCM9_30)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM9_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM9_30.tiff')), overwrite=TRUE)
  suitability_current9 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MIROC5_30_8.5_random.tif')))
  suitability_current9[MESS < 0] <- NA
  suitability_current9 <- crop(suitability_current9, extent(MESS))
  writeRaster(suitability_current9, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM9.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM9.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 MIROC5 (GCM9)'))
  dev.off()
}

GCM10_30 <- list.files('CHELSA_RASTERS_1km/2030_correct/RCP85/MPIESMMR', '\\.tif$', full.names=TRUE)
GCM10_30  <- stack(GCM10_30)

names(GCM10_30)    <- c("biol_01",
                        "biol_02",
                        "biol_03",
                        "biol_04",
                        "biol_05",
                        "biol_06",
                        "biol_07",
                        "biol_08",
                        "biol_09",
                        "biol_10",
                        "biol_11",
                        "biol_12",
                        "biol_13",
                        "biol_14",
                        "biol_15",
                        "biol_16",
                        "biol_17",
                        "biol_18",
                        "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM10_30[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM10_30.tiff')), overwrite=TRUE)
  suitability_current10 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MPIESMMR_30_8.5_random.tif')))
  suitability_current10[MESS < 0] <- NA
  suitability_current10 <- crop(suitability_current10, extent(MESS))
  writeRaster(suitability_current10, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2030_GCM10.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2030_GCM10.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2030 MPIESMMR (GCM10)'))
  dev.off()
}




GCM1_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/ACCESS10', '\\.tif$', full.names=TRUE)
GCM1_50  <- stack(GCM1_50)

names(GCM1_50)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM1_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM1_50.tiff')), overwrite=TRUE)
  suitability_current11 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_ACCESS_50_8.5_random.tif')))
  suitability_current11[MESS < 0] <- NA
  suitability_current11 <- crop(suitability_current11, extent(MESS))
  writeRaster(suitability_current11, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM1.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM1.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 ACCESS (GCM1)'))
  dev.off()
}

GCM2_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/CESM1BGC', '\\.tif$', full.names=TRUE)
GCM2_50  <- stack(GCM2_50)

names(GCM2_50)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM2_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM2_50.tiff')), overwrite=TRUE)
  suitability_current12 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_CESM1BGC_50_8.5_random.tif')))
  suitability_current12[MESS < 0] <- NA
  suitability_current12 <- crop(suitability_current12, extent(MESS))
  writeRaster(suitability_current12, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM2.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM2.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 CESM1BGC (GCM2)'))
  dev.off()
}


GCM3_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/CESM1CAM5', '\\.tif$', full.names=TRUE)
GCM3_50  <- stack(GCM3_50)

names(GCM3_50)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM3_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM3_50.tiff')), overwrite=TRUE)
  suitability_current13 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_CESM1CAM5_50_8.5_random.tif')))
  suitability_current13[MESS < 0] <- NA
  suitability_current13 <- crop(suitability_current13, extent(MESS))
  writeRaster(suitability_current13, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM3.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM3.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 CESM1CAM5'))
  dev.off()
}



GCM4_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/CMCCCM', '\\.tif$', full.names=TRUE)
GCM4_50  <- stack(GCM4_50)

names(GCM4_50)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM4_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM4_50.tiff')), overwrite=TRUE)
  suitability_current14 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_CMCCCM_50_8.5_random.tif')))
  suitability_current14[MESS < 0] <- NA
  suitability_current14 <- crop(suitability_current14, extent(MESS))
  writeRaster(suitability_current14, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM4.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM4.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 CMCCCM (GCM4)'))
  dev.off()
}

GCM5_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/FIOESM', '\\.tif$', full.names=TRUE)
GCM5_50  <- stack(GCM5_50)

names(GCM5_50)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM5_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM5_50.tiff')), overwrite=TRUE)
  suitability_current15 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FIOESM_50_8.5_random.tif')))
  suitability_current15[MESS < 0] <- NA
  suitability_current15 <- crop(suitability_current15, extent(MESS))
  writeRaster(suitability_current15, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM5.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM5.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 FIOESM (GCM5)'))
  dev.off()
}

GCM6_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/GISSE2H', '\\.tif$', full.names=TRUE)
GCM6_50  <- stack(GCM6_50)

names(GCM6_50)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM6_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM6_50.tiff')), overwrite=TRUE)
  suitability_current16 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_GISSE2H_50_8.5_random.tif')))
  suitability_current16[MESS < 0] <- NA
  suitability_current16 <- crop(suitability_current16, extent(MESS))
  writeRaster(suitability_current16, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM6.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM6.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 GISSE2H (GCM6)'))
  dev.off()
}


GCM7_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/Inmcm4', '\\.tif$', full.names=TRUE)
GCM7_50  <- stack(GCM7_50)

names(GCM7_50)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM7_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM7_50.tiff')), overwrite=TRUE)
  suitability_current17 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_Inmcm4_50_8.5_random.tif')))
  suitability_current17[MESS < 0] <- NA
  suitability_current17 <- crop(suitability_current17, extent(MESS))
  writeRaster(suitability_current17, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM7.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM7.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 Inmcm4 (GCM7)'))
  dev.off()
}

GCM8_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/IPSLCM5AMR', '\\.tif$', full.names=TRUE)
GCM8_50  <- stack(GCM8_50)

names(GCM8_50)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM8_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM8_50.tiff')), overwrite=TRUE)
  suitability_current18 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_IPSLCM5AMR_50_8.5_random.tif')))
  suitability_current18[MESS < 0] <- NA
  suitability_current18 <- crop(suitability_current18, extent(MESS))
  writeRaster(suitability_current18, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM8.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM8.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 IPSLCM5AMR (GCM8)'))
  dev.off()
}

GCM9_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/MIROC5', '\\.tif$', full.names=TRUE)
GCM9_50  <- stack(GCM9_50)

names(GCM9_50)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM9_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM9_50.tiff')), overwrite=TRUE)
  suitability_current19 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MIROC5_50_8.5_random.tif')))
  suitability_current19[MESS < 0] <- NA
  suitability_current19 <- crop(suitability_current19, extent(MESS))
  writeRaster(suitability_current19, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM9.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM9.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 MIROC5 (GCM9)'))
  dev.off()
}

GCM10_50 <- list.files('CHELSA_RASTERS_1km/2050/RCP85/MPIESMMR', '\\.tif$', full.names=TRUE)
GCM10_50  <- stack(GCM10_50)

names(GCM10_50)    <- c("biol_01",
                        "biol_02",
                        "biol_03",
                        "biol_04",
                        "biol_05",
                        "biol_06",
                        "biol_07",
                        "biol_08",
                        "biol_09",
                        "biol_10",
                        "biol_11",
                        "biol_12",
                        "biol_13",
                        "biol_14",
                        "biol_15",
                        "biol_16",
                        "biol_17",
                        "biol_18",
                        "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM10_50[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM10_50.tiff')), overwrite=TRUE)
  suitability_current20 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MPIESMMR_50_8.5_random.tif')))
  suitability_current20[MESS < 0] <- NA
  suitability_current20 <- crop(suitability_current20, extent(MESS))
  writeRaster(suitability_current20, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2050_GCM10.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2050_GCM10.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2050 MPIESMMR (GCM10)'))
  dev.off()
}






GCM1_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/ACCESS10', '\\.tif$', full.names=TRUE)
GCM1_70  <- stack(GCM1_70)

names(GCM1_70)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM1_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM1_70.tiff')), overwrite=TRUE)
  suitability_current21 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_ACCESS_70_8.5_random.tif')))
  suitability_current21[MESS < 0] <- NA
  suitability_current21 <- crop(suitability_current21, extent(MESS))
  writeRaster(suitability_current21, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM1.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM1.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 ACCESS (GCM1)'))
  dev.off()
}

GCM2_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/CESM1BGC', '\\.tif$', full.names=TRUE)
GCM2_70  <- stack(GCM2_70)

names(GCM2_70)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM2_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM2_70.tiff')), overwrite=TRUE)
  suitability_current22 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_CESM1BGC_70_8.5_random.tif')))
  suitability_current22[MESS < 0] <- NA
  suitability_current22 <- crop(suitability_current22, extent(MESS))
  writeRaster(suitability_current22, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM2.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM2.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 CESM1BGC (GCM2)'))
  dev.off()
}


GCM3_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/CESM1CAM5', '\\.tif$', full.names=TRUE)
GCM3_70  <- stack(GCM3_70)

names(GCM3_70)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM3_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM3_70.tiff')), overwrite=TRUE)
  suitability_current23 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_CESM1CAM5_70_8.5_random.tif')))
  suitability_current23[MESS < 0] <- NA
  suitability_current23 <- crop(suitability_current23, extent(MESS))
  writeRaster(suitability_current23, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM3.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM3.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 CESM1CAM5'))
  dev.off()
}



GCM4_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/CMCCCM', '\\.tif$', full.names=TRUE)
GCM4_70  <- stack(GCM4_70)

names(GCM4_70)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM4_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM4_70.tiff')), overwrite=TRUE)
  suitability_current24 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_CMCCCM_70_8.5_random.tif')))
  suitability_current24[MESS < 0] <- NA
  suitability_current24 <- crop(suitability_current24, extent(MESS))
  writeRaster(suitability_current24, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM4.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM4.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 CMCCCM (GCM4)'))
  dev.off()
}

GCM5_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/FIOESM', '\\.tif$', full.names=TRUE)
GCM5_70  <- stack(GCM5_70)

names(GCM5_70)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM5_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM5_70.tiff')), overwrite=TRUE)
  suitability_current25 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FIOESM_70_8.5_random.tif')))
  suitability_current25[MESS < 0] <- NA
  suitability_current25 <- crop(suitability_current25, extent(MESS))
  writeRaster(suitability_current25, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM5.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM5.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 FIOESM (GCM5)'))
  dev.off()
}

GCM6_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/GISSE2H', '\\.tif$', full.names=TRUE)
GCM6_70  <- stack(GCM6_70)

names(GCM6_70)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM6_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM6_70.tiff')), overwrite=TRUE)
  suitability_current26 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_GISSE2H_70_8.5_random.tif')))
  suitability_current26[MESS < 0] <- NA
  suitability_current26 <- crop(suitability_current26, extent(MESS))
  writeRaster(suitability_current26, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM6.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM6.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 GISSE2H (GCM6)'))
  dev.off()
}


GCM7_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/Inmcm4', '\\.tif$', full.names=TRUE)
GCM7_70  <- stack(GCM7_70)

names(GCM7_70)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM7_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM7_70.tiff')), overwrite=TRUE)
  suitability_current27 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_Inmcm4_70_8.5_random.tif')))
  suitability_current27[MESS < 0] <- NA
  suitability_current27 <- crop(suitability_current27, extent(MESS))
  writeRaster(suitability_current27, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM7.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM7.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 Inmcm4 (GCM7)'))
  dev.off()
}

GCM8_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/IPSLCM5AMR', '\\.tif$', full.names=TRUE)
GCM8_70  <- stack(GCM8_70)

names(GCM8_70)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM8_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM8_70.tiff')), overwrite=TRUE)
  suitability_current28 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_IPSLCM5AMR_70_8.5_random.tif')))
  suitability_current28[MESS < 0] <- NA
  suitability_current28 <- crop(suitability_current28, extent(MESS))
  writeRaster(suitability_current28, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM8.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM8.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 IPSLCM5AMR (GCM8)'))
  dev.off()
}

GCM9_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/MIROC5', '\\.tif$', full.names=TRUE)
GCM9_70  <- stack(GCM9_70)

names(GCM9_70)    <- c("biol_01",
                       "biol_02",
                       "biol_03",
                       "biol_04",
                       "biol_05",
                       "biol_06",
                       "biol_07",
                       "biol_08",
                       "biol_09",
                       "biol_10",
                       "biol_11",
                       "biol_12",
                       "biol_13",
                       "biol_14",
                       "biol_15",
                       "biol_16",
                       "biol_17",
                       "biol_18",
                       "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM9_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM9_70.tiff')), overwrite=TRUE)
  suitability_current29 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MIROC5_70_8.5_random.tif')))
  suitability_current29[MESS < 0] <- NA
  suitability_current29 <- crop(suitability_current29, extent(MESS))
  writeRaster(suitability_current29, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM9.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM9.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 MIROC5 (GCM9)'))
  dev.off()
}

GCM10_70 <- list.files('CHELSA_RASTERS_1km/2070/RCP85/MPIESMMR', '\\.tif$', full.names=TRUE)
GCM10_70  <- stack(GCM10_70)

names(GCM10_70)    <- c("biol_01",
                        "biol_02",
                        "biol_03",
                        "biol_04",
                        "biol_05",
                        "biol_06",
                        "biol_07",
                        "biol_08",
                        "biol_09",
                        "biol_10",
                        "biol_11",
                        "biol_12",
                        "biol_13",
                        "biol_14",
                        "biol_15",
                        "biol_16",
                        "biol_17",
                        "biol_18",
                        "biol_19")


for (i in 1:length(maxent_fitted)) {
  mx_fit <- readRDS(maxent_fitted[[i]])
  clim1 <- c("biol_01","biol_04","biol_05","biol_12","biol_14","biol_15" )
  MESS <- rmaxent::similarity(GCM10_70[[(clim1)]], mx_fit$swd)$similarity_min
  writeRaster(MESS, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_GCM10_70.tiff')), overwrite=TRUE)
  suitability_current30 <- raster(file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MPIESMMR_70_8.5_random.tif')))
  suitability_current30[MESS < 0] <- NA
  suitability_current30 <- crop(suitability_current30, extent(MESS))
  writeRaster(suitability_current30, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_FINAL_2070_GCM10.tiff')), overwrite=TRUE)
  
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_MESS_2070_GCM10.png'))
  png(f, 9, 6, units='in', res=300)
  plot(MESS, main=paste0(sp.nm[[i]], ' - Multivariate environmental similarity surfaces (MESS) 2070 MPIESMMR (GCM10)'))
  dev.off()
}




###################----------------------------------###################
###################--->   new Ensemble           <---###################
###################----------------------------------###################

spp_occ=  list.files('__OUTDIR__/plants_clim1/', 'occ\\.shp', recursive=TRUE, full=TRUE)
sp <- gsub('__OUTDIR__/plants_clim1//', '', spp_occ)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)
path <- '__OUTDIR__/plants_clim1/'


for (i in 1:length(spp_occ)) {
  GCM1_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_FINAL_2030_GCM1.tif'))))
  GCM2_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2030_GCM2.tif'))))
  GCM3_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2030_GCM3.tif'))))
  GCM4_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2030_GCM4.tif'))))
  GCM5_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2030_GCM5.tif'))))
  GCM6_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2030_GCM6.tif'))))
  GCM7_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2030_GCM7.tif'))))
  GCM8_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2030_GCM8.tif'))))
  GCM9_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2030_GCM9.tif'))))
  GCM10_30 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_FINAL_2030_GCM10.tif'))))
  
  Ens_10GCMs_30 <- (GCM1_30 + GCM2_30 + GCM3_30 + GCM4_30 + GCM5_30 + GCM6_30 + GCM7_30 + GCM8_30 + GCM9_30 +  GCM10_30)/10
  writeRaster(Ens_10GCMs_30, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_Ensem_30.tiff')), overwrite=TRUE)
  
}


for (i in 1:length(spp_occ)) {
  GCM1_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2050_GCM1.tif'))))
  GCM2_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2050_GCM2.tif'))))
  GCM3_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2050_GCM3.tif'))))
  GCM4_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2050_GCM4.tif'))))
  GCM5_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2050_GCM5.tif'))))
  GCM6_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2050_GCM6.tif'))))
  GCM7_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2050_GCM7.tif'))))
  GCM8_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2050_GCM8.tif'))))
  GCM9_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2050_GCM9.tif'))))
  GCM10_50 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_FINAL_2050_GCM10.tif'))))
  
  Ens_10GCMs_50 <- (GCM1_50 + GCM2_50 + GCM3_50 + GCM4_50 + GCM5_50 + GCM6_50 + GCM7_50 + GCM8_50 + GCM9_50 +  GCM10_50)/10
  writeRaster(Ens_10GCMs_50, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_Ensem_50.tiff')), overwrite=TRUE)
  
}


for (i in 1:length(spp_occ)) {
  GCM1_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2070_GCM1.tif'))))
  GCM2_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2070_GCM2.tif'))))
  GCM3_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2070_GCM3.tif'))))
  GCM4_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2070_GCM4.tif'))))
  GCM5_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2070_GCM5.tif'))))
  GCM6_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2070_GCM6.tif'))))
  GCM7_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2070_GCM7.tif'))))
  GCM8_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2070_GCM8.tif'))))
  GCM9_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]],'_FINAL_2070_GCM9.tif'))))
  GCM10_70 <- raster(file.path(path, paste0(sp[i], paste0("/", sp[[i]], '_FINAL_2070_GCM10.tif'))))
  
  Ens_10GCMs_70 <- (GCM1_70 + GCM2_70 + GCM3_70 + GCM4_70 + GCM5_70 + GCM6_70 + GCM7_70 + GCM8_70 + GCM9_70 +  GCM10_70)/10
  writeRaster(Ens_10GCMs_70, file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_Ensem_70.tiff')), overwrite=TRUE)
  
}



###################----------------------------------###################
###################--->        Plotting          <---###################
###################----------------------------------###################

# create folders for pngs if they don't exist
if(!dir.exists('__OUTDIR__/plants_clim1/2000')) dir.create('__OUTDIR__/plants_clim1/2000')
if(!dir.exists('__OUTDIR__/plants_clim1/2030')) dir.create('__OUTDIR__/plants_clim1/2030')
if(!dir.exists('__OUTDIR__/plants_clim1/2050')) dir.create('__OUTDIR__/plants_clim1/2050')
if(!dir.exists('__OUTDIR__/plants_clim1/2070')) dir.create('__OUTDIR__/plants_clim1/2070')



colfunc <- colorRampPalette(c("#dededb", "#FFFFB3","#8CFF19","#003300"))
files <- list.files('__OUTDIR__/plants_clim1/','_FINAL_CURRENT\\.tif$|_Ensem_[357]0\\.tif$', recursive=TRUE, full.names=TRUE)
#ff <- list.files('C:/Users/shab0021/Dropbox/1600SP/Rony_Data/Rony_data/Data/OUTPUT_10GCM_10SP_DATA/plants_clim1/','_current_aus_random_Alber\\.tif$|_Ensem_[57]0\\.tif$', recursive=TRUE, full.names=TRUE)
lapply(files, function(f) {
  r <- raster(f)
  if (is.na(raster::maxValue(r))) {
    show(paste(f, 'has no data to plot'))
    return()
  }
  r <- r/1000
  sp <- basename(dirname(f))
  yr <- sub('.*_([^_]+)\\.tif$', '\\1', f)
  yr <- if(yr=='CURRENT') 2000 else as.numeric(paste0(20, yr))
  f_out <- sprintf('%s/%s/%s', dirname(dirname(f)), yr, sub('\\.tif$', '.png', basename(f)))
  # option 1 (raster plot)
  #plot(r, axes=FALSE, xlab='', ylab='', main=sprintf('%s (%s)', sp, yr), col=colfunc(30))
  
  # option 2 (rasterVis plot)
  tryCatch({
	  p <- rasterVis::levelplot(r, col.regions=colfunc, at=seq(0, 1, length=30), scales=list(draw=FALSE),
				    colorkey=list(height=0.6), main=sprintf('%s (%s)', sp, yr),
				    margin=FALSE)
	  png(f_out, 13, 13, res=200, units='cm')
	  print(p)
	  dev.off()
          f_out
  }, error = function(e) {
  	message(e)
	NA
  })
})



###################----------------------------------###################
###################--->    #plot occ_for AUS
###################----------------------------------###################
library(ggplot2)
library(rgdal)

spp_occ= list.files('__OUTDIR__/plants_clim1/', 'occ\\.shp', recursive=TRUE, full=TRUE)

sp <- gsub('__OUTDIR__/plants_clim1//', '', spp_occ)
sp <- gsub('/occ.shp', '', sp)
sp.nm <- gsub('_', ' ', sp)

aus_shp <- shapefile("additional_data_for_modeling/Australia/Aust_states.shp")
aus_df <- fortify(aus_shp, region = 'State')
for (i in 1:length(spp_occ)) {
  sp_oc <- readOGR(spp_occ[[i]])
  sp_oc <- spTransform(sp_oc, CRS(proj4string(aus_shp)))
  sp_oc_aus <- raster::crop(sp_oc, aus_shp)
  p <- ggplot() +
    geom_polygon(data = aus_df , aes(x = long, y = lat, group = group), fill = 'white', color = 'black')
  if (is.null(sp_oc_aus)) {
    show(paste('No occ aus data to plot for', sp[i]))
  } else {
    p <- p + geom_point(data= data.frame(sp_oc_aus@coords), aes(x=coords.x1, y=coords.x2), color="#0A8064", size = 0.55)
  }
  p <- p + ggtitle(paste0(sp.nm[[i]], ' occurrence points in Australia'))
  f <- file.path(paste0('__OUTDIR__/plants_clim1/', sp[i]), paste0(sp[[i]], '_occ_AUS.png'))
  png(f, 6, 5, units='in', res=200)
  print(p)
  dev.off()
}


###################----------------------------------###################
###################---> zonal extract               <---################
###################----------------------------------###################

postcode_aus <- shapefile('postcode_file/POA_2016_AUST.shp') 
postcode_aus <- spTransform(postcode_aus, CRS="+proj=aea +lat_0=0 +lon_0=132 +lat_1=-18 +lat_2=-36 +x_0=0 +y_0=0 +ellps=GRS80 +units=m +no_defs")

files <- list.files('__OUTDIR__/plants_clim1', '\\.tif$', full=TRUE, recursive=TRUE) %>% 
    grep("Ensem|FINAL_CURRENT", ., value=TRUE, invert=F) %>%
    split(sub('__OUTDIR__/plants_clim1/([^/]+)/.*', '\\1', .))

list_s <- lapply(files, stack)

dir.create(file.path('__OUTDIR__/plants_clim1', 'zonal'), recursive=TRUE)

for (i in 1:length(list_s)) {
    s <- list_s[[i]]
    # name the cells appropriately
    cells <- cellFromPolygon(s, postcode_aus) %>%
        setNames(postcode_aus$POA_CODE16)
    # make sure to use na.rm=TRUE for mean call
    vals <- as.data.frame(s[unlist(cells)]) %>%
        mutate(postcode=rep(names(cells), lengths(cells))) %>%
        group_by(postcode) %>%
        summarise_all(mean, na.rm=TRUE)
    # write csv
    for (j in 2:ncol(vals)) {
        f <- file.path('__OUTDIR__/plants_clim1', 'zonal', paste0(names(vals)[j], '.csv'))
        cols <- c("postcode", names(vals)[j])
        write.table(vals %>% select(all_of(cols)), file=f, sep=",", row.names=FALSE, col.names=FALSE)
    }
}
