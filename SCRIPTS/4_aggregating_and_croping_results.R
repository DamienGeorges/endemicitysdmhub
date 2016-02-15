##' ------------------------------------------------------------------
##' @title Reshaping the modelling outputs to be able to use them in
##'   the next step mechanistic model
##' @author Damien G
##' @date 09/02/2016
##' 
##' @description
##' In this script we will project sdm produced via the script 2_species_niche_modelling.R
##' 
##' Species list :
##' - Campanula pulla
##' - Dianthus alpinus
##' - Festuca peudodura
##' - Primula clusiana
##'
##' Formal predictions:
##'   - 100m resolution for the full alpine arc (ETRS89)
##'   
##' Reshaped predictions:
##'   - 250m resolution for the area c(264089.6, 450820.4, 1870309, 1969299)
##'     proj = +proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs
##'   
##'   
##' @licence GPL-2
##' ----------------------------------------------------------------------------

setwd("/nfs_scratch2/dgeorges/Endemicity_SDM/WORKDIR/")
library(raster)
library(dplyr)

sp.list <- c("Camppull", "Dianalpi", "Festpseu", "Primclus")
scen.dir <- c("EOBS_1970_2005", "rcp26_ICHEC-EC-EARTH_SMHI-RCA4", "rcp45_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4", "rcp85_ICHEC-EC-EARTH_DMI-HIRHAM5")
scen.year <- seq(2020, 2090, 10)
mod.id <- c("EMcaByTSS", "EMmedianByTSS", "EMwmeanByTSS")

projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

new.proj <- "+proj=aea +lat_1=43 +lat_2=62 +lat_0=30 +lon_0=10 +x_0=0 +y_0=0 +ellps=intl +units=m +no_defs"
new.ext <- extent(c(264089.6, 450820.4, 1870309, 1969299))
new.ext.etrs89.coord <- coordinates(spTransform(SpatialPoints(data.frame(x = c(new.ext@xmin, new.ext@xmax), y = c(new.ext@ymin, new.ext@ymax)),
                                proj4string = CRS(new.proj)), CRS(projETRS89)))
new.ext.etrs89 <- extent(c(new.ext.etrs89.coord[1,1], new.ext.etrs89.coord[2,1], new.ext.etrs89.coord[1,2], new.ext.etrs89.coord[2,2]))

new.ras.ref <- raster(new.ext, crs = CRS(new.proj), resolution = 250)

output.dir <- "/nfs_scratch2/dgeorges/Endemicity_SDM/RESULTS/crop_and_scaled_outputs"
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)

## define which tiles intersect with the new area
tiles.table <- read.table("../DATA/raster_tiles_100000m_45t.txt", sep = "\t", header = TRUE)  
tiles.table <- tiles.table %>% mutate(is.in.new.ext = (xu >= new.ext.etrs89@xmin & 
                                                         xl <= new.ext.etrs89@xmax &
                                                         yu >= new.ext.etrs89@ymin & 
                                                         yl <= new.ext.etrs89@ymax))
tiles.inter.id <- which(tiles.table$is.in.new.ext)

# sp_ <- sp.list[1]
# scen_ <- scen.dir[1]
# year_ <- scen.year[1]
# mod_ <- mod.id[1]

for(sp_ in sp.list){
  cat("\n\n\n>", sp_)
  for(scen_ in scen.dir){
    cat("\n\t>", scen_)
    for(mod_ in mod.id){
      cat("\n\t\t>", mod_)
      ## deal with current conditions
      if(scen_ == "EOBS_1970_2005"){
        ## load the formal raster
        ras.form <- raster(paste0(sp_, "/proj_", scen_, "_100m",  
                                  "/individual_projections/", sp_, "_", mod_, "_mergedAlgo_mergedRun_mergedData.grd"))
        ## reproject and save the raster
        ras.new.name <- file.path(output.dir, paste0(sp_, "_", mod_, "_", scen_, ".img"))
        ras.new <- projectRaster(ras.form, new.ras.ref, filename = ras.new.name, overwrite = TRUE)
      } else {
        ## deal with future conditions
        for(year_ in scen.year){
          cat("\t", year_)
          ## get the list of raster files
          ras.form.files <- paste0(sp_, "/proj_", scen_, "_", year_, "_", tiles.inter.id,  
                                   "/individual_projections/", sp_, "_", mod_, "_mergedAlgo_mergedRun_mergedData.grd") 
          ## load raster tiles
          ras.form.list <- lapply(ras.form.files, raster)
          
          ## merge raster tiles
          ras.form <- do.call(merge, ras.form.list)
          
          ## reproject and save the raster
          ras.new.name <- file.path(output.dir, paste0(sp_, "_", mod_, "_", scen_, "_", year_, ".img"))
          ras.new <- projectRaster(ras.form, new.ras.ref, filename = ras.new.name, overwrite = TRUE)
        } ## end loop over years
      }
    } ## end loop over models
  } ## end loop over scenar
}## end loop over species


## extracting ensemble models scores and thresholds and store them into a table
library(biomod2)
# sp_ <- sp.list[1]

bm.eval.list <- vector(mode = 'list', length(sp.list))
names(bm.eval.list) <- sp.list

for(sp_ in sp.list){
  bm.em <- get(load(paste0(sp_, "/", sp_, ".endemicensemble.models.out")))
  bm.eval <- get_evaluations(bm.em, as.data.frame = TRUE)
  bm.eval <- bm.eval %>% mutate(sp.name = sp_, mod.name.short = sub("_mergedAlgo_mergedRun_mergedData", "", Model.name))
  bm.eval.list[[sp_]] <- bm.eval
}

bm.eval <- do.call(rbind, bm.eval.list)
write.table(bm.eval, file = file.path(output.dir, "ensemble_models_scores.txt"), sep = "\t", 
            quote = FALSE, col.names = TRUE, row.names = FALSE)


