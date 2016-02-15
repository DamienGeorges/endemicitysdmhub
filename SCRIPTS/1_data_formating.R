## create the suitable version of environmental raster for Endemicity study

setwd("~/Work/PROJECTS/Endemicity_SDM/WORKDIR/")
library(raster)
library(foreach)
library(doParallel)


projETRS89 <- "+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs"

input.path <- "/media/georgeda/equipes/emabio/GIS_DATA/Alpes/PHYSIQUE/CLIMATE/Alp_dwnscl_100m_15yr_newprec"
ref.dir <- "EOBS_1970_2005"
scen.dir <- c("rcp26_ICHEC-EC-EARTH_SMHI-RCA4", "rcp45_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4", "rcp85_ICHEC-EC-EARTH_DMI-HIRHAM5")
scen.year <- seq(2020, 2090, 10)

output.dir <- "~/Work/PROJECTS/Endemicity_SDM/DATA/well_formated_data/"

## we take bio_1 in current condition as reference
ref.grid <- raster(file.path(input.path, ref.dir, "bio", "bio_1.asc"), crs = CRS(projETRS89))
ref.grid <- reclassify(ref.grid, c(-Inf, Inf, 1),  overwrite = TRUE)
writeRaster(ref.grid, "refgrid.grd") 

env.var <- c("bio_1", "bio_12")

## save current bioclim
for(env.var_ in env.var){
  cat("\n> reshaping", env.var_)
  output.dir_ <- file.path(output.dir, ref.dir)
  dir.create(output.dir_, showWarnings = FALSE, recursive = TRUE)
  bio.tmp <- raster(file.path(input.path, ref.dir, "bio", paste0(env.var_, ".asc")), crs = CRS(projETRS89))
  bio.tmp <- projectRaster(bio.tmp, ref.grid)
  bio.tmp <- bio.tmp * ref.grid 
  writeRaster(bio.tmp, filename = file.path(output.dir_, paste0(env.var_, ".grd")), overwrite = TRUE)
}


# scen.dir_ <- scen.dir[1]
# scen.year_ <- scen.year[1]
# env.var_ <- env.var[1]
## save future bioclim
cl <- makeCluster(8)
registerDoParallel(cl)

for(scen.dir_ in scen.dir){
  cat(">", scen.dir_)
  p.out <- foreach(scen.year_ = scen.year, 
                   .export = c("input.path", "output.dir", "scen.dir_", "env.var"),
                   .packages = c("raster")) %dopar% {
    for(env.var_ in env.var){
      cat("\n> reshaping", env.var_, scen.dir_, scen.year_)
      output.dir_ <- file.path(output.dir, scen.dir_, scen.year_)
      dir.create(output.dir_, showWarnings = FALSE, recursive = TRUE)
      bio.tmp <- raster(file.path(input.path, scen.dir_, "bio", paste0("bio", scen.year_, sub("bio", "", env.var_))), crs = CRS(projETRS89))
      bio.tmp <- projectRaster(bio.tmp, ref.grid)
      bio.tmp <- bio.tmp * ref.grid 
      writeRaster(bio.tmp, filename = file.path(output.dir_, paste0(env.var_, ".grd")), overwrite = TRUE)
    }
  }
}

stopCluster(cl)

## save carbonate raster
carbon <- raster("/media/georgeda/equipes/emabio/GIS_DATA/Alpes/PHYSIQUE/TOPO/SOIL/carbon_whole_alps.img")
carbon[carbon[] > 100] <- NA
carbon[carbon[] < 0] <- NA
carbon <- projectRaster(carbon, ref.grid)
carbon <- carbon * ref.grid
writeRaster(carbon, filename = file.path(output.dir, "carbon.grd"), overwrite = TRUE)

