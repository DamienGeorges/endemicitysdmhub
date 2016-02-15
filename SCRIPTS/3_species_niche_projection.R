##' ------------------------------------------------------------------
##' @title Projection of 4 austiria endemics species niches
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
##' Input data characteristics :
##' - Studied area : 
##'     - Austria for calibration (releves at 250m resolution / climate 100m)
##'     - Alpine arc (100m) for projection
##' - explanatory variables : 
##'     - soil (% of carbonate in soil)
##'     - temperture (bio1 == Annual Mean Temperature)
##'     - precipitation (bio12 == Annual Precipitation)
##'     
##' Projection characteristics:
##'   - resolution: 100m
##'   - dates: 2020 to 2090 every 10 years
##'   - scenar: "rcp26_ICHEC-EC-EARTH_SMHI-RCA4", "rcp45_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4", 
##'             "rcp85_ICHEC-EC-EARTH_DMI-HIRHAM5"
##' 
##' @log
##'   - 10/02/2016 (damien): with the 100m resolution median ensemble models seems to
##'     fail quite often (too much memory consumtion?), we tried to rerun the campain 
##'     omitting median ensemble models
##' 
##' @licence GPL-2
##' ----------------------------------------------------------------------------

rm(list=ls())

## getting parameters ----------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
sp.name <- as.character(args[1])
scen <- as.character(args[2])
date <- as.character(args[3])
tile <- as.numeric(args[4])

# ## test
# sp.name <- "Camppull"
# scen <- "rcp26_ICHEC-EC-EARTH_SMHI-RCA4"
# date <- "2020"
# tile <- 1

check.computed <- TRUE ## not recompute already computed projections

cat("\n> sp.name:", sp.name, "\t scen:", scen, "\t date:", date, "\t tile:", tile, "\n")

## define the correct workdir --------------------------------------------------
## on Luke
setwd("/nfs_scratch2/dgeorges/Endemicity_SDM/WORKDIR/")

## library loading -------------------------------------------------------------
.libPaths('/nfs_scratch2/emabio/R_PKG_LUKE')
library(biomod2)
rasterOptions(tmpdir = "/nfs_scratch2/dgeorges/R_raster_georges/",
              tmptime = 24, ## time after which raster tmp files will be deleted
              chunksize = 5e+08, ## size of blocks that will be written on hardrive (for I/O optimisation)
              maxmemory = 1e+09, ## max number of cell loaded in the memory (for I/O optimisation)
              overwrite = TRUE)


## define data paths -----------------------------------------------------------
data_dir <- "../DATA"
output_dir <- "../RESULTS"

## load explanatory varaiables -------------------------------------------------
expl.100m <- raster::stack( c( bio1 = file.path("../DATA/well_formated_data", scen, date, "bio_1.grd"),
                              bio12 = file.path("../DATA/well_formated_data", scen, date, "bio_12.grd"),
                              carbon = "../DATA/well_formated_data/carbon.grd"))
## get the extent of the tile
tiles.table <- read.table("../DATA/raster_tiles_100000m_45t.txt", sep = "\t", header = TRUE)
tile.extent <- raster::extent(tiles.table[tile, "xl"], tiles.table[tile, "xu"],
                              tiles.table[tile, "yl"], tiles.table[tile, "yu"])

## crop the environmental stack
expl.100m <- raster::stack(raster::crop(expl.100m, tile.extent))

## load models and ensemble models
bm.mod <- get(load(paste0(sp.name, "/", sp.name,".endemic.models.out")))
bm.em <- get(load(paste0(sp.name, "/", sp.name, ".endemicensemble.models.out", sep="")))

if(check.computed){
  ## generate the path to the projection and ensemble projection files
  bm.proj.100m.file <- file.path(getwd(), sp.name, paste0("proj_", scen, "_", date, '_', tile), paste0(sp.name, ".", scen, "_", date, '_', tile, ".projection.out"))
  bm.ensproj.100m.file <- file.path(getwd(), sp.name, paste0("proj_", scen, "_", date, '_', tile), paste0(sp.name, ".", scen, "_", date, '_', tile, ".ensemble.projection.out"))
  ## check file existance
  bm.proj.100m.exists <- file.exists(bm.proj.100m.file)
  bm.ensproj.100m.exists <- file.exists(bm.ensproj.100m.file)
  ## load existing files
  if(bm.proj.100m.exists) bm.proj.100m <- get(load(bm.proj.100m.file))
  if(bm.ensproj.100m.exists) bm.ensproj.100m <- get(load(bm.ensproj.100m.file))  
} 

# do projections ---------------------------------------------------------------
if(!exists('bm.proj.100m')){
  bm.proj.100m <- BIOMOD_Projection(modeling.output = bm.mod,
                                    new.env = expl.100m,
                                    proj.name = paste0(scen, '_', date, '_', tile),
                                    selected.models = 'all',
                                    compress = 'xz',
                                    do.stack = FALSE,
                                    build.clamping.mask = FALSE)
}

if(!exists('bm.ensproj.100m')){
  bm.ensproj.100m <- BIOMOD_EnsembleForecasting(projection.output = bm.proj.100m,
                                                EM.output = bm.em,
#                                                 selected.models = grep('EMmedian', get_built_models(bm.em), value = TRUE, invert = TRUE),
                                                binary.meth = 'TSS',
                                                do.stack = FALSE,
                                                build.clamping.mask = FALSE)
}

q("no")

## end of script ---------------------------------------------------------------

# ## generate parms files --------------------------------------------------------
# sp.name <- c("Camppull", "Dianalpi", "Festpseu", "Primclus") 
# scen <- c("rcp26_ICHEC-EC-EARTH_SMHI-RCA4", "rcp45_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4", "rcp85_ICHEC-EC-EARTH_DMI-HIRHAM5")
# date <- seq(2020, 2090, 10)
# tile <- 1:45
# 
# params <- expand.grid(sp.name = sp.name, 
#                       scen = scen, 
#                       date = date,
#                       tile = tile)
# ## note: some tiles didn't have any values so they will failed but it is not a problem
# tiles.no.vals <- c(4, 5, 6, 7, 8, 9, 13, 15, 16, 17, 18, 27, 37, 38, 39)
# params <- params %>% filter(!is.element(tile, tiles.no.vals))
# 
# # ## debug
# # params <- params[1:5, ]
# # params <- params[ !proj.file.list$proj.file.ok, ]
# 
# write.table(params, file = "/nfs_scratch2/dgeorges/Endemicity_SDM/SCRIPTS/3_species_niche_projection.params", sep = "\t", 
#             quote = FALSE, col.names = FALSE, row.names = FALSE, append = FALSE)
# # ## end generate parms files ----------------------------------------------------

# ## check campain state ---------------------------------------------------------
# library(dplyr)
# params <- read.table("/nfs_scratch2/dgeorges/Endemicity_SDM/SCRIPTS/3_species_niche_projection.params",
#                      header = FALSE, sep = "\t", stringsAsFactors = FALSE, 
#                      col.names = c("sp.name", "scen", "date", "tile"))
# ## note: some tiles didn't have any values so they will failed but it is not a problem
# tiles.no.vals <- c(4, 5, 6, 7, 8, 9, 13, 15, 16, 17, 18, 27, 37, 38, 39)
# params <- params %>% filter(!is.element(tile, tiles.no.vals))
# 
# out.dir <- "/nfs_scratch2/dgeorges/Endemicity_SDM/WORKDIR"
# proj.file.list <- params %>% mutate(proj.file = file.path(out.dir, sp.name, paste0("proj_", scen, "_", date, "_", tile), paste0(sp.name, ".", scen, "_", date, "_", tile, ".projection.out")),
#                                     ens.proj.file = file.path(out.dir, sp.name, paste0("proj_", scen, "_", date, "_", tile), paste0(sp.name, ".", scen, "_", date, "_", tile, ".ensemble.projection.out")),
#                                     proj.file.ok = file.exists(proj.file),
#                                     ens.proj.file.ok = file.exists(ens.proj.file))
# mean(proj.file.list$proj.file.ok)
# mean(proj.file.list$ens.proj.file.ok)
# sum(!proj.file.list$ens.proj.file.ok)
# 
# proj.file.list %>% filter(!proj.file.ok)
# proj.file.list %>% filter(!ens.proj.file.ok)
# 
# ## end check campain state -----------------------------------------------------