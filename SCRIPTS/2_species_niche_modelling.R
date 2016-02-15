##' ------------------------------------------------------------------
##' @title Modeling of 4 austiria endemics species niches
##' @author Damien G
##' @date 09/02/2016
##' 
##' @description
##' In this script we will model the current climatic niches
##' of a list of 4 autrian endemic species.
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
##' Modeling characteristics :
##' - 4 models : SRE, GLM (quadratic), GAM (mgcv), RF
##' - models evaluation by TSS
##' - Presences/absences dataset (No pseudo-absences)
##' - 5 Cross Validation at 70%
##' - 3 variable importance calculation
##' - 3 ensemble models (median, weighted mean, commitee averaging)
##' 
##' @log
##'   - 09/02/2016 (damien)
##'     We redo the full modelling/projection procedure for our 4 endemic species.
##'     Changes compare to the previous version are:
##'       - explanatory var have been changed (new: bio1, bio12, carbon 
##'         old: bio2,bio14, carbon)
##'       - we are now using the 100m res version of downscaled bioclim
##'       - we are considering different scenario of climatic change
##' 
##' @licence GPL-2
##' ----------------------------------------------------------------------------


rm(list=ls())

## getting parameters ----------------------------------------------------------
args <- commandArgs(trailingOnly=TRUE)
sp.name <- as.character(args[1])

## define the correct workdir --------------------------------------------------
## on Luke
setwd("/nfs_scratch2/dgeorges/Endemicity_SDM/WORKDIR/")

## library loading -------------------------------------------------------------
.libPaths('/nfs_scratch2/emabio/R_PKG_LUKE')
library(biomod2)
rasterOptions(overwrite = TRUE, 
              tmpdir = "/nfs_scratch2/dgeorges/R_raster_georges/")

## define data paths -----------------------------------------------------------
data_dir <- "../DATA"
output_dir <- "../RESULTS"

##' 0. load data ---------------------------------------------------------------

## load the table containing species occ for our 4 endemics
load(file.path(data_dir,"species_occ.RData"))

## get species pres/abs and associated coordiantes
sp.dat <- species_occ[, sp.name]
sp.coords <- species_occ[, c("x", "y")]

## load current environment
curr100m <- raster::stack( c( bio1 = '../DATA/well_formated_data/EOBS_1970_2005/bio_1.grd',
                              bio12 = '../DATA/well_formated_data/EOBS_1970_2005/bio_12.grd',
                              carbon = '../DATA/well_formated_data/carbon.grd'))

# 1. data initialisation -------------------------------------------------------
bm.init <- BIOMOD_FormatingData( resp.var = sp.dat, 
                                 expl.var = curr100m,
                                 resp.xy = sp.coords,
                                 resp.name = sp.name)
  
# 2. modeling options def -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
bm.opt <- BIOMOD_ModelingOptions( GLM = list(type = 'quadratic', interaction.level = 0),
                                  GAM = list(algo = 'GAM_mgcv', k = 3))

# 3. do modeling ---------------------------------------------------------------
bm.mod <- BIOMOD_Modeling(bm.init, 
                          models = c('GLM','GAM','RF','SRE'),
                          models.options = bm.opt,
                          NbRunEval = 5, 
                          DataSplit = 70, 
#                           Prevalence = 0.5, 
                          VarImport = 3, 
                          models.eval.meth = c('TSS'),
                          SaveObj = TRUE,
                          rescal.all.models = FALSE,
                          do.full.models = FALSE,
                          modeling.id = 'endemic')
  
# 3. do ensemble modeling ------------------------------------------------------
bm.em <- BIOMOD_EnsembleModeling( modeling.output = bm.mod,
                                  chosen.models = 'all',
                                  eval.metric = c('TSS'),
                                  eval.metric.quality.threshold = c(0.3),
                                  prob.mean = F,
                                  prob.cv = F,
                                  prob.ci = F,
                                  prob.ci.alpha = 0.05,
                                  prob.median = T,
                                  committee.averaging = T,
                                  prob.mean.weight = T,
                                  prob.mean.weight.decay = 'proportional',
                                  em.by = 'all',
                                  VarImport = F) ## realy time consuming

# 4. make projections ----------------------------------------------------------
bm.proj.100m <- BIOMOD_Projection(modeling.output = bm.mod,
                                  new.env = curr100m,
                                  proj.name = "EOBS_1970_2005_100m",
                                  selected.models = 'all',
                                  compress = 'xz',
                                  do.stack = FALSE,
                                  build.clamping.mask = FALSE)

bm.ensproj.100m <- BIOMOD_EnsembleForecasting(projection.output = bm.proj.100m,
                                              EM.output = bm.em,
                                              binary.meth = 'TSS',
                                              do.stack = FALSE,
                                              build.clamping.mask = FALSE)

q("no")
