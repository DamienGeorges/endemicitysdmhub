# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# 	Modeling of 4 Austria endemics species niches
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# Authors -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
Damien G. for modeling part
Julien R. for data formating part
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# Date -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
14/02/25 - formal version
16/02/11 - new version (see log)
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# Log -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-= #
- 16/02/11 (damien g.)
  Previous version of the models was done with bio2 and 
  bio14 as climatic prediction. Because they are complex 
  to explain and because they are not so much changing
  in the future, we decided to replace them by bio1 and 
  bio12.
  We also change climatic data use to the latest one we 
  have that are at 100m resolution on the whole alpine 
  arc.
  Projections are done for current and future conditions
  ("rcp26_ICHEC-EC-EARTH_SMHI-RCA4", 
  "rcp45_CNRM-CERFACS-CNRM-CM5_SMHI-RCA4" and 
  "rcp85_ICHEC-EC-EARTH_DMI-HIRHAM5") from 2020 to 2090 
  every 10 years. Because of the resolution, for memory
  consuming purpose, the alpine arc was splitted into 
  45 tilles (10km x 10km). This tilles can be merged in
  R with raster::merge() function.

- 16/02/15 (damien g.)
  Reproject ensemble models outputs to fit a 250m resolution
  grid.
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# Summary =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
We are interested here in modeling niche of 4 endemics
species of Austria. We want to build quite easy interpretable
models so we use only few explanatory variables.
The SDM modeling procedure is done with biomod2 package.

** Species **
 - Campanula pulla
 - Dianthus alpinus
 - Festuca peudodura
 - Primula clusiana
Data are "true" presences/absences data

** Study area **
 - Austria for calibration (because of endemic)
 - Alpine arc for projections (some part are exclude because
	soil data (% of carbonate) was not computed yet... We 
	will be able to re-project it on the whole arc latter 
	if needed)

** Explanatory variables **
	- soil (% of carbonate in soil)
	- temperature (bio1 == Mean Anual Temperature)
	- precipitation (bio12 == Annual Precipitation)
This 3 variables are almost independent one with another

** Data resolution and projection system **
 - for calibration : 100m resolution data on a regular grid
 - proj : ETRS89

** Modeling parameters **
 - 4 models : SRE, GLM (quadratic), GAM (mgcv), RF
 - models evaluation by TSS
 - Presences/absences dataset (No pseudo-absences)
 - 5 Cross Validation at 70%
 - 3 variable importance calculation
 - 3 ensemble models (median, weighted mean, committee averaging)
 - projection at 100m resolution over the whole alpine arc

# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# Files description =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
The whole simulation is composed of 4 directories :

** DATA **
The data used for this simulation.
 -  species_occ.RData, an Robject (a 
		data.frame with xy coordinates, the presences/absences
		data for our 4 species (250m resolution data) 
 - raster_tiles_1000000m_45.txt a table containing the extent 
   of every alpine arc tile considered at projection step.
 - well_formated_data is the directory containing the explanatory
		data used to make projections. Raster are available at
    100m  for current and 3 scenario and are coherent in term 
    of covered area.

** RESULTS **
This file contains the main results of the SDM modeling.
Each species has 2 .csv files (should be .zip compressed files).
This tables contain prediction of all sdm models (25) over
all cells of study area at 250m and 1km resolution. Predictions
are on a 0-1000 scale corresponding to 0-1 probabilities * 1000. 
Note : This files are quite big and the same information stored
	model by model is available in WORDIR/sp_name/pro_current_xx/
	under raster version. This files are maybe easier to manage.
Note2 : Another file called "modeling_characteristics.txt"
	is a file where all modeling architecture specification
	are stored (R version, packages version... ). It should be useful
	in case we want to rerun something.
Note3 : Explanatory variables for both resolution (exactly same format)
	are also stored within this directory ("Explanatory_variables_xx").

** SCRIPT **
 - 1_data_formatting.R :  the script to transform raw data into
		well formated ones.
 - 2_species_niche_modelling.R : the main script use to produce
		SDMs. This file is associated to 2_species_niche_modelling.oar
    and 2_species_niche_modelling.params that are respectivelly the 
    script to lunch a oar grid campain and the associated parameters
  - 3_species_niche_projection.R : the script where all the future
    projections are computed. This file is associated to 3_species_niche_projection.oar
    and 3_species_niche_projection.params that are respectivelly the 
    script to lunch a oar grid campain and the associated parameters
  - 4_aggregating_and_cropung_results.R: This script will reproject formmal
    100m ETRS89 resolution ensemble projections rasters into 250m ones
    on a reduced area.


** WORKDIR **
The outputs of biomod2 modeling (a directory by species).
Each directory is composed of :
	- models : the directory where all biomod2 models objects
		are stored
	- proj_xx : the dir where
		models projections are stored. In this dir you have :
			- individual_projections : all individual models
				projections under rater format
			- xx.currentyy.projection.out and xx.currentyy.ensembleprojection.out
				the biomod2 projection objects
	- xx.models.out & xx.ensemble.models.out : the biomod2
		modelling output objects.
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

# Notes =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #
# =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=- #

