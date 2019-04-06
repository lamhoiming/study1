# R script to read in RIOPS or GIOPS netCDF files and export snow and ice layers
# as separate geoTiffs
#
# Source: 
# Modelled ice and snow values
# http://dd.weather.gc.ca/model_giops/doc/README_GIOPS.txt
# http://navigator.oceansdata.ca/thredds/catalog.html
# RIOPS daily averages, http://navigator.oceansdata.ca
#
# Data from 
# http://navigator.oceansdata.ca/thredds/dodsC/riops/daily/20161226_2D.nc.html
# are projected to:
# CRS("+proj=stere +a=6371229.0 +b=6371229.0 +lat_0=90 +lat_ts=90 +lon_0=-100 +k=0.9330124 +x_0=4245000.0 +y_0=5295000.0 +units=m +no_defs") 
# which is projected to the Caandian Polar Stereographic projection:
# CRS("+proj=stere +lat_0=90 +lat_ts=90 +lon_0=-100 +k=0.994 +x_0=2000000 +y_0=2000000 +datum=WGS84 +units=m +no_defs") 
#
# The following type of warning can be ignored:
# 	Warning messages:
# 	1: In .getCRSfromGridMap4(atts) : cannot process these parts of the CRS:
# 	earth_radius=6371229
#
####################################################################################

library(rgdal)
library(raster)

# INPUT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Set RIOPS/GIOPS file
RIOPS_dir <- "C:/Users/Torsten/Documents/CIS_DA/RIOPS/2016/"
setwd(RIOPS_dir)

# Projection of the input RIOPS netCDF files from navigator.oceandata.ca
RIOPS_pol_ster <- CRS("+proj=stere +a=6371229.0 +b=6371229.0 +lat_0=90 +lat_ts=90 +lon_0=-100 +k=0.9330124 +x_0=4245000.0 +y_0=5295000.0 +units=m +no_defs") 

# Projection of output files (Canadian Polar Stereographic)
Can_pol_ster <- CRS("+proj=stere +lat_0=90 +lat_ts=90 +lon_0=-100 +k=0.994 +x_0=2000000 +y_0=2000000 +datum=WGS84 +units=m +no_defs") 

# Set the desired output resolution (m)
Out_resolution <- 12500

# Set the output resolution text; this is used in the filesnames
Out_res_text <- "12p5km"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> End of input

# Loop through RIOPS netCDF files in the directory
RIOPS_files <- dir(pattern = '*.nc')

for(g in 1:length(RIOPS_files)){
	RIOPS_file <- RIOPS_files[g]
	
	# Get RIOPS file name
	RIOPS_base <- sub("^([^.]*).*", "\\1", basename(RIOPS_file)) 
	
	# Snow image
	RIOPS_snow <- raster(paste(RIOPS_dir,RIOPS_file,sep=""),varname="isnowvol")	
	names(RIOPS_snow) <- "Snow_cm"
	crs(RIOPS_snow) <- RIOPS_pol_ster
	RIOPS_snow_pr <- projectRaster(RIOPS_snow, res=Out_resolution, crs=Can_pol_ster, method="bilinear")
	writeRaster(RIOPS_snow_pr, filename=paste(RIOPS_base,"_",Out_res_text,"_snow.tif",sep=""), format="GTiff", overwrite=TRUE)

	# Ice image
	RIOPS_ice <- raster(paste(RIOPS_dir,RIOPS_file,sep=""),varname="iicevol")	
	names(RIOPS_ice) <- "Ice_m"
	crs(RIOPS_ice) <- RIOPS_pol_ster
	RIOPS_ice_pr <- projectRaster(RIOPS_ice, res=Out_resolution, crs=Can_pol_ster, method="bilinear")
	writeRaster(RIOPS_ice_pr, filename=paste(RIOPS_base,"_",Out_res_text,"_ice.tif",sep=""), format="GTiff", overwrite=TRUE)
}

# Clear variables
# rm(list=ls(all=TRUE))
gc()
print("Done")
# END