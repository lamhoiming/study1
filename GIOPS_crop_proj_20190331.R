# Open raw GIOPS or RIOPS file from MRD, crop to Arctic, 
# rasterise to a desired resolution and project to a 
# Canadian Polar Stereographic projection, and save as .tif
#
# Input files are e.g. SDV_20160413_00000000_0_meters_above_sea_level.txt
# these files use the WGS84 datum
#
# Run separately for snow and ice files by commenting out the lines you don't want running
#
# The following type of warning can be ignored:
# 	Warning messages:
# 	1: In .getCRSfromGridMap4(atts) : cannot process these parts of the CRS:
# 	earth_radius=6371229
#
####################################################################################

library(sp)
library(spdep)
library(rgdal)
library(ncdf4)
library(raster)

# INPUT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Set working directory
wdir <- "D:/Data/CryoSat-2/study1/"
# Assumes the MRD text files do not have the "." in the middle of the filename: e.g.
# "SDV_20160302_00000000_0_meters_above_sea_level", instead of "SDV_20160302_00000000_0.0_meters_above_sea_level"

setwd(wdir)

# Projections
Can_pol_ster <- CRS("+proj=stere +lat_0=90 +lat_ts=90 +lon_0=-100 +k=0.994 +x_0=2000000 +y_0=2000000 +datum=WGS84 +units=m +no_defs") 
WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") 

# Set the desired output resolution (m)
Out_resolution <- 15000

# Set the output resolution text; this is used in the filesnames
Out_res_text <- "15km"

# Set RIOPS file from http://navigator.oceansdata.ca to use as the extent; this file will be emptied, so it doesn't matter
# which file you indicate here
RIOPS_dir <- "D:/Data/CryoSat-2/study1/"
RIOPS_file <- "riops_201612_2D.nc"

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> End of input

# Projection used by RIOPS in http://navigator.oceansdata.ca 
RIOPS_pol_ster <- CRS("+proj=stere +a=6371229.0 +b=6371229.0 +lat_0=90 +lat_ts=90 +lon_0=-100 +k=0.9330124 +x_0=4245000.0 +y_0=5295000.0 +units=m +no_defs") 
RIOPS_snow <- raster(paste(RIOPS_dir,RIOPS_file,sep=""),varname="isnowvol")	
names(RIOPS_snow) <- "Snow_cm"     #rename
crs(RIOPS_snow) <- RIOPS_pol_ster #define coordinates system
RIOPS_snow_pr <- projectRaster(RIOPS_snow, res=Out_resolution, crs=Can_pol_ster, method="bilinear") #project to grid

# Loop through GIOPS/RIOPS MRD files in the working directory
RIOPS_files <- dir(pattern = '*level')

for(g in 1:length(RIOPS_files)){
	RIOPS_file <- RIOPS_files[g]

	# Get RIOPS file name
	# RIOPS_base <- sub("^([^.]*).*", "\\1", basename(RIOPS_file)) 
	RIOPS_base <- sub("^SDV_([^.]*)_00000000_0_meters_above_sea_level", "\\1", basename(RIOPS_file))   # snow
	# RIOPS_base <- sub("^GE_([^.]*)_00000000_0_meters_above_sea_level", "\\1", basename(RIOPS_file))  # ice
	# RIOPS_base <- sub("^TMI_([^.]*)_00000000_0_meters_above_sea_level", "\\1", basename(RIOPS_file))  # surf temps
	
	G_txt <- read.table(RIOPS_file)
	names(G_txt) <- c("SDV","lat","lon","Use")	# snow thickness
	# names(G_txt) <- c("GE","lat","lon","Use")		# ice thickness
	# names(G_txt) <- c("TMI","lat","lon","Use")		# ice/snow surface temperature

	# Crop to north of 45 degrees and omit land values (zeros) => ocean and >45N
	G_North <- subset(G_txt, G_txt$lat >=  45 & G_txt$Use ==  1)

	# Get spatial coordinates
	sp_points <- cbind(G_North$lon, G_North$lat) 

	# Create spatial data frame
	G_North_ll <- SpatialPointsDataFrame(G_North, coords=sp_points, proj4string=WGS84) 

	# Project to Can Pol Ster
	G_North_PS <- spTransform(G_North_ll, Can_pol_ster) 

		# writeOGR(G_North_PS, wdir, "GIOPS_CanPolSter_20160302_TMI", driver="ESRI Shapefile",
				# morphToESRI=TRUE, overwrite_layer=TRUE)	

	# Empty Navigator RIOPS file to rasterise onto
	RIOPS_snow_pr_empty <- RIOPS_snow_pr * NA

	G_North_PS_rf <- rasterize(G_North_PS, RIOPS_snow_pr_empty, G_North_PS$SDV, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	# snow
	# G_North_PS_rf <- rasterize(G_North_PS, RIOPS_snow_pr_empty, G_North_PS$GE, fun=mean, background=NA,
		# mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)		# ice
	# G_North_PS_rf <- rasterize(G_North_PS, RIOPS_snow_pr_empty, G_North_PS$TMI, fun=mean, background=NA,
		# mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)		# surface temps
		
	writeRaster(G_North_PS_rf, filename=paste("GIOPS_snow_",RIOPS_base,"_",Out_res_text,".tif",sep=""), format="GTiff", overwrite=TRUE)
	# writeRaster(G_North_PS_rf, filename=paste("GIOPS_ice_",RIOPS_base,"_",Out_res_text,".tif",sep=""), format="GTiff", overwrite=TRUE)
	# writeRaster(G_North_PS_rf, filename=paste("GIOPS_Tsurf_",RIOPS_base,"_",Out_res_text,".tif",sep=""), format="GTiff", overwrite=TRUE)

}
# End
