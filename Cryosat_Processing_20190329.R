# R script to read in CryoSat-2 data from CPOM,
# add RIOPS/GIOPS snow and ice thickness fields,
# add an OSI-SAF ice type field, 
# and calculate three new ice thickness products:
# 	Ice thickness with a Delta-s correction for FYI as a function of snow thickness, based on Nandan et al. (2018)
# 	Ice thickness with the Warren et al. (1999) snow estimate replaced by RIOPS, or GIOPS, snow thickness estimates
# 	Ice thickness with the Delta-s correction using RIOPS, or GIOPS, snow thickness estimates
# The average for each field is calculated for a selected grid size (based on the resolution of the input RIOPS/GIOPS file),
# which is set by the processing in the R script "RIOPS_nc_processing_20190331.R".
# It also reads in CryoVEx data, which are also averaged over the grid cell.
#
# The output files are: 
#	Daily files (if indicated) for the Individual Passes (IP) points, aggregated to the grid
#		these are tab-delimited text files that have had their projection
#		removed and long/lat coordinates have been added (uses the WGS84 datum), e.g.:
#		20140317_C2_IP_ge7n_12p5km.txt 
#		where ge7n stands for >= 7 samples, and 12p5km = 12.5 km
#		file naming is automatic based on input files and values set in Input
#	A shapefile of IP aggregated points (based on grid), e.g.:
#		C2_IP_12p5km_20140325_20140401_IP_only.shp
# 	A shapefile of all aggregated points (based on grid), e.g.:
# 		C2_IP_12p5km_20140315_20140401.shp
#	A tab-delimited text file of all aggregated points (based on grid), 
#		this file has had its projection removed and long/lat coordinates
# 		have been added (uses the WGS84 datum), e.g.:
#		C2_IP_12p5km_20140315_20140401.txt
#
####################################################################################

library("sp")
library("raster")  
library("ncdf4")
library("rgdal") 
library("rgeos") 
library("maptools")
library("spatialEco")

# INPUT <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Set the working directory 
# wdir <- "C:/Users/Torsten/Documents/CIS_DA/Cryosat/TorstenG-MarApr2014/1km/CryoVEx"
wdir <- "C:/Users/Torsten/Documents/CIS_DA/Cryosat/CPOM_Experimental/Ind_pass"

setwd(wdir)

# Get list of CPOM CryoSat-2 retrievals; these must all be in the wdir
C2ip_files <- dir(pattern = "\\.thk$") 	# individual passes
C2fb_files <- dir(pattern = "\\.fb$") 	# radar freeboard before assumptions
C2_files <- dir(pattern = "\\.map$") 	# 2-day 1 km CPOM points

# Set analysis projection
Can_pol_ster <- CRS("+proj=stere +lat_0=90 +lat_ts=90 +lon_0=-100 +k=0.994 +x_0=2000000 +y_0=2000000 +datum=WGS84 +units=m +no_defs") 

# Set CryoVEx data; data are already combined and projected to Can_pol_ster
Cvex_dir <- "C:/Users/Torsten/Documents/CIS_DA/Cryosat/CryoVex/2014/York/D4AEMThicknessData" #switch for ASIRAS
Cvex_file <- paste(Cvex_dir,"/","201403_18_19_21_allfinal_CanPolStr.shp",sep="")
# Cvex_file <- paste(Cvex_dir,"/","201403_26_29_30_01_allfinal_CanPolStr.shp",sep="")
# Cvex_file <- paste(Cvex_dir,"/","201403_18to01_allfinal_CanPolStr.shp",sep="")
#
# For Cvex comparisons, use 2-day CPOM data from 15,16,18,20,   25,27,29,31 (these are for date and next day)
# and use Ind pass and Fr data from 15,16,17,18,19,20,21,  25,26,27,28,29,30,31,01

Use_Cvex <- 0 # 1=yes, 0=no

# Set OSI-SAF ice type dir and file (only used for 2-day 1 km products)
# From http://osisaf.met.no/p/ice/#type ftp://osisaf.met.no/archive/ice/type/2014/03/
# OSI-SAF projection for netCDF files
osi_pol_ster <- CRS("+proj=stere +a=6378273 +b=6356889.44891 +lat_0=90 +lat_ts=70 +lon_0=-45")
osi_dir <- "C:/Users/Torsten/Documents/CIS_DA/Cryosat/OSI_SAF_icetype/"
# osi_file <- "ice_type_nh_polstere-100_multi_201612261200.nc"
osi_file <- "ice_type_nh_polstere-100_multi_201403211200.nc"

# Set RIOPS/GIOPS file (always include these - sometimes they are just used for the grid)
# RIOPS daily mean data; already projected to CanPolSter during cropping and rasterising in script "RIOPS_nc_processing_20190331.R"
# The RIOPS data were originally 2D netCDF files from http://navigator.oceansdata.ca/thredds/catalog.html and 
# had a projection of CRS("+proj=stere +a=6371229.0 +b=6371229.0 +lat_0=90 +lat_ts=90 +lon_0=-100 +k=0.9330124 +x_0=4245000.0 +y_0=5295000.0 +units=m +no_defs") 
# Additional documentation at http://dd.weather.gc.ca/model_giops/doc/README_GIOPS.txt
RIOPS_dir <- "C:/Users/Torsten/Documents/CIS_DA/RIOPS/2016/"
# Select resolution of analysis grid using the relevant files below
# RIOPS_snow_file <- "riops_201612_2D_12p5km_snow.tif"
# RIOPS_ice_file <- "riops_201612_2D_12p5km_ice.tif"
# RIOPS_snow_file <- "riops_201612_2D15km_snow.tif"
# RIOPS_ice_file <- "riops_201612_2D15km_ice.tif"
RIOPS_snow_file <- "riops_201612_2D40km_snow.tif"
RIOPS_ice_file <- "riops_201612_2D40km_ice.tif"
#
# GIOPS data have also already been projected to CanPolSter during cropping and rasterising in script "GIOPS_crop_proj_20190331.R"
# Original GIOPS datum is CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
GIOPS_Snow_dir <- "C:/Users/Torsten/Documents/CIS_DA/GIOPS/2014/Snow/"
GIOPS_Snow_file <- "GIOPS_snow_20140319_12p5km.tif"	# projected to CanPolSter during cropping and rasterising
# GIOPS_Snow_file <- "GIOPS_snow_20140326_12p5km.tif"	# projected to CanPolSter during cropping and rasterising
# GIOPS_Snow_file <- "GIOPS_snow_20140319_5km.tif"
# GIOPS_Snow_file <- "GIOPS_snow_20140326_5km.tif"
#
GIOPS_Ice_dir <- "C:/Users/Torsten/Documents/CIS_DA/GIOPS/2014/Ice/"
GIOPS_Ice_file <- "GIOPS_ice_20140319_12p5km.tif"
# GIOPS_Ice_file <- "GIOPS_ice_20140326_12p5km.tif"
# GIOPS_Ice_file <- "GIOPS_ice_20140319_5km.tif"
# GIOPS_Ice_file <- "GIOPS_ice_20140326_5km.tif"

# Set which to use, RIOPS or GIOPS
RIOPS_GIOPS <- 0	# 1=RIOPS, 0=GIOPS

# Set point spacing text as per RIOPS or GIOPS files - used in filename for .txt output
spacing <- "15km"

# Set if daily output files should be produced
dailies <- 0	# 1=yes, 0=no

# For the daily files, set the number of samples that must satisfy >= # ;
# e.g. number of individual retrievals within each grid
# cell (set by RIOPS/GIOPS file spacing above)
NumSamples <- 1

# Set MYI snow thickness for month in question (in m) - see table below (values from Andy Ridout at CPOM)
#Hs_Tilling <- 0.2332	# December
 Hs_Tilling <- 0.3160	# March

#Ps_Tilling <- 288.2128/1000	# December
Ps_Tilling <- 315.4449/1000		# March

# Monthly Arctic-wide snow thickness (m) and density (kg/m^3) used in CPOM retrievals
# 01    0.2668  287.5529
# 02    0.2888  305.2287
# 03    0.3160  315.4449
# 04    0.3132  309.8278
# 05    0.3259  320.9877
# 06    0.2771  332.0151
# 07    0.0338  358.8291
# 08    0.0251  238.2517
# 09    0.0841  235.9338
# 10    0.1716  268.1507
# 11    0.2136  279.7318
# 12    0.2332  288.2128

# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> End of input

# CPOM 1 km retrievals
if(length(C2_files) != 0) {
	
	# Dummy line for collector
	CPOMs <- c(NA, NA, NA, NA, NA, NA, NA)
	# Header
	names(CPOMs) <- c("Lat", "Long", "SeaIceThick",
					"StDev", "n_samples", "CoG", "C2date")
	for(f in 1:length(C2_files)){
		file1 <- C2_files[f]
		
		CPOM <- read.table(file1, header=FALSE)
		names(CPOM) <- c("Lat", "Long", "SeaIceThick",
						"StDev", "n_samples", "CoG")
						
		# Get file name (e.g. thk_20140320_1km.map)
		CPOM_base <- sub("^([^.]*).*", "\\1", basename(file1)) 
		CPOM$C2date <- sub("^thk_(.*)_1km$", "\\1", CPOM_base, perl = TRUE)
	
		CPOMs <- rbind(CPOMs, CPOM)
	}
	
	# Delete dummy line
	CPOMs <- CPOMs[-1,]		

	CPOMs <- CPOMs[!duplicated(CPOMs[,c('Lat','Long')]),]	# removes duplicates which occur if consecutive day files are used

	CPOMs$C2ID <- row.names(CPOMs) # assign IDs
	
	CPOMs <- subset(CPOMs, CPOMs$n_samples >= 15)	# only use 1 km points if they contain 15 samples (15 used due to extrapolation uncertainty)

	# Put CS-2 points into a spatial data frame and projects to Can_pol_ster
	# Data are originally in long/lat using WGS84
	cxy <- CPOMs[,c(2,1)]
	CPOMs_ll <- SpatialPointsDataFrame(coords = cxy, data = CPOMs,
				proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
	CPOMs_CanPolSter <- spTransform(CPOMs_ll, Can_pol_ster)
		# writeOGR(CPOMs_CanPolSter, wdir, "CPOMs_CanPolSter_20140315_0321", driver="ESRI Shapefile",
				# morphToESRI=TRUE, overwrite_layer=TRUE)
	
	C2_passes <- CPOMs_CanPolSter
}

# CPOM Individual Passes retrievals
if(length(C2ip_files) != 0) {
	
	Ind_passes <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA)
	names(Ind_passes) <- c("Lat", "Long", "Valid", "ipSeaIceThick", "ipIceType",
						"Conf_IceType", "SeaIce_Conc", "ipC2date", "ipC2time") 
						
	for(ip in 1:length(C2ip_files)){
		file_ip <- C2ip_files[ip]

		Ind_pass <- read.table(file_ip, header=FALSE)
		names(Ind_pass) <- c("PacketID", "Block", "Time", "Lat", "Long",
								"Valid", "ipSeaIceThick", "ipIceType", "Conf_IceType",
								"SeaIce_Conc", "CoGq")
		Ind_pass <- subset(Ind_pass, select = Lat:SeaIce_Conc)				
		Ind_pass <- subset(Ind_pass, Valid == 1) #may switch off
		
		if(nrow(Ind_pass) != 0) {
			# Get file name (e.g. cry_NO_20140317T000913.thk)
			CPOMip_base <- sub("^([^.]*).*", "\\1", basename(file_ip)) 
			Ind_pass$ipC2date <- sub("^cry_NO_(.*)T.*$", "\\1", CPOMip_base, perl = TRUE)
			Ind_pass$ipC2time <- sub("^cry_NO_.*T(.*)$", "\\1", CPOMip_base, perl = TRUE)

			Ind_passes <- rbind(Ind_passes, Ind_pass)
		}
			# Format of .thk files
			# --------------------
			# Col 1 : Source Packet ID
			# Col 2 : Block number (0->19)
			# Col 3 : Time
			# Col 4 : Latitude
			# Col 5 : Longitude
			# Col 6 : Valid (1=Yes, 0=No) [See note 1]
			# Col 7 : Sea Ice thickness (m)
			# Col 8 : Sea Ice Type [See note 2]
				# (0 -> Unset)
				# (1 -> no ice or very open ice)
				# (2 -> relatively young ice)
				# (3 -> ice that survived a summer melt)
				# (4 -> ambiguous ice type)
			# Col 9 : Confidence in Sea Ice Type
				# (0 -> Not Processed)
				# (1 -> Computation Failed)
				# (2 -> Use With Care)
				# (3 -> Acceptable)
				# (4 -> Good)
				# (5 -> Excellent)
			# Col 10 : Sea Ice Conc. (%) [Required for volume]
			# Col 11 : Spare

			# NOTE 1: A valid freeboard results in an invalid thickness when the
			# freeboard corrected for the snow layer is outside a certain
			# range, or the sea ice type is 0, 1 or 4.

			# NOTE 2: Sea ice concentration and sea ice type come from different sources.
			# This can result in a valid freeboard having a sea ice type of 0 when
			# near the coast, very occasionally 1, and 4 near the FYI/MYI boundary.
			# As noted above the thickness record will be invalid.
	}
	
	Ind_passes <- Ind_passes[-1,]
	Ind_passes$ipID <- row.names(Ind_passes)
	
	ipxy <- Ind_passes[,c(2,1)]
	C2ip_ll <- SpatialPointsDataFrame(coords = ipxy, data = Ind_passes,
							   proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
	C2ip_pr <- spTransform(C2ip_ll, Can_pol_ster)
		# writeOGR(C2ip_pr, wdir, "C2ip_CanPolSter_20140315_0321", driver="ESRI Shapefile",
				# morphToESRI=TRUE, overwrite_layer=TRUE) #for testing without for loop
		# writeOGR(gg, wdir, "C2ip_CanPolSter_20161231", driver="ESRI Shapefile",
				# morphToESRI=TRUE, overwrite_layer=TRUE)

}

# CPOM Radar freeboard files (before assumptions added)
if(length(C2fb_files) != 0) {
	
	Radar_fbs <- c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)
	names(Radar_fbs) <- c("Lat", "Long", "R_freeboard", "SLA", "Backscatter",
							"SeaIce_Conc", "FbIceType", "Conf_IceType", "FbSnow_flag",
							"FbSnow_depth", "FbSnow_Dens", "fbC2date", "fbC2time")	
							
	for(frbd in 1:length(C2fb_files)){
		file_fr <- C2fb_files[frbd]

		Radar_fb <- read.table(file_fr, header=FALSE)
		names(Radar_fb) <- c("PacketID", "Block", "Time", "Lat", "Long",
								"R_freeboard", "SLA", "Backscatter", "FbSeaIce_Conc",
								"FbIceType", "Conf_IceType", "FbSnow_flag", "FbSnow_depth",
								"FbSnow_dens", "CoGq")	
								
		Radar_fb <- subset(Radar_fb, select = Lat:FbSnow_dens)
		Radar_fb <- subset(Radar_fb, FbSnow_flag < 3)
		
		if(nrow(Radar_fb) != 0) {
			# Get file name (e.g. cry_NO_20140317T000913.fb)
			CPOMfb_base <- sub("^([^.]*).*", "\\1", basename(file_fr)) 
			Radar_fb$fbC2date <- sub("^cry_NO_(.*)T.*$", "\\1", CPOMfb_base, perl = TRUE)
			Radar_fb$fbC2time <- sub("^cry_NO_.*T(.*)$", "\\1", CPOMfb_base, perl = TRUE)
			
			# Radar_fb <- subset(Radar_fb, Conf_IceType > 3)
			Radar_fbs <- rbind(Radar_fbs, Radar_fb)
		} 
			# Format of .fb Files
			# -------------------
			# Column 1 : Source Packet ID
			# Column 2 : Block Number (0->19)
			# Column 3 : Time
			# Column 4 : Latitude
			# Column 5 : Longitude
			# Column 6 : Freeboard [this is radar freeboard]
			# Column 7 : SLA interpolated at this location
			# Column 8 : Backscatter #not in this file
			# Column 9 : Sea Ice Concentration (%)
			# Column 10 : Sea Ice Type (0:Unset, 1:Open Water, 2:FYI, 3:MYI, 4:Ambigous)
			# Column 11 : Confidence in Sea Ice Type (0:Not Processed, 1:Computation Failed,
							# 2:Use With Care, 3:Acceptable, 4:Good, 5:Excellent)
			# Column 12 : Snow flag (0:OK, 1:Not in Polygon, 2: Lat<60N, 3:Daft Values)
			# Column 13 : Snow depth (Adjusted for ice type) (m)
			# Column 14 : Snow Density (kg m-3)
			# Column 15 : Spare

			# NOTE: Sea ice concentration and sea ice type come from different sources.
			# This can result in a valid freeboard having a sea ice type of 0 when
			# near the coast, very occasionally 1, and 4 near the FYI/MYI boundary.
	}
	
	Radar_fbs <- Radar_fbs[-1,]
	Radar_fbs$fbID <- row.names(Radar_fbs)
		
	fbxy <- Radar_fbs[,c(2,1)]
	C2fb_ll <- SpatialPointsDataFrame(coords = fbxy, data = Radar_fbs,
							   proj4string = crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
	C2fb_pr <- spTransform(C2fb_ll, Can_pol_ster)
	# writeOGR(C2fb_pr, wdir, "C2fb_CanPolSter_20140315_0321", driver="ESRI Shapefile",
			# morphToESRI=TRUE, overwrite_layer=TRUE)
		
}

# Get modelled snow and ice values
RIOPS_snow <- raster(paste(RIOPS_dir,RIOPS_snow_file,sep=""))	 #use a dummy riops file
names(RIOPS_snow) <- "Snow_cm"
crs(RIOPS_snow) <- Can_pol_ster
RIOPS_snow_pr <- RIOPS_snow

RIOPS_ice <- raster(paste(RIOPS_dir,RIOPS_ice_file,sep=""))	
names(RIOPS_ice) <- "Ice_m"
crs(RIOPS_ice) <- Can_pol_ster # affirm projection
RIOPS_ice_pr <- RIOPS_ice

GIOPS_snow <- raster(paste(GIOPS_Snow_dir,GIOPS_Snow_file,sep=""),varname="SDV")
crs(GIOPS_snow) <- Can_pol_ster # affirm projection
names(GIOPS_snow) <- "Snow_cm"

GIOPS_ice <- raster(paste(GIOPS_Ice_dir,GIOPS_Ice_file,sep=""),varname="GE")
crs(GIOPS_ice) <- Can_pol_ster # affirm projection
names(GIOPS_ice) <- "Ice_m"


# Read in CryoVEx files; already combined (as .shp) and projected on can pol ster
#switch for ASIRAS fields; now for York
if(Use_Cvex == 1) {
	Cvex_shp <- shapefile(Cvex_file)
	crs(Cvex_shp) <- Can_pol_ster	# affirm projection
	names(Cvex_shp) <- c("Year","Month","Day","GPSseconds","Fid",
							"Lat","Long","Distance","Thickness","Height",
							"layer","path")
	Cvex_shp$CvexID <- seq.int(nrow(Cvex_shp))	# assign IDs
		# writeOGR(Cvex_shp, wdir, "Cvex_shp_crs_set", driver="ESRI Shapefile",
				# morphToESRI=TRUE, overwrite_layer=TRUE)
	Cvex_shp$CvexDate <- as.numeric(paste(Cvex_shp$Year,"0",Cvex_shp$Month,Cvex_shp$Day,sep=""))
}

# For 2-day 1 km product add OSI-SAF ice type column
if(length(C2_files) != 0) {
	# OSI-SAF ice type (in their projection)
	osi_saf <- raster(paste(osi_dir,osi_file,sep=""),varname="ice_type")
	crs(osi_saf) <- osi_pol_ster
	names(osi_saf) <- "osi_saf_ice_type"
	extent(osi_saf) <- c(-3850000, 3750000, -5350000, 5850000)
	osi_saf_pr <- projectRaster(osi_saf, crs=Can_pol_ster, method="ngb")

	C2_passes_sp_osi <- extract(osi_saf_pr, C2_passes, sp=TRUE)
	C2_passes_sp_pr <- C2_passes_sp_osi 
		# writeOGR(C2_passes_sp_pr, wdir, "C2_passes_sp_pr_201403_15to01", driver="ESRI Shapefile",
				# morphToESRI=TRUE, overwrite_layer=TRUE)
			
	# Project OSI-SAF data for viewing; optional
	# osi_saf_dis <- disaggregate(osi_saf, fact=2)
	# osi_saf_pr <- projectRaster(osi_saf_dis, crs=Can_pol_ster, method="ngb")
	# writeRaster(osi_saf_pr, filename="osi_saf_pr_201403211200_fact2.tif", format="GTiff", overwrite=TRUE)
}

# Combine retrieval points with RIOPS or GIOPS data 
if(RIOPS_GIOPS == 1) {	# RIOPS
	if(length(C2_files) != 0) {
		C2_passes_sp_pr <- extract(RIOPS_ice, C2_passes_sp_pr, sp=TRUE)
		C2_passes_sp_pr <- extract(RIOPS_snow, C2_passes_sp_pr, sp=TRUE)
	}
	if(length(C2ip_files) != 0) {
		C2ip_pr <- extract(RIOPS_ice, C2ip_pr, sp=TRUE)
		C2ip_pr <- extract(RIOPS_snow, C2ip_pr, sp=TRUE)
	}
	if(length(C2fb_files) != 0) {
		C2fb_pr <- extract(RIOPS_ice, C2fb_pr, sp=TRUE)
		C2fb_pr <- extract(RIOPS_snow, C2fb_pr, sp=TRUE)
	}
	if(Use_Cvex == 1) {
		Cvex_shp <- extract(RIOPS_ice, Cvex_shp, sp=TRUE)
		Cvex_shp <- extract(RIOPS_snow, Cvex_shp, sp=TRUE)
	}
} else {  # GIOPS
	if(length(C2_files) != 0) {
		C2_passes_sp_pr <- extract(GIOPS_ice, C2_passes_sp_pr, sp=TRUE)
		C2_passes_sp_pr <- extract(GIOPS_snow, C2_passes_sp_pr, sp=TRUE)
	}
	if(length(C2ip_files) != 0) {
		C2ip_pr <- extract(GIOPS_ice, C2ip_pr, sp=TRUE)
		C2ip_pr <- extract(GIOPS_snow, C2ip_pr, sp=TRUE)
	}
	if(length(C2fb_files) != 0) {
		C2fb_pr <- extract(GIOPS_ice, C2fb_pr, sp=TRUE)
		C2fb_pr <- extract(GIOPS_snow, C2fb_pr, sp=TRUE)
	}
	if(Use_Cvex == 1) {
		Cvex_shp <- extract(GIOPS_ice, Cvex_shp, sp=TRUE)
		Cvex_shp <- extract(GIOPS_snow, Cvex_shp, sp=TRUE)
	}
}

# Separate MYI and FYI	
if(length(C2_files) != 0) {
	C2_passes_sp_pr_MYI <- subset(C2_passes_sp_pr, C2_passes_sp_pr$osi_saf_ice_type == 3)
	C2_passes_sp_pr_FYI <- subset(C2_passes_sp_pr, C2_passes_sp_pr$osi_saf_ice_type == 2)
}
if(length(C2ip_files) != 0) {
	C2ip_pr_MYI <- subset(C2ip_pr, C2ip_pr$ipIceType == 3) #from osisaf
	C2ip_pr_FYI <- subset(C2ip_pr, C2ip_pr$ipIceType == 2)
}
if(length(C2fb_files) != 0) {
	C2fb_pr_MYI <- subset(C2fb_pr, C2fb_pr$FbIceType == 3)
	C2fb_pr_FYI <- subset(C2fb_pr, C2fb_pr$FbIceType == 2)
}

# Adjustments for MYI ---------------------
Pw <- 1.025
Ps <- Ps_Tilling
Pice <- 0.882
Hs_MYI <- Hs_Tilling
#
# Reverse CPOM assumptions to find original radar freeboard
if(length(C2_files) != 0) {
	C2_passes_sp_pr_MYI$Fi <- ((C2_passes_sp_pr_MYI$SeaIceThick * (Pw-Pice) - Hs_MYI*Ps) #check for validity
							/ Pw) - 0.25*Hs_MYI #adjustment for lower speed in snow , check with vishnu
}
if(length(C2ip_files) != 0) {
	C2ip_pr_MYI$ipFi <- ((C2ip_pr_MYI$ipSeaIceThick * (Pw-Pice) - Hs_MYI*Ps)
							/ Pw) - 0.25*Hs_MYI
}
#							
# Calculate ice thickness based on new assumptions: Replacing W99 with riops/giops
if(length(C2_files) != 0) {
	Hs_RIOPS_MYI_cm <- C2_passes_sp_pr_MYI$Snow_cm
	Hs_RIOPS_MYI_m <- Hs_RIOPS_MYI_cm/100
	C2_passes_sp_pr_MYI$Ti_RIOPS <- ((C2_passes_sp_pr_MYI$Fi + 0.25*Hs_RIOPS_MYI_m) * Pw
								+ Hs_RIOPS_MYI_m*Ps)/(Pw-Pice)
	C2_passes_sp_pr_MYI$Ti_DeltaS <- NA
	C2_passes_sp_pr_MYI$Ti_RIOPS_DeltaS <- NA 
}
if(length(C2ip_files) != 0) {
	ipHs_RIOPS_MYI_cm <- C2ip_pr_MYI$Snow_cm #using reversed engineered fb
	ipHs_RIOPS_MYI_m <- ipHs_RIOPS_MYI_cm/100
	C2ip_pr_MYI$ipTi_RIOPS <- ((C2ip_pr_MYI$ipFi + 0.25*ipHs_RIOPS_MYI_m) * Pw
								+ ipHs_RIOPS_MYI_m*Ps)/(Pw-Pice) #hydrostatic using riops snow to calculate ice thickness
	C2ip_pr_MYI$ipTi_DeltaS <- NA
	C2ip_pr_MYI$ipTi_RIOPS_DeltaS <- NA #no salinity adjustment for myi

}
if(length(C2fb_files) != 0) {
	fbHs_RIOPS_MYI_cm <- C2fb_pr_MYI$Snow_cm #using cpom fb
	fbHs_RIOPS_MYI_m <- fbHs_RIOPS_MYI_cm/100
	C2fb_pr_MYI$fbTi <- ((C2fb_pr_MYI$R_freeboard + 0.25*Hs_MYI) * Pw
								+ Hs_MYI*Ps)/(Pw-Pice)
	C2fb_pr_MYI$fbTi_RIOPS <- ((C2fb_pr_MYI$R_freeboard + 0.25*fbHs_RIOPS_MYI_m) * Pw
								+ fbHs_RIOPS_MYI_m*Ps)/(Pw-Pice)
	C2fb_pr_MYI$fbTi_DeltaS <- NA
	C2fb_pr_MYI$fbTi_RIOPS_DeltaS <- NA
}

								
# Adjustments for FYI ---------------------------
Pw <- 1.025
Ps <- Ps_Tilling
Pice <- 0.917
Hs_FYI <- Hs_Tilling/2
Hs_cm <- Hs_FYI*100
DeltaS <- 1.4022229 + 0.9114689*Hs_cm - 0.0437265*Hs_cm^2 + 0.00061*Hs_cm^3	# Nandan et al. (2017)
#
# Reverse assumptions to find original radar freeboard
if(length(C2_files) != 0) {
	C2_passes_sp_pr_FYI$Fi <- ((C2_passes_sp_pr_FYI$SeaIceThick * (Pw-Pice) - Hs_FYI*Ps)
								/ Pw) - 0.25*Hs_FYI
}
if(length(C2ip_files) != 0) {
	C2ip_pr_FYI$ipFi <- ((C2ip_pr_FYI$ipSeaIceThick * (Pw-Pice) - Hs_FYI*Ps)
							/ Pw) - 0.25*Hs_FYI
}
#
# Calculate ice thickness based on new assumptions
if(length(C2_files) != 0) {
	C2_passes_sp_pr_FYI$Ti_DeltaS <- (((C2_passes_sp_pr_FYI$Fi + 0.25*Hs_FYI) - DeltaS/100)
									* Pw + Hs_FYI*Ps)/(Pw-Pice)
	Hs_RIOPS_cm <- C2_passes_sp_pr_FYI$Snow_cm
	Hs_RIOPS_m <- Hs_RIOPS_cm/100
	Hs_RIOPS_DeltaS <- 1.4022229 + 0.9114689*Hs_RIOPS_cm - 0.0437265*Hs_RIOPS_cm^2 + 0.00061*Hs_RIOPS_cm^3	# Nandan et al. (2017)
	C2_passes_sp_pr_FYI$Ti_RIOPS <- ((C2_passes_sp_pr_FYI$Fi + 0.25*Hs_RIOPS_m)
									* Pw + Hs_RIOPS_m*Ps)/(Pw-Pice)
	C2_passes_sp_pr_FYI$Ti_RIOPS_DeltaS <- (((C2_passes_sp_pr_FYI$Fi + 0.25*Hs_RIOPS_m)
										- Hs_RIOPS_DeltaS/100)*Pw + Hs_RIOPS_m*Ps)/(Pw-Pice)
	# Combine MYI and FYI following adjustments
	C2_passes_sp_pr_comb <- rbind(C2_passes_sp_pr_FYI,C2_passes_sp_pr_MYI)
}
if(length(C2ip_files) != 0) {
	C2ip_pr_FYI$ipTi_DeltaS <- (((C2ip_pr_FYI$ipFi + 0.25*Hs_FYI) - DeltaS/100)
								* Pw + Hs_FYI*Ps)/(Pw-Pice)
	ipHs_RIOPS_FYI_cm <- C2ip_pr_FYI$Snow_cm
	ipHs_RIOPS_FYI_m <- ipHs_RIOPS_FYI_cm/100
	ipHs_RIOPS_DeltaS <- 1.4022229 + 0.9114689*ipHs_RIOPS_FYI_cm - 0.0437265*ipHs_RIOPS_FYI_cm^2 + 0.00061*ipHs_RIOPS_FYI_cm^3	# Nandan et al. (2017)
	C2ip_pr_FYI$ipTi_RIOPS <- ((C2ip_pr_FYI$ipFi + 0.25*ipHs_RIOPS_FYI_m)
								* Pw + ipHs_RIOPS_FYI_m*Ps)/(Pw-Pice)
	C2ip_pr_FYI$ipTi_RIOPS_DeltaS <- (((C2ip_pr_FYI$ipFi + 0.25*ipHs_RIOPS_FYI_m)
										- ipHs_RIOPS_DeltaS/100)*Pw + ipHs_RIOPS_FYI_m*Ps)/(Pw-Pice)
	C2ip_pr_comb <- rbind(C2ip_pr_FYI,C2ip_pr_MYI)
	# plot(C2ip_pr_comb$ipSeaIceThick,C2ip_pr_comb$ipTi_RIOPS_DeltaS)	
}
if(length(C2fb_files) != 0) {
	C2fb_pr_FYI$fbTi <- (((C2fb_pr_FYI$R_freeboard + 0.25*Hs_FYI))
										* Pw + Hs_FYI*Ps)/(Pw-Pice)
	C2fb_pr_FYI$fbTi_DeltaS <- (((C2fb_pr_FYI$R_freeboard + 0.25*Hs_FYI) - DeltaS/100)
										* Pw + Hs_FYI*Ps)/(Pw-Pice)
	fbHs_RIOPS_FYI_cm <- C2fb_pr_FYI$Snow_cm
	fbHs_RIOPS_FYI_m <- fbHs_RIOPS_FYI_cm/100
	fbHs_RIOPS_DeltaS <- 1.4022229 + 0.9114689*fbHs_RIOPS_FYI_cm - 0.0437265*fbHs_RIOPS_FYI_cm^2 + 0.00061*fbHs_RIOPS_FYI_cm^3	# Nandan et al. (2017)
	C2fb_pr_FYI$fbTi_RIOPS <- ((C2fb_pr_FYI$R_freeboard + 0.25*fbHs_RIOPS_FYI_m)
								* Pw + fbHs_RIOPS_FYI_m*Ps)/(Pw-Pice)
	C2fb_pr_FYI$fbTi_RIOPS_DeltaS <- (((C2fb_pr_FYI$R_freeboard + 0.25*fbHs_RIOPS_FYI_m)
										- fbHs_RIOPS_DeltaS/100)*Pw + fbHs_RIOPS_FYI_m*Ps)/(Pw-Pice)
	C2fb_pr_comb <- rbind(C2fb_pr_FYI,C2fb_pr_MYI)
}

# Average each field based on the RIOPS/GIOPS spacing, this includes
# the date and time fields, so these are a function of the date range of
# the files in the wdir

# Create empty grids based on RIOPS/GIOPS spacing
	if(RIOPS_GIOPS == 1) {	# RIOPS
		RIOPS_empty_grid <- RIOPS_snow_pr * NA
	} else {  				#GIOPS
		RIOPS_empty_grid <- GIOPS_snow * NA
	}

Grid_list <- c("C2_date_grid","C2_ice_grid","RIOPS_ice_grid","OSISAF_icetype_grid",
					"Fr_reverse_grid","Ti_DeltaS_grid","Ti_RIOPS_grid",
					"Ti_RIOPS_DeltaS_grid","RIOPS_snow_grid","C2_count_grid",
					"ipC2_date_grid","ipC2_time_grid","ipC2_ice_grid","ipRIOPS_ice_grid",
					"ipOSISAF_icetype_grid","ipFr_reverse_grid","ipTi_DeltaS_grid",
					"ipTi_RIOPS_grid","ipTi_RIOPS_DeltaS_grid","ipRIOPS_snow_grid",
					"ipIceType_grid","ipConf_IceType_grid","ipSIC_grid","ipC2_count_grid",
					"fbC2_date_grid","fbC2_time_grid","fbC2_ice_grid",
					"fbRIOPS_ice_grid","fbOSISAF_icetype_grid",
					"fbFb_grid","fbTi_DeltaS_grid","fbTi_RIOPS_grid",
					"fbTi_RIOPS_DeltaS_grid","fbRIOPS_snow_grid","fbC2_count_grid",
					"Cvex_date_grid","Cvex_ice_grid","CvexRIOPS_ice_grid","CvexRIOPS_snow_grid","Cvex_count_grid")	

# Produce daily files if indicated
if(dailies == 1) {
	# For each day extract the mean values for each grid cell (based on RIOPS or GIOPS file spacing)
	# and write to a .txt file in long/lat WGS84
	for(dd in min(C2ip_pr_comb$ipC2date, na.rm=TRUE):max(C2ip_pr_comb$ipC2date, na.rm=TRUE)) {
		daily <- subset(C2ip_pr_comb, C2ip_pr_comb$ipC2date == dd)
		
		# Empty grid for current day
		for(gr in 1:length(Grid_list)) {
			Grid <- paste(as.symbol(Grid_list[gr]),sep="")
			assign(Grid, RIOPS_empty_grid)
		}

		# Individual passes
		if(length(C2ip_files) != 0) {
		ipC2_date_grid <- rasterize(daily, RIOPS_empty_grid, 
				as.numeric(daily$ipC2date), fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
		ipC2_time_grid <- rasterize(daily, RIOPS_empty_grid, 
				as.numeric(daily$ipC2time), fun=median, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
		ipC2_ice_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$ipSeaIceThick, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
		ipRIOPS_ice_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$Ice_m, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
		ipOSISAF_icetype_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$ipIceType, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
		ipFr_reverse_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$ipFi, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
		ipTi_DeltaS_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$ipTi_DeltaS, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
		ipTi_RIOPS_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$ipTi_RIOPS, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)			
		ipTi_RIOPS_DeltaS_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$ipTi_RIOPS_DeltaS, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)		
		ipRIOPS_snow_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$Snow_cm, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
		ipIceType_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$ipIceType, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
		ipConf_IceType_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$Conf_IceType, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)		
		ipSIC_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$SeaIce_Conc, fun=mean, background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
		ipC2_count_grid <- rasterize(daily, RIOPS_empty_grid, 
				daily$ipSeaIceThick, fun='count', background=NA,
				mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
		}

		C2_stack <- stack(ipC2_date_grid, ipC2_time_grid, ipC2_ice_grid, ipRIOPS_ice_grid, 
						ipOSISAF_icetype_grid,
						ipFr_reverse_grid, ipTi_DeltaS_grid, ipTi_RIOPS_grid,
						ipTi_RIOPS_DeltaS_grid, ipRIOPS_snow_grid, 
						ipIceType_grid, ipConf_IceType_grid, ipSIC_grid, ipC2_count_grid)
		names(C2_stack) <- c("ip_date", "ip_time", "ip_Ti", "ipR_Ti", "ipOSAFIt",
						"ipFr_Rev", "ipTiDs", "ipTiRs", "ipTiRsDs", "ipRsnow",
						"ipOS_IT", "ipConfIT", "ipSIC", "ipPntCnt")
						
		# Convert raster cells to points
		GridPoints <- rasterToPoints(C2_stack, spatial=TRUE)
		names(GridPoints) <- c("ip_date", "ip_time", "ip_Ti", "ipR_Ti", "ipOSAFIt",
						"ipFr_Rev", "ipTiDs", "ipTiRs", "ipTiRsDs", "ipRsnow",
						"ipOS_IT", "ipConfIT", "ipSIC", "ipPntCnt")

		# On first day, start accumulator for daily grids
		if(dd == min(C2ip_pr_comb$ipC2date, na.rm=TRUE)) {
			GridPoints_many <- GridPoints
		} else {
			GridPoints_many <- rbind(GridPoints_many, GridPoints)
		}
		
		# Only use points that have a sufficient number of samples; set at Input
		GridPoints <- subset(GridPoints, GridPoints$ipPntCnt >= NumSamples)

		# Subset to only the IP points for operational file output
		IPgrid <- GridPoints  
		IPgrid <- subset(IPgrid, select=c(ip_date:ipPntCnt))

		# For each day, produce an output file of the IP data
		if(nrow(daily) > 0) {
			# Create an unprojected long/lat data set using WGS84 datum
			daily_ll <- spTransform(IPgrid, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
			daily_ll$long <- coordinates(daily_ll)[,1]
			daily_ll$lat <- coordinates(daily_ll)[,2]
			daily_ll <- daily_ll[c("long","lat","ip_date","ip_time",
							"ip_Ti","ipR_Ti","ipOSAFIt","ipFr_Rev",
							"ipTiDs","ipTiRs","ipTiRsDs","ipRsnow",
							"ipOS_IT","ipConfIT","ipSIC","ipPntCnt")] 
			write.table(daily_ll, paste(dd,"_C2_IP_ge",NumSamples,"n_",spacing,".txt", sep=""),
						sep="\t", row.names=FALSE, append = FALSE)
		}
		
	} # End daily loop
} # End dailies

# Rasterize all points onto a grid; the mean/median/counts are for all of the
# points falling in the date/time range of the input files in wdir
#
# Rasterize C2 points onto a grid
# CPOM 1 km 15n
if(length(C2_files) != 0) {
C2_date_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		as.numeric(C2_passes_sp_pr_comb$C2date), fun=median, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
C2_ice_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		C2_passes_sp_pr_comb$SeaIceThick, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
RIOPS_ice_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		C2_passes_sp_pr_comb$Ice_m, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
OSISAF_icetype_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		C2_passes_sp_pr_comb$osi_saf_ice_type, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
Fr_reverse_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		C2_passes_sp_pr_comb$Fi, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
Ti_DeltaS_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		C2_passes_sp_pr_comb$Ti_DeltaS, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
Ti_RIOPS_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		C2_passes_sp_pr_comb$Ti_RIOPS, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)			
Ti_RIOPS_DeltaS_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		C2_passes_sp_pr_comb$Ti_RIOPS_DeltaS, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)		
RIOPS_snow_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		C2_passes_sp_pr_comb$Snow_cm, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
C2_count_grid <- rasterize(C2_passes_sp_pr_comb, RIOPS_empty_grid, 
		C2_passes_sp_pr_comb$SeaIceThick, fun='count', background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
}

# Individual passes
if(length(C2ip_files) != 0) {
ipC2_date_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		as.numeric(C2ip_pr_comb$ipC2date), fun=median, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
ipC2_time_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		as.numeric(C2ip_pr_comb$ipC2time), fun=median, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
ipC2_ice_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$ipSeaIceThick, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
ipRIOPS_ice_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$Ice_m, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
ipOSISAF_icetype_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$ipIceType, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
ipFr_reverse_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$ipFi, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
ipTi_DeltaS_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$ipTi_DeltaS, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
ipTi_RIOPS_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$ipTi_RIOPS, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)			
ipTi_RIOPS_DeltaS_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$ipTi_RIOPS_DeltaS, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)		
ipRIOPS_snow_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$Snow_cm, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
ipIceType_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$ipIceType, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
ipConf_IceType_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$Conf_IceType, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)		
ipSIC_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$SeaIce_Conc, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
ipC2_count_grid <- rasterize(C2ip_pr_comb, RIOPS_empty_grid, 
		C2ip_pr_comb$ipSeaIceThick, fun='count', background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
}

# Radar freeboard
if(length(C2fb_files) != 0) {
fbC2_date_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		as.numeric(C2fb_pr_comb$fbC2date), fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
fbC2_time_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		as.numeric(C2fb_pr_comb$fbC2time), fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
fbC2_ice_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		C2fb_pr_comb$fbTi, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
fbRIOPS_ice_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		C2fb_pr_comb$Ice_m, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
fbOSISAF_icetype_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		C2fb_pr_comb$FbIceType, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
fbFb_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		C2fb_pr_comb$R_freeboard, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
fbTi_DeltaS_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		C2fb_pr_comb$fbTi_DeltaS, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
fbTi_RIOPS_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		C2fb_pr_comb$fbTi_RIOPS, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)			
fbTi_RIOPS_DeltaS_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		C2fb_pr_comb$fbTi_RIOPS_DeltaS, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)		
fbRIOPS_snow_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		C2fb_pr_comb$Snow_cm, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
fbC2_count_grid <- rasterize(C2fb_pr_comb, RIOPS_empty_grid, 
		C2fb_pr_comb$fbTi, fun='count', background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
}
	
# CryoVEx data
if(Use_Cvex == 1) {
Cvex_date_grid <- rasterize(Cvex_shp, RIOPS_empty_grid, 
		Cvex_shp$CvexDate, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
Cvex_ice_grid <- rasterize(Cvex_shp, RIOPS_empty_grid, 
		Cvex_shp$Thickness, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
CvexRIOPS_ice_grid <- rasterize(Cvex_shp, RIOPS_empty_grid, 
		Cvex_shp$Ice_m, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
CvexRIOPS_snow_grid <- rasterize(Cvex_shp, RIOPS_empty_grid, 
		Cvex_shp$Snow_cm, fun=mean, background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)	
Cvex_count_grid <- rasterize(Cvex_shp, RIOPS_empty_grid, 
		Cvex_shp$Thickness, fun='count', background=NA,
		mask=FALSE, update=TRUE, updateValue='all', na.rm=TRUE)
}
	
C2_stack <- stack(C2_date_grid, C2_ice_grid, RIOPS_ice_grid, OSISAF_icetype_grid, 
				Fr_reverse_grid, Ti_DeltaS_grid, Ti_RIOPS_grid, 
				Ti_RIOPS_DeltaS_grid, RIOPS_snow_grid, C2_count_grid,
				ipC2_date_grid, ipC2_time_grid, ipC2_ice_grid, ipRIOPS_ice_grid, 
				ipOSISAF_icetype_grid,
				ipFr_reverse_grid, ipTi_DeltaS_grid, ipTi_RIOPS_grid,
				ipTi_RIOPS_DeltaS_grid, ipRIOPS_snow_grid, 
				ipIceType_grid, ipConf_IceType_grid, ipSIC_grid, ipC2_count_grid,
				fbC2_date_grid, fbC2_time_grid, fbC2_ice_grid, 
				fbRIOPS_ice_grid, fbOSISAF_icetype_grid,
				fbFb_grid, fbTi_DeltaS_grid, fbTi_RIOPS_grid,
				fbTi_RIOPS_DeltaS_grid, fbRIOPS_snow_grid, fbC2_count_grid, 
				Cvex_date_grid, Cvex_ice_grid, CvexRIOPS_ice_grid, CvexRIOPS_snow_grid,
				Cvex_count_grid)
names(C2_stack) <- c("C2date", "C2_Ti", "RIOPS_Ti", "OSAF_It", "C2Fr_rev",
				"C2TiDs", "C2TiRs", "C2TiDsRs", "RIOPS_Hs", "C2PntCnt",
				"ip_date", "ip_time", "ip_Ti", "ipR_Ti", "ipOSAFIt",
				"ipFr_Rev", "ipTiDs", "ipTiRs", "ipTiRsDs", "ipRsnow",
				"ipOS_IT", "ipConfIT", "ipSIC", "ipPntCnt", 
				"fb_date", "fb_time", "fb_Ti", "fbRIOPTi", "fbOSAFIt",
				"fb_Fb", "fbTiDs", "fbTiRs", "fbTiRsDs", "fbRsnow", "fbPntCnt", 
				"CvexDate", "Cvex_TiHs", "CvexRTi", "CvexRsno", "CvexPCnt")
	
GridPoints_all <- rasterToPoints(C2_stack, spatial=TRUE) #place values at centre of each grid cell
names(GridPoints_all) <- c("C2date", "C2_Ti", "RIOPS_Ti", "OSAF_It", "C2Fr_rev",
				"C2TiDs", "C2TiRs", "C2TiDsRs", "RIOPS_Hs", "C2PntCnt",
				"ip_date", "ip_time", "ip_Ti", "ipR_Ti", "ipOSAFIt",
				"ipFr_Rev", "ipTiDs", "ipTiRs", "ipTiRsDs", "ipRsnow",
				"ipOS_IT", "ipConfIT", "ipSIC", "ipPntCnt", 
				"fb_date", "fb_time", "fb_Ti", "fbRIOPTi", "fbOSAFIt",
				"fb_Fb", "fbTiDs", "fbTiRs", "fbTiRsDs", "fbRsnow", "fbPntCnt", 
				"CvexDate", "Cvex_TiHs", "CvexRTi", "CvexRsno", "CvexPCnt")

# Output file naming
accum_out <- paste("C2_IP_",spacing,"_",
			min(GridPoints_all$ip_date, na.rm=TRUE),"_",
			max(GridPoints_all$ip_date, na.rm=TRUE), sep="")

# Write only the IP points to a shapefile
writeOGR(GridPoints_many, wdir, paste(accum_out,"_IP_only",sep=""),
		driver="ESRI Shapefile", morphToESRI=TRUE, overwrite_layer=TRUE)

# Write the accumulated data (all data in wdir) to a shapefile and a tab-delimited text file
writeOGR(GridPoints_all, wdir, accum_out,
		driver="ESRI Shapefile", morphToESRI=TRUE, overwrite_layer=TRUE)
#
accum_ll <- spTransform(GridPoints_all, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
write.table(accum_ll, paste(accum_out,".txt", sep=""), sep="\t",
		row.names=FALSE, append = FALSE)


gc()
print("Done")
# END
