#######################
### Updated: Febuary 13, 2024
#' smaller grid size to fit the SLV
#' 
#' Description: 
#' All Sectors regridded and unit conversions:
#' Regridded from 4km grid cells in Lambert Conical Conformal
#' Unit conversion from absolute emissions (mT hr-1) to flux emissions ( umol s-1 m-2)
#' 
#' File combinations 
### By: Haley Humble
#######################

rm(list=ls()) #Clear environment to avoid issues

#Install Packages:
#Load Packages
library(ncdf4) # file reading/writing
library(RNetCDF) #work with NetCDF
library(tidyverse) #work with data frames
library(lubridate)
library(dplyr)
library(raster)
library(abind)
library(geosphere)

###Define variable of interest in the NC file
var_of_interest = 'HC01' #methane is HC01 -- hydrocarbon with 1 carbon. 
date <- "2020_01"
name <- "CH4" #saving name

#############  CH4 On Road Gasoline ################# 
print('Analyzing On Road Gas')
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OnRoad_Gasoline/weekdy"))
folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OnRoad_Gasoline/weekdy")
files <- list.files(path = folder, pattern = "*.nc" , full.names = FALSE)

# empty arrays to use later for the SLC cropped area in 0.01o
onroad_Gas_0to23z <- array(numeric(),dim=c(100,60,0))
output <- array(numeric(),dim=c(100,60,0))
 
for(file in files){
  setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OnRoad_Gasoline/weekdy"))
  nc <- nc_open(file)
  lon <- (ncvar_get(nc,'XLONG'))
  lat <- (ncvar_get(nc,"XLAT"))
  CH4 <- (ncvar_get(nc,var_of_interest)) #focus on CH4
  dim(CH4)
  nc_close(nc) #close file helps with data space
  
  # Create loop to read through the [,,i] for the hourly slices
  
    for(i in 1:12){
    # Take the First slice out of 12, One slice for each hour
    # known for the NOAA CSL emission inventory files
    CH4_slice <- (CH4[,,i]) 
    
    ####
    # Convert Units of the absolute emissions:
    # Start in mole hr^-1 
    CH4_moles <- (CH4_slice/3600) #hours to seconds
    CH4_moles <- (CH4_moles*1000000) #moles to micro moles
    # Now umol s^-1 
    # Convert to flux emissions later in code after we rasterize values into the LCC grid

    #### FOR MOLES AND METRIC TONNES
    #     # Convert Units of the absolute emissions:
    #     # Start in mT hr^-1 
    #     NOX_slice <- (NOX_slice*1000000) #metric tons to grams
    #     NOX_moles <- (NOX_slice/46.005) #grams to moles
    #     NOX_moles <- (NOX_moles/3600) #hours to seconds
    #     NOX_moles <- (NOX_moles*1000000) #moles to micro moles
    #     # Now umol s^-1 
    #     # Convert to flux emissions later in code after we rasterize values into the LCC grid
    
    # Sum absolute emissions of the US for sanity check in umol s-1:
    old_emiss.US <- sum(CH4_moles) # umol s-1 per a grid cell
    
    # create the LCC coordinates into the raster of CH4_moles
    nrow <- 1332
    ncol <- 1008
    dat <- raster(
    xmn = 1, xmx = ncol,
    ymn = 1, ymx = nrow,
    nrows = nrow, ncols = ncol,
    vals = (CH4_moles)) 
    
    crs.lcc <- CRS("+proj=lcc +lat_0=39.34594 +lon_0=-97 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs")
    
    # get extent of the domain
    df <- data.frame(x= as.vector(lon), y= as.vector(lat))
    ext <- extent(df[, c('x', 'y')])
    
    # create a raster for extent with nrow and ncol
    rast <- raster(
      ext, nrow = nrow, ncol = ncol,
      crs = crs.lcc)
    
    # rasterize the data values
    rast <- rasterize(
    df[, c('x', 'y')], rast, CH4_moles)
    
    # resample Projection, dummy raster 0.01o gridding 
    blank <- raster(xmn=-138, xmx = -59, ymn= 19, ymx = 58, ncol = 7900, nrow= 3900) 
    crs(blank) <- "+proj=longlat"
    
    # Need to convert emissions to emissions per m2 before regridding. 
    # This will convert absolute emissions into flux emissions
    # To compute multiply by 4000x4000 to go from 4km2 (the grid cell size) to m2.
    emiss.lcc <-  rast/(4000*4000) #units: umol s-1 per a grid cell >>>> umol s-1 m-2
  
    # For bilinear resample of the original LCC grid:
    # from LCC to long lat:
    CH4_emis.bi <- resample(emiss.lcc, blank, method = "bilinear") # umol s-1 m-2
    #image(log10(CH4_emis.bi))
    
    # SANITY CHECK:
    # Convert back to umol s-1 PER grid cell for bug checking after re sample:
    # Due to the area function the return of the area checking will be in km2 because of the lat.lon CRS
    new_area.US <- raster::area(CH4_emis.bi)
    new_area.US <- new_area.US * (1000*1000) # from km2 to m2 for the absolute emissions check
    # check it in absolute emissions umole s-1
    check_emiss.bi <- CH4_emis.bi * new_area.US
    new_emiss.US <- cellStats(check_emiss.bi, stat = "sum") # umole s-1
    percent_diff <- (abs(old_emiss.US- new_emiss.US)/old_emiss.US)*100 #Percent difference <= 4%
    print(paste0("Percent difference = ",round(percent_diff,2),"%"))
   
    
    # If the Percent difference is less than 4% continue to crop for the SLV:
    # Define the extent with longitude and latitude bounds
    # check SLV emission Prior to Regridding: 
    lon_domain_bounds <- c(-112.5, -111.5)
    lat_domain_bounds <- c(40.4, 41.0)
    lonlat_domain_bounds <- c(lon_domain_bounds, lat_domain_bounds)
    
    domain_extent <- as(extent(lonlat_domain_bounds), 'SpatialPolygons')
    crs(domain_extent) <- "+proj=longlat"
    SLC.zoom <- crop(CH4_emis.bi, domain_extent) 
    # units of umol s-1 m-2 
    # saved into RDS file then places in reconfigured .nc
  
#Combine Slices of the pizza:
    if(i==1){
      CH4_output <- SLC.zoom}
    else{
      CH4_output<- stack(CH4_output, SLC.zoom)
    }
  } #end of slices addition loop for each file


  #bind CH4_output for both files 00 to 24 Z
  if(file == "onroad_00to12Z.nc" ){
    CH4_onroad_Gas_0to23z <- CH4_output
      }
  else{
    CH4_onroad_Gas_0to23z <- stack(CH4_onroad_Gas_0to23z, CH4_output)
  }
  
 } # Out of lat long determination & SLV Crop loop!
#test:
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
#setwd("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/RDS_CH4")
saveRDS(CH4_onroad_Gas_0to23z, file = paste0("CH4_",date,"_onroad_Gas_0to_23z_WKD_SLV.RDS"))


############# CH4  On Road Diesel ################# 
print('Analyzing On Road Diesel')
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OnRoad_Diesel/weekdy"))
folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OnRoad_Diesel/weekdy")

files <- list.files(path = folder, pattern = "*.nc" , full.names = FALSE)

#empty arrays to stuff later:
CH4_onroad_Diesel_0to23z <- array(numeric(),dim=c(100,60,0))
CH4_output <- array(numeric(),dim=c(100,60,0)) #checked an correct

for(file in files){
  
  setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OnRoad_Diesel/weekdy"))
  nc <- nc_open(file)
  lon <- (ncvar_get(nc,'XLONG'))
  lat <- (ncvar_get(nc,"XLAT"))
  CH4 <- (ncvar_get(nc,var_of_interest)) #focus on CH4
  dim(CH4)
  nc_close(nc) #close file helps with data space
  
  # Create loop to read through the [,,i] for the hourly slices
  
  for(i in 1:12){
    # Take the First slice out of 12, One slice for each hour
    # known for the NOAA CSL emission inventory files
    CH4_slice <- (CH4[,,i]) 
    
    ####
    # Convert Units of the absolute emissions:
    # Start in mole hr^-1 
    CH4_moles <- (CH4_slice/3600) #hours to seconds
    CH4_moles <- (CH4_moles*1000000) #moles to micro moles
    # Now umol s^-1 
    # Convert to flux emissions later in code after we rasterize values into the LCC grid
    
    # Sum absolute emissions of the US for sanity check in umol s-1:
    old_emiss.US <- sum(CH4_moles) # umol s-1 per a grid cell
    
    # create the LCC coordinates into the raster of CH4_moles
    nrow <- 1332
    ncol <- 1008
    dat <- raster(
      xmn = 1, xmx = ncol,
      ymn = 1, ymx = nrow,
      nrows = nrow, ncols = ncol,
      vals = (CH4_moles)) 
    
    crs.lcc <- CRS("+proj=lcc +lat_0=39.34594 +lon_0=-97 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs")
    
    # get extent of the domain
    df <- data.frame(x= as.vector(lon), y= as.vector(lat))
    ext <- extent(df[, c('x', 'y')])
    
    # create a raster for extent with nrow and ncol
    rast <- raster(
      ext, nrow = nrow, ncol = ncol,
      crs = crs.lcc)
    
    # rasterize the data values
    rast <- rasterize(
      df[, c('x', 'y')], rast, CH4_moles)
    
    # resample Projection, dummy raster 0.01o gridding 
    blank <- raster(xmn=-138, xmx = -59, ymn= 19, ymx = 58, ncol = 7900, nrow= 3900) 
    crs(blank) <- "+proj=longlat"
    
    # Need to convert emissions to emissions per m2 before regridding. 
    # This will convert absolute emissions into flux emissions
    # To compute multiply by 4000x4000 to go from 4km2 (the grid cell size) to m2.
    emiss.lcc <-  rast/(4000*4000) #units: umol s-1 per a grid cell >>>> umol s-1 m-2
    
    # For bilinear resample of the original LCC grid:
    CH4_emis.bi <- resample(emiss.lcc, blank, method = "bilinear") # umol s-1 m-2
    #image(log10(CH4_emis.bi))
    
    # SANITY CHECK:
    # Convert back to umol s-1 PER grid cell for bug checking after re sample:
    # Due to the area function the return of the area checking will be in km2 because of the lat.lon CRS
    new_area.US <- raster::area(CH4_emis.bi)
    new_area.US <- new_area.US * (1000*1000) # from km2 to m2 for the absolute emissions check
    # check it in absolute emissions umole s-1
    check_emiss.bi <- CH4_emis.bi * new_area.US
    new_emiss.US <- cellStats(check_emiss.bi, stat = "sum") # umole s-1
    percent_diff <- (abs(old_emiss.US- new_emiss.US)/old_emiss.US)*100 #Percent difference <= 4%
    print(paste0("Percent difference = ",round(percent_diff,2),"%"))
    
    # Define the extent with longitude and latitude bounds
    lon_domain_bounds <- c(-112.5, -111.5)
    lat_domain_bounds <- c(40.4, 41.0)
    lonlat_domain_bounds <- c(lon_domain_bounds, lat_domain_bounds)
    
    domain_extent <- as(extent(lonlat_domain_bounds), 'SpatialPolygons')
    crs(domain_extent) <- "+proj=longlat"
    SLC.zoom <- crop(CH4_emis.bi, domain_extent) 
    # units of umol s-1 m-2 
    # saved into RDS file then places in reconfigured .nc
    
    #Combine Slices of the pizza:
    
    if(i==1){
      CH4_output <- SLC.zoom}
    else{
      CH4_output<- stack(CH4_output, SLC.zoom)
    }
    
  } #end of slices addition loop for each file
  
  #bind CH4_output for both files 00 to 24 Z
  if(file == "onroad_00to12Z.nc" ){
    CH4_onroad_Diesel_0to23z <- CH4_output
  }
  else{
    CH4_onroad_Diesel_0to23z <- stack(CH4_onroad_Diesel_0to23z, CH4_output)
  }
  
} # Out of lat long determination & SLV Crop loop!

#setwd("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/RDS_CH4")
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
saveRDS(CH4_onroad_Diesel_0to23z, file = paste0("CH4_",date,"_onroad_Diesel_0to_23z_WKD_SLV.RDS"))


############# CH4 Off Road ################# 
print('Analyzing Off Road')
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OffRoad/weekdy"))
folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OffRoad/weekdy")

files <- list.files(path = folder, pattern = "*.nc" , full.names = FALSE)

#empty arrays to stuff later:
CH4_offroad_0to23z <- array(numeric(),dim=c(100,60,0))
CH4_output <- array(numeric(),dim=c(100,60,0)) #checked an correct

for(file in files){
  setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OffRoad/weekdy"))
  nc <- nc_open(file)
  lon <- (ncvar_get(nc,'XLONG'))
  lat <- (ncvar_get(nc,"XLAT"))
  CH4 <- (ncvar_get(nc,var_of_interest)) #focus on CH4
  dim(CH4)
  nc_close(nc) #close file helps with data space
  
  # Create loop to read through the [,,i] for the hourly slices
  
  for(i in 1:12){
    # Take the First slice out of 12, One slice for each hour
    # known for the NOAA CSL emission inventory files
    CH4_slice <- (CH4[,,i]) 
    
    ####
    # Convert Units of the absolute emissions:
    # Start in moles hr^-1 
    CH4_moles <- (CH4_slice/3600) #hours to seconds
    CH4_moles <- (CH4_moles*1000000) #moles to micro moles
    # Now umol s^-1 
    # Convert to flux emissions later in code after we rasterize values into the LCC grid
    
    # Sum absolute emissions of the US for sanity check in umol s-1:
    old_emiss.US <- sum(CH4_moles) # umol s-1 per a grid cell
    
    # create the LCC coordinates into the raster of CH4_moles
    nrow <- 1332
    ncol <- 1008
    dat <- raster(
      xmn = 1, xmx = ncol,
      ymn = 1, ymx = nrow,
      nrows = nrow, ncols = ncol,
      vals = (CH4_moles)) 
    
    crs.lcc <- CRS("+proj=lcc +lat_0=39.34594 +lon_0=-97 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs")
    
    # get extent of the domain
    df <- data.frame(x= as.vector(lon), y= as.vector(lat))
    ext <- extent(df[, c('x', 'y')])
    
    # create a raster for extent with nrow and ncol
    rast <- raster(
      ext, nrow = nrow, ncol = ncol,
      crs = crs.lcc)
    
    # rasterize the data values
    rast <- rasterize(
      df[, c('x', 'y')], rast, CH4_moles)
    
    # resample Projection, dummy raster 0.01o gridding 
    blank <- raster(xmn=-138, xmx = -59, ymn= 19, ymx = 58, ncol = 7900, nrow= 3900) 
    crs(blank) <- "+proj=longlat"
    
    # Need to convert emissions to emissions per m2 before regridding. 
    # This will convert absolute emissions into flux emissions
    # To compute multiply by 4000x4000 to go from 4km2 (the grid cell size) to m2.
    emiss.lcc <-  rast/(4000*4000) #units: umol s-1 per a grid cell >>>> umol s-1 m-2
    
    # For bilinear resample of the original LCC grid:
    CH4_emis.bi <- resample(emiss.lcc, blank, method = "bilinear") # umol s-1 m-2
    #image(log10(CH4_emis.bi))
    
    # SANITY CHECK:
    # Convert back to umol s-1 PER grid cell for bug checking after re sample:
    # Due to the area function the return of the area checking will be in km2 because of the lat.lon CRS
    new_area.US <- raster::area(CH4_emis.bi)
    new_area.US <- new_area.US * (1000*1000) # from km2 to m2 for the absolute emissions check
    # check it in absolute emissions umole s-1
    check_emiss.bi <- CH4_emis.bi * new_area.US
    new_emiss.US <- cellStats(check_emiss.bi, stat = "sum") # umole s-1
    percent_diff <- (abs(old_emiss.US- new_emiss.US)/old_emiss.US)*100 #Percent difference <= 4%
    print(paste0("Percent difference = ",round(percent_diff,2),"%"))
    
    
    # If the Percent difference is less than 4% continue to crop for the SLV:
    
    # Define the extent with longitude and latitude bounds
    lon_domain_bounds <- c(-112.5, -111.5)
    lat_domain_bounds <- c(40.4, 41.0)
    lonlat_domain_bounds <- c(lon_domain_bounds, lat_domain_bounds)
    
    domain_extent <- as(extent(lonlat_domain_bounds), 'SpatialPolygons')
    crs(domain_extent) <- "+proj=longlat"
    SLC.zoom <- crop(CH4_emis.bi, domain_extent) 
    # units of umol s-1 m-2 
    # saved into RDS file then places in reconfigured .nc
    #Combine Slices of the pizza:
    
    if(i==1){
      CH4_output <- SLC.zoom}
    else{
      CH4_output<- stack(CH4_output, SLC.zoom)
    }
    
  } #end of slices addition loop for each file
  
  #bind CH4_output for both files 00 to 24 Z
  if(file == "offroad_00to12Z.nc" ){
    CH4_offroad_0to23z <- CH4_output
  }
  else{
    CH4_offroad_0to23z <- stack(CH4_offroad_0to23z,CH4_output)
  }

#Hour_slices <- raster()
} # Out of lat long determination & SLV Crop loop!

#setwd("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/RDS_CH4")
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
saveRDS(CH4_offroad_0to23z, file =paste0("CH4_",date,"_offroad_0to_23z_WKD_SLV.RDS"))


#############  CH4 Oil and Gas ################# 
print('Analyzing Oil and Gas')
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OG/weekdy"))
folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OG/weekdy/")

files <- list.files(path = folder, pattern = "*.nc" , full.names = FALSE)

#empty arrays to stuff later:
CH4_OG_0to23z <- array(numeric(),dim=c(100,60,0))
CH4_output <- array(numeric(),dim=c(100,60,0)) #checked an correct

for(file in files){
  setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Area_OG/weekdy"))
  
  nc <- nc_open(file)
  lon <- (ncvar_get(nc,'XLONG'))
  lat <- (ncvar_get(nc,"XLAT"))
  CH4 <- (ncvar_get(nc,var_of_interest)) #focus on CH4
  dim(CH4)
  nc_close(nc) #close file helps with data space
  
  # Create loop to read through the [,,i] for the hourly slices
  
  for(i in 1:12){
    # Take the First slice out of 12, One slice for each hour
    # known for the NOAA CSL emission inventory files
    CH4_slice <- (CH4[,,i]) 
    
    ####
    # Convert Units of the absolute emissions:
    # Start in mole hr^-1 
    CH4_moles <- (CH4_slice/3600) #hours to seconds
    CH4_moles <- (CH4_moles*1000000) #moles to micro moles
    # Now umol s^-1 
    # Convert to flux emissions later in code after we rasterize values into the LCC grid
    
    # Sum absolute emissions of the US for sanity check in umol s-1:
    old_emiss.US <- sum(CH4_moles) # umol s-1 per a grid cell
    
    # create the LCC coordinates into the raster of CH4_moles
    nrow <- 1332
    ncol <- 1008
    dat <- raster(
      xmn = 1, xmx = ncol,
      ymn = 1, ymx = nrow,
      nrows = nrow, ncols = ncol,
      vals = (CH4_moles)) 
    
    crs.lcc <- CRS("+proj=lcc +lat_0=39.34594 +lon_0=-97 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs")
    
    # get extent of the domain
    df <- data.frame(x= as.vector(lon), y= as.vector(lat))
    ext <- extent(df[, c('x', 'y')])
    
    # create a raster for extent with nrow and ncol
    rast <- raster(
      ext, nrow = nrow, ncol = ncol,
      crs = crs.lcc)
    
    # rasterize the data values
    rast <- rasterize(
      df[, c('x', 'y')], rast, CH4_moles)
    
    # resample Projection, dummy raster 0.01o gridding 
    blank <- raster(xmn=-138, xmx = -59, ymn= 19, ymx = 58, ncol = 7900, nrow= 3900) 
    crs(blank) <- "+proj=longlat"
    
    # Need to convert emissions to emissions per m2 before regridding. 
    # This will convert absolute emissions into flux emissions
    # To compute multiply by 4000x4000 to go from 4km2 (the grid cell size) to m2.
    emiss.lcc <-  rast/(4000*4000) #units: umol s-1 per a grid cell >>>> umol s-1 m-2
    
    # For bilinear resample of the original LCC grid:
    CH4_emis.bi <- resample(emiss.lcc, blank, method = "bilinear") # umol s-1 m-2
    #image(log10(CH4_emis.bi))
    
    # SANITY CHECK:
    
    # Convert back to umol s-1 PER grid cell for bug checking after re sample:
    # Due to the area function the return of the area checking will be in km2 because of the lat.lon CRS
    new_area.US <- raster::area(CH4_emis.bi)
    new_area.US <- new_area.US * (1000*1000) # from km2 to m2 for the absolute emissions check
    # check it in absolute emissions umole s-1
    check_emiss.bi <- CH4_emis.bi * new_area.US
    new_emiss.US <- cellStats(check_emiss.bi, stat = "sum") # umole s-1
    percent_diff <- (abs(old_emiss.US- new_emiss.US)/old_emiss.US)*100 #Percent difference <= 4%
    print(paste0("Percent difference = ",round(percent_diff,2),"%"))

    # If the Percent difference is less than 4% continue to crop for the SLV:
    
    # Define the extent with longitude and latitude bounds
    lon_domain_bounds <- c(-112.5, -111.5)
    lat_domain_bounds <- c(40.4, 41.0)
    lonlat_domain_bounds <- c(lon_domain_bounds, lat_domain_bounds)
    
    domain_extent <- as(extent(lonlat_domain_bounds), 'SpatialPolygons')
    crs(domain_extent) <- "+proj=longlat"
    SLC.zoom <- crop(CH4_emis.bi, domain_extent) 
    # units of umol s-1 m-2 
    # saved into RDS file then places in reconfigured .nc
    #Combine Slices of the pizza:
    
    if(i==1){
      CH4_output <- SLC.zoom}
    else{
      CH4_output<- stack(CH4_output, SLC.zoom)
    }
  
} # Out of lat long determination & SLV Crop loop!
  
#bind CH4_output for both files:
if(file == "AreaFOGnei17_00to12Z.nc" ){
  CH4_OG_0to23z <- CH4_output
}
else{
  CH4_OG_0to23z <- stack(CH4_OG_0to23z, CH4_output)
}
  
} # Out of lat long determination & SLV Crop loop!

#setwd("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/RDS_CH4")
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
saveRDS(CH4_OG_0to23z, file = paste0("CH4_",date,"_OG_0to_23z_WKD_SLV.RDS"))


#############  CH4 Area AG ################# 
print('Analyzing Area AG')
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/AreaAG/weekdy"))
folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/AreaAG/weekdy/")
files <- list.files(path = folder, pattern = "*.nc" ,full.names = FALSE)

# empty arrays to use later for the SLC cropped area in 0.01o
CH4_AG_0to23z <- array(numeric(),dim=c(100,60,0))
CH4_output <- array(numeric(),dim=c(100,60,0))

 for(file in files){
   setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/AreaAG/weekdy"))
   nc <- nc_open(file)#
   lon <- (ncvar_get(nc,'XLONG'))
   lat <- (ncvar_get(nc,"XLAT"))
   CH4 <- (ncvar_get(nc,var_of_interest)) #focus on CH4
   dim(CH4)


   # Create loop to read through the [,,i] for the hourly slices

  for(i in 1:12){
     # Take the First slice out of 12, One slice for each hour
     # known for the NOAA CSL emission inventory files
     CH4_slice <- (CH4[,,i])
     print(sum(CH4_slice))
     
     ####
     # Convert Units of the absolute emissions:
     # Start in mole hr^-1 
     CH4_moles <- (CH4_slice/3600) #hours to seconds
     CH4_moles <- (CH4_moles*1000000) #moles to micro moles
     # Now umol s^-1 
     # Convert to flux emissions later in code after we rasterize values into the LCC grid

     # Sum absolute emissions of the US for sanity check in umol s-1:
      old_emiss.US <- sum(CH4_moles) # umol s-1 per a grid cell

     # create the LCC coordinates into the raster of CH4_moles
     nrow <- 1332
     ncol <- 1008
     dat <- raster(
     xmn = 1, xmx = ncol,
     ymn = 1, ymx = nrow,
    nrows = nrow, ncols = ncol,
       vals = (CH4_moles))

     crs.lcc <- CRS("+proj=lcc +lat_0=39.34594 +lon_0=-97 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs")

    # get extent of the domain
    df <- data.frame(x= as.vector(lon), y= as.vector(lat))
    ext <- extent(df[, c('x', 'y')])
    
    # create a raster for extent with nrow and ncol
    rast <- raster(
      ext, nrow = nrow, ncol = ncol,
      crs = crs.lcc)

     # rasterize the data values
     rast <- rasterize(
       df[, c('x', 'y')], rast, CH4_moles)

     # resample Projection, dummy raster 0.01o gridding
     blank <- raster(xmn=-138, xmx = -59, ymn= 19, ymx = 58, ncol = 7900, nrow= 3900)
     crs(blank) <- "+proj=longlat"

     # Need to convert emissions to emissions per m2 before regridding.
     # This will convert absolute emissions into flux emissions
     # To compute multiply by 4000x4000 to go from 4km2 (the grid cell size) to m2.
     emiss.lcc <-  rast/(4000*4000) #units: umol s-1 per a grid cell >>>> umol s-1 m-2

     # For bilinear resample of the original LCC grid:
     # from LCC to long lat:
     CH4_emis.bi <- resample(emiss.lcc, blank, method = "bilinear") # umol s-1 m-2
     #image(log10(CH4_emis.bi))

     # SANITY CHECK:
     # Convert back to umol s-1 PER grid cell for bug checking after re sample:
     # Due to the area function the return of the area checking will be in km2 because of the lat.lon CRS
     new_area.US <- raster::area(CH4_emis.bi)
     new_area.US <- new_area.US * (1000*1000) # from km2 to m2 for the absolute emissions check
     # check it in absolute emissions umole s-1
     check_emiss.bi <- CH4_emis.bi * new_area.US
     new_emiss.US <- cellStats(check_emiss.bi, stat = "sum") # umole s-1
     percent_diff <- (abs(old_emiss.US- new_emiss.US)/old_emiss.US)*100 #Percent difference <= 4%
     print(paste0("Percent difference = ",round(percent_diff,2),"%"))


     # If the Percent difference is less than 4% continue to crop for the SLV:
     # Define the extent with longitude and latitude bounds
     # check SLV emission Prior to Regridding:
     lon_domain_bounds <- c(-112.5, -111.5)
     lat_domain_bounds <- c(40.4, 41.0)
     lonlat_domain_bounds <- c(lon_domain_bounds, lat_domain_bounds)

     domain_extent <- as(extent(lonlat_domain_bounds), 'SpatialPolygons')
     crs(domain_extent) <- "+proj=longlat"
     SLC.zoom <- crop(CH4_emis.bi, domain_extent)
     # units of umol s-1 m-2
     # saved into RDS file then places in reconfigured .nc

     #Combine Slices of the pizza:
     if(i==1){
       CH4_output <- SLC.zoom}
     else{
       CH4_output<- stack(CH4_output, SLC.zoom)
     }
   } #end of slices addition loop for each file


   #bind CH4_output for both files 00 to 24 Z
   if(file == "AreaAG_00to12Z.nc" ){
     CH4_AG_0to23z <- CH4_output
   }
   else{
     CH4_AG_0to23z <- stack(CH4_AG_0to23z, CH4_output)
     }
   nc_close(nc) #close file helps with data space
 } # Out of lat long determination & SLV Crop loop!
 #test:
CH4_AG_0to23z <- array(numeric(0),dim=c(100,60,0))
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
saveRDS(CH4_AG_0to23z, file = paste0("CH4_",date,"_AG_NO.OG_0to_23z_WKD_SLV.RDS"))

# AREA INDUSTRY
print('Analyzing Area Industry')
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/AreaIndustry/weekdy"))
folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/AreaIndustry/weekdy")

files <- list.files(path = folder, pattern = "*.nc" , full.names = FALSE)

#empty arrays to stuff later:
CH4_Ind_0to23z <- array(numeric(),dim=c(100,60,0))
CH4_output <- array(numeric(),dim=c(100,60,0)) #checked an correct

for(file in files){
  
  setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/AreaIndustry/weekdy"))
  nc <- nc_open(file)
  lon <- (ncvar_get(nc,'XLONG'))
  lat <- (ncvar_get(nc,"XLAT"))
  CH4 <- (ncvar_get(nc,var_of_interest)) #focus on CH4
  dim(CH4)
  
  
  # Create loop to read through the [,,i] for the hourly slices
  
  for(i in 1:12){
    # Take the First slice out of 12, One slice for each hour
    # known for the NOAA CSL emission inventory files
    CH4_slice <- (CH4[,,i]) 
    print(sum(CH4_slice))
    
    ####
    # Convert Units of the absolute emissions:
    # Start in mole hr^-1 
    CH4_moles <- (CH4_slice/3600) #hours to seconds
    CH4_moles <- (CH4_moles*1000000) #moles to micro moles
    # Now umol s^-1 
    # Convert to flux emissions later in code after we rasterize values into the LCC grid
    
    # Sum absolute emissions of the US for sanity check in umol s-1:
    old_emiss.US <- sum(CH4_moles) # umol s-1 per a grid cell
    
    # create the LCC coordinates into the raster of CH4_moles
    nrow <- 1332
    ncol <- 1008
    dat <- raster(
      xmn = 1, xmx = ncol,
      ymn = 1, ymx = nrow,
      nrows = nrow, ncols = ncol,
      vals = (CH4_moles)) 
    
    crs.lcc <- CRS("+proj=lcc +lat_0=39.34594 +lon_0=-97 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs")
    
    # get extent of the domain
    df <- data.frame(x= as.vector(lon), y= as.vector(lat))
    ext <- extent(df[, c('x', 'y')])
    
    # create a raster for extent with nrow and ncol
    rast <- raster(
      ext, nrow = nrow, ncol = ncol,
      crs = crs.lcc)
    
    # rasterize the data values
    rast <- rasterize(
      df[, c('x', 'y')], rast, CH4_moles)
    
    # resample Projection, dummy raster 0.01o gridding 
    blank <- raster(xmn=-138, xmx = -59, ymn= 19, ymx = 58, ncol = 7900, nrow= 3900) 
    crs(blank) <- "+proj=longlat"
    
    # Need to convert emissions to emissions per m2 before regridding. 
    # This will convert absolute emissions into flux emissions
    # To compute multiply by 4000x4000 to go from 4km2 (the grid cell size) to m2.
    emiss.lcc <-  rast/(4000*4000) #units: umol s-1 per a grid cell >>>> umol s-1 m-2
    
    # For bilinear resample of the original LCC grid:
    CH4_emis.bi <- resample(emiss.lcc, blank, method = "bilinear") # umol s-1 m-2
    #image(log10(CH4_emis.bi))
    
    # SANITY CHECK:
    # Convert back to umol s-1 PER grid cell for bug checking after re sample:
    # Due to the area function the return of the area checking will be in km2 because of the lat.lon CRS
    new_area.US <- raster::area(CH4_emis.bi)
    new_area.US <- new_area.US * (1000*1000) # from km2 to m2 for the absolute emissions check
    # check it in absolute emissions umole s-1
    check_emiss.bi <- CH4_emis.bi * new_area.US
    new_emiss.US <- cellStats(check_emiss.bi, stat = "sum") # umole s-1
    percent_diff <- (abs(old_emiss.US- new_emiss.US)/old_emiss.US)*100 #Percent difference <= 4%
    print(paste0("Percent difference = ",round(percent_diff,2),"%"))
    
    # Define the extent with longitude and latitude bounds
    lon_domain_bounds <- c(-112.5, -111.5)
    lat_domain_bounds <- c(40.4, 41.0)
    lonlat_domain_bounds <- c(lon_domain_bounds, lat_domain_bounds)
    
    domain_extent <- as(extent(lonlat_domain_bounds), 'SpatialPolygons')
    crs(domain_extent) <- "+proj=longlat"
    SLC.zoom <- crop(CH4_emis.bi, domain_extent) 
    # units of umol s-1 m-2 
    # saved into RDS file then places in reconfigured .nc
    
    #Combine Slices of the pizza:
    
    if(i==1){
      CH4_output <- SLC.zoom}
    else{
      CH4_output<- stack(CH4_output, SLC.zoom)
    }
    
  } #end of slices addition loop for each file
  
  #bind CH4_output for both files 00 to 24 Z
  if(file == "AreaIndustry_00to12Z.nc" ){
    CH4_Ind_0to23z <- CH4_output
  }
  else{
    CH4_Ind_0to23z <- stack(CH4_Ind_0to23z, CH4_output)
  }
  nc_close(nc) #close file helps with data space
} # Out of lat long determination & SLV Crop loop!
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
saveRDS(CH4_Ind_0to23z, file =  paste0("CH4_",date,"_Industry_NO.OG_0to_23z_WKD_SLV.RDS"))

############# CH4 VCP ################# 
print('Analyzing VCP')
 setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/AreaVCP/weekdy"))
 folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/AreaVCP/weekdy/")
 
 files <- list.files(path = folder, pattern = "*.nc" , full.names = FALSE)
 
 #empty arrays to stuff later:
 CH4_VCP_0to23z <- array(numeric(),dim=c(100,60,0))
 CH4_output <- array(numeric(),dim=c(100,60,0)) #checked an correct
 
 for(file in files){
   setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/AreaVCP/weekdy"))
   nc <- nc_open(file)
   lon <- (ncvar_get(nc,'XLONG'))
   lat <- (ncvar_get(nc,"XLAT"))
   CH4 <- (ncvar_get(nc,var_of_interest)) #focus on CH4
   dim(CH4)
   
   # Create loop to read through the [,,i] for the hourly slices
   
   for(i in 1:12){
      #Take the First slice out of 12, One slice for each hour
     # known for the NOAA CSL emission inventory files
     CH4_slice <- (CH4[,,i]) 
     
     ####
     # Convert Units of the absolute emissions:
     # Start in mole hr^-1 
     CH4_moles <- (CH4_slice/3600) #hours to seconds
     CH4_moles <- (CH4_moles*1000000) #moles to micro moles
     # Now umol s^-1 
     # Convert to flux emissions later in code after we rasterize values into the LCC grid
     
     # Sum absolute emissions of the US for sanity check in umol s-1:
     old_emiss.US <- sum(CH4_moles) # umol s-1 per a grid cell
     
    # create the LCC coordinates into the raster of CH4_moles
     nrow <- 1332
     ncol <- 1008
     dat <- raster(
       xmn = 1, xmx = ncol,
       ymn = 1, ymx = nrow,
       nrows = nrow, ncols = ncol,
       vals = (CH4_moles)) 
     
     crs.lcc <- CRS("+proj=lcc +lat_0=39.34594 +lon_0=-97 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs")
     
     # get extent of the domain
     df <- data.frame(x= as.vector(lon), y= as.vector(lat))
     ext <- extent(df[, c('x', 'y')])
     
     # create a raster for extent with nrow and ncol
     rast <- raster(
       ext, nrow = nrow, ncol = ncol,
       crs = crs.lcc)
     
     # rasterize the data values
     rast <- rasterize(
       df[, c('x', 'y')], rast, CH4_moles)
     
     # resample Projection, dummy raster 0.01o gridding 
     blank <- raster(xmn=-138, xmx = -59, ymn= 19, ymx = 58, ncol = 7900, nrow= 3900) 
     crs(blank) <- "+proj=longlat"
     
     # Need to convert emissions to emissions per m2 before regridding. 
     # This will convert absolute emissions into flux emissions
     # To compute multiply by 4000x4000 to go from 4km2 (the grid cell size) to m2.
     emiss.lcc <-  rast/(4000*4000) #units: umol s-1 per a grid cell >>>> umol s-1 m-2
     
     # For bilinear resample of the original LCC grid:
     CH4_emis.bi <- resample(emiss.lcc, blank, method = "bilinear") # umol s-1 m-2
     #image(log10(CH4_emis.bi))
     
     # SANITY CHECK:
     # Convert back to umol s-1 PER grid cell for bug checking after re sample:
     # Due to the area function the return of the area checking will be in km2 because of the lat.lon CRS
     new_area.US <- raster::area(CH4_emis.bi)
     new_area.US <- new_area.US * (1000*1000) # from km2 to m2 for the absolute emissions check
     # check it in absolute emissions umole s-1
      check_emiss.bi <- CH4_emis.bi * new_area.US
     new_emiss.US <- cellStats(check_emiss.bi, stat = "sum") # umole s-1
     percent_diff <- (abs(old_emiss.US- new_emiss.US)/old_emiss.US)*100 #Percent difference <= 4%
     print(paste0("Percent difference = ",round(percent_diff,2),"%"))
     
     
     #If the Percent difference is less than 4% continue to crop for the SLV:
    
    # Define the extent with longitude and latitude bounds
     lon_domain_bounds <- c(-112.5, -111.5)
     lat_domain_bounds <- c(40.4, 41.0)
     lonlat_domain_bounds <- c(lon_domain_bounds, lat_domain_bounds)
     
     domain_extent <- as(extent(lonlat_domain_bounds), 'SpatialPolygons')
     crs(domain_extent) <- "+proj=longlat"
     SLC.zoom <- crop(CH4_emis.bi, domain_extent) 
     # units of umol s-1 m-2 
     # saved into RDS file then places in reconfigured .nc
     #Combine Slices of the pizza:
     
     if(i==1){
       CH4_output <- SLC.zoom}
     else{
       CH4_output<- stack(CH4_output, SLC.zoom)
     }
     
   } #end of slices addition loop for each file
   
   #bind CH4_output for both files 00 to 24 Z
   if(file == "AreaVCP_00to12Z.nc" ){
     CH4_VCP_0to23z <- CH4_output
   }
   else{
     CH4_VCP_0to23z <- stack(CH4_VCP_0to23z,CH4_output)
   }
   nc_close(nc) #close file helps with data space
   #Hour_slices <- raster()
 } # Out of lat long determination & SLV Crop loop!

 setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
saveRDS(CH4_VCP_0to23z, file = paste0("CH4_",date,"_VCP_NO.OG_0to_23z_WKD_SLV.RDS"))

########## OTHER NO OIL AND GAS ###########
print('Analyzing Other No Oil and Gas')
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/OtherArea/weekdy"))
folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Area/OtherArea/weekdy/")

files <- list.files(path = folder, pattern = "*.nc" , full.names = FALSE)

#empty arrays to stuff later:
CH4_Other_0to23z <- array(numeric(),dim=c(100,60,0))
CH4_output <- array(numeric(),dim=c(100,60,0)) #checked an correct

for(file in files){
  
  nc <- nc_open(file)
  lon <- (ncvar_get(nc,'XLONG'))
  lat <- (ncvar_get(nc,"XLAT"))
  CH4 <- (ncvar_get(nc,var_of_interest)) #focus on CH4
  dim(CH4)
  
  # Create loop to read through the [,,i] for the hourly slices
  
  for(i in 1:12){
    # Take the First slice out of 12, One slice for each hour
    # known for the NOAA CSL emission inventory files
    CH4_slice <- (CH4[,,i]) 
    print(sum(CH4_slice))
    ####
    # Convert Units of the absolute emissions:
    # Start in mole hr^-1 
    CH4_moles <- (CH4_slice/3600) #hours to seconds
    CH4_moles <- (CH4_moles*1000000) #moles to micro moles
    # Now umol s^-1 
    # Convert to flux emissions later in code after we rasterize values into the LCC grid
    
    # Sum absolute emissions of the US for sanity check in umol s-1:
    old_emiss.US <- sum(CH4_moles) # umol s-1 per a grid cell
    
    # create the LCC coordinates into the raster of CH4_moles
    nrow <- 1332
    ncol <- 1008
    dat <- raster(
      xmn = 1, xmx = ncol,
      ymn = 1, ymx = nrow,
      nrows = nrow, ncols = ncol,
      vals = (CH4_moles)) 
    
    crs.lcc <- CRS("+proj=lcc +lat_0=39.34594 +lon_0=-97 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs")
    
    # get extent of the domain
    df <- data.frame(x= as.vector(lon), y= as.vector(lat))
    ext <- extent(df[, c('x', 'y')])
    
    # create a raster for extent with nrow and ncol
    rast <- raster(
      ext, nrow = nrow, ncol = ncol,
      crs = crs.lcc)
    
    # rasterize the data values
    rast <- rasterize(
      df[, c('x', 'y')], rast, CH4_moles)
    
    # resample Projection, dummy raster 0.01o gridding 
    blank <- raster(xmn=-138, xmx = -59, ymn= 19, ymx = 58, ncol = 7900, nrow= 3900) 
    crs(blank) <- "+proj=longlat"
    
    # Need to convert emissions to emissions per m2 before regridding. 
    # This will convert absolute emissions into flux emissions
    # To compute multiply by 4000x4000 to go from 4km2 (the grid cell size) to m2.
    emiss.lcc <-  rast/(4000*4000) #units: umol s-1 per a grid cell >>>> umol s-1 m-2
    
    # For bilinear resample of the original LCC grid:
    CH4_emis.bi <- resample(emiss.lcc, blank, method = "bilinear") # umol s-1 m-2
    #image(log10(CH4_emis.bi))
    
    # SANITY CHECK:
    
    # Convert back to umol s-1 PER grid cell for bug checking after re sample:
    # Due to the area function the return of the area checking will be in km2 because of the lat.lon CRS
    new_area.US <- raster::area(CH4_emis.bi)
    new_area.US <- new_area.US * (1000*1000) # from km2 to m2 for the absolute emissions check
    # check it in absolute emissions umole s-1
    check_emiss.bi <- CH4_emis.bi * new_area.US
    new_emiss.US <- cellStats(check_emiss.bi, stat = "sum") # umole s-1
    percent_diff <- (abs(old_emiss.US- new_emiss.US)/old_emiss.US)*100 #Percent difference <= 4%
    print(paste0("Percent difference = ",round(percent_diff,2),"%"))
    
    # If the Percent difference is less than 4% continue to crop for the SLV:
    
    # Define the extent with longitude and latitude bounds
    lon_domain_bounds <- c(-112.5, -111.5)
    lat_domain_bounds <- c(40.4, 41.0)
    lonlat_domain_bounds <- c(lon_domain_bounds, lat_domain_bounds)
    
    domain_extent <- as(extent(lonlat_domain_bounds), 'SpatialPolygons')
    crs(domain_extent) <- "+proj=longlat"
    SLC.zoom <- crop(CH4_emis.bi, domain_extent) 
    # units of umol s-1 m-2 
    # saved into RDS file then places in reconfigured .nc
    #Combine Slices of the pizza:
    
    if(i==1){
      CH4_output <- SLC.zoom}
    else{
      CH4_output<- stack(CH4_output, SLC.zoom)
    }
    
  } # Out of lat long determination & SLV Crop loop!
  
  #bind CH4_output for both files:
  if(file == "OtherArea_00to12Z.nc" ){
    CH4_Other_0to23z <- CH4_output
  }
  else{
    CH4_Other_0to23z <- stack(CH4_Other_0to23z, CH4_output)
  }
  nc_close(nc) #close file helps with data space
} # Out of lat long determination & SLV Crop loop!

setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
saveRDS(CH4_Other_0to23z, file = paste0("CH4_",date,"_Other_NO.OG_0to_23z_WKD_SLV.RDS"))

############# Point Sources 00 to 12Z ################# 
############ CH4 ##############
print('Analyzing Point Sources 00 to 12z')

setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Sources/WKD_00to12z/"))
folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Sources/WKD_00to12z/")
files_00to12Z <- list.files(path = folder, pattern = "*00to12Z.nc" , full.names = FALSE)

#Place over all outputs out of the Loops:
X_output <- as.array(0, nrow = 60, ncol = 100)

for(file in files_00to12Z){ # Loops through all files in list
  
  setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Sources/WKD_00to12z/"))
  ncin <- nc_open(file)
  lon <- (ncvar_get(ncin,'XLONG'))
  lon <- as.numeric(lon)
  lat <- (ncvar_get(ncin,"XLAT"))
  lat <- as.numeric(lat)
  X.alls <- ncvar_get(ncin,var_of_interest) #focus on CH4
  print(sum(X.alls))
  nc_close(ncin) #important!!! Close the connection before doing other things
  
  for(i in 1:12){   #Loops through Slices in files
    
    print(i)
    X.slice <- (X.alls[,i]) 
    
    # This is the rows and columns used for the area files:
    lat_min <- 40.4
    lat_max <- 41.0
    lon_min <- -112.5
    lon_max <- -111.5
    lat_extent <- c(lat_min, lat_max)
    lon_extent <- c(lon_min, lon_max)
    # data with specified latitudes and longitudes from the ncdf file
    # Full data set for the USA, subset for the SLV
    X.slice <- as.data.frame(X.slice)
    df <- data.frame(lon,lat,X.slice)
    coords <- cbind(df$lon, df$lat)
    
    subset_df.a <- subset(df, df$lon >= lon_min & df$lon <= lon_max)
    subset_df <- subset(subset_df.a, subset_df.a$lat >= lat_min & subset_df.a$lat <= lat_max)
    data <- subset_df 
    raster_template <- raster(extent(lon_extent, lat_extent), resolution = 0.01)
    grid_matrix <- rasterize(data[, c("lon", "lat")], raster_template, data$X.slice)
    crs(grid_matrix) <- "+proj=longlat"
    
    #X <- (grid_matrix*1000000) #metric tons to grams
    #X_moles <- (X/44.01) #grams to moles
    #CH4 is in mol hr-1
    X_moles <- grid_matrix
    X_moles <- (X_moles/3600) #hours to seconds
    X_moles <- (X_moles*1000000) #moles to micro moles
    
    area.SLV <- raster::area(X_moles)*(1000*1000) # from km2 to m2 for the absolute emissions check
    grid <- X_moles / area.SLV # this gives us umol m-2 s-1
    
    #Combine slices of arrays:
    if(i==1){
      X_output <- grid}
    else{
      X_output <- stack(X_output, grid)
    }
  } # Loop for slices
  
  setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
  output_file <- paste0("CH4_",date,"_00.12_WKD_SLV", file, ".rds")
  # Save the processed data as RDS
  saveRDS(X_output, file = output_file) 
  
} # Loop for files 

############# Point Sources 12 to 23Z ################# 
############ CH4 ##############
print('Analyzing Point Sources 12 to 23Z')
setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Sources/WKD_12to23z/"))
folder <- paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Sources/WKD_12to23z")
files_12to23Z <- list.files(path = folder, pattern = "*12to24Z.nc" , full.names = FALSE)

#Place over all outputs out of the Loops:
X_output <- as.array(0, nrow = 60, ncol = 100)

for(file in files_12to23Z){ # Loops through all files in list
  setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/EI_",date,"/Point_Sources/WKD_12to23z/"))
  ncin <- nc_open(file)
  lon <- (ncvar_get(ncin,'XLONG'))
  lon <- as.numeric(lon)
  lat <- (ncvar_get(ncin,"XLAT"))
  lat <- as.numeric(lat)
  X.alls <- ncvar_get(ncin,var_of_interest) #focus on CH4
  print(sum(X.alls))
  nc_close(ncin) #important!!! Close the connection before doing other things
  
  for(i in 1:12){   #Loops through Slices in files
    
    print(i)
    X.slice <- (X.alls[,i]) 
    
    # This is the rows and columns used for the area files:
    lat_min <- 40.4
    lat_max <- 41.0
    lon_min <- -112.5
    lon_max <- -111.5
    lat_extent <- c(lat_min, lat_max)
    lon_extent <- c(lon_min, lon_max)
    # data with specified latitudes and longitudes from the ncdf file
    # Full data set for the USA, subset for the SLV
    X.slice <- as.data.frame(X.slice)
    df <- data.frame(lon,lat,X.slice)
    coords <- cbind(df$lon, df$lat)
    
    subset_df.a <- subset(df, df$lon >= lon_min & df$lon <= lon_max)
    subset_df <- subset(subset_df.a, subset_df.a$lat >= lat_min & subset_df.a$lat <= lat_max)
    data <- subset_df 
    raster_template <- raster(extent(lon_extent, lat_extent), resolution = 0.01)
    grid_matrix <- rasterize(data[, c("lon", "lat")], raster_template, data$X.slice)
    crs(grid_matrix) <- "+proj=longlat"
    
    #X <- (grid_matrix*1000000) #metric tons to grams
    #CH4 already in mol hr-1
    X_moles <- grid_matrix
    X_moles <- (X_moles/3600) #hours to seconds
    X_moles <- (X_moles*1000000) #moles to micro moles
    
    area.SLV <- raster::area(X_moles)*(1000*1000) # from km2 to m2 for the absolute emissions check
    grid <- X_moles / area.SLV # this gives us umol m-2 s-1
    
    #Combine slices of arrays:
    if(i==1){
      X_output <- grid}
    else{
      X_output <- stack(X_output, grid)
    }
    
  } # Loop for slices
  
  setwd(paste0("/uufs/chpc.utah.edu/common/home/lin-group10/HKH/CH4_RDS_test/",date))
  output_file <- paste0("CH4_",date,"_12.24_WKD_SLV_", file, ".rds")
  # Save the processed data as RDS
  saveRDS(X_output, file = output_file)
} # Loop for files


