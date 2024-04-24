#######################
### Updated: May, 2024
#' smaller grid size to fit the SLV
#' 
#' Description: 
#' Function definitions for working with NOAA CSL emissions inventories
#' 
#' All Sectors regridded and unit conversions:
#' Regridded from 4km grid cells in Lambert Conical Conformal
#' Unit conversion from absolute emissions (mT hr-1) to flux emissions ( umol s-1 m-2)
#' 
#' File combinations 
### By: Haley Humble, Aaron Meyer
#######################

#Load Packages
library(ncdf4) # file reading/writing
library(RNetCDF) #work with NetCDF
library(tidyverse) #work with data frames
library(lubridate)
library(dplyr)
library(raster)
library(abind)
library(geosphere)

#' Loads a variable from a netcdf file
#' 
#' @param full_file_path The full path to the netcdf file
#' @param var_name The variable of interest in the CSL inventories
#' @returns List containing the Latitude array, Longitude array, and variable of interest array
load_nc_var <- function(full_file_path,var_name){
    nc <- nc_open(full_file_path)
    lon <- (ncvar_get(nc,'XLONG'))
    lat <- (ncvar_get(nc,"XLAT"))
    var_arr <- (ncvar_get(nc,var_name)) 
    nc_close(nc) #close file helps with data space
    output <- list(lat,lon,var_arr)
    names(output) <- c('lat','lon',var_name)
    return(output)
}

raw_data_path <- "/uufs/chpc.utah.edu/common/home/lin-group10/HKH/Public/NOAA_CSL_Raw_Data/"
year_month <- '2020_01'
source_ID <- 'Area_OffRoad'
day_type <- 'weekdy' #'weekdy' for weekday, 'satdy' for saturday, 'sundy' for sunday
var_name <- 'HC01'

nc_path <- file.path(raw_data_path,paste0('EI_',year_month),source_ID,day_type)
nc_files <- list.files(path = nc_path, pattern = "*.nc" , full.names = TRUE)

full_nc_path <- nc_files[[1]]


nc_data = nc_open(full_nc_path)
lon <- (ncvar_get(nc_data,"XLONG"))
lat <- (ncvar_get(nc_data,"XLAT"))
var_arr <- (ncvar_get(nc_data,var_name)) 
nc_close(nc_data) #close file helps with data space


dat <- raster(
xmn = 1, xmx = 1008,
ymn = 1, ymx = 1332,
nrows = nrow, ncols = ncol,
vals = CH4_slice)

# get extent of the domain
df <- data.frame(x= as.vector(lon), y= as.vector(lat))
ext <- extent(df[, c('x', 'y')])

crs.lcc <- CRS("+proj=lcc +lat_0=39.34594 +lon_0=-97 +lat_1=33 +lat_2=45 +x_0=0 +y_0=0 +R=6370000 +units=m +no_defs")

nrow <- 1332
ncol <- 1008

# create a raster for extent with nrow and ncol
rast <- raster(
  ext, nrow = nrow, ncol = ncol,
  crs = crs.lcc)


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