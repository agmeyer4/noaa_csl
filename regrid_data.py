'''
Module: regrid_data.py
Author: Aaron G. Meyer (agmeyer4@gmail.com)
Last Updated: May 2, 2024
Description:
A script to regrid base data from LCC into lat/lon coordinates 
'''

#Import modules
import xarray as xr
import pandas as pd
import os
import noaa_csl_funcs as ncf
import pyproj
import time
import numpy as np
import xesmf as xe

# Define Classes
class RegridInputs():
    '''A class to hold the inputs for the CSL regridder. These are the defaults.'''

    def __init__(self,**kwargs):
        self.grid_out = {  
                        'lat': np.arange(19, 58, 0.1), #Center Point Spacing Lat
                        'lon': np.arange(-138, -59, 0.1), #Center Point Spacing Lon
                        'lat_b': np.arange(18.95, 58.05, 0.1), # Boundary Spacing Lat 
                        'lon_b': np.arange(-138.05, -58.95, 0.1), # Boundary Spacing Lon
                        }
        self.method = 'conservative'
        self.input_dims=('south_north','west_east')
        self.weights_path = '/uufs/chpc.utah.edu/common/home/u0890904/NOAA_CSL/noaa_csl/regridding/saved_weights'
        self.weights_file = 'lcc_to_latlon_onroad_gasoline_2019_1_weekdy.nc'
        self.regridded_path =  '/uufs/chpc.utah.edu/common/home/lin-group9/agm/NOAA_CSL_Data/regridded'
        for k,v in kwargs.items():
            if k in self.__dict__.keys():
                setattr(self,k,v)

def make_regrid_path(bau_or_covid,regridded_path,sector,year,month,day_type):
    '''Create a path for the regridded data to be saved to mirroring the base data structure
    
    Args:
    bau_or_covid (str) : "COVID" default, or "BAU business as usual
    regridded_path (str) : path to where the regridded data will be saved
    sector (str) : the sector ex (area_OG)
    year (int) : integer year
    month (int) : integer month
    day_type (str) : 'satdy', 'sundy', or 'weekdy'
    
    Returns:
    full path (str) : string to the day path where the nc file will go an dthe path that was created
    
    Raises:
    FileExistsError: This is a check to make sure you don't overwrite data that already exists 
    '''
    yr_str = ncf.yr_to_yrstr(sector,year,bau_or_covid) #get the year string
    month_str = ncf.month_int_to_str(month) #get the month string 
    full_path = os.path.join(regridded_path,sector,yr_str,month_str,day_type) #define the full path
    try: 
        os.makedirs(full_path) #try creating the path
        return full_path #and return it 
    except FileExistsError: #if it exists, raise error
        raise FileExistsError('Folder already exists. You may end up overwriting data if you continue.')
    
def regrid_and_save(BCH,unit_converter,csl_regridder,sector,year,month,day_type,save_fname = 'default',sanity_check_specs = None,grid_area_path='./regridding/grid_area/grid_out_area.nc'):
    '''Main function to regrid and save a single base .nc file given the parameters input
    
    Args:
    BCH (ncf.Base_CSL_Handler object) : an object which handles the formatting of base data
    unit_converter (ncf.CSL_Unit_converter object) : an object which converts units in the xarray datasets
    csl_regridder (ncf.CSL_Regridder object) : an object defining the regridding parameters
    sector (str) : the sector to regrid (ex "area_OG")
    year (int) : year integer
    month (int) : month integer
    day_type (str) : 'satdy', 'sundy', or 'weekdy'
    save_fname (str, optional) : if you want a specific filename when saving, it can be input here. otherise the default name will be {sector}_regridded.nc
    sanity_check_specs (list, optional) : the species (xarray data_vars) that you want to do a sanity check on -- compares sums of regrid vs sum of base data on these species
    grid_area_path (str, optional) : path to a file containing the grid cell area .nc file for the regridded data. defaults to the save one in this git created using cdo 
    '''

    day_regrid_path = make_regrid_path(BCH.bau_or_covid,csl_regridder.inputs.regridded_path,sector,year,month,day_type) #define the path to the 
    base_ds = BCH.load_fullday_nc(sector,year,month,day_type) #load the day's netcdf
    ds_flx = unit_converter.absolute_to_flux(base_ds) #convert from units of per gridcell in the LCC projection to per m2 (flux) prior to regridding
    regridded_ds = csl_regridder.regrid_ds(ds_flx) #regrid the dataset 
    regridded_ds = ncf.add_githash_to_ds(regridded_ds) #add the current git hash to track how the data was processed

    if sanity_check_specs is not None: #if we ant to sanity check 
        grid_area = xr.open_dataset(grid_area_path) #load the grid area ds
        for species in sanity_check_specs: #only do the species identified
            perc_diff = ncf.sanity_check(base_ds,regridded_ds,grid_area,species) #do the sanity check
            if perc_diff > 2.0: #if its more than 2 percent, alert
                print(f'*********Warning*********')
                print(f'{round(perc_diff,1)}% difference between base and regridded sums for {species} {sector} {year} {month} {day_type}')
                print('**************************')

    if save_fname == 'default': #if we're using the default save filename
        save_fname = f'{sector}_regridded.nc' #define it here
    print('Saving regridded dataset')
    regridded_ds.to_netcdf(os.path.join(day_regrid_path,save_fname)) #and save the regridded_ds to an nc file. 

#Define main function
def main():
    '''The main function to regrid data. This can be done ad-hoc but is currently set up to do a bulk regrid (all of the base data from 01-2019 to 08-2021). 
    I believe regridding the entire dataset at once is likely beneficial, so that there are not inconsistancies between datasets that may have been processed
    in different ways. To combat this, I have added the git hash to each regridded nc file as an attribute, so it should be reproducable. 
    
    Running python regrid_data.py will:
    Define the parameters of the regrid
    Create the path for the new .nc file to be saved
    Create the regridder based on the regrid inputs
    Load the base data
    Convert from absolute to flux units
    Regrid the data
    Sanity check on the defined species
    Save the .nc file 
    Alert the user of times to track how long things are taking
    '''

    print('Regridding NOAA data from base (Lambert Conical Conformal) to regridded (Lat/Lon)')

    t1 = time.time()
    print(f'Git hash of this regrid = {ncf.get_githash()}')
    print(f'Starting regrid at {t1}')

    base_path = '/uufs/chpc.utah.edu/common/home/lin-group9/agm/NOAA_CSL_Data/base' #where the data downloaded using data_download.py lives
    print(f'Saving regridded .nc files to {base_path}')
    bau_or_covid = 'COVID'
    BCH = ncf.Base_CSL_Handler(base_path,bau_or_covid) #setup the base data handler
    unit_converter = ncf.CSL_Unit_Converter() #setup the unit converter

    weights_file = 'lcc_to_latlon_onroad_gasoline_2019_1_weekdy.nc' #or 'create', for this I am using a weights file created from the onroad gas
    inputs = RegridInputs(weights_file = weights_file)#Define the inputs. these are the defaults 
    csl_regridder = ncf.CSL_Regridder(inputs) #create the regridder class

    sanity_check_specs = ['CO2','CO','HC01','NOX'] #the species we are most interested in
    area_sectors = [s for s in ncf.listdir_visible(base_path) if s.startswith('area')] #all of the area sectors in the base path
    years = [2019,2020,2021] #all of the 
    months = list(range(1,13)) #all of the months 
    day_types = ['weekdy','satdy','sundy'] #all of the day types
    for sector in area_sectors: #loop through the sectors
        for year in years: #loop through years
            for month in months: #loop through months
                if (year==2021) & (month>8): #if we're in 2021 theres only data up to august
                    continue  #so just get to the end of the array
                for day_type in day_types: #loop through the day types
                    print(f'Regridding {sector} {year} {month} {bau_or_covid} {day_type}')
                    try: #put it in a try loop so we can log time
                        regrid_and_save(BCH,unit_converter,csl_regridder,sector,year,month,day_type,sanity_check_specs=sanity_check_specs) #the main regrid function
                    except Exception as e: #grab exceptions
                        print(f'Error at {time.time()}') #print the time
                        raise Exception(e) #still give the excption
                    print('') #print a blank line between .nc files

    #get the business as usual data regridded for 2020, only applies to traffic data
    bau_or_covid = 'BAU'
    BCH = ncf.Base_CSL_Handler(base_path,bau_or_covid) #setup the base data handler
    bau_sectors = ['area_onroad_gasoline','area_onroad_diesel','area_offroad'] #the traffic data with bau data
    years = [2020] #only 2020
    months = list(range(1,13)) #all of the months 
    day_types = ['weekdy','satdy','sundy'] #all of the day types
    for sector in bau_sectors: #loop through the sectors
        for year in years: #loop through years
            for month in months: #loop through months
                for day_type in day_types: #loop through the day types
                    print(f'Regridding {sector} {year} {month} {bau_or_covid} {day_type}')
                    try: #put it in a try loop so we can log time
                        regrid_and_save(BCH,unit_converter,csl_regridder,sector,year,month,day_type,sanity_check_specs=sanity_check_specs) #the main regrid function
                    except Exception as e: #grab exceptions
                        print(f'Error at {time.time()}') #print the time
                        raise Exception(e) #still give the excption
                    print('') #print a blank line between .nc files

    t2 = time.time()
    print(f'ended at {t2}')
    print(f'total runtime (seconds) = {t2-t1}')

if __name__ == "__main__":
    main()
