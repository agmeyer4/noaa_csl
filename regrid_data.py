'''
Module: regrid_data.py
Author: Aaron G. Meyer (agmeyer4@gmail.com)
Last Updated: May 2, 2024
Description:
A script to regrid base data into lat/lon coordinates 
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
    yr_str = ncf.yr_to_yrstr(sector,year,bau_or_covid)
    month_str = ncf.month_int_to_str(month)
    full_path = os.path.join(regridded_path,sector,yr_str,month_str,day_type)
    try:
        os.makedirs(full_path)
        return full_path
    except FileExistsError:
        raise FileExistsError('Folder already exists. You may end up overwriting data if you continue.')
    
def regrid_and_save(BCH,unit_converter,cr,sector,year,month,day_type,save_fname = 'default'):
    day_regrid_path = make_regrid_path(BCH.bau_or_covid,cr.inputs.regridded_path,sector,year,month,day_type)
    base_ds = BCH.load_fullday_nc(sector,year,month,day_type)
    ds_flx = unit_converter.absolute_to_flux(base_ds)
    regridded_ds = cr.regrid_ds(ds_flx)
    if save_fname == 'default':
        save_fname = f'{sector}_regridded.nc'
    print('Saving')
    regridded_ds.to_netcdf(os.path.join(day_regrid_path,save_fname))

#Define main function
def main():
    t1 = time.time()

    # Regrid here
    base_path = '/uufs/chpc.utah.edu/common/home/lin-group9/agm/NOAA_CSL_Data/base' #where the data downloaded using data_download.py lives
    bau_or_covid = 'COVID'
    BCH = ncf.Base_CSL_Handler(base_path,bau_or_covid)
    unit_converter = ncf.CSL_Unit_Converter() #setup the unit converter

    weights_file = 'lcc_to_latlon_onroad_gasoline_2019_1_weekdy.nc' #or 'create'
    inputs = RegridInputs(weights_file = weights_file)#Define the inputs. these are the defaults
    cr = ncf.CSL_Regridder(inputs) #create the regridder class

    sectors = ['area_onroad_gasoline']
    years = [2019]
    months = [1]
    day_types = ['satdy','sundy','weekdy']
    for sector in sectors:
        for year in years:
            for month in months:
                for day_type in day_types:
                    print(f'Regridding {sector} {year} {month} {day_type}')
                    regrid_and_save(BCH,unit_converter,cr,sector,year,month,day_type)

    t2 = time.time()
    print(f'total runtime (seconds) = {t2-t1}')

if __name__ == "__main__":
    main()
