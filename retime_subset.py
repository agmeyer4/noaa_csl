#Import modules
import xarray as xr
import pandas as pd
import os
import pyproj
import numpy as np
import xesmf as xe
import calendar
import datetime
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.img_tiles as cimgt
import matplotlib.pyplot as plt
import sys
import noaa_csl_funcs as ncf

def rework_datetime(ds):
    combined_ds_list = []
    for yr_mo in ds.yr_mo.values:
        year, month = map(int,yr_mo.split('-'))
        month_ds = ds.sel(yr_mo=yr_mo)

        dates_in_month = pd.date_range(start=f'{year}-{month:02}-01', end=f'{year}-{month:02}-{calendar.monthrange(year, month)[1]}')
        dates_by_day_type = {
            'weekdy':[date for date in dates_in_month if date.weekday() < 5],
            'satdy':[date for date in dates_in_month if date.weekday() == 5],
            'sundy':[date for date in dates_in_month if date.weekday() == 6]
        }

        new_month_ds_list = []
        for day_type,dates in dates_by_day_type.items():
            subds = month_ds.sel(day_type=day_type).drop_vars(['day_type','yr_mo'])
            subds = subds.assign_coords({'date':dates})
            datetimes = [pd.Timestamp(date) + pd.Timedelta(hours=int(hour)) for date in subds.coords['date'].values for hour in subds.coords['utc_hour'].values]
            datetime_index = pd.DatetimeIndex(datetimes)

            subds = subds.stack({'datetime':('date','utc_hour')})#.assign_coords({'datetime':datetime_index})
            subds = subds.drop_vars(['date','utc_hour','datetime'])
            subds = subds.assign_coords({'datetime':datetime_index})
            new_month_ds_list.append(subds)

        new_month_ds = xr.concat(new_month_ds_list,dim='datetime').sortby('datetime')
        combined_ds_list.append(new_month_ds)

    combined_ds = xr.concat(combined_ds_list,dim='datetime').sortby('datetime')
    return combined_ds


#Define main function
def main():
    # Load from the regridded files for a single year
    regrid_id = '_0.025deg'
    regridded_path = f'/uufs/chpc.utah.edu/common/home/lin-group9/agm/NOAA_CSL_Data/regridded{regrid_id}'
    RCH = ncf.Regridded_CSL_Handler(regridded_path)
    dt1  = pd.to_datetime('2019-01-01 00') 
    dt2 = pd.to_datetime('2019-12-31 23') 
    day_types = ['weekdy','satdy','sundy'] #a list with any or all of 'weekdy','satdy','sundy'
    species = ['CO2','CO','HC01','HC02','HC14','NH3','NOX','SO2']
    dataset_extent = {'lon_min':-112.1,
                    'lon_max':-111.6,
                    'lat_min':40.3,
                    'lat_max':41.1} 
    sector_types = ['area','point']

    combined_dss = {}
    for sector_type in sector_types:
        print('Loading data for sector type:',sector_type)
        #Get the paths to the files that match the criteria
        days_paths = RCH.get_days_in_range(dt1,dt2,day_types,sector_type) 
        files = RCH.get_files_in_days(days_paths)

        load_vars = [] #the variables we want to load
        if sector_type == 'area': #for area sources, just load the species defined above
            load_vars = species.copy() 
        if sector_type == 'point': #for point, often useful to have the stack/type information 
            load_vars = species.copy()
            load_vars.extend(['ITYPE','STKht','STKdiam','STKtemp','STKve','STKflw','FUGht']) #so add it to the actual species

        #Load the files with xarray, preprocessing them so they can be combined by coordinates
        ds_list = [] #initialize the list of datasets
        for file in files:
            ds = RCH.preprocess_regridded(xr.open_dataset(file,chunks = {'utc_hour':1}),dataset_extent)[load_vars] #prepreprocess the file, open with dask chunking, and only keep the species of interest
            ds_list.append(ds)  
        ds_combined = xr.combine_by_coords(ds_list,combine_attrs='drop_conflicts') #this is the combined dataset!
        combined_dss[sector_type] = ds_combined

    dt_dss = {}
    for sector_type in sector_types:
        print('Reworking datetime for sector type:',sector_type)
        dt_dss[sector_type] = rework_datetime(combined_dss[sector_type])

    #Save to nc
    save_path = '/uufs/chpc.utah.edu/common/home/u0890904/LAIR_1/Data/NC'
    save_prefix = 'slv_2019_dt'
    for sector_type in sector_types:
        print('Saving data for sector type:',sector_type)
        dt_dss[sector_type].to_netcdf(os.path.join(save_path,f"{save_prefix}_{sector_type}.nc"))

if __name__ == "__main__":
    main()
