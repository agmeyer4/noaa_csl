'''
Module: noaa_csl_funcs.py
Author: Aaron G. Meyer (agmeyer4@gmail.com)
Last Updated: May 2, 2024
Description:
Useful functions for dealing with noaa csl data

More information can be found here: https://csl.noaa.gov/groups/csl7/measurements/2020covid-aqs/emissions/
'''

##################################################################################################################################
#Import Packages
import os
import subprocess
import shutil
import glob
import xarray as xr
import pyproj
import sys
import numpy as np
import xesmf as xe
import pandas as pd
import git
import calendar

##################################################################################################################################
# Define Functions
def get_githash():
    '''Gets the git hash of the current code to save as metadata in regridded/changed objects
    
    Returns:
    githash (str) : the hash of the current git commit
    '''
    repo = git.Repo(search_parent_directories=True)
    sha = repo.head.object.hexsha
    return sha

def add_githash_to_ds(ds):
    '''Adds the current git hash as an attribute to a dataset
    
    Args:
    ds (xarray.DataSet): the dataset to add the hash to 
    
    Returns:
    ds (xarray.DataSet) : the dataset with the current git hash added as "sha" attribute 
    '''

    git_sha = get_githash()
    ds.attrs['git_sha'] = git_sha
    return ds


def replace_all_strs(text,dic):
    '''Replaces all strings contained in text matching the keys of dic with the values of dic
    
    Args:
    text (str) : the string we want to replace stuff in
    dic (dict) : a dictionary with keys as what should be replaced in text and values with what the key should be replaced with
    
    Returns :
    text (str): the new text with keys replaced
    '''

    for key,value in dic.items():
        text = text.replace(key,value)
    return text

def month_int_to_str(month_int):
        '''Converts a month integer to a character string like "MonthXX"
        
        Args: 
        month_int (int) : the integer of the month (1,2,12,etc)
        
        Returns (str) : a string of form "MonthXX" where XX is the two digit month (1 --> Month01)'''

        return f'Month{month_int:02d}'

def yr_to_yrstr(sector,year,bau_or_covid):
        '''Converts a year to a string, and adds "COVID" or "BAU" to onroad/offroad data for 2020
        
        Args:
        sector (str) : the sector string (example = area_OG)
        year (int) : the year integer
        
        Returns:
        year_str (str) : the year in string format and bau or covid if needed
        '''

        if sector not in ['area_onroad_gasoline','area_onroad_diesel','area_offroad']: #if it's not an onroad/offroad source, no bau or covid
            return str(year) #so just return the string of the int
        if year == 2020: #if it's 2020, there is a covid option for onroad/offroad
            year_str = str(year)+bau_or_covid #use the input bau_or_covid variable to get correct string of form 2020COVID or 2020BAU
        else:#if it's not 2020
            year_str = str(year) #just use the year string 
        return year_str

def check_space(path,excep_thresh='8Tb'):
    '''Checks the amount of space on a filesystem given a path, and raise an error if the amount of space is below the threshold
    
    Args:
    path (str) : path to check how much room there would be based on the filesystem
    excep_thresh (str) : threshold below which to raise an exception, defaults '8T'
    
    Returns:
    usage_details (tuple) : output of shutil.disk_usage(path)
    
    Raises:
    MemoryError : if the amount of space on the filesystem containing the path is below the excep_thresh
    '''

    disk_usage = shutil.disk_usage(path)
    total_bytes,used_bytes,free_bytes = int(disk_usage.total), int(disk_usage.used), int(disk_usage.free)
    thresh_bytes = human_to_bytes(excep_thresh)
    if free_bytes < thresh_bytes:
        print(f'Path: {path}')
        print(f'Total space: {bytes_to_human(total_bytes)}')
        print(f'Used space: {bytes_to_human(used_bytes)}')
        print(f'Free space: {bytes_to_human(free_bytes)}')
        print(f'Required space (a user defined property): {excep_thresh}\n')
        raise MemoryError(f'There is less than {excep_thresh} in the destination filesystem. Ensure there is enough room for your data download, or change the threshold')
    else:
        print(f'Sufficient space detected in {path}')
        print(f'Found {bytes_to_human(free_bytes)} free space')
        return disk_usage

def human_to_bytes(human_readable, units=['b', 'Kb', 'Mb', 'Gb', 'Tb', 'Pb']):
    '''Makes a human readable storage size into a bytes integer
    #TODO Different definitions for this and the other direction (1000 bytes per kb vs 1024)
    
    Args:
    human_readable (str) : a human readable string representing space (something like 2.1Gb, or 5Tb)
    units (list) : list of strings defining the units 
    
    Returns:
    bytes (int) : integer representing how many bytes the human readable string is
    '''

    # Split the input string into value and unit
    value_str = ''.join([ch for ch in human_readable if ch.isdigit() or ch == '.'])
    unit_str = ''.join([ch for ch in human_readable if not ch.isdigit() and ch != '.']).strip()
    value = float(value_str)# Convert value string to a float
    
    # Find the unit in the provided units list
    try:
        unit_index = units.index(unit_str)
    except ValueError:
        raise ValueError(f"Invalid unit '{unit_str}' in input string.")
    
    multiplier = 1024 ** unit_index # Calculate the multiplier based on the unit index
    bytes = int(value * multiplier)   # Convert the human-readable value to bytes    
    return bytes

def bytes_to_human(bytes,units=['b','Kb','Mb','Gb','Tb','Pb']):
    '''Converts an integer bytes into a human readable string notation
    #TODO Different definitions for this and the other direction (1000 bytes per kb vs 1024)
    
    Args: 
    bytes (int) : integer bytes
    units (list) : unit convention for factors of 10 
    
    Returns (str) : the human readable string
    '''

    return str(bytes) + units[0] if bytes < 1024 else bytes_to_human(bytes>>10, units[1:])

def listdir_visible(path,add_path = False):
    '''Function to list only "visible" files/folders (not starting with a period)
    
    Args:
    path (str) : the filepath to list elements of
    
    Returns:
    vis_list (list) : list of files within the input path not starting with a period
    '''
    
    vis_list = [f for f in os.listdir(path) if not f.startswith('.')]
    if add_path:
        vis_list = [os.path.join(path,f) for f in vis_list]
    return vis_list

def sanity_check(og_ds,regridded_ds,regridded_ds_cellarea,species):
    '''A sanity check using daily sums of a single species to see if the total before and after regridding is close
    
    Args:
    og_ds (xr.Dataset) : the original dataset loaded from "base" noaa_csl data in LCC 
    regridded_ds (xr.Dataset) : the regridded dataset in WGS
    regridded_ds_cellarea (xr.Dataset) : a dataset with a cell_area data var giving the cell area in m2 for each cell in the regridded ds. from a loaded nc created by cdo gridarea 
    species (str) : the species to check
    
    Returns:
    perc_diff (float) : the difference (in percent) between the daily sum of all grid cells in the original versus the regridded dataset
    '''
    
    original_sum = float(og_ds[species].sum().values)
    #print('original_sum = ', original_sum)
    regrid_daysum = regridded_ds[species].sum(dim='utc_hour')
    regrid_daysum_absolute = regrid_daysum * regridded_ds_cellarea['cell_area']
    regrid_sum = float(regrid_daysum_absolute.sum().values)
    #print('regrid_sum = ', regrid_sum)
    try:
        perc_diff = abs(regrid_sum-original_sum)/((regrid_sum+original_sum)/2)*100
    except ZeroDivisionError:
        print(f'{species}: Sums are 0')
        return 0
    print(f'{species} sum diff = {round(perc_diff,3)}%')
    return perc_diff

def slice_extent(ds,extent):
    '''Slice a dataset to a bounding box defined by the extent argument
    
    Args:
    ds (xr.DataSet) : the dataset to process
    extent (dict) : a dictionary with 4 elements defining the bounding box -- must include 'lat_min', 'lat_max', 'lon_min', 'lon_max' -- inclusive.

    Returns:
    ds (xr.DataSet) : the dataset clipped to within the extent range
    '''

    try:
        grid_type = ds.attrs['grid_type']
    except:
        grid_type = 'area'
    if grid_type == 'area': #if it's an area grid
        ds = ds.sel(lat=slice(extent['lat_min'],
                              extent['lat_max']),
                    lon=slice(extent['lon_min'],
                              extent['lon_max'])) #we can slice on the lat lon coordinates
    elif grid_type == 'point': #if it's a point grid
        ds = ds.where(((ds.lat>=extent['lat_min'])&
                       (ds.lat<=extent['lat_max'])&
                       (ds.lon>=extent['lon_min'])&
                       (ds.lon<=extent['lon_max'])).compute(),drop=True) #we have to select the individual points using .where, and compute them. The drop term gets rid of NA's not in the extent
    else: 
        raise TypeError('Did not recognize the grid type, unsure how to slice')
    
    return ds 

def ncount_satsunwkd(year,month):
    '''Gets the number of saturdays, sundays and weekdays in a given month/year
    
    Args:
    year (int) : the year
    month (int) : the month, as an integer
    
    Returns:
    sat_count (int) : number of saturdays
    sun_count (int) : number of sundays
    weekd_count (int) : number of weekdays
    '''

    num_days_in_month = calendar.monthrange(year,month)[1]
    month_str = f'{month:02d}'

    dow_ints = list(pd.date_range(start=f'{year}-{month_str}-01',end=f'{year}-{month_str}-{num_days_in_month}').weekday)
    sat_count = len([ dow for dow in dow_ints if dow == 5 ])
    sun_count = len([ dow for dow in dow_ints if dow == 6 ])
    weekd_count = len([ dow for dow in dow_ints if dow < 5 ])
    return sat_count,sun_count,weekd_count

def get_daytype(dow):
    '''Finds the day type from the numeric day of week
    
    Args:
    dow (int) : day integer of the week Monday = 0, Sunday = 6

    Returns:
    day_type (str) : as satdy, sundy, or weekdy
    '''
    if dow == 5:
        day_type = 'satdy'
    elif dow == 6:
        day_type = 'sundy'
    elif dow <5:
        day_type = 'weekdy'
    else:
        raise ValueError(f'Invalid day of week integer {dow} (should be between 0 and 6)')
    return day_type

def get_point_sum_in_gc(row,point_df,species):
    '''Calculates the sum of all of the point sources within a gridcell, according to a row with lat min/max and lon min/max repping that gridcell
    
    Args:
    row (pd.DataFrame row) : a row or dictionary containing lat_min, lat_max, lon_min, lon_max -- the extent in which we want to sum
    point_df (pd.DataFrame) : the point_df we want to get the sum within the extent of 
    species (str) : species, should be one of the columns of point_df
    
    Returns:
    (float) : the sum of all of the point sources for the in put species within the gridcell defined by the row
    '''

    subpointdf = point_df.loc[(point_df['lat'] >= row['lat_min'])&
                 (point_df['lat'] <= row['lat_max'])&
                 (point_df['lon'] >= row['lon_min'])&
                 (point_df['lon'] <= row['lon_max'])]
    return subpointdf[species].sum()

def get_cellsize(ds):
    '''Calculates the cell size for lat and lon by subtracting the second element from the first (assumes equal grid spacing throughout)
    
    Args:
    ds (xarray.Dataset) : dtaset with lat and lon coordinates
    
    Returns:
    lat_size : the size of a grid cell in latitude (degrees)
    lon_size : the size of a grid cell in longitude (degrees)
    '''

    lat_size = float(ds['lat'][1]-ds['lat'][0]) #subtract the second element from the first
    lon_size = float(ds['lon'][1]-ds['lon'][0])
    return lat_size,lon_size

def pointdf_to_ds(point_df,area_ds,species):
    '''Transforms a point source dataframe into a dataset on the same grid spacing as an input area dataset
    
    Args:
    point_df (pd.DataFrame) : dataframe with lat, lon and species that you wnat to transform into a dataset
    area_ds (xarray.Dataset) : the form of the dataset you want to transform into -- will use this as the grid spacing
    
    Returns 
    point_ds (xarray.Dataset) : a dataset on the same grid spacing as the area dataset with all point sources in each gridcell SUMMED
    '''
     
    lat_size,lon_size = get_cellsize(area_ds) #get the grid spacing
    if len(point_df.columns) == 3: #find the species, this assumes there are only 3 columns TODO deal with more columns
        for col in point_df.columns:
            if (col=='lat')|(col=='lon'):
                continue
            else:
                species = col
    df = area_ds.to_dataframe(name=species).reset_index()[['lat','lon']] #transoform the template ds into a df
    df['lat_min'] = df.apply(lambda row: row['lat']-lat_size/2,axis = 1) # each min/max is calculated using the cell sizes
    df['lat_max'] = df.apply(lambda row: row['lat']+lat_size/2,axis = 1)
    df['lon_min'] = df.apply(lambda row: row['lon']-lon_size/2,axis = 1)
    df['lon_max'] = df.apply(lambda row: row['lon']+lon_size/2,axis = 1)

    df[species] = df.apply(lambda row: get_point_sum_in_gc(row,point_df,species),axis = 1) #add a column by summing all points in the gridcell
    point_ds = df[['lat','lon',species]].set_index(['lat','lon']).to_xarray() #make it a dataset
    return point_ds

class Base_CSL_Handler:
    '''This class is built to handle the file storage and naming conventions for the "base" NOAA CSL inventory data, as downloaded and 
    orgainized using data_download.py. These are the raw netcdf files obtained from the NOAA servers'''
    
    def __init__(self,base_path,bau_or_covid):
        self.base_path = base_path #location of the base data directory
        if bau_or_covid not in ['BAU','COVID']: #define if you want to use "business as usual" or "covid" for 2020 traffic data
            raise ValueError('The bau_or_covid input must be one of "BAU" or "COVID"')
        self.bau_or_covid = bau_or_covid

    def load_fullday_nc(self,full_sector,year,month,day_type,chunks = {}):
        '''Loads a full day of data from two 12-hour nc files in a day folder
        
        Args:
        full_sector (str) : name of the sector (like area_OG)
        year (int) : int or string of year
        month (int) : integer month
        day_type (str) : either "satdy", "sundy", or weekdy
        chunks (dict) : a dictionary defining chunks to use in dask. defaults to no chunks, example chunks = {'Time':1}
        
        Returns:
        full_ds (xarray.DataSet) : xarray dataset concatenated on the "utc_hour" dimension for both 00 and 12 day halves
        '''

        ds_00 = self.load_fmt_single_nc(full_sector,year,month,day_type,'00',chunks=chunks) #load the first half
        ds_12 = self.load_fmt_single_nc(full_sector,year,month,day_type,'12',chunks=chunks) #load the second half
        full_ds = xr.concat([ds_00,ds_12],dim = 'utc_hour',data_vars = 'minimal') #combine them
        return full_ds

    def load_fmt_single_nc(self,full_sector,year,month,day_type,hour_start,chunks={}):
        '''Loads and formats a single nc file into an xarray dataset
        
        Args:
        full_sector (str) : name of the sector (like area_OG)
        year (int) : int or string of year
        month (int) : integer month
        day_type (str) : either "satdy", "sundy", or weekdy
        hour_start (str) : either "00" for file starting at 0utc or "12" for file starting at 12utc
        chunks (dict) : a dictionary defining chunks to use in dask. defaults to no chunks, example chunks = {'Time':1}

        Returns:
        ds (xarray.DataSet) : the loaded and formatted dataset for the single nc file        
        '''

        nc_fpath = self.get_full_fname(full_sector,year,month,day_type,hour_start) #get the full path of the nc file
        ds = xr.open_dataset(nc_fpath,chunks=chunks) #open the dataset with one time chunk

        #Define attributes from the loaded file that aren't natively in the nc metatdata
        if full_sector.startswith('area'):
            grid_type = 'area'
        elif full_sector.startswith('point'):
            grid_type = 'point'
        else:
            raise ValueError(f'The input {full_sector} for full_sector does not start with "area" or "point". Something wrong here.')
        sector_id = '_'.join(full_sector.split('_')[1:]) #the full_sector without the "point" or "area" prepend
        
        #set the attributes
        ds.attrs['grid_type'] = grid_type
        ds.attrs['sector_id'] = sector_id
        ds.attrs['year'] = year
        ds.attrs['month'] = month
        ds.attrs['day_type'] = day_type
        ds.attrs['nc_fpath'] = nc_fpath
        if (year == 2020) & (sector_id in ['onroad_gasoline','onroad_diesel','offroad']): #if it's covid changable, add the attr
            ds.attrs['bau_or_covid'] = self.bau_or_covid

        #Deal with time dimension
        new_timedim_name = 'utc_hour' #new name of time dimension for clarity
        ds = self.rename_timedim(ds,newname=new_timedim_name) #rename the time dimension
        if hour_start == '12': #if it's a "12z" file, need to add 12 to the time dimension as when loaded, they are [0:12] (indicating hours since 12). 
            ds[new_timedim_name] = ds[new_timedim_name]+12
        else: #if its a '00' file, the indicies are correct and 0 = 00 hours , 1 = 01 hours etc. 
            ds[new_timedim_name] = ds[new_timedim_name]

        #In some of the sectors, there is an unkown datavariable called "Times" which is confusing. I just delete it if it's there. 
        if 'Times' in list(ds.data_vars):
            ds = ds.drop_vars('Times')

        for data_var in list(ds.data_vars.keys()): #for all of the data variables
            if data_var == 'NOX': #nox is a special case with a weird unit
                ds[data_var].attrs['units'] = 'metric_Ton(NO2equiv) hr^-1' #get rid of the space between NO2 and equiv so the unit conveter works
            if grid_type == 'area':
                ds[data_var].attrs['units'] = ds[data_var].attrs['units'] + ' gridcell^-1' #add ' gridcell^-1' to the units for clarity. all NOAA CSL base data is absolute unites i.e. "per gridcell"

        return ds

    def get_full_fname(self,full_sector,year,month,day_type,hour_start):
        '''Gets the full path and filename for a given base data file
        
        Args:
        full_sector (str) : name of the sector (like area_OG)
        year (int) : int or string of year
        month (int) : integer month
        day_type (str) : either "satdy", "sundy", or weekdy
        hour_start (str) : either "00" for file starting at 0utc or "12" for file starting at 12utc
        
        Returns:
        day_fullpath (str) : full path to the data file defined by the inputs. If no files, returns None
        '''

        year_str = yr_to_yrstr(full_sector,year,self.bau_or_covid) #make the year a string
        month_str = month_int_to_str(month) #make the month a string of the form "MonthXX"
        day_folder = os.path.join(self.base_path,full_sector,year_str,month_str,day_type) #get the day folder
        day_files = self.list_day_files(day_folder) #list all the files in the day folder
        for day_fname in day_files: #loop through files
            day_fullpath = os.path.join(day_folder,day_fname) #get the full path
            if self.hour_matches(day_fname,hour_start): #if the file matches the hour, that's the one we want
                return day_fullpath #so return it
            else: #otherwise
                continue #keep looping through the files in the day folder
        print("No matching files found") #if it gets here, there were no files that matched all of the input parameters
        return None #so we return none
    
    def hour_matches(self,fname,hour_start):
        '''Check that the hour start matches in the filename
        
        Args:
        fname (str) : the name of the file (not the full path)
        hour_start (str) : the hour start to check if it is in the file
        
        Returns:
        (bool) : True if the hour start is found in the fname, false if not
        '''

        split_fname = fname.split('_') #split the fname by underscore
        fname_hr = split_fname[1][0:2] #the first two characters of the split filename's second element is the hour in the filename
        if fname_hr not in ['00','12']: #it should be one of these two
            print(fname,fname_hr) #if it isnt, print it and raise an error
            raise ValueError('The "start hour" in the filename should be "00" or "12" -- something went wrong here')
        if hour_start == fname_hr: #ir the hour start matches
            return True #return true
        else: #otherwise
            return False #return false
        
    def list_day_files(self,day_folder):
        '''Lists the files in a day folder and checks that there are the right number (2)
        
        Args: 
        day_folder (str) : full filepath to the "day" folder (ex ~/data/area_OG/2019/Month01/satdy)
        
        Returns:
        files (list) : list of files in the day folder
        
        Raises:
        Exception : this happens when the number of files in this folder is not 2 -- something may be wrong with the data organizatino
        '''
        
        files = os.listdir(day_folder) #lists the files in the folder
        if len(files)!=2: #check that there are exactly two
            raise Exception(f'There are more than 2 nc files in {day_folder}. Something might be wrong here')
        return files
    
    def rename_timedim(self,ds,newname = 'utc_hour'):
        '''Renames the time dimension
        
        Args:
        ds (xarray.DataSet) : an xarray dataset with a time dimension
        newname (str, optional) : the new name of the time dimension. defaults to 'utc_hour' as defined in base data
        
        Returns:
        ds (xarray.DataSet) : the xarray dataset with the changed time dimension
        '''

        og_time_dim = self.get_og_timedim_name(ds) #get the original name of the time dimension in the dataset (they can differ so need a function)
        ds = ds.rename_dims({og_time_dim:newname}) #rename the time dimension from the original to the new
        return ds

    def get_og_timedim_name(self,ds):
        '''Gets the origial name for the time dimension in a dataset. In different data sets there are different labels for time, so we want to find what it is
        
        Args:
        ds (xarray.DataSet) : the dataset to find the time dimension in 
        
        Returns:
        time_dim (str) : the string name of the time dimension in the input dataset
        
        Rasies:
        ValueError : if there are multiple time dimensions or no time dimensions detected
        '''
         
        dims = list(ds.dims) #list the dims of the dataset
        time_dims = [] #initialize the time dim list 
        for dim in dims: #loop through the dims
            if dim.lower() in ['time','times']: #it could be a combination of Time, Times, time or times
                time_dims.append(dim) #if it fits the criteria, add it to the list
        if len(time_dims) > 1: 
            raise ValueError("There appears to be more than one labeled time dimension in the original dataset")
        elif len(time_dims) == 0:
            raise ValueError("Could not detect a time dimension in the input dataset")
        else:
            time_dim = time_dims[0] #should be the only one
        return time_dim
    
class CSL_Unit_Converter:
    '''
    Class to handle unit conversino in the CSL datasets. Should work on both "base" data (prior to regridding), and "regridded" data. Since
    we need to convert absolute units (per gridcell) to standard m^-2 before regridding, this is used prior to regridding. It can also be used
    after regridding to convert mass and time units
    '''

    def __init__(self):
        #self.species_details = self.load_species_details(species_details_fullpath)
        pass

    def absolute_to_flux(self,ds,species_list = 'all',grid_cell_area = 16000000,drop_species = True):
        '''Converts from absolute uints (per grid cell) to flux units (per meter^2)
        
        Args:
        ds (xr.DataSet): an xarray dataset with absolute units (per grid cell)
        species_list (str, list) : a list of species (data_vars) in the ds that you want to convert. defaults to 'all' 
        grid_cell_area (int) : the area in meters^2 of a grid cell. For NOAA CSL area emissions, in Lambert Conformal, the grid cells are 
                                approximately 4km x 4km = 4000m x 4000m = 16000000m^2
        drop_species (bool) : True if you want only the input species_list to be returned in the ds, if false, will return all species with only the species_list ones converted
        
        Returns:
        ds (xr.DataSet): an xarray dataset with flux units (meter^-2)
        '''

        if species_list == 'all': # if we want to convert all the variables
            species_list = self.get_data_vars(ds) #set the species list to the entire list of data_vars
        if species_list == None: #if there is no species, it's probably a data array instead of a dataset. 
            raise TypeError('Input was likely a data array not a dataset. Need a dataset.')
        
        out_ds = ds.copy() #use the original dataset as the starting point
        for species in species_list: #loop through all of the species to convert
            old_unit_list = self.unit_finder(ds[species]) #find the old units of the species
            if 'gridcell^-1' not in old_unit_list: #if there arene units of gridcell^-1, this unit conversion shouldn't be done!
                raise ValueError(f'The dataset did not have units in gridcell^-1 for species {species}')
            
            out_ds[species] = out_ds[species]/grid_cell_area #key step, divide each value by the grid cell area in m^2 to get to m^-2

            out_ds[species].attrs = ds[species].attrs #keep the attributes
            new_units = self.unit_replacer(old_unit_list,'gridcell^-1','meters^-2') #find the new units 
            out_ds[species].attrs['units'] = new_units #replace the units attribute with the new units
        if drop_species: #if you wanted to only get the species you converted
            out_ds = out_ds[species_list] #subset to them
        return out_ds
    
    def mass_conversion(self,da,from_unit='metric_Ton',to_unit='grams'):
        '''Converts between units of mass
        
        Args:
        ds (xr.DataSet): an xarray dataset with an attribute "units" that has a mass element
        from_unit (str): the unit to convert from
        to_unit (str): the unit to convert to
        
        Returns:
        out_ds (xr.DataSet): an xarray dataset with units of grams
        '''

        mass_unit_to_grams= {
            'metric_Ton':1E6,
            'metric_Ton(NO2equiv)':1E6,
            'kilograms':1E3,
            'grams':1
        }
        
        unit_list = self.unit_finder(da)
        if from_unit not in unit_list:
            raise ValueError(f'Trying to convert metric_Ton to grams, but no {from_unit} in units attribute')
        if from_unit not in list(mass_unit_to_grams.keys()):
            raise ValueError(f'Could not find {from_unit} in conversion dictionary -- options are {list(mass_unit_to_grams.keys)}')
        if to_unit not in list(mass_unit_to_grams.keys()):
            raise ValueError(f'Could not find to_unit={to_unit} in the conversion dictionary -- options are {list(mass_unit_to_grams.keys)})')
        

        gram_converter = mass_unit_to_grams[from_unit]
        if to_unit == 'grams':
            out_converter = gram_converter
        else:
            out_converter = gram_converter/mass_unit_to_grams[to_unit]
        
        out_da = da * gram_converter
        out_da.attrs = da.attrs #keep the attributes
        new_units = self.unit_replacer(unit_list,from_unit,to_unit) #find the new units 
        out_da.attrs['units'] = new_units #replace the units attribute with the new units
        return out_da
    
    def perhr_to_persec(self,ds):
        '''Converts from uints per hour to units per second
        
        Args:
        ds (xr.DataSet): an xarray dataset with units in unit/hour
        
        Returns:
        ds (xr.DataSet): an xarray dataset with unit/second
        '''
        
        s_per_h = 3600
        return ds/s_per_h
    
    def unit_finder(self,da):
        '''Finds the units of a data array
        
        Args:
        da (xarray.DataArray) : the data array to find the units of 
        
        Returns:
        unit_list (list) : the units in list form, each unit is its own element
        '''

        unit_string = da.attrs['units'] #get the units from the attributes
        unit_list = unit_string.split() #split it on spaces
        return unit_list
    
    def unit_replacer(self,unit_list,from_unit,to_unit):
        '''Replace a unit in the unit_list with desired new unit
        
        Args:
        unit_list (list) : units in list form
        from_unit (str) : the unit to be replaced
        to_unit (str) : the new unit string to be put in
        
        Returns:
        new_units (str) : the new units, in string form
        '''

        if from_unit not in unit_list:
            raise ValueError(f'from_unit ({from_unit}) is not in the unit list ({unit_list})')

        new_unit_list = [] #initialize the new units list
        for unit in unit_list: #loop through units
            if unit == from_unit: #if the unit matches the from unit
                new_unit_list.append(to_unit) #replace it with the to_unit
            else: #otherwise it isnt the from_unit
                new_unit_list.append(unit) #so keep it the same
        new_units = ' '.join(new_unit_list) #join them on space to get the string
        return new_units
    
    def get_data_vars(self,ds):
        '''Gets a list of the data variables in a dataset 
        
        Args: 
        ds (xarray.DataSet): xarray dataset we want to find the data variables of 
        
        Returns:
        data_var_list (list) : a list of strings representing the data_vars of ds
        '''

        data_var_list = list(ds.data_vars.keys())
        return data_var_list

class CSL_Regridder:
    '''A class to regrid the noaa_csl base data from lambert conformal coordinates into latitude/longitude coordinates'''

    def __init__(self,regrid_inputs):
        ''' 
        Args: 
        regrid_inputs (inst RegridInputs) : an instance of the RegridInputs class
        '''

        self.inputs = regrid_inputs

    def regrid_ds(self,ds,**kwargs):
        '''Regrid the input dataset. Creates a regridder if there is not one present
        
        Args:
        ds (xr.Dataset) : the xarray dataset to regrid, should be loaded and formated nc file from base noaa csl
        **kwargs : arguments that can be passed to xe.Regridder like grid_out, method, weights_path, or input_dims

        Returns
        regridded_ds (xr.Dataset) : dataset regridded in lat/lon
        '''

        try: #if a regridder exists in this class, just use that
            self.regridder
        except AttributeError: #if there isn't one yet, create it
            self.regridder = self.create_regridder(ds,**kwargs)
        
        #print('Regridding')
        regridded_ds = self.regridder(ds,keep_attrs = True)#regrid the dataset
        old_attrs = list(regridded_ds.attrs.keys())
        kept_attrs = {} #select the attributes to keep (don't want to keep attributes from the old grid)
        for attr in old_attrs:
            if attr not in ['grid_type','sector_id','year','month','day_type','nc_fpath']: #these are the attributes that go with the dataset, not the grid
                del regridded_ds.attrs[attr] #delete the old attributes
            else:
                continue

        return regridded_ds


    def create_regridder(self,ds,**kwargs):
        '''Creates the regridder object that can be used to regrid datasets
        
        Args:
        ds (xarray.DataSet) : the input dataset to use for the "in_grid". this should be "base" data
        weightsfile_or_create (str, optional) : if 'create', the weights in the regridder will be created using the grid_in and grid_out. If this variable
                                                is a full filepath, it will use that filepath as the saved weights instead of creating anew
        **kwargs : can override using self.inputs by adding a keyword arg matching the xe.Regridder (grid_out,method,input_dims)
        '''
        print('Creating regridder')
        grid_in = self.create_ingrid(ds) #create the input grid from the dataset
        grid_out = self.kwargs_or_inputs(kwargs,'grid_out') #get the grid out from kwargs first then inputs
        method = self.kwargs_or_inputs(kwargs,'method') #get the method from kwargs first then inputs
        input_dims = self.kwargs_or_inputs(kwargs,'input_dims') #get the input dims from kwargs firs then inputs
        weightsfile_or_create = self.kwargs_or_inputs(kwargs,'weights_file')

        if weightsfile_or_create == 'create': #if we're creating it from the in and out grids 
            regridder = xe.Regridder(grid_in,grid_out,method=method,input_dims=input_dims) #don't reuse weights, and just create it
        else: #if not, we're expecting a filepath 
            weights_path = self.kwargs_or_inputs(kwargs,'weights_path') #path to the weights 
            weights_full_filepath = os.path.join(weights_path,weightsfile_or_create) #full path to the weights file
            if not os.path.isfile(weights_full_filepath): #if there isn't a file there, something went wrong
                raise FileExistsError(f'No weights file found at {weights_full_filepath}')
            
            #otherwise, create the regridder with the saved weights file
            regridder = xe.Regridder(grid_in,grid_out,method=method,input_dims=input_dims,reuse_weights=True,filename=weights_full_filepath)

        regridder.ds_attrs = ds.attrs #add the dataset attributes to the regridder attributes to keep track of things
        return regridder

    def create_ingrid(self,ds):
        '''Creates the grid that will be input into the regridder from a "base" noaa_csl dataset. Adapted from Colin Harkins at NOAA
        
        Args:
        ds (xarray.DataSet) : the xarray dataset to use as the grid_in
        
        Returns:
        grid_in (dict) : a grid_in type dictionary that can be fed to xe.Regridder
        '''

        #create the projection transformers. 'wgs_to_lcc' converts EPSG4326 to lambert conformal, 'lcc_to_wgs' does the opposite
        wgs_to_lcc, lcc_to_wgs = self.create_transformers(ds)  

        # Calculate the easting and northings of the domain center point
        e,n = wgs_to_lcc.transform(ds.CEN_LON, ds.CEN_LAT) #use the attribributes to transform lat lons defined in the ds as center to lcc

        # Grid parameters from the dataset
        dx, dy = ds.DX, ds.DY 
        nx, ny = ds.dims['west_east'], ds.dims['south_north']

        # bottom left corner of the domain
        x0 = -(nx-1) / 2. * dx + e
        y0 = -(ny-1) / 2. * dy + n

        # Calculating the boundary X-Y Coordinates
        x_b, y_b = np.meshgrid(np.arange(nx+1) * dx + x0 -dx/2, np.arange(ny+1) * dy + y0 -dy/2)
        x_bc, y_bc = lcc_to_wgs.transform(x_b, y_b)

        #define the input grid
        grid_in = {
            'lat':ds['XLAT'].values,
            'lon':ds['XLONG'].values,
            'lat_b':y_bc,
            'lon_b':x_bc
        }

        return grid_in

    def create_transformers(self,ds):
        '''Create transformers for going from lambert conformal to WGS and vice versa
        
        Args:
        ds (xarray.dataset) : use the dataset to define the lcc projection, as noaa_csl "base" data is in LCC
        
        Returns:
        wgs_to_lcc (pyproj.Transformer) : a transformer object to transform from WGS coordinates to LCC coordinates
        lcc_to_wgs (pyproj.Transformer) : a transformer object to transform from LCC coordinates to WGS coordinates
        '''

        proj4_str = self.proj4_from_ds(ds) #get the proj4 string from the dataset
        lcc_crs = pyproj.CRS.from_proj4(proj4_str) #define the lcc coordinate reference system using the proj4 string
        wgs_crs = pyproj.CRS.from_epsg(4326) #define wgs coordinates as espg 4326
        wgs_to_lcc = pyproj.Transformer.from_crs(wgs_crs,lcc_crs,always_xy=True) #create one transformer
        lcc_to_wgs = pyproj.Transformer.from_crs(lcc_crs,wgs_crs,always_xy=True) #create the other transformer

        return wgs_to_lcc, lcc_to_wgs

    def proj4_from_ds(self,ds,map_proj=None,earth_rep="sphere"):
        '''Get a proj4 string defining a projection from a noaa_csl "base" dataset 
        
        Args:
        ds (xarray.Dataset) : a dataset to get the projection from
        map_proj (proj,optional) :the pyproj.Proj proj argument. allow the user to directly define the map proj, defaults to getting it from the dataset
        earth_rep (str) : change the earth representation of the input projection TODO, only can do "sphere" right now
        '''
         
        if map_proj == None: #if the user didn't define a map projection
            try: #try it so we can catch bad ones
                if ds.MAP_PROJ_CHAR == 'Lambert Conformal': #teh datasets from noaa_csl "base" data should have the "Lambert Conformal" attribute, so confirm this
                    map_proj = 'lcc' #if so, the map_proj code is lcc
                else: #if not, something unexpected happend
                    raise ValueError(f"Unknown map projection in ds.MAP_PROJ_CHAR: {ds.MAP_PROJ_CHAR}")
            except: #if it failed, something weird happened
                raise Exception('No map projection found')
            
        if earth_rep == 'sphere': #can only do spheres for now. this is how the noaa_csl "base" data is represented (like WRF) see references
            my_proj = pyproj.Proj(proj=map_proj, # projection type: Lambert Conformal Conic
                                lat_1=ds.TRUELAT1, lat_2=ds.TRUELAT2, # Cone intersects with the sphere
                                lat_0=ds.MOAD_CEN_LAT, lon_0=ds.STAND_LON, # Center point
                                a=6370000, b=6370000) # The Earth is a perfect sphere
        else:
            raise Exception("Haven't dealt with non spheres yet")
        
        return my_proj.to_proj4()
    
    def kwargs_or_inputs(self,kwargs,key):
        '''Get a value from either a dict of keyword arguments as preference, or if not use self.inputs
        
        Args:
        kwargs (dict) : keyword arguments
        key (str) : the key to extract the value of
        
        Returns:
        value : the value prioritizing kwargs, if not self.inputs
        '''
        
        if key in list(kwargs.keys()):
            return kwargs[key]
        else:
            return getattr(self.inputs,key)

    def save_regrid_weights(self,regridder,save_path=None,fname=None):
        '''Saves the regridded weights to a file
        
        Args:
        regridder (xe.Regridder) : an instance of the regridder class to save
        save_path (str, optional) : the path in which to store the weights. defaults to the inputs.weights_path
        fname (str, optional) : the name of the saved weights file. defaults to a filename based on attributes of the dataset
        '''

        if save_path is None:
            save_path = self.inputs.weights_path
        
        if fname is None:
            fname = f"lcc_to_latlon_{regridder.ds_attrs['sector_id']}_{regridder.ds_attrs['year']}_{regridder.ds_attrs['month']}_{regridder.ds_attrs['day_type']}.nc"
        
        regridder.to_netcdf(os.path.join(save_path,fname))

class Regridded_CSL_Handler:
    '''A class to handle NOAA CSL inventory data that has been regridded and organized by regrid_data.py'''

    def __init__(self,regridded_path,bau_or_covid='COVID'):
        '''Everything revolves around the "regridded_path", which determines sectors via the filenames contained within'''

        self.regridded_path = regridded_path
        self.sectors = self.get_sectors()
        self.bau_or_covid = bau_or_covid

    def get_sectors(self):
        '''Lists the sectors in the regridded data storage path'''

        sector_list = listdir_visible(self.regridded_path)
        sectors = {'area':[],'point':[]}
        for sector in sector_list:
            if 'area' in sector:
                sectors['area'].append(sector)
            elif 'point' in sector:
                sectors['point'].append(sector)
            else:
                raise ValueError(f"Unexpected sector type {sector}, not point or area.")
        return sectors

    def get_sector_subset_list(self,sector_subset):
        '''Gets a subset of the sectors which could be one, all, or some
        
        Args:
        sector_subset (str,list): "all" will return all sectors,  'point' will return point sectors, 'area' will return area sectors, a list will just return that list
        
        Returns:
        sector_subset_list (list) : list of sectors in the subset. 
        '''

        if sector_subset == 'all':
            sector_subset_list = []
            for k,v in self.sectors.items():
                sector_subset_list.extend(v)
            return sector_subset_list
        elif type(sector_subset)==str:
            return self.sectors[sector_subset]
        else:
            return sector_subset

    def get_days_in_range(self,dt1,dt2,day_types,sector_subset = 'all',add_path=True):
        '''Gets all filepaths to the day_type level that are within a datetime range
        
        Args:
        dt1 (datetime.date) : a date, datetime, etc to start the range (will only use year and month)
        dt2 (datetime.date) : a date, datetime, etc to end the range (will only use year and month)
        sectors (list) : list of sectors to include in the list
        day_types (list) : list of day types to include in the list
        add_path (bool, optional) : if true (default) it will add the regridded path to each element

        Returns:
        days_in_range (list) : list of paths to files that are within the date range and sector, day_types, etc. 
        '''

        dates_list = pd.date_range(dt1,dt2,freq = 'MS') #get a list of all the months between the dts
        sector_subset_list = self.get_sector_subset_list(sector_subset)
        days_in_range = []
        for date in dates_list:
            for sector in sector_subset_list:
                for day_type in day_types:
                    day_path = f'{sector}/{yr_to_yrstr(sector,date.year,self.bau_or_covid)}/{month_int_to_str(date.month)}/{day_type}'
                    if add_path:
                        days_in_range.append(os.path.join(self.regridded_path,day_path))
                    else:
                        days_in_range.append(day_path)
        return days_in_range
    
    def get_files_in_days(self,days_paths):
        '''Gets the files that exist in the paths
        
        Args:
        days_path (list) : list of paths to the days folders
        
        Returns:
        files (list) : list of files in those days' paths'''

        files = []
        for day_path in days_paths:
            files.extend(listdir_visible(day_path,add_path=True))
        return files

    def preprocess_regridded(self,ds,extent=None):
        '''Preprocesses the regridded dataset when loaded to add attributes needed for concatenation
        
        Args:
        ds (xr.DataSet) : the dataset to process
        extent (dict) : a dictionary with 4 elements defining the bounding box -- must include 'lat_min', 'lat_max', 'lon_min', 'lon_max'. These are inclusive. Defaults to "None" which will return the whole ds
        
        Returns 
        ds (xr.DataSet) : the dataset, with added coordinates taken from the attributes, sliced to the input extent
        '''

        grid_type = ds.attrs['grid_type']
        if grid_type == 'area':
            ds = ds.assign_coords(sector = 'area_'+ ds.attrs['sector_id']) #add back the area, was cut off in attributes for some reason
        elif grid_type == 'point':
            ds = ds.assign_coords(sector = 'point_'+ ds.attrs['sector_id']) #add back the point, was cut off in attributes for some reason

        ds = ds.assign_coords(day_type = ds.attrs['day_type']) #assign the day_type coordinate
        ds = ds.assign_coords(yr_mo=f"{ds.attrs['year']}-{ds.attrs['month']}") #assign the year/month coordinate
        ds = ds.expand_dims(dim=['sector','day_type','yr_mo']) #make the coordinates into dimensions
        if extent is not None:
            ds = slice_extent(ds,extent) #slice to the bounding box extent

        try:
            del ds.attrs['nc_fpath'] #we don't need this attribute -- it points to an old nc path. 
        except:
            pass

        return ds

def main():
    print(get_githash())

if __name__ == "__main__":
    main()
