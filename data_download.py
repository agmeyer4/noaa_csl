'''
Module: data_download.py
Author: Aaron G. Meyer (agmeyer4@gmail.com)
Last Updated: May 2, 2024
Description:
This module allows a user to easily download inventory data from the NOAA servers based on sector, year, and month. The download class
will find the correct url for getting the tar file, extract the tar file, format the filenames and directories in a consistent way, delete
or organize data specific to each sector, and remove the superfluous tar and extracted downloaded files.

More information can be found here: https://csl.noaa.gov/groups/csl7/measurements/2020covid-aqs/emissions/

**** 
A warning about storage space: One month of data from all sectors is approximately 131 Gigabytes. There is a function here that will alert
the user if there is less than a defined amount (default = 6Tb) of space on the filesystem where the data will end up (base_data_storage_path). 
The dataset from 2019 Month01 to 2021 Month08 is approximately 4.5 Terabytes. Downloading this entire dataset will take ~13 hours. 
****


TODO Don't download if the data is already there
'''

##################################################################################################################################
#Import Packages
import os
import subprocess
import shutil
import glob
import noaa_csl_funcs as ncf

##################################################################################################################################
# Define Classes
class NOAA_CSL_download:
    '''This class provides all functionality to retrieve, extract and organize NOAA CSL emission inventories'''

    base_noaa_url = 'https://csl.noaa.gov/groups/csl7/measurements/2020covid-aqs/data/emissions' #the base directory on the NOAA server for the CSL

    #The sector_details dictionary defines specifics about the different types of sectors including the name, the folder on the NOAA server, how the 
    #filenames are stored on the server, what years exist, the type of grid it is (area or point), and the path prefix template (to the level of the)
    #folder with satdy, sundy, weekdy) once the tar is extracted
    sector_details = {
        'area_onroad_gasoline':{
            'url_folder':'FIVE/onroad_gasoline',
            'url_fname_temp':'FIVE_onroad-gas_year_MonthXX.tar',
            'years':['2017','2018','2019','2020BAU','2020COVID','2021'],
            'extracted_preday_prefix':'MonthXX',
            'grid_type':'area'
        },
        'area_onroad_diesel':{
            'url_folder':'FIVE/onroad_diesel',
            'url_fname_temp':'FIVE_onroad-dsl_year_MonthXX.tar',
            'years':['2017','2018','2019','2020BAU','2020COVID','2021'],
            'extracted_preday_prefix':'MonthXX',
            'grid_type':'area'
        },
        'area_offroad':{
            'url_folder':'FIVE/offroad',
            'url_fname_temp':'FIVE_offroad_year_MonthXX.tar.gz',
            'years':['2018','2019','2020BAU','2020COVID','2021'],
            'extracted_preday_prefix':'MonthXX',
            'grid_type':'area'
        },
        'area_OG':{
            'url_folder':'NRT_area-point_emissions/FOG-NEI17_Oil_and_Gas',
            'url_fname_temp':'areaOG_emis_year_MonthXX.tar.gz',
            'years':['2019','2020','2021'],
            'extracted_preday_prefix':'wrk/d2/charkins/CSD_Emissions/conus4k/AREAyr2d_OG/MonthXX/AreaFOGnei17',
            'grid_type':'area'
        },
        'area_AG':{
            'url_folder':'NRT_area-point_emissions/area',
            'url_fname_temp':'areaAG_emis_year_MonthXX.tar.gz',
            'years':['2019','2020','2021'],
            'extracted_preday_prefix':'AreaAG',
            'grid_type':'area'
        },
        'area_Industry':{
            'url_folder':'NRT_area-point_emissions/area',
            'url_fname_temp':'areaIndustry_emis_year_MonthXX.tar.gz',
            'years':['2019','2020','2021'],
            'extracted_preday_prefix':'AreaIndustry',
            'grid_type':'area'
        },
        'area_Other':{
            'url_folder':'NRT_area-point_emissions/area',
            'url_fname_temp':'areaOtherArea_emis_year_MonthXX.tar.gz',
            'years':['2019','2020','2021'],
            'extracted_preday_prefix':'OtherArea',
            'grid_type':'area'
        },
        'area_VCP':{
            'url_folder':'NRT_area-point_emissions/area',
            'url_fname_temp':'areaVCP_emis_year_MonthXX.tar.gz',
            'years':['2019','2020','2021'],
            'extracted_preday_prefix':'AreaVCP',
            'grid_type':'area'
        },        
        'point':{
            'url_folder':'NRT_area-point_emissions/point',
            'url_fname_temp':'point_emis_year_MonthXX.tar.gz',
            'years':['2019','2020','2021'],
            'extracted_preday_prefix':'MonthXX',
            'grid_type':'point'
        }     
    }

    #Below defines the change in notation for point sources, keys are the new (local) names, values are the extracted point source sector names
    point_extract_structure = {
        'point_Other':'OtherPoint',
        'point_EGU':'PtEGU',
        'point_Industry':'PtIndustry',
        'point_OG':'PtOG',
        'point_VCP':'PtVCP'
    }
    
    def __init__(self,base_data_storage_path,bau_or_covid='COVID'):
        '''
        Args:
        base_data_storage_path (str) : path to where the sector data will be stored
        bau_or_covid (str) : for 2020 FIVE vehicle data, which version to use. Options are "BAU"-business as usual or "COVID"-covid
        '''

        self.base_data_storage_path = base_data_storage_path
        self.bau_or_covid = bau_or_covid

    def retrieve_format_data(self,sector,year,month):
        '''Main function to download NOAA CSL data for a specific sector, year, and month and format the directories nicely
        
        Args:
        sector (str) : sector to download, from the self.sector_details.keys() dictionary
        year (int) : integer year to download
        month (int) : integer month to download
        '''
        
        month_str = ncf.month_int_to_str(month) #convert the integer month into the month string used by NOAA
        if (year == 2020) & (sector in ['area_onroad_gasoline','area_onroad_diesel','area_offroad']): #in 2020, there are two versios of the FIVE inventory, choose the one defined in the class instatiation
            year_str = str(year)+self.bau_or_covid
        else:
            year_str = str(year) #get the year in string form
        print(f'Downloading, extracting, and organizing data for {sector} {year_str} {month_str}')
        if year_str not in self.sector_details[sector]['years']: #check that the data should be available
            raise ValueError(f'The year {year_str} is not included in the sector details for {sector}, it likely does not exist on the server')
        full_url = self.full_url_create(sector,year_str,month_str) #create the URL to get the tar file from
        tar_fname = os.path.split(full_url)[1] #get the name of the tar file
        os.makedirs(os.path.join(self.base_data_storage_path,'.temp'),exist_ok = True) #create a temporary path for downloading in the base data dir if necessary
        download_path = os.path.join(self.base_data_storage_path,'.temp') #define the download path as a temp path
        self.wget_file(full_url,download_path) #download the file
        self.extract_tar(os.path.join(download_path,tar_fname)) #extract it
        dir_dict = self.dir_dict_setup(sector,year_str,month_str,download_path) #define the directory organization for moving files 
        self.mv_extracted(dir_dict) #move the files to the correct locations
        self.sector_extra_org(sector,dir_dict) #do any extra sector-specific cleanup
        self.clear_folder(download_path) #remove everything from the temporary folder

    def full_url_create(self,sector,year_str,month_str):
        '''Creates the full url from which to get the file off NOAA's servers
        
        Args:
        sector (str) : the sector key
        year_str (str) : string of year
        month_str (str) : string of month

        Returns:
        full_url (str) : the full url of the tar file to be downloaded from NOAA
        '''

        sec_det = self.sector_details[sector] #get the details about the specific sector
        url_fname = self.url_fname_create(sec_det['url_fname_temp'],year_str,month_str) #create the filename
        full_url = os.path.join(self.base_noaa_url,sec_det['url_folder'],year_str,url_fname) #join them all
        if (sector in ['area_onroad_gasoline','area_onroad_diesel']) & (year_str=='2021'): #there is a .gz on the FIVE 2021 data
            full_url = full_url+'.gz'
        return full_url
    
    def wget_file(self,full_url,download_path):
        '''Download a file using wget to a specific path
        
        Args:
        full_url (str) : the full input to the wget function (file to get)
        download_path (str) : where the download should go
        '''

        print(f'Downloading {full_url} to {download_path}')
        command = ['wget',full_url,'-P',download_path] #create the wget command 
        proc = subprocess.Popen(command,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) #call the command
        proc.communicate() #show the stdout and sterr 

    def extract_tar(self,full_fname):
        '''Extract a tar file to the path it is already in 
        
        Args:
        full_fname (str) : full filepath to the tar file to be exctracted
        '''

        path,fname = os.path.split(full_fname) #get the path and the name
        print(f'Extracting {fname} into {path}')
        command = ['tar','xvf',os.path.join(path,fname),'-C',path] #define the tar extract command
        proc = subprocess.Popen(command,stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT) #call the command
        proc.communicate() #output

    def dir_dict_setup(self,sector,year_str,month_str,download_path):
        '''Sets up the directory dictionary to define where files should be moved from and two, handling specific prefixes
        
        Args:
        sector (str) : the sector
        year_str (str) : the string of the year
        month_str (str) : the formatted month string (MonthXX)
        download_path (str) : full path to where the tar files are downloaded
        
        Returns:
        dir_dict (dict) : dictionary containing the sectors and the "from" and "to" paths defining where files should come from and go
        '''
        yr2d = year_str[2:4]
        prefix_temp = self.sector_details[sector]['extracted_preday_prefix'] #get the prefix template for the sector
        prefix = ncf.replace_all_strs(prefix_temp,{'MonthXX':month_str,'yr2d':yr2d}) #get the real prefix by replacing the month (if necessary)
        dir_dict = {} #initialize the data dict
        if sector == 'point': #for the "point" sector, we must define different "subsectors" that will rise to top level sectors when moving
            for new_subsec, old_subsec in self.point_extract_structure.items():
                to_path = os.path.join(self.base_data_storage_path,new_subsec,year_str,month_str)
                try: 
                    os.makedirs(to_path)
                except:
                    print(f'Path exists at {to_path}')
                dir_dict[new_subsec] = {
                    'from_path':os.path.join(download_path,prefix,old_subsec),
                    'to_path':to_path
                }
        else: #if not point, then area
            to_path = os.path.join(self.base_data_storage_path,sector,year_str,month_str) #the to path is defined here
            try: #try to make the directory
                os.makedirs(to_path)
            except FileExistsError: #if it exists
                print(f'Path exists at {to_path}') #just tell us and continue
            dir_dict[sector] = { 
                'from_path':os.path.join(download_path,prefix), #from path is the download path plus the found prefix
                'to_path':to_path
            }#define the dir_dict
        return dir_dict

    def mv_extracted(self,dir_dict):
        '''Moves files from where they were extracted from the tar file into the place they should be stored, provided by the dir_dict
        
        Args:
        dir_dict (dict) : a dictionary with keys as the sector id and values as another dictionary, indicating where files should be moved from and moved to
        '''

        for sector_id,path_details in dir_dict.items(): #loop through the items in the dir dict
            for f in os.listdir(path_details['from_path']): #for all of the files in the "from path"
                mv_from = os.path.join(path_details['from_path'],f) #get the full from path
                mv_to = os.path.join(path_details['to_path']) #get the full to path
                shutil.move(mv_from,mv_to) #move the files

    def sector_extra_org(self,sector,dir_dict):
        '''Sector specific extra organization/file moving or deletion -- each sector has specific things that may need to be formated
        
        Args:
        sector (str) : the top level sector
        dir_dict (dict) : dictionary with the directory path details of the files
        '''

        if self.sector_details[sector]['grid_type'] == 'area': #if it's an area source
            path_with_days = dir_dict[sector]['to_path'] #the "path with days" or path with satdy, sundy, weekdy should just be the "to_path"
        if sector in ['area_AG']: #if the sector is area_AG
            to_remove = glob.glob(f'{path_with_days}/*/*.rds') #list unnecessary symlinked rds files
            to_remove.extend(glob.glob(f'{path_with_days}/*/.RData')) #add unnecessary symlinked .RData files to list
            for file in to_remove: #loop through files to remove
                os.remove(file) #remove the files
        if sector == 'point':
            to_remove = []
            for to_from in dir_dict.values():
                to_path = to_from['to_path']
                to_remove.extend(glob.glob(f'{to_path}/*/*gridded*'))
            for file in to_remove:
                os.remove(file)

    def clear_folder(self,path):
        '''Clears all files and directories from an input folder **danger, deletes files!** 
        
        Args:
        path (str) : path to a directory which is to be emptied (directory itself will remain)
        '''

        for f in os.listdir(path): #get all the files and folders in the path
            fullpath = os.path.join(path,f) #join the paths
            if os.path.isfile(fullpath): #if it is a file
                self.delete_file(fullpath) #delete the file
            if os.path.isdir(fullpath): #if its a directory
                self.delete_tree(fullpath) #delete the full tree

    def delete_file(self,full_path):
        '''Deletes a file **danger deletes files!**
        
        Args:
        full_path (str) : full path to a file to be deleted
        '''

        print(f'Deleting {full_path}')
        os.remove(full_path)
    
    def delete_tree(self,tree_path):
        '''Deletes a directory and all of its subcomponents **danger deletes files!**
        
        Args:
        tree_path (str) : the path to the directory which is to be completely deleted
        '''

        print(f'Deleting tree {tree_path}')
        shutil.rmtree(tree_path)
    
    def url_fname_create(self,url_fname_temp,year_str,month_str):
        '''Creates the filename in the NOAA server for a specific year/month given the template
        
        Args:
        url_fname_temp (str) : a string that has the template for the filename in the NOAA server
        year_str (str) : string of the year to replace
        month_str (str) : string of the month to replace
        
        Returns (str) : a string with the year and month replaced representing the filename in the NOAA server'''

        dic = {'year':year_str,'MonthXX':month_str}
        return ncf.replace_all_strs(url_fname_temp,dic)



##################################################################################################################################
# Main Function
def main():
    import time
    t1 = time.time()
    base_data_storage_path = '/uufs/chpc.utah.edu/common/home/lin-group9/agm/NOAA_CSL_Data/base' #where the data will be stored
    ncf.check_space(base_data_storage_path,excep_thresh='6Tb') #ensure there is enough space in the director
    ncd = NOAA_CSL_download(base_data_storage_path,bau_or_covid='BAU') #setup the downloader

    # Single sector, single month
    #ncd.retrieve_format_data('area_OG',2020,1)



    # All sectors, single month
    # for sector in ncd.sector_details.keys(): #loop through sectors
    #     if sector in ['area_onroad_diesel','area_onroad_gasoline','area_offroad','area_OG']:
    #         continue
    #     ncd.retrieve_format_data(sector,2021,5)
    

    ## All sectors, specific range
    year = 2020
    for month in range(2,13):
        for sector in ncd.sector_details.keys(): #loop through sectors
            if sector in ['area_onroad_diesel','area_onroad_gasoline','area_offroad']:
                ncd.retrieve_format_data(sector,year,month)

    # Full Download
    # for year in [2019,2020,2021]:
    #     for month in range(1,13):
    #         if (year==2021)&(month>8):
    #             continue
    #         for sector in ncd.sector_details.keys(): #loop through sectors
    #             ncd.retrieve_format_data(sector,year,month)

    t2 = time.time()
    print(f'total runtime (seconds) = {t2-t1}')

if __name__ == "__main__":
    main()