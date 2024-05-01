import os
import subprocess
import shutil
import glob

class NOAA_CSL_download:
    base_noaa_url = 'https://csl.noaa.gov/groups/csl7/measurements/2020covid-aqs/data/emissions'
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
            'extracted_preday_prefix':'wrk/d2/charkins/CSD_Emissions/conus4k/AREA19_OG/MonthXX/AreaFOGnei17',
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

    point_extract_structure = {

    }
    
    def __init__(self,base_data_storage_path):
        self.base_data_storage_path = base_data_storage_path

    def retrieve_format_data(self,sector,year,month):
        month_str = self.month_int_to_str(month)
        year_str = str(year)
        if year_str not in self.sector_details[sector]['years']:
            raise ValueError(f'The year {year_str} is not included in the sector details for {sector}, it likely does not exist on the server')
        full_url = self.full_url_create(sector,year_str,month_str)
        tar_fname = os.path.split(full_url)[1]
        download_path = os.path.join(self.base_data_storage_path,'.temp')#,sector,year_str,month_str)
        self.download_tar(full_url,download_path)
        self.extract_tar(os.path.join(download_path,tar_fname))
        store_path = self.sector_dir_setup(sector,year_str,month_str)
        self.mv_extracted(sector,month_str,download_path,store_path)
        self.sector_extra_org(sector,store_path)
        self.delete_tar(os.path.join(download_path,tar_fname))

    def sector_dir_setup(self,sector,*args):
        if self.sector_details[sector]['grid_type'] == 'area':
            store_path = os.path.join(self.base_data_storage_path,sector,*args)
            try:
                os.makedirs(store_path)
            except FileExistsError:
                print(f'Path to {store_path} already exists')
                return store_path
        elif self.sector_details[sector]['grid_type'] == 'point':
            pass
        return store_path


    def sector_extra_org(self,sector,path_with_days):
        if sector in ['area_AG']:
            to_remove = glob.glob(f'{path_with_days}/*/*.rds')
            to_remove.extend(glob.glob(f'{path_with_days}/*/.RData'))
            for file in to_remove:
                os.remove(file)


    def mv_extracted(self,sector,month_str,download_path,store_path):
        print(f'Moving extracted day files to {store_path}')
        prefix_temp = self.sector_details[sector]['extracted_preday_prefix']
        prefix = replace_all_strs(prefix_temp,{'MonthXX':month_str})
        extracted_path = os.path.join(download_path,prefix)

        for f in os.listdir(extracted_path):
            shutil.move(os.path.join(extracted_path,f),store_path)

        top_level_prefix = prefix.split('/')[0]
        rm_tree = os.path.join(download_path,top_level_prefix)
        print(f'Removing tree {rm_tree}')
        shutil.rmtree(rm_tree)

    def delete_tar(self,full_tar_path):
        print(f'Deleting tar {full_tar_path}')
        command = ['rm',full_tar_path]
        proc = subprocess.Popen(command)
        proc.communicate()

    def extract_tar(self,full_fname):
        '''Extract a tar file to the path it is already in 
        
        Args:
        full_fname (str) : full filepath to the tar file to be exctracted
        '''

        path,fname = os.path.split(full_fname) #get the path and the name
        print(f'Extracting {fname} into {path}')
        command = ['tar','xvf',os.path.join(path,fname),'-C',path] #define the tar extract command
        proc = subprocess.Popen(command) #call the command
        proc.communicate() #output

    def download_tar(self,full_url,download_path):
        '''Download a tar file using wget to a specific path
        
        Args:
        full_url (str) : the full input to the wget function (file to get)
        download_path (str) : where the download should go
        '''

        print(f'Downloading {full_url} to {download_path}')
        command = ['wget',full_url,'-P',download_path] #create the wget command 
        proc = subprocess.Popen(command) #call the command
        proc.communicate() #show the stdout and sterr 

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
        return full_url
    
    def month_int_to_str(self,month_int):
        '''Converts a month integer to a two character string
        
        Args: 
        month_int (int) : the integer of the month (1,2,12,etc)
        
        Returns (str) : a string of form "MonthXX" where XX is the two digit month (1 --> Month01)'''

        return f'Month{month_int:02d}'
    
    def url_fname_create(self,url_fname_temp,year_str,month_str):
        '''Creates the filename in the NOAA server for a specific year/month given the template
        
        Args:
        url_fname_temp (str) : a string that has the template for the filename in the NOAA server
        year_str (str) : string of the year to replace
        month_str (str) : string of the month to replace
        
        Returns (str) : a string with the year and month replaced representing the filename in the NOAA server'''

        dic = {'year':year_str,'MonthXX':month_str}
        return replace_all_strs(url_fname_temp,dic)
    
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

if __name__ == "__main__":
    ncd = NOAA_CSL_download('/uufs/chpc.utah.edu/common/home/lin-group9/agm/NOAA_CSL')
    ncd.retrieve_format_data('area_AG',2019,1)
