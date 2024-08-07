# NOAA CSL Combined Inventory
This python package is designed to deal with data from the NOAA CSL emissions inventory found at https://csl.noaa.gov/groups/csl7/measurements/2020covid-aqs/emissions/.   

To run the code, you will need an appropriate conda environment. This can be created using the environment file by running (within the main git folder): 
```
> conda env create --file=environment.yml
```

Key elements of this package are currently:     
1. data_download.py : Download all or a subset of the data from the link above, organize, and standardize the structure of the dataset. Note that downloading the entire dataset (through Aug 2021) will take ~20hrs and use ~5Tb of space. More details in download_data.py.

2. regrid_data.py : Regrid the downloaded data from Lambert Conformal Coordinates to lat/long (WGS) coordinates, at a default of 0.1&deg; grid spacing. 

3. noaa_csl_funcs.py : Useful functions to handle paths and organization, deal with base data inconsistencies, deal with units/attributes, sanity checks, regridding library, etc. 

4. regridding : Folder containing .nc files that are used in the regrid. 

5. slurm : Folder containing scripts and logs for running the download and regrid scripts on slurm

6. ipynbs : folder containing Jupyter notebooks used for illustrating examples, walking through code, debugging, and development. **A good place to start if you are familiar with using Jupyter is ipynbs/regridded_handling.ipynb. This will allow you to load data and play around** ***even if you aren't python savvy, take a look at these files directly in github to see what is possible***

For those on Utah CHPC, the data archives can be found at /uufs/chpc.utah.edu/common/home/lin-group9/agm/NOAA_CSL_Data. 

## Folder Structure of Data
```
├── NOAA_CSL_Data
│   ├── base
│   │   ├── sectors
│   │   │   ├── year 
│   │   │   │   ├── month 
│   │   │   │   │   ├── daytype 
│   │   │   │   │   │   ├── sector_00to12z.nc 
│   │   │   │   │   │   ├── sector_12to24z.nc 
│   ├── regridded
│   │   ├── sectors
│   │   │   ├── year 
│   │   │   │   ├── month 
│   │   │   │   │   ├── daytype 
│   │   │   │   │   │   ├── sector_regridded.nc 
```
Where "base" data is the downloaded and organized data from the CSL servers, "regridded" data is are the regridded (or standardized in the case of point sectors) files, and:
```
sectors: 
    area_ag = agriculture area sources (NEI2017)
    area_Industry = industry area sources (NEI2017)
    area_offroad = offroad vehicles area sources (FIVE)
    area_OG = oil and gas area sources (FOG and NEI2017)
    area_onroad_diesel = onroad diesel area sources (FIVE)
    area_onroad_gasoline  = onroad gasoline area sources (FIVE)
    area_Other  = other area sources (NEI2017)
    area_VCP = volatile chemical product area sources (NEI2017)
    point_EGU = energy generation units(NEI2017)
    point_Industry = industry (NEI2017)
    point_OG = oil and gas point sources (NEI2017)
    point_Other = other point sources (NEI2017)
    point_VCP = vcp point sources (NEI2017)

year:   
    2019
    2020 : all but onroad/offroad sectors
    2020BAU : for onroad and offroad sectors, business as usual
    2020COVID: for onroad and offroad sectors, covid
    2021

month: 
    format = Month01

daytype: 
    satdy = saturday
    sundy = sunday
    weekday = monday through friday
```

*A note about point vs area sectors: The point data is simply lat/lon pairs and emisisions, so no "regrid" needs to be done. To retain the most information (not lose spatial resolution), I simply combined the 12hr base data and standardized attributes, naming convention, etc in the "regridded" folder path. Part of my ongoing work is building a streamlined way to aggregate all of the regridded area and point sources.

The regridded area data is at its inherent resolution of approximately 4kmx4km. However, the point source data is precise. Thus, in a representation of emissions, we have sub-grid scale point data. Keeping it as such instead of incorporating into the grid allows for more sophisticated analysis on grid scaling applications. 

## Ongoing Work
The following ongoing development can be found in ./ipynbs/noaa_csl_dev.ipynb
* Streamline loading of different sets of regridded data (spatial boundaries, datetime ranges, etc). 
* Continue work on a "unit conversion" class to easily convert units for any of the species listed (especially HCXX variables). Tune to SAPRC molecular weights. 
* Improve visualization capabilities
# Contact and Acknowledgements
This work is being carried out under the direction of Dr. John C. Lin in the Land-Atmosphere Interactions Research (LAIR) group at the University of Utah, Atmospheric Sciences Department. 

Please contact me directly with any questions or comments at agmeyer4@gmail.com. 

Produced by Aaron G. Meyer, 2024
