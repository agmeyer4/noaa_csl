# NOAA CSL Combined Inventory
This python package is designed to deal with data from the NOAA CSL emissions inventory found at https://csl.noaa.gov/groups/csl7/measurements/2020covid-aqs/emissions/.   

To run the code, you will need an appropriate conda environment. This can be created using the environment file by running (within the main git folder): 
```
> conda env create --file=environment.yml
```

Key capabilities of this package are currently:     
1. Download all or a subset of the data from the link above, organize, and standardize the structure of the dataset. Note that downloading the entire dataset (through Aug 2021) will take ~20hrs and use ~5Tb of space. More details in download_data.py.
2. Regrid the downloaded data from Lambert Conformal Coordinates to lat/long (WGS) coordinates. Functions to deal with base data inconsistencies, deal with units/attributes, etc. including sanity checks. More information in regrid_base_data.py and noaa_csl_dev.ipynb. 
3. Begin to visualize and analyze regridded data. 
4. Develop useful functions for handling this large data archive, primarily using xarray. 

Walkthroughs can be found in noaa_csl_dev.ipynb, which has examples and descriptions for how the base data handler and regridders work. 

For those on Utah CHPC, the data archives can be found at /uufs/chpc.utah.edu/common/home/lin-group9/agm/NOAA_CSL_Data. 

# Contact and Acknowledgements
This work is being carried out under the direction of Dr. John C. Lin in the Land-Atmosphere Interactions Research (LAIR) group at the University of Utah, Atmospheric Sciences Department. 

Please contact me directly with any questions or comments at agmeyer4@gmail.com. 

Produced by Aaron G. Meyer, 2024
