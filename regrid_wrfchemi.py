# Uses _____ env
# from pyproj import Proj, transform
import numpy as np
import xarray as xr
import xesmf as xe
# import ESMF
import pyproj


class inputs():
    # Define all variables related to input and output and regions
    # File names to read in
    fn_base1 = 'wrfchemi_00z_d01'
    fn_base2 = 'wrfchemi_12z_d01'
    fn_ext = '.nc'
    
    # file name for save out
    fn_out = 'wrfchemi_regridded.nc'
    
    input_dir = '/data/colharkins/EI_noaa/VCP/VCP_EI/'
    out_dir = '/data/colharkins/EI_noaa/VCP/VCP_EI/'
    base_dir = '/data/colharkins/EI_noaa/VCP/VCP_EI/'
    
    days = ['weekdy','satdy','sundy']
    years = [2011]
    months = {'Month01':1,'Month02':2,
        'Month03':3,'Month04':4,
        'Month05':5,'Month06':6,
        'Month07':7,'Month08':8,
        'Month09':9,'Month10':10,
        'Month11':11,'Month12':12}

        # Converting Everything into kg/m^2/hour 
    MW = {'E_ALD': 44.05, 'E_CH4': 16,'E_CO' : 28,'E_CO2' : 44,'E_CSL' : 106.16,'E_ETH' : 28.05,
             'E_HC3': 44.1,'E_HC5' : 58.12,'E_HC8' : 58.12,'E_HCHO' : 30.031, 'E_HONO': 47, 'E_ISO': 68.12,
             'E_KET': 58.08, 'E_NH3': 17, 'E_NO': 30.01, 'E_NO2': 46.0055, 'E_OL2': 28.05, 'E_OLI': 42.08,
             'E_OLT' : 42.08,'E_ORA2': 74.08, 'E_SO2': 64, 'E_TERP': 136.23,'E_TOL': 92.14, 'E_UNID': 128.18,'E_XYL': 106.16}
    
    Aerosols = ['E_EC', 'E_NO3', 'E_ORG', 'E_PM25', 'E_PM10', 'E_SO4']
    
    reuse_weight_file = True
# End inputs class

def makeregridder(ds):
    #Calculating Projection from WRF LCC to WRF lat long space 
    wrf_proj = pyproj.Proj(proj='lcc', # projection type: Lambert Conformal Conic
                           lat_1=ds.TRUELAT1, lat_2=ds.TRUELAT2, # Cone intersects with the sphere
                           lat_0=ds.MOAD_CEN_LAT, lon_0=ds.STAND_LON, # Center point
                           a=6370000, b=6370000) # This is it! The Earth is a perfect sphere
    
    # Calculating the Projection from WRF lat long to WGS lat long
    # Easting and Northings of the domains center point
    wgs_proj = pyproj.Proj(proj='latlong', datum='WGS84')
    e, n = pyproj.transform(wgs_proj, wrf_proj, ds.CEN_LON, ds.CEN_LAT)
    
    # Grid parameters
    dx, dy = ds.DX, ds.DY
    nx, ny = ds.dims['west_east'], ds.dims['south_north']
    # Down left corner of the domain
    x0 = -(nx-1) / 2. * dx + e
    y0 = -(ny-1) / 2. * dy + n

    # 2d grid
    xx, yy = np.meshgrid(np.arange(nx) * dx + x0, np.arange(ny) * dy + y0)
    
    # Calculating the boundary X-Y Coordinates
    x_b, y_b = np.meshgrid(np.arange(nx+1) * dx + x0 -dx/2, np.arange(ny+1) * dy + y0 -dy/2)

    #Transformation of Center X-Y to Center Lat-Lon
    xc, yc = pyproj.transform(wrf_proj, wgs_proj, xx, yy)
    
    #Transformation of Boundary X-Y to Center Lat_Lon
    x_bc, y_bc = pyproj.transform(wrf_proj, wgs_proj, x_b, y_b)
    
    #Pulling out 1 variable to calc transformation with
    dr = ds.E_ISO[0,0,:,:]
    dr_int = dr.to_dataset(name='E_ISO')
    
    #Grid Spacing For prior grid
    in_coords = {'lat': dr_int['XLAT'].values, #Center Point Spacing Lat
                        'lon': dr_int['XLONG'].values, #Center Point Spacing Lon
                        'lat_b': y_bc, # Boundary Spacing Lat 
                        'lon_b': x_bc, # Boundary Spacing Lon
                       }
    
    #Grid Spacing For New Grid
    out_coords = {'lat': np.arange(23, 54, 0.1), #Center Point Spacing Lat
                        'lon': np.arange(-126, -63, 0.1), #Center Point Spacing Lon
                        'lat_b': np.arange(22.95, 54.05, 0.1), # Boundary Spacing Lat 
                        'lon_b': np.arange(-126.05, -62.95, 0.1), # Boundary Spacing Lon
                       }
    
    # Constructing Regridder
    regridder = xe.Regridder(in_coords, out_coords, method='conservative',reuse_weights=inputs.reuse_weight_file)
    return regridder, xc, yc

#end make regridder 

def reformat_ds(ds, lonc, latc):
    # Renaming Files to proper Names, setting coordinates
    ds = ds.rename({'south_north':'y'})
    ds = ds.rename({'west_east':'x'})
    ds = ds.rename({'emissions_zdim':'z'})
    ds = ds.assign({'lat': (['y', 'x'], latc),'lon': (['y', 'x'], lonc)})
    ds = ds.set_coords('lat')
    ds = ds.set_coords('lon')
    ds = ds.set_coords('Times')
    
    # Dropping Some Variables 
    
    # ds = ds.drop('E_ORGI_BB') 
    # ds = ds.drop('E_ORGJ_BB') 
    
    # Combining and renaming some variables
    
    # ds['E_NAAI'] = ds['E_NAAI'] + ds['E_NAAJ'] # Combine 
    # ds = ds.rename({'E_NAAI':'E_NAA'})
    # ds = ds.drop('E_NAAJ')
    
    ds['E_NO3I'] = ds['E_NO3I'] + ds['E_NO3J'] # Combine 
    ds = ds.rename({'E_NO3I':'E_NO3'})
    ds = ds.drop('E_NO3J')
    
    ds['E_ECI'] = ds['E_ECI'] + ds['E_ECJ'] # Combine 
    ds = ds.rename({'E_ECI':'E_EC'})
    ds = ds.drop('E_ECJ')
    
    ds['E_PM25I'] = ds['E_PM25I'] + ds['E_PM25J'] # Combine 
    ds = ds.rename({'E_PM25I':'E_PM25'})
    ds = ds.drop('E_PM25J')
    
    ds['E_SO4I'] = ds['E_SO4I'] + ds['E_SO4J'] # Combine 
    ds = ds.rename({'E_SO4I':'E_SO4'})
    ds = ds.drop('E_SO4J')
    
    ds['E_ORGI'] = ds['E_ORGI'] + ds['E_ORGJ'] # Combine 
    ds = ds.rename({'E_ORGI':'E_ORG'})
    ds = ds.drop('E_ORGJ')
    
    # ds['E_ORGI_A'] = ds['E_ORGI_A'] + ds['E_ORGJ_A'] # Combine 
    # ds = ds.rename({'E_ORGI_A':'E_ORG_A'})
    # ds = ds.drop('E_ORGJ_A')
    
    return ds
# end reformat_ds

def convert_units(ds_out):
    ds_out_convert = ds_out
    # Convert other units
    # originally in moles/km^2/hour, converting to kg/m^2/s
    for key in inputs.MW:
        ds_out_convert[key] = ds_out[key]*inputs.MW[key]/3600/1000000/1000 # 1hour/3600s 1m2/1e6km2 1kg/1000g
        
    # convert aerosol units
    for word in inputs.Aerosols:
        ds_out_convert[word] = ds_out[word]/1000000000 # convert from ug/m^3 m/s to kg/m^2/s by * 1e9
    return ds_out_convert
    
# End convert_units

def reformat_ds_out(ds_out_convert,day_half):
    #More formatting to COARDS
    # # Renaming Dimentions 
    ds_out_convert = ds_out_convert.rename({'z':'lev'})
    ds_out_convert = ds_out_convert.rename({'Time':'time'})
    ds_out_convert = ds_out_convert.rename({'Times':'time'})# rename coordinate of time
    
    
    ## Changing time to be what it should be for 7/4/2018
    #ds_out_convert['time'] = [162216, 162217, 162218, 162219, 162220, 162221, 162222, 162223, 162224, 162225, 162226, 162227]
    if day_half == '00z':
        ds_out_convert['time'] = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    elif day_half == '12z':
        ds_out_convert['time'] = [12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23]
        
    # # Setting time Coordinate Attributes
    #ds_out_convert['time'].attrs={'long_name': 'time', 'units': 'hours since 2000-01-01 00:00:00', 'calendar' : 'standard', 'axis' : 'T'}
    ds_out_convert['time'].attrs={'long_name': ('time'), 'units': ('hours since 2010-07-01 00:00:00 GMT'), 'calendar' : ('standard'), 'axis' : ('T')}
    #ds_out_convert['time'].attrs={'long_name': ('time'), 'units': ('secs since 1970-01-01 00:00:00'), 'calendar' : ('standard'), 'axis' : ('T')}

    # Setting lev Coordinate Attributes
    ds_out_convert['lev'].attrs={'long_name': ('GEOS-Chem level'), 'units': ('level'), 'positive' : ('up'), 'axis' : ('Z')}
    
    # Setting Lat Coordinate Attributes
    ds_out_convert['lat'].attrs={'long_name': ('Latitude'), 'units': ('degrees_north'),'axis' : ('Y')}
    
    # Setting Lon Coordinate Attributes
    ds_out_convert['lon'].attrs={'long_name': ('Longitude'), 'units': ('degrees_east'),'axis' : ('X')}
    
    #ds_out_convert['lev'] = [1, 3, 5, 7, 9, 11, 12, 13, 15, 16, 18, 19, 21, 22, 26, 28, 30, 33, 35, 39]
    ds_out_convert['lev'] = [1, 2, 4, 6, 8, 10, 11, 13, 14, 16, 17, 19, 20, 22, 25, 28, 30, 33, 35, 38]
    
    new_index = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47]
    ds_out_convert = ds_out_convert.reindex({"lev": new_index})
    ds_out_convert = ds_out_convert.fillna(0)
    
    # Adding Variable Attributes to all
    ds_out_convert['E_ALD'].attrs={'long_name': ('Aldehydes'), 'units': ('kg/m2/s')}
    ds_out_convert['E_CH4'].attrs={'long_name': ('Methane'), 'units': ('kg/m2/s')}
    ds_out_convert['E_CO'].attrs={'long_name': ('CO'), 'units': ('kg/m2/s')}
    ds_out_convert['E_CO2'].attrs={'long_name': ('CO2'), 'units': ('kg/m2/s')}
    ds_out_convert['E_CSL'].attrs={'long_name': ('Cresols Phenols'), 'units': ('kg/m2/s')}
    ds_out_convert['E_EC'].attrs={'long_name': ('Elemental Carbon'), 'units': ('kg/m2/s')}
    ds_out_convert['E_ETH'].attrs={'long_name': ('Ethane'), 'units': ('kg/m2/s')}
    ds_out_convert['E_HC3'].attrs={'long_name': ('HC3'), 'units': ('kg/m2/s')}
    ds_out_convert['E_HC5'].attrs={'long_name': ('HC5'), 'units': ('kg/m2/s')}
    ds_out_convert['E_HC8'].attrs={'long_name': ('HC8'), 'units': ('kg/m2/s')}
    ds_out_convert['E_HCHO'].attrs={'long_name': ('Formaldehyde'), 'units': ('kg/m2/s')}
    ds_out_convert['E_HONO'].attrs={'long_name': ('Nitrous Acid'), 'units': ('kg/m2/s')}
    ds_out_convert['E_ISO'].attrs={'long_name': ('Isoprene'), 'units': ('kg/m2/s')}
    ds_out_convert['E_UNID'].attrs={'long_name': ('Unidentified VOCs'), 'units': ('kg/m2/s')}
    ds_out_convert['E_KET'].attrs={'long_name': ('Ketones'), 'units': ('kg/m2/s')}
    #ds_out_convert['E_NAA'].attrs={'long_name': ('Naphthaleneacetic acid'), 'units': ('kg/m2/s')}
    ds_out_convert['E_NH3'].attrs={'long_name': ('Ammonia'), 'units': ('kg/m2/s')}
    ds_out_convert['E_NO'].attrs={'long_name': ('Nitric Oxide'), 'units': ('kg/m2/s')}
    ds_out_convert['E_NO2'].attrs={'long_name': ('Nitrogen Dioxide'), 'units': ('kg/m2/s')}
    ds_out_convert['E_NO3'].attrs={'long_name': ('Nitrate'), 'units': ('kg/m2/s')}
    ds_out_convert['E_OL2'].attrs={'long_name': ('Ethylene'), 'units': ('kg/m2/s')}
    ds_out_convert['E_OLI'].attrs={'long_name': ('Internal Alkenes'), 'units': ('kg/m2/s')}
    ds_out_convert['E_OLT'].attrs={'long_name': ('Terminal Alkenes'), 'units': ('kg/m2/s')}
    ds_out_convert['E_ORA2'].attrs={'long_name': ('C2+ Organic Acids'), 'units': ('kg/m2/s')}
    ds_out_convert['E_ORG'].attrs={'long_name': ('Organic Carbon'), 'units': ('kg/m2/s')}
    #ds_out_convert['E_ORG_A'].attrs={'long_name': ('Aircraft Organic Carbon'), 'units': ('kg/m2/s')}
    ds_out_convert['E_PM25'].attrs={'long_name': ('Unspeciated PM2.5'), 'units': ('kg/m2/s')}
    ds_out_convert['E_PM10'].attrs={'long_name': ('PM10'), 'units': ('kg/m2/s')}
    ds_out_convert['E_SO2'].attrs={'long_name': ('Sulfur Dioxide'), 'units': ('kg/m2/s')}
    ds_out_convert['E_SO4'].attrs={'long_name': ('Sulfate'), 'units': ('kg/m2/s')}
    ds_out_convert['E_TERP'].attrs={'long_name': ('Terpenes'), 'units': ('kg/m2/s')}
    ds_out_convert['E_TOL'].attrs={'long_name': ('Toluene and less reactive aromatics'), 'units': ('kg/m2/s')}
    ds_out_convert['E_XYL'].attrs={'long_name': ('Xylenes and more reactive aromatics'), 'units': ('kg/m2/s')}

    return ds_out_convert
# end reformat_ds_out

def combine_dayhalves(ds_00z, ds_12z):
    merger = xr.concat([ds_00z, ds_12z], dim='time')
    # merger['time'].astype(int)
    merger['time'].attrs={'long_name': ('time'), 'units': ('hours since 2010-07-01 00:00:00 GMT'), 'calendar' : ('standard'), 'axis' : ('T')}
    return merger
# End Combine_dayhalves

def set_output_attrs(ds_out):
    for ii in ds_out.data_vars:
        ds_out[ii].encoding={'dtype': 'float32', 'chunksizes': (1,1,310,630),'zlib': True, 'complevel': 1}
        #ds_out[ii].attrs={'description': (ds[ii].attrs['description'])}
        #print(ds_out[ii].encoding)
    # setting global attributes
    ds_out.attrs['Title']='Emission inventory Regridded from WRF to Lat Lon grid'
    ds_out.attrs['Conventions']='COARDS'
    ds_out.attrs['History']=''
    #del ds_out_convert.attrs['global_attributes']
    return ds_out
# End set_output_attrs
    
def main():

    ds_example = xr.open_dataset(inputs.base_dir  + str(inputs.years[0]) + '/' + list(inputs.months)[0] + '/' + inputs.days[0] + '/' +inputs.fn_base1
                                 , chunks={'Time': 1,'emissions_zdim': 1}) 
    
    regridder, lonc, latc = makeregridder(ds_example)
    
    for year in inputs.years:
        for month in inputs.months.keys():
            for day in inputs.days:
                dir_name = inputs.base_dir  + str(year) + '/' + month + '/' + day + '/'
                fn_00z = dir_name + inputs.fn_base1
                fn_12z = dir_name + inputs.fn_base2
                
                print('Reading Files:   ')
                print(fn_00z)
                print(fn_12z)
                
                ds_in_00z = xr.open_dataset(fn_00z,chunks={'Time': 1,'emissions_zdim': 1})
                ds_in_12z = xr.open_dataset(fn_12z,chunks={'Time': 1,'emissions_zdim': 1})
    
                ds_in_00z = reformat_ds(ds_in_00z,lonc,latc)
                ds_in_12z = reformat_ds(ds_in_12z,lonc,latc)
                
                # Call Regridder on the full Dataset
                ds_out_00z = regridder(ds_in_00z)
                ds_out_12z = regridder(ds_in_12z)
                
                ds_out_00z = convert_units(ds_out_00z)
                ds_out_12z = convert_units(ds_out_12z)
                
                ds_out_00z = reformat_ds_out(ds_out_00z,'00z')
                ds_out_12z = reformat_ds_out(ds_out_12z,'12z')
                
                ds_out = combine_dayhalves(ds_out_00z, ds_out_12z)
                
                ds_out = set_output_attrs(ds_out)    
                
                out_dir = inputs.out_dir  + str(year) + '/' + month + '/' + day + '/'
                
                # Save out the File to NetCDF
                ds_out.to_netcdf(out_dir + inputs.fn_out,format='netCDF4',engine='netcdf4')
#End main

if __name__ == "__main__":
    main()


