
# coding: utf-8

# Read in and download MEPS data

# In[ ]:


import sys
sys.path.append('/uio/kant/geo-metos-u7/franzihe/Documents/Thesis/Python')
import time
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
#import fill_values as fv
#import calc_station_properties as cs

import createFolder as cF
from scipy.integrate import simps
import fcts_read_stat as rs


# In[ ]:


thredds      = 'http://thredds.met.no/thredds/dodsC/meps25epsarchive'

stn_name     = 'Haukeliseter'
stn_lat      = 59.8
stn_lon      = 7.2
#stn_name     = 'Stavanger'
#stn_lat      = 58.87
#stn_lon      = 5.67

#month        = 12
#day          = 24
forecasttime = '00'
all_em       = 0      # 1==yes, 0==no
em0          = 1

if all_em == 1:
#    met_file = 'meps_subset_2_5km_'
    layers = 'sfc_hybrid5_allEM'
elif em0 == 1:
 #   met_file = 'meps_mbr0_full_backup_2_5km_'
  #  met_file2 = 'meps_mbr0_vc_2_5km_'
    layers = 'hybrid65_EM0'


# In[ ]:


#year         = 2016

### Upslope
# Oct 2016
#month = 10
#t     = [28]

# Nov 2016
#month = 11
#t = np.arange(11,19)
#t = np.append(t,[28,29])



# Dez 2016
#month = 12
#t = np.arange(21,28)

#year         = 2017
# Jan 2017
#month = '01'
#t = np.arange(2,14)
#t = np.append(t,[28,29])

# Feb 2017
#month = '02'
#t = np.arange(1,5)
#m = ['11', '12', '01', '02', '03']
#m = ['01', '02', '03']
m = ['12', '01', '02', '03']




# In[ ]:


main_dir       = '../../Data/MEPS'



# In[ ]:


def read_for_station(thredds,year,month,day,forecasttime,stn_lat,stn_lon,dirnc):
    
    if all_em == 1:
        met_files = ['meps_subset_2_5km_' ]
        memb = np.arange(0,10)
        k = 0
    elif em0 == 1:
        met_files = ['meps_subset_2_5km_' , 'meps_mbr0_full_backup_2_5km_'#, 'meps_mbr0_vc_2_5km_'
                ]
        memb = np.arange(0,1)
        k = 1
    
    fn = dict()
    for i in range(0,np.shape(met_files)[0]):
 #       fn[i] = netCDF4.Dataset('%s/%s/%s/%s/%s%s%s%sT%sZ.nc' %(thredds,year,month,day,met_files[i],year,month,day,forecasttime),
  #                               'r')
        try:
            fn[i] = netCDF4.Dataset('%s/%s/%s/%s/%s%s%s%sT%sZ.nc' %(thredds,year,month,day,met_files[i],year,month,day,forecasttime),
                                 'r')
        except OSError:
            print('no file found: %s/%s/%s/%s/%s%s%s%sT%sZ.nc' %(thredds,year,month,day,met_files[i],year,month,day,forecasttime))
            return

## Latitudes
## [y = 949][x = 739]
    latitude = fn[0].variables['latitude']

## Longitudes 
## [y = 949][x = 739]
    longitude = fn[0].variables['longitude']

# Now find the absolute value of the difference between the  station's lat/lon with every point in the grid. 
    x,y = rs.find_station_yx(latitude, longitude, stn_lat, stn_lon)
    y = y[0]
    x = x[0]

    

######## with Vertical Levels ( hybrid )#################################
## hybrid levels atmosphere_hybrid_sigma_pressure_coordinate
# formula: p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)
# positive: down


### variables to calculate pressure
#p0 = fn.variables['p0']    ## p0: p0
    try:
      ap = fn[k].variables['ap']    ## ap: ap
    except KeyError:
      ap = fn[k].variables['ap0']
    try:
      b = fn[k].variables['b']      ## b: b
    except KeyError:
      b = fn[k].variables['b0']

    surface_air_pressure = fn[0].variables['surface_air_pressure']
    air_temperature_0m   = fn[0].variables['air_temperature_0m']
    try:
      air_temperature_ml   = fn[k].variables['air_temperature_ml']
    except KeyError:
      print('no air_temperature_ml found --> no file saved!!! %s/%s%s%s_%s.nc' %(dirnc,year,month,day,forecasttime))
      return
      

### mask arrays
    surface_air_pressure, dtype_sap  = rs.mask_array(surface_air_pressure, #ens_memb, 
                                                  y, x,EM=surface_air_pressure.shape[2])
    air_temperature_0m,   dtype_at0m = rs.mask_array(air_temperature_0m,   #ens_memb, 
                                                  y, x,EM=air_temperature_0m.shape[2])
    air_temperature_ml,   dtype_atml = rs.mask_array(air_temperature_ml,   #ens_memb, 
                                                  y, x,EM=air_temperature_ml.shape[2])

    for ens_memb in memb:
### Transfer pressure coordinate
        p_interface = ap[:]+b[:]*surface_air_pressure[:,:,ens_memb]

### 1) Connect interface values and surface values for pressure
### Pressure
        p_interfaces2 = np.concatenate((p_interface[:,:],surface_air_pressure[:,:,ens_memb]),axis=1)

# transform hybrid sigma pressure coordinates at interface levels to pressure at model levels
        p_ml = np.empty([p_interfaces2.shape[0],p_interfaces2.shape[1]-1])
        for i in range(0,p_interfaces2.shape[1]-1):
            p_ml[:,i] = 1/2 * (p_interfaces2[:,i] + p_interfaces2[:,i+1])
        p_ml = np.concatenate((p_ml[:,:],surface_air_pressure[:,:,ens_memb]),axis=1)
    
### --> Now the pressure is calculated at each model level

### 2) Connect model levels and surface values for temperature
### Temperature
        temperature_ml = np.concatenate((air_temperature_ml[:air_temperature_ml.shape[0],:,ens_memb],
                                        air_temperature_0m[:air_temperature_ml.shape[0],:,ens_memb]),axis=1)

    

        dz, dgeop = rs.get_thickness(p_ml, temperature_ml)

        if ens_memb == 0:
            pressure_in_modellev = p_ml[:,0:-1]
            thickness_m   = dz
            thickness_phi = dgeop  
        else:
            pressure_in_modellev = np.dstack((pressure_in_modellev,p_ml[:,0:5]))
            thickness_m   = np.dstack((thickness_m,dz))
            thickness_phi = np.dstack((thickness_phi,dgeop))

    pressure_in_modellev = np.ma.array(pressure_in_modellev,mask=np.ma.is_masked(pressure_in_modellev), 
                                       fill_value = np.nan)
    thickness_m = np.ma.array(thickness_m,mask=np.ma.is_masked(thickness_m), 
                                       fill_value = np.nan)
    thickness_phi = np.ma.array(thickness_phi,mask=np.ma.is_masked(thickness_phi), 
                                       fill_value = np.nan)


# Read in all values needed to present the microphysics
## Time
    time_arr = fn[k].variables['time']
    ensemble_member_arr = fn[k].variables['ensemble_member']

## heights
    hybrid_arr = fn[k].variables['hybrid']

    if all_em == 1:
        height0_arr = fn[k].variables['height0']
        height1_arr = fn[k].variables['height1']
        height3_arr = fn[k].variables['height3']
        height_above_msl_arr = fn[k].variables['height_above_msl']

    
######## with Vertical Levels ( height0 ) #################################
        air_temperature_0m = rs.get_value_at_station(fn[k], 'air_temperature_0m', y,x)
        #graupelfall_amount = rs.get_value_at_station(fn[k], 'graupelfall_amount', y,x)
        liquid_water_content_of_surface_snow = rs.get_value_at_station(fn[k], 'liquid_water_content_of_surface_snow',y,x)
        precipitation_amount_acc = rs.get_value_at_station(fn[k], 'precipitation_amount_acc',y,x)
        snowfall_amount_acc = rs.get_value_at_station(fn[k], 'snowfall_amount_acc',y,x)
        #rainfall_amount = rs.get_value_at_station(fn[k], 'rainfall_amount',y,x)
        #snowfall_amount = rs.get_value_at_station(fn[k], 'snowfall_amount',y,x)
        surface_air_pressure = rs.get_value_at_station(fn[k], 'surface_air_pressure',y,x)
        surface_geopotential = rs.get_value_at_station(fn[k], 'surface_geopotential',y,x)

######## with Vertical Levels ( height1 )#################################
        air_temperature_2m = rs.get_value_at_station(fn[k],'air_temperature_2m',y,x)
        #specific_humidity_2m = rs.get_value_at_station(fn[k],'specific_humidity_2m',y,x)

######## with Vertical Levels ( height3 )#################################
        x_wind_10m = rs.get_value_at_station(fn[k],'x_wind_10m',y,x)
        y_wind_10m = rs.get_value_at_station(fn[k],'y_wind_10m',y,x)

######## with Vertical Levels ( height_above_msl )#################################
        air_pressure_at_sea_level = rs.get_value_at_station(fn[k],'air_pressure_at_sea_level',y,x)

######## with Vertical Levels ( hybrid )#################################
## hybrid levels atmosphere_hybrid_sigma_pressure_coordinate
# formula: p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)
# positive: down

## values in fn file
    air_temperature_ml = rs.get_value_at_station(fn[k],'air_temperature_ml',y,x)
    specific_humidity_ml = rs.get_value_at_station(fn[k],'specific_humidity_ml',y,x)
    x_wind_ml = rs.get_value_at_station(fn[k],'x_wind_ml',y,x)
    y_wind_ml = rs.get_value_at_station(fn[k],'y_wind_ml',y,x)

    if em0 == 1:
        atmosphere_cloud_condensed_water_content_ml= rs.get_value_at_station(fn[k],'atmosphere_cloud_condensed_water_content_ml',y,x)
        pressure_departure = rs.get_value_at_station(fn[k],'pressure_departure',y,x)
        atmosphere_cloud_ice_content_ml = rs.get_value_at_station(fn[k],'atmosphere_cloud_ice_content_ml',y,x)
        upward_air_velocity_ml = rs.get_value_at_station(fn[k],'upward_air_velocity_ml',y,x)
    
## values in 2nd file
#       graupelfall_amount_ml = rs.get_value_at_station(fn[2],'graupelfall_amount_ml',y,x)
#       rainfall_amount_ml = rs.get_value_at_station(fn[2],'rainfall_amount_ml',y,x)
#       snowfall_amount_ml = rs.get_value_at_station(fn[2],'snowfall_amount_ml',y,x)

### write netCDF file
    f = netCDF4.Dataset('%s/%s%s%s_%s.nc' %(dirnc,year,month,day,forecasttime), 'w')

### create dimensions
    f.createDimension('time', time_arr.shape[0])
    f.createDimension('hybrid', hybrid_arr.shape[0])
    f.createDimension('ensemble_member', ensemble_member_arr.shape[0])

    if all_em == 1:
        f.createDimension('height0', height0_arr.shape[0])
        f.createDimension('height1', height1_arr.shape[0])
        f.createDimension('height3', height3_arr.shape[0])
        f.createDimension('height_above_msl', height_above_msl_arr.shape[0])

    t = f.createVariable('time', time_arr.dtype,'time',zlib = True)
    t[:] = time_arr[:]


    if all_em == 1:
    
######## with Vertical Levels ( height0 ) #################################
        h = f.createVariable('height0', height0_arr.dtype, 'height0', zlib=True)
        h[:] = height0_arr[:]
        dim = ('time', 'height0', 'ensemble_member')

        at_0m = rs.get_netCDF_variable(f,'air_temperature_0m', air_temperature_0m,dim)
        #ga_0m = rs.get_netCDF_variable(f,'graupelfall_amount', graupelfall_amount,dim)
        lwc_0m = rs.get_netCDF_variable(f,'liquid_water_content_of_surface_snow', liquid_water_content_of_surface_snow,dim)
        pr_0m = rs.get_netCDF_variable(f,'precipitation_amount_acc',precipitation_amount_acc,dim)
        sn_0m = rs.get_netCDF_variable(f,'snowfall_amount_acc',snowfall_amount_acc,dim)
        #ra_0m = rs.get_netCDF_variable(f,'rainfall_amount',rainfall_amount,dim)
        #sa_0m = rs.get_netCDF_variable(f,'snowfall_amount',snowfall_amount,dim)
        ps = rs.get_netCDF_variable(f,'surface_air_pressure',surface_air_pressure,dim)
        geop = rs.get_netCDF_variable(f,'surface_geopotential',surface_geopotential,dim)

######## with Vertical Levels ( height1 )#################################
        h1 = f.createVariable('height1', height1_arr.dtype, 'height1', zlib=True)
        h1[:] = height1_arr[:]
        dim = ('time', 'height1', 'ensemble_member')

        at_2m = rs.get_netCDF_variable(f,'air_temperature_2m', air_temperature_2m,dim)
        #sh_2m = cs.rs.get_netCDF_variable(f,'specific_humidity_2m',specific_humidity_2m,dim)

######## with Vertical Levels( height3 )#################################
        h3 = f.createVariable('height3', height3_arr.dtype, 'height3', zlib=True)
        h3[:] = height3_arr[:]
        dim = ('time', 'height3', 'ensemble_member')
        
        xwind_10m = rs.get_netCDF_variable(f,'x_wind_10m', x_wind_10m,dim)
        ywind_10m = rs.get_netCDF_variable(f,'y_wind_10m', y_wind_10m,dim)

######## with Vertical Levels ( height_above_msl )#################################
        h_asl = f.createVariable('height_above_msl', height_above_msl_arr.dtype, 'height_above_msl', zlib=True)
        h_asl[:] = height_above_msl_arr[:]
        dim = ('time', 'height_above_msl', 'ensemble_member')

        pressure_sea_level = rs.get_netCDF_variable(f,'air_pressure_at_sea_level', air_pressure_at_sea_level,dim)

######## with Vertical Levels ( hybrid )#################################
    hyb = f.createVariable('hybrid', hybrid_arr.dtype, 'hybrid', zlib=True)
        
    hyb[:] = hybrid_arr[:]

    if all_em == 1:
        dim = ('time','hybrid', 'ensemble_member')
    elif em0 == 1:
        dim = ('time','hybrid')

    at_ml = rs.get_netCDF_variable(f,'air_temperature_ml',air_temperature_ml,dim)
    sh_ml = rs.get_netCDF_variable(f,'specific_humidity_ml',specific_humidity_ml,dim)
    xwind_ml = rs.get_netCDF_variable(f,'x_wind_ml',x_wind_ml,dim)
    ywind_ml = rs.get_netCDF_variable(f,'y_wind_ml',y_wind_ml,dim)

    if em0 == 1:
        ccw_ml = rs.get_netCDF_variable(f,'atmosphere_cloud_condensed_water_content_ml',
                                     atmosphere_cloud_condensed_water_content_ml,dim)
        cic_ml = rs.get_netCDF_variable(f,'atmosphere_cloud_ice_content_ml',atmosphere_cloud_ice_content_ml,dim)
        pres_dep_ml = rs.get_netCDF_variable(f,'pressure_departure',pressure_departure,dim)
        up_air_velo_ml = rs.get_netCDF_variable(f,'upward_air_velocity_ml', upward_air_velocity_ml, dim)
#       sf_ml = rs.get_netCDF_variable(f,'snowfall_amount_ml',snowfall_amount_ml,dim)
#       rf_ml = rs.get_netCDF_variable(f,'rainfall_amount_ml',rainfall_amount_ml,dim)
#       gf_ml = rs.get_netCDF_variable(f,'graupelfall_amount_ml',graupelfall_amount_ml,dim)

    pres_ml = rs.get_netCDF_variable(f,'pressure_ml',pressure_in_modellev[:time_arr.shape[0],:],dim)
    dz_ml = rs.get_netCDF_variable(f,'layer_thickness',thickness_m[:time_arr.shape[0],:],dim)
    dgeop_ml = rs.get_netCDF_variable(f,'geop_layer_thickness',thickness_phi[:time_arr.shape[0],:],dim)

    f.close()

    for i in range(0,np.shape(met_files)[0]):
        fn[i].close()
    print('file written: %s/%s%s%s_%s.nc' %(dirnc,year,month,day,forecasttime))


# In[ ]:


#%%time
for month in m:
  if month == '11':
    t = np.arange(8,31)
  if month == '12': #or month == '03':
    t = np.arange(16,32)
  if month == '01' or month == '03':
    t = np.arange(1,32)
#  if month == '12':
 #   t = np.arange(19,32)
  #if month == '01' or month == '03':
   # t = np.arange(1,32)
  if month == '02':
    t = np.arange(1,29)
  if month == '11' or month == '12':
    year = '2016'
  if month == '01' or month == '02' or month == '03':
    year = '2017'
  dirnc          = '%s/%s/%s%s/%s_%s' %(main_dir,stn_name,year,month,layers,forecasttime)
  for day in t:
    if day < 10:
        day = '0%s' %(day)
        
    start_time = time.time()
    
    ### direction where files should be saved
    cF.createFolder('%s' %(dirnc))
    read_for_station(thredds,year,month,day,forecasttime,stn_lat,stn_lon,dirnc)
 #   print('file written: %s/%s%s%s_%s.nc' %(dirnc,year,month,day,forecasttime))
    
    print("--- %s seconds ---" % round(time.time() - start_time, 2))

