#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys

sys.path.append('/home/franzihe/Documents/Python/Thesis/')
#sys.path.append('/uio/kant/geo-metos-u7/franzihe/Documents/Thesis/Python')
import time
import netCDF4
import numpy as np
# import fill_values as fv
# import calc_station_properties as cs

import createFolder as cF
import fcts_read_stat as rs
import gc
# %%

# %%

# In[ ]:

thredds = 'http://thredds.met.no/thredds/dodsC/meps25epsarchive'

stn_name = 'Kiruna'
stn_lat = 67.84
stn_lon = 20.41

# month        = 12
# day          = 24
forecasttime = '00'
# m = ['09', '10', '11', '12', '01', '02', '03','04']
# m = ['10', '11',
m = ['03', '04']

# In[ ]:


main_dir = '../../Data/MEPS'


# In[ ]:


def read_for_station(thredds, year, month, day, forecasttime, stn_lat, stn_lon, dirnc):
    met_files = ['meps_full_2_5km_']

    try:
        fn = netCDF4.Dataset(
            '%s/%s/%s/%s/%s%s%s%sT%sZ.nc' % (thredds, year, month, day, met_files[0], year, month, day, forecasttime),
            'r')
    except OSError:
        print('no file found: %s/%s/%s/%s/%s%s%s%sT%sZ.nc' % (
        thredds, year, month, day, met_files[0], year, month, day, forecasttime))
        return

    ## Latitudes
    ## [y = 949][x = 739]
    latitude = fn.variables['latitude']

    ## Longitudes
    ## [y = 949][x = 739]
    longitude = fn.variables['longitude']

    # Now find the absolute value of the difference between the  station's lat/lon with every point in the grid.
    x, y = rs.find_station_yx(latitude, longitude, stn_lat, stn_lon)
    y = y[0]
    x = x[0]

    del latitude, longitude

    ######## with Vertical Levels ( hybrid )#################################
    ## hybrid levels atmosphere_hybrid_sigma_pressure_coordinate
    # formula: p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)
    # positive: down

    ### variables to calculate pressure
    # p0 = fn.variables['p0']    ## p0: p0
    ap = fn.variables['ap']  ## ap: ap
    b = fn.variables['b']  ## b: b

    surface_air_pressure = fn.variables['surface_air_pressure']
    air_temperature_0m = fn.variables['air_temperature_0m']
    air_temperature_ml = fn.variables['air_temperature_ml']

    ### mask arrays
    surface_air_pressure, dtype_sap = rs.mask_array(surface_air_pressure,
                                                    y, x, EM=surface_air_pressure.shape[2])
    air_temperature_0m, dtype_at0m = rs.mask_array(air_temperature_0m,
                                                   y, x, EM=air_temperature_0m.shape[2])
    air_temperature_ml, dtype_atml = rs.mask_array(air_temperature_ml,
                                                   y, x, EM=air_temperature_ml.shape[2])

    for ens_memb in np.arange(0, 10):
        ### Transfer pressure coordinate
        p_interface = ap[:] + b[:] * surface_air_pressure[:, :, ens_memb]

        ### 1) Connect interface values and surface values for pressure
        ### Pressure
        p_interfaces2 = np.concatenate((p_interface[:, :], surface_air_pressure[:, :, ens_memb]), axis=1)

        # transform hybrid sigma pressure coordinates at interface levels to pressure at model levels
        p_ml = np.empty([p_interfaces2.shape[0], p_interfaces2.shape[1] - 1])
        for i in range(0, p_interfaces2.shape[1] - 1):
            p_ml[:, i] = 1 / 2 * (p_interfaces2[:, i] + p_interfaces2[:, i + 1])
        p_ml = np.concatenate((p_ml[:, :], surface_air_pressure[:, :, ens_memb]), axis=1)

        ### --> Now the pressure is calculated at each model level

        ### 2) Connect model levels and surface values for temperature
        ### Temperature
        temperature_ml = np.concatenate((air_temperature_ml[:, :, ens_memb],
                                         air_temperature_0m[:, :, ens_memb]), axis=1)

        dz, dgeop = rs.get_thickness(p_ml, temperature_ml)

        if ens_memb == 0:
            pressure_in_modellev = p_ml[:, 0:-1]
            thickness_m = dz
            thickness_phi = dgeop
        else:
            pressure_in_modellev = np.dstack((pressure_in_modellev, p_ml[:, 0:-1]))
            thickness_m = np.dstack((thickness_m, dz))
            thickness_phi = np.dstack((thickness_phi, dgeop))

    pressure_in_modellev = np.ma.array(pressure_in_modellev, mask=np.ma.is_masked(pressure_in_modellev),
                                       fill_value=np.nan)
    thickness_m = np.ma.array(thickness_m, mask=np.ma.is_masked(thickness_m),
                              fill_value=np.nan)
    thickness_phi = np.ma.array(thickness_phi, mask=np.ma.is_masked(thickness_phi),
                                fill_value=np.nan)

    del ap, b, p_ml, dz, dgeop

    # Read in all values needed to present the microphysics
    ## Time
    time_arr = fn.variables['time'][:]
    ensemble_member_arr = fn.variables['ensemble_member'][:]

    ## heights
    height_above_msl_arr = fn.variables['height_above_msl'][:]
    height0_arr = fn.variables['height0'][:]
    height1_arr = fn.variables['height1'][:]
    height7_arr = fn.variables['height7'][:]
    hybrid_arr = fn.variables['hybrid'][:]
    atm_as_single_layer_arr = fn.variables['atmosphere_as_single_layer'][:]

    ######## with Vertical Levels ( height0 ) #################################
    liquid_water_content_of_surface_snow = rs.get_value_at_station(fn, 'liquid_water_content_of_surface_snow', y, x)
    air_temperature_0m = rs.get_value_at_station(fn, 'air_temperature_0m', y, x)
    rainfall_amount = rs.get_value_at_station(fn, 'rainfall_amount', y, x)
    snowfall_amount = rs.get_value_at_station(fn, 'snowfall_amount', y, x)
    graupelfall_amount = rs.get_value_at_station(fn, 'graupelfall_amount', y, x)
    surface_air_pressure = rs.get_value_at_station(fn, 'surface_air_pressure', y, x)
    surface_geopotential = rs.get_value_at_station(fn, 'surface_geopotential', y, x)
    precipitation_amount_acc = rs.get_value_at_station(fn, 'precipitation_amount_acc', y, x)
    precipitation_type = rs.get_value_at_station(fn, 'precipitation_type', y, x)

    ######## with Vertical Levels ( height1 )#################################
    air_temperature_2m = rs.get_value_at_station(fn, 'air_temperature_2m', y, x)
    relative_humidity_2m = rs.get_value_at_station(fn, 'relative_humidity_2m', y, x)
    specific_humidity_2m = rs.get_value_at_station(fn, 'specific_humidity_2m', y, x)

    ######## with Vertical Levels ( height7 )#################################
    x_wind_10m = rs.get_value_at_station(fn, 'x_wind_10m', y, x)
    y_wind_10m = rs.get_value_at_station(fn, 'y_wind_10m', y, x)

    ######## with Vertical Levels ( height_above_msl )#################################
    air_pressure_at_sea_level = rs.get_value_at_station(fn, 'air_pressure_at_sea_level', y, x)

    ######## with Vertical Levels ( hybrid )#################################
    ## hybrid levels atmosphere_hybrid_sigma_pressure_coordinate
    # formula: p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)
    # positive: down

    specific_humidity_ml = rs.get_value_at_station(fn, 'specific_humidity_ml', y, x)
    atmosphere_cloud_condensed_water_content_ml = rs.get_value_at_station(fn,
                                                                          'atmosphere_cloud_condensed_water_content_ml',
                                                                          y, x)
    atmosphere_cloud_ice_content_ml = rs.get_value_at_station(fn, 'atmosphere_cloud_ice_content_ml', y, x)
    snowfall_amount_ml = rs.get_value_at_station(fn, 'snowfall_amount_ml', y, x)
    rainfall_amount_ml = rs.get_value_at_station(fn, 'rainfall_amount_ml', y, x)
    graupelfall_amount_ml = rs.get_value_at_station(fn, 'graupelfall_amount_ml', y, x)
    pressure_departure = rs.get_value_at_station(fn, 'pressure_departure', y, x)
    air_temperature_ml = rs.get_value_at_station(fn, 'air_temperature_ml', y, x)
    x_wind_ml = rs.get_value_at_station(fn, 'x_wind_ml', y, x)
    y_wind_ml = rs.get_value_at_station(fn, 'y_wind_ml', y, x)

    ######## with Vertical Levels ( atmosphere as single layer )#################################
    cloud_base_altitude = rs.get_value_at_station(fn, 'cloud_base_altitude', y, x)
    cloud_top_altitude = rs.get_value_at_station(fn, 'cloud_top_altitude', y, x)
    integral_graupel_wrt_height = rs.get_value_at_station(fn, 'integral_graupel_wrt_height', y, x)

    ### write netCDF file
    f = netCDF4.Dataset('%s/%s%s%s_%s.nc' % (dirnc, year, month, day, forecasttime), 'w')

    ### create dimensions
    f.createDimension('time', time_arr.shape[0])
    f.createDimension('ensemble_member', ensemble_member_arr.shape[0])
    f.createDimension('height_above_msl', height_above_msl_arr.shape[0])
    f.createDimension('height0', height0_arr.shape[0])
    f.createDimension('height1', height1_arr.shape[0])
    f.createDimension('height7', height7_arr.shape[0])
    f.createDimension('hybrid', hybrid_arr.shape[0])
    f.createDimension('atmosphere_as_single_layer', atm_as_single_layer_arr.shape[0])

    t = f.createVariable('time', time_arr.dtype, 'time', zlib=True)
    t[:] = time_arr[:]

    ######## with Vertical Levels ( height0 ) #################################
    h = f.createVariable('height0', height0_arr.dtype, 'height0', zlib=True)
    h[:] = height0_arr[:]
    dim = ('time', 'height0', 'ensemble_member')

    lwc_0m = rs.get_netCDF_variable(f, 'liquid_water_content_of_surface_snow', liquid_water_content_of_surface_snow,
                                    dim)
    at_0m = rs.get_netCDF_variable(f, 'air_temperature_0m', air_temperature_0m, dim)
    ra_0m = rs.get_netCDF_variable(f, 'rainfall_amount', rainfall_amount, dim)
    sn_0m = rs.get_netCDF_variable(f, 'snowfall_amount', snowfall_amount, dim)
    ga_0m = rs.get_netCDF_variable(f, 'graupelfall_amount', graupelfall_amount, dim)
    ps_0m = rs.get_netCDF_variable(f, 'surface_air_pressure', surface_air_pressure, dim)
    geop_0m = rs.get_netCDF_variable(f, 'surface_geopotential', surface_geopotential, dim)
    pr_0m = rs.get_netCDF_variable(f, 'precipitation_amount_acc', precipitation_amount_acc, dim)
    precip_type_0m = rs.get_netCDF_variable(f, 'precipitation_type', precipitation_type, dim)

    ######## with Vertical Levels ( height1 )#################################
    h1 = f.createVariable('height1', height1_arr.dtype, 'height1', zlib=True)
    h1[:] = height1_arr[:]
    dim = ('time', 'height1', 'ensemble_member')

    at_2m = rs.get_netCDF_variable(f, 'air_temperature_2m', air_temperature_2m, dim)
    rh_2m = rs.get_netCDF_variable(f, 'relative_humidity_2m', relative_humidity_2m, dim)
    sh_2m = rs.get_netCDF_variable(f, 'specific_humidity_2m', specific_humidity_2m, dim)

    ######## with Vertical Levels( height7 )#################################
    h7 = f.createVariable('height7', height7_arr.dtype, 'height7', zlib=True)
    h7[:] = height7_arr[:]
    dim = ('time', 'height7', 'ensemble_member')

    xwind_10m = rs.get_netCDF_variable(f, 'x_wind_10m', x_wind_10m, dim)
    ywind_10m = rs.get_netCDF_variable(f, 'y_wind_10m', y_wind_10m, dim)

    ######## with Vertical Levels ( height_above_msl )#################################
    h_asl = f.createVariable('height_above_msl', height_above_msl_arr.dtype, 'height_above_msl', zlib=True)
    h_asl[:] = height_above_msl_arr[:]
    dim = ('time', 'height_above_msl', 'ensemble_member')

    pressure_sea_level = rs.get_netCDF_variable(f, 'air_pressure_at_sea_level', air_pressure_at_sea_level, dim)

    ######## with Vertical Levels ( hybrid )#################################
    hyb = f.createVariable('hybrid', hybrid_arr.dtype, 'hybrid', zlib=True)
    hyb[:] = hybrid_arr[:]
    dim = ('time', 'hybrid', 'ensemble_member')

    sh_ml = rs.get_netCDF_variable(f, 'specific_humidity_ml', specific_humidity_ml, dim)
    ccw_ml = rs.get_netCDF_variable(f, 'atmosphere_cloud_condensed_water_content_ml',
                                    atmosphere_cloud_condensed_water_content_ml, dim)
    cic_ml = rs.get_netCDF_variable(f, 'atmosphere_cloud_ice_content_ml', atmosphere_cloud_ice_content_ml, dim)
    sf_ml = rs.get_netCDF_variable(f, 'snowfall_amount_ml', snowfall_amount_ml, dim)
    rf_ml = rs.get_netCDF_variable(f, 'rainfall_amount_ml', rainfall_amount_ml, dim)
    gf_ml = rs.get_netCDF_variable(f, 'graupelfall_amount_ml', graupelfall_amount_ml, dim)
    pres_dep_ml = rs.get_netCDF_variable(f, 'pressure_departure', pressure_departure, dim)
    at_ml = rs.get_netCDF_variable(f, 'air_temperature_ml', air_temperature_ml, dim)

    xwind_ml = rs.get_netCDF_variable(f, 'x_wind_ml', x_wind_ml, dim)
    ywind_ml = rs.get_netCDF_variable(f, 'y_wind_ml', y_wind_ml, dim)

    pres_ml = rs.get_netCDF_variable(f, 'pressure_ml', pressure_in_modellev, dim)
    dz_ml = rs.get_netCDF_variable(f, 'layer_thickness', thickness_m, dim)
    dgeop_ml = rs.get_netCDF_variable(f, 'geop_layer_thickness', thickness_phi, dim)

    ######## with Vertical Levels ( atmosphere as single layer )#################################
    asl = f.createVariable('atmosphere_as_single_layer', atm_as_single_layer_arr.dtype, 'atmosphere_as_single_layer',
                           zlib=True)
    asl[:] = atm_as_single_layer_arr[:]
    dim = ('time', 'atmosphere_as_single_layer', 'ensemble_member')

    cb_alt = rs.get_netCDF_variable(f, 'cloud_base_altitude', cloud_base_altitude, dim)
    ct_alt = rs.get_netCDF_variable(f, 'cloud_top_altitude', cloud_top_altitude, dim)
    int_graupel = rs.get_netCDF_variable(f, 'integral_graupel_wrt_height', integral_graupel_wrt_height, dim)

    f.close()
    fn.close()

    del time_arr, ensemble_member_arr, height_above_msl_arr, height0_arr, height1_arr, height7_arr, hybrid_arr, atm_as_single_layer_arr
    gc.collect()

    print('file written: %s/%s%s%s_%s.nc' % (dirnc, year, month, day, forecasttime))


# In[ ]:


# %%time
for month in m:
    if month == '10':  # or
        t = np.arange(28, 32)
    if month == '12' or month == '03':
        t = np.arange(1, 32)
    if month == '01':
        t = np.arange(31, 32)
    if month == '09' or month == '11' or month == '04':
        t = np.arange(1, 31)
    if month == '02':
        t = np.arange(28, 29)
    if month == '09' or month == '10' or month == '11' or month == '12':
        year = '2017'
    if month == '01' or month == '02' or month == '03' or month == '04':
        year = '2018'

    dirnc = '%s/%s/%s%s/%s' % (main_dir, stn_name, year, month, forecasttime)
    for day in t:
        if day < 10:
            day = '0%s' % (day)

        start_time = time.time()

        ### direction where files should be saved
        cF.createFolder('%s' % (dirnc))
        read_for_station(thredds, year, month, day, forecasttime, stn_lat, stn_lon, dirnc)

        print("--- %s seconds ---" % round(time.time() - start_time, 2))
