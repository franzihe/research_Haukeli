# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:light
#     text_representation:
#       extension: .py
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.3.0
#   kernelspec:
#     display_name: Python 3
#     language: python
#     name: python3
# ---

# +
import sys

#sys.path.append('/home/franzihe/Documents/Python/Thesis')
sys.path.append('/uio/kant/geo-metos-u1/franzihe/Documents/Thesis/Python')
import time
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
# import fill_values as fv
# import calc_station_properties as cs

import createFolder as cF

import fcts_read_stat as rs
import gc


# +
thredds = 'http://thredds.met.no/thredds/dodsC/metusers/bjorgjke-3mnd_ws/'
#run = 'XCCR'
run = 'CTRL'

stn_name = 'Haukeliseter'
stn_lat = 59.81
stn_lon = 7.21

forecasttime = '00'
m = [#'12',
'01','02']
# -

main_dir = '../../../Data/MEPS'
dirnc = '%s/%s/%s' % (main_dir, stn_name, run)


def mask_array(variable, y, x):
    if np.ma.is_masked(variable[:, :, y, x]):
        mask = np.ma.getmaskarray(variable[:, :, y, x])
        fill_value = np.nan
        marr = np.ma.array(variable[:, :, y, x],
                           mask=mask,
                           fill_value=fill_value)
        dtype = marr.filled().dtype
        filled = marr.filled()
    else:
        fill_value = np.nan
        marr = variable[:, :, y, x]
        filled = marr
        dtype = marr.dtype

    return (filled, dtype)


def get_value_at_station(fn, variable, y, x):
    variable = fn.variables[variable]
    variable, dtype = mask_array(variable, y, x)
    variable = np.fliplr(variable)
    variable = np.ma.masked_where(np.isnan(variable), variable)
    return (variable)


def read_for_station(thredds, year, month, day, forecasttime, stn_lat, stn_lon, dirnc):
    try:
        fn = netCDF4.Dataset('%s/%s/fc%s%s%s%s.nc' % (thredds, run, year, month, day, forecasttime), 'r')
    except OSError:
        print('no file found: %s/%s/fc%s%s%s%s.nc' % (thredds, run, year, month, day, forecasttime))
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
    surface_air_pressure, dtype_sap = mask_array(surface_air_pressure, y, x)
    air_temperature_0m, dtype_at0m = mask_array(air_temperature_0m, y, x)
    air_temperature_ml, dtype_atml = mask_array(air_temperature_ml, y, x)

    ### Transfer pressure coordinate
    p_interface = ap[:] + b[:] * surface_air_pressure[:, :]

    ### 1) Connect interface values and surface values for pressure
    ### Pressure
    p_interfaces2 = np.concatenate((p_interface[:, :], surface_air_pressure[:, :]), axis=1)

    # transform hybrid sigma pressure coordinates at interface levels to pressure at model levels
    p_ml = np.empty([p_interfaces2.shape[0], p_interfaces2.shape[1] - 1])
    for i in range(0, p_interfaces2.shape[1] - 1):
        p_ml[:, i] = 1 / 2 * (p_interfaces2[:, i] + p_interfaces2[:, i + 1])
    p_ml = np.concatenate((p_ml[:, :], surface_air_pressure[:, :]), axis=1)

    ### --> Now the pressure is calculated at each model level

    ### 2) Connect model levels and surface values for temperature
    ### Temperature
    temperature_ml = np.concatenate((air_temperature_ml[:, :], air_temperature_0m[:, :]), axis=1)

    dz, dgeop = rs.get_thickness(p_ml, temperature_ml)

    pressure_in_modellev = p_ml[:, 0:-1]  # does not include the surface air pressure
    thickness_m = dz
    thickness_phi = dgeop

    pressure_in_modellev = np.ma.array(pressure_in_modellev,
                                       mask=np.ma.is_masked(pressure_in_modellev),
                                       fill_value=np.nan)
    thickness_m = np.ma.array(thickness_m,
                              mask=np.ma.is_masked(thickness_m),
                              fill_value=np.nan)
    thickness_phi = np.ma.array(thickness_phi,
                                mask=np.ma.is_masked(thickness_phi),
                                fill_value=np.nan)

    # Read in all values needed to present the microphysics
    ## Time
    time_arr = fn.variables['time']
    # ensemble_member_arr = fn.variables['ensemble_member']

    ## heights
    hybrid_arr = fn.variables['hybrid']
    height0_arr = fn.variables['height0']
    height1_arr = fn.variables['height1']
    height3_arr = fn.variables['height3']
    height_above_msl_arr = fn.variables['height_above_msl']

    ######## with Vertical Levels ( height0 ) #################################
    air_temperature_0m = get_value_at_station(fn, 'air_temperature_0m', y, x)
    liquid_water_content_of_surface_snow = get_value_at_station(fn, 'liquid_water_content_of_surface_snow', y, x)
    rainfall_amount = get_value_at_station(fn, 'rainfall_amount', y, x)
    snowfall_amount = get_value_at_station(fn, 'snowfall_amount', y, x)
    graupelfall_amount = get_value_at_station(fn, 'graupelfall_amount', y, x)
    surface_air_pressure = get_value_at_station(fn, 'surface_air_pressure', y, x)
    surface_geopotential = get_value_at_station(fn, 'surface_geopotential', y, x)
    precipitation_amount_acc = get_value_at_station(fn, 'precipitation_amount_acc', y, x)
    integral_of_snowfall_amount_wrt_time = get_value_at_station(fn, 'integral_of_snowfall_amount_wrt_time', y, x)
    integral_of_rainfall_amount_wrt_time = get_value_at_station(fn, 'integral_of_rainfall_amount_wrt_time', y, x)
    integral_of_graupelfall_amount_wrt_time = get_value_at_station(fn, 'integral_of_graupelfall_amount_wrt_time', y, x)
    surface_snow_sublimation_amount_acc = get_value_at_station(fn, 'surface_snow_sublimation_amount_acc', y, x)

    ######## with Vertical Levels ( height1 )#################################
    air_temperature_2m = get_value_at_station(fn, 'air_temperature_2m', y, x)
    relative_humidity_2m = get_value_at_station(fn, 'relative_humidity_2m', y, x)
    specific_humidity_2m = get_value_at_station(fn, 'specific_humidity_2m', y, x)

    ######## with Vertical Levels ( height3 )#################################
    x_wind_10m = get_value_at_station(fn, 'x_wind_10m', y, x)
    y_wind_10m = get_value_at_station(fn, 'y_wind_10m', y, x)

    ######## with Vertical Levels ( height_above_msl )#################################
    air_pressure_at_sea_level = get_value_at_station(fn, 'air_pressure_at_sea_level', y, x)

    ######## with Vertical Levels ( hybrid )#################################
    ## hybrid levels atmosphere_hybrid_sigma_pressure_coordinate
    # formula: p(n,k,j,i) = ap(k) + b(k)*ps(n,j,i)
    # positive: down
    specific_humidity_ml = get_value_at_station(fn, 'specific_humidity_ml', y, x)
    atmosphere_cloud_condensed_water_content_ml = get_value_at_station(fn,
                                                                       'mass_fraction_of_cloud_condensed_water_in_air_ml',
                                                                       y, x)
    atmosphere_cloud_ice_content_ml = get_value_at_station(fn, 'mass_fraction_of_cloud_ice_in_air_ml', y, x)
    atmosphere_cloud_snow_content_ml = get_value_at_station(fn, 'mass_fraction_of_snow_in_air_ml', y, x)
    atmosphere_cloud_rain_content_ml = get_value_at_station(fn, 'mass_fraction_of_rain_in_air_ml', y, x)
    atmosphere_cloud_graupel_content_ml = get_value_at_station(fn, 'mass_fraction_of_graupel_in_air_ml', y, x)
    pressure_departure = get_value_at_station(fn, 'pressure_departure', y, x)
    air_temperature_ml = get_value_at_station(fn, 'air_temperature_ml', y, x)
    x_wind_ml = get_value_at_station(fn, 'x_wind_ml', y, x)
    y_wind_ml = get_value_at_station(fn, 'y_wind_ml', y, x)

    ### write netCDF file

    f = netCDF4.Dataset('%s/%s%s%s_%s.nc' % (dirnc, year, month, day, forecasttime), 'w')

    ### create dimensions
    f.createDimension('time', time_arr.shape[0])
    f.createDimension('hybrid', hybrid_arr.shape[0])
    f.createDimension('height0', height0_arr.shape[0])
    f.createDimension('height1', height1_arr.shape[0])
    f.createDimension('height3', height3_arr.shape[0])
    f.createDimension('height_above_msl', height_above_msl_arr.shape[0])

    t = f.createVariable('time', time_arr.dtype, 'time', zlib=True)
    t[:] = time_arr[:]

    ######## with Vertical Levels ( height0 ) #################################
    h = f.createVariable('height0', height0_arr.dtype, 'height0', zlib=True)
    h[:] = height0_arr[:]
    dim = ('time', 'height0',)

    at_0m = rs.get_netCDF_variable(f, 'air_temperature_0m', air_temperature_0m, dim)
    lwc_0m = rs.get_netCDF_variable(f, 'liquid_water_content_of_surface_snow', liquid_water_content_of_surface_snow,
                                    dim)
    ra_0m = rs.get_netCDF_variable(f, 'rainfall_amount', rainfall_amount, dim)
    sa_0m = rs.get_netCDF_variable(f, 'snowfall_amount', snowfall_amount, dim)
    ga_0m = rs.get_netCDF_variable(f, 'graupelfall_amount', graupelfall_amount, dim)
    ps = rs.get_netCDF_variable(f, 'surface_air_pressure', surface_air_pressure, dim)
    geop = rs.get_netCDF_variable(f, 'surface_geopotential', surface_geopotential, dim)
    pr_0m = rs.get_netCDF_variable(f, 'precipitation_amount_acc', precipitation_amount_acc, dim)
    int_snow_wrt_time = rs.get_netCDF_variable(f, 'integral_of_snowfall_amount_wrt_time',
                                               integral_of_snowfall_amount_wrt_time, dim)
    int_rain_wrt_time = rs.get_netCDF_variable(f, 'integral_of_rainfall_amount_wrt_time',
                                               integral_of_rainfall_amount_wrt_time, dim)
    int_grauple_wrt_time = rs.get_netCDF_variable(f, 'integral_of_graupelfall_amount_wrt_time',
                                                  integral_of_graupelfall_amount_wrt_time, dim)
    sfc_snow_sub = rs.get_netCDF_variable(f, 'surface_snow_sublimation_amount_acc', surface_snow_sublimation_amount_acc,
                                          dim)

    ###### with Vertical Levels ( height1 )#################################
    h1 = f.createVariable('height1', height1_arr.dtype, 'height1', zlib=True)
    h1[:] = height1_arr[:]
    dim = ('time', 'height1',)

    at_2m = rs.get_netCDF_variable(f, 'air_temperature_2m', air_temperature_2m, dim)
    rel_2m = rs.get_netCDF_variable(f, 'relative_humidity_2m', relative_humidity_2m, dim)
    sh_2m = rs.get_netCDF_variable(f, 'specific_humidity_2m', specific_humidity_2m, dim)

    ######## with Vertical Levels( height3 )#################################
    h3 = f.createVariable('height3', height3_arr.dtype, 'height3', zlib=True)
    h3[:] = height3_arr[:]
    dim = ('time', 'height3',)

    xwind_10m = rs.get_netCDF_variable(f, 'x_wind_10m', x_wind_10m, dim)
    ywind_10m = rs.get_netCDF_variable(f, 'y_wind_10m', y_wind_10m, dim)

    ######## with Vertical Levels ( height_above_msl )#################################
    h_asl = f.createVariable('height_above_msl', height_above_msl_arr.dtype, 'height_above_msl', zlib=True)
    h_asl[:] = height_above_msl_arr[:]
    dim = ('time', 'height_above_msl',)

    pressure_sea_level = rs.get_netCDF_variable(f, 'air_pressure_at_sea_level', air_pressure_at_sea_level, dim)

    ######## with Vertical Levels ( hybrid )#################################
    hyb = f.createVariable('hybrid', hybrid_arr.dtype, 'hybrid', zlib=True)
    hyb[:] = hybrid_arr[:]
    dim = ('time', 'hybrid',)

    sh_ml = rs.get_netCDF_variable(f, 'specific_humidity_ml', specific_humidity_ml, dim)
    ccw_ml = rs.get_netCDF_variable(f, 'atmosphere_cloud_condensed_water_content_ml',
                                    atmosphere_cloud_condensed_water_content_ml, dim)
    cic_ml = rs.get_netCDF_variable(f, 'atmosphere_cloud_ice_content_ml', atmosphere_cloud_ice_content_ml, dim)
    csc_ml = rs.get_netCDF_variable(f, 'atmosphere_cloud_snow_content_ml', atmosphere_cloud_snow_content_ml, dim)
    crc_ml = rs.get_netCDF_variable(f, 'atmosphere_cloud_rain_content_ml', atmosphere_cloud_rain_content_ml, dim)
    cgc_ml = rs.get_netCDF_variable(f, 'atmosphere_cloud_graupel_content_ml', atmosphere_cloud_graupel_content_ml, dim)
    pres_dep_ml = rs.get_netCDF_variable(f, 'pressure_departure', pressure_departure, dim)
    at_ml = rs.get_netCDF_variable(f, 'air_temperature_ml', air_temperature_ml, dim)
    xwind_ml = rs.get_netCDF_variable(f, 'x_wind_ml', x_wind_ml, dim)
    ywind_ml = rs.get_netCDF_variable(f, 'y_wind_ml', y_wind_ml, dim)

    pres_ml = rs.get_netCDF_variable(f, 'pressure_ml', pressure_in_modellev, dim)
    dz_ml = rs.get_netCDF_variable(f, 'layer_thickness', thickness_m, dim)
    dgeop_ml = rs.get_netCDF_variable(f, 'geop_layer_thickness', thickness_phi, dim)

    f.close()
    print('file written: %s/%s%s%s_%s.nc' % (dirnc, year, month, day, forecasttime))


for month in m:
    if month == '12' or month == '01':
        t = np.arange(30, 32)
    if month == '02':
        t = np.arange(1, 29)
    if month == '12':
        year = '2016'
    if month == '01' or month == '02':
        year = '2017'
    for day in t:
        if day < 10:
            day = '0%s' % day

        start_time = time.time()

        ### direction where files should be saved
        cF.createFolder(dirnc)
        read_for_station(thredds, year, month, day, forecasttime, stn_lat, stn_lon, dirnc)
        print('--- %s seconds ---' % round(time.time() - start_time, 2))
