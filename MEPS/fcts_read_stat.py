
# coding: utf-8

# Read in and download MEPS data

# In[1]:


import sys
sys.path.append('/uio/kant/geo-metos-u7/franzihe/Documents/Thesis/Python')
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import datetime
import pandas as pd
#import fill_values as fv
#import calc_station_properties as cs

import createFolder as cF
from scipy.integrate import simps



# In[2]:


def find_station_yx(latitude, longitude, stn_lat, stn_lon):
# find the absolute value of the difference between the  station's lat/lon with every point in the grid. 
# This tells us how close a point is to the particular latitude and longitude.
    abslat = np.abs(latitude[:,:]-stn_lat)
    abslon = np.abs(longitude[:,:]-stn_lon)

# Now we need to combine these two results. We will use numpy.maximum, which takes two arrays and finds the local 
# maximum.
    c = np.maximum(abslon, abslat)

# If you don't like flattened arrays, you can also get the row/column index like this
    y, x = np.where(c == np.min(c))
    return(x,y);


# In[3]:


def mask_array(variable, #ens_memb, 
               y, x, EM):
    if EM == 10:
        if np.ma.is_masked(variable[:,:,:,#ens_memb,
                                    y,x]):
            mask = np.ma.getmaskarray(variable[:,:,:,#ens_memb,
                                               y,x])  
            fill_value = np.nan
            marr = np.ma.array(variable[:,:,:,#ens_memb,
                                        y,x], mask = mask, fill_value = fill_value)
            dtype = marr.filled().dtype
            filled = marr.filled()
        else:
            fill_value = np.nan
            marr = variable[:,:,:,#ens_memb,
                        y,x]
            filled = marr
            dtype = marr.dtype
    elif EM == 1:
        if np.ma.is_masked(variable[:,:,:,y,x]):
            mask = np.ma.getmaskarray(variable[:,:,:,y,x])  
            fill_value = np.nan
            marr = np.ma.array(variable[:,:,:,y,x], mask = mask, fill_value = fill_value)
            dtype = marr.filled().dtype
            filled = marr.filled()
        else:
            fill_value = np.nan
            marr = variable[:,:,:,y,x]
            filled = marr
            dtype = marr.dtype
    return(filled, dtype)
    

def mask_array2(var):
    var = np.ma.masked_array(var#[:,:,:]
                             ,mask = np.ma.getmaskarray(var),fill_value=np.nan)
    return(var)
# In[4]:


def get_thickness(pressure, temperature):
    '''Calculuate the thickness of each layer in m or geop height'''

### 3) to convert pressure-levels into actual heights use the hypsometric equation --> Temperature and pressure
# are needed. After J. E. Martin: Mid-Latitude Atmospheric Dynamics Eq. 3.6
    Rd = 287.    # gas constant for dry air [J kg^-1 K^-1]
    g = 9.81     # Standard gravity [m s^-2]
    
### calculate the pressure-weighted, column averaged temperature as from J. E. Martin: Mid-Latitude Atmospheric 
# Dynamics Book, Eq. below Eq. 3.6
    temp_mean = []
    for i in range(0, temperature.shape[1]-1):
        numT = simps(y=temperature[:,i:(i+2)], x=np.log(pressure[:,i:(i+2)]), 
                     dx = np.log(pressure[:,i:(i+2)]),even='last')
        denomT = simps(y=np.ones(temperature[:,i:(i+2)].shape), x = np.log(pressure[:,i:(i+2)]), 
                       dx = np.log(pressure[:,i:(i+2)]),even='last')
        t_mean = numT/denomT
        temp_mean.append(t_mean)

# get temperature and pressure, and value so that array zero contains low levels (transpose or flip)
    temp_mean = (np.fliplr(np.transpose(temp_mean)))
    pressure = np.fliplr(pressure)
    temperature = np.fliplr(temperature)

# calculate layer thickness
    thickness = []
    geop_th = []
    for i in range(0, temp_mean.shape[1]):
   # if (i+1) == pres.shape[1]:
    #    continue
        p1 = pressure[:,i]
        p2 = pressure[:,(i+1)]   
        dz = (Rd * temp_mean[:,i])/g * np.log((p1/p2))    # thickness in [m]
        dgeop = (Rd * temp_mean[:,i])* np.log((p1/p2))    # thickness in [J/kg]
        thickness.append(dz)
        geop_th.append(dgeop)

# transpose array, so that array starts at lower levels
    thickness = np.transpose(thickness)
    geop_th = np.transpose(geop_th)
    
    return(thickness, geop_th)


# In[5]:


def get_value_at_station(fn, variable, y, x):
    variable = fn.variables[variable]
    variable, dtype = mask_array(variable,y,x,variable.shape[2])
    variable = np.fliplr(variable)
    variable = np.ma.masked_where(np.isnan(variable), variable)
    return(variable)


# In[6]:


def get_netCDF_variable(f, var_name, var, dim):
    v_0m = f.createVariable(varname=var_name, datatype=var.dtype, dimensions=dim,
                           fill_value=var.fill_value, zlib=True)
    v_0m[:] = var[:]
    return(v_0m)

# def fill_nan(var):
    mask = np.ma.getmaskarray(var)
    marr = np.ma.array(var, mask = mask, fill_value = np.nan)
    x = marr.filled(np.nan)
    return(x)
# In[7]:




