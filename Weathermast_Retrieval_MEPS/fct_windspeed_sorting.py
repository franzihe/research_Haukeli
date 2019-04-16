#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import sys
sys.path.append('/uio/kant/geo-metos-u7/franzihe/Documents/Thesis/Python')

import numpy as np
import netCDF4
import pandas as pd
import matplotlib.pyplot as plt

import createFolder as cF
np.warnings.filterwarnings('ignore')


# In[ ]:
m = ['11', '12', '01', '02', '03']

def get_array_of_dictionary(obs_precip):
    vals_west = []
    for v in obs_precip.items():
        vals_west.append(v)
    vals_west = np.array(vals_west)
    return vals_west


# In[ ]:


def get_array_of_dictionary_MEPS(meps_precip):
    vals_west = []
    for v in meps_precip.values():
        vals_west.append(v)
    vals_west = np.array(vals_west)

    keys_west = []
    for k in meps_precip.keys():
        keys_west.append(k)
    keys_west = np.array(keys_west)

    arr = np.empty(shape=(vals_west.shape[0], vals_west.shape[1]+1))
    arr[:] = np.nan
    arr[:,0] = np.array(keys_west)
    for ens_memb in range(10):
        arr[:,ens_memb+1] = vals_west[:,ens_memb]
        
    return arr


# In[ ]:


def get_precip_amount_for_wind_speed(obs_precip_diff, obs_wind_dir, obs_west_idx, obs_east_idx, obs_wind_speed):
    
    west_WS_0_4 = dict()
    east_WS_0_4 = dict()
    west_WS_4_8 = dict()
    east_WS_4_8 = dict()
    west_WS_8_12 = dict()
    east_WS_8_12 = dict()
    west_WS_12_16 = dict()
    east_WS_12_16 = dict()
    west_WS_16_20 = dict()
    east_WS_16_20 = dict()
    west_WS_20 = dict()
    east_WS_20 = dict()

    precip_west_0_4 = dict()
    precip_east_0_4 = dict()
    precip_west_4_8 = dict()
    precip_east_4_8 = dict()
    precip_west_8_12 = dict()
    precip_east_8_12 = dict()
    precip_west_12_16 = dict()
    precip_east_12_16 = dict()
    precip_west_16_20 = dict()
    precip_east_16_20 = dict()
    precip_west_20 = dict()
    precip_east_20 = dict()

    for month in m:
        if month == '11':
    #        t = np.arange(8,31)
            t = [9, 12, 13, 14, 15, 16, 17, 23, 24, 26,29,30]
        if month == '12': #or  or 
     #       t = np.arange(1,32)
            t = [7, 9, 10, 11, 12, 13, 15, 21, 22, 23, 24, 25, 26, 29, 30, 31]
        if month == '01':
            t = [2, 3, 6, 8, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 23, 24, 25, 28, 29]
        if month == '02':
      #      t = np.arange(1,29)
            t = [2, 3, 4, 6]
        if month == '03':
            t = [9, 14, 16, 17, 19, 21, 24, 30, 31]
        if month == '11' or month == '12':
            year = '2016'
        if month == '01' or month == '02' or month == '03':
            year = '2017'

        for day in t:

            if day < 10:
                day = '0%s' %day
            else:
                day = '%s' %day

            Date = int(year+month+day)

            if np.count_nonzero(~np.isnan(obs_wind_dir[Date][:])) != 1440 and                np.count_nonzero(~np.isnan(obs_wind_dir[Date][:])) != 24:
            ######################
    #            print([Date], np.count_nonzero(~np.isnan(obs_wind_dir[Date][:])))
                west_WS_0_4[Date] = np.empty(shape=(1,))
                west_WS_0_4[Date][:] = np.nan    ## west 0-4
                east_WS_0_4[Date] = np.empty(shape=(1,))
                east_WS_0_4[Date][:] = np.nan    ## east 0-4


                precip_west_0_4[Date] = np.nan
                precip_east_0_4[Date] = np.nan
             ############################   
                west_WS_4_8[Date] = np.empty(shape=(1,))
                west_WS_4_8[Date][:] = np.nan    ## west 4-8
                east_WS_4_8[Date] = np.empty(shape=(1,))
                east_WS_4_8[Date][:] = np.nan    ## east 4-8


                precip_west_4_8[Date] = np.nan
                precip_east_4_8[Date] = np.nan
            ##############################    
                west_WS_8_12[Date] = np.empty(shape=(1,))
                west_WS_8_12[Date][:] = np.nan    ## west 8-12
                east_WS_8_12[Date] = np.empty(shape=(1,))
                east_WS_8_12[Date][:] = np.nan    ## east 8-12

                precip_west_8_12[Date] = np.nan
                precip_east_8_12[Date] = np.nan
             ################################   
                west_WS_12_16[Date] = np.empty(shape=(1,))
                west_WS_12_16[Date][:] = np.nan    ## west 12-16
                east_WS_12_16[Date] = np.empty(shape=(1,))
                east_WS_12_16[Date][:] = np.nan    ## east 12-16

                precip_west_12_16[Date] = np.nan
                precip_east_12_16[Date] = np.nan
              #####################################  
                west_WS_16_20[Date] = np.empty(shape=(1,))
                west_WS_16_20[Date][:] = np.nan    ## west 16-20
                east_WS_16_20[Date] = np.empty(shape=(1,))
                east_WS_16_20[Date][:] = np.nan    ## east 16-20

                precip_west_16_20[Date] = np.nan
                precip_east_16_20[Date] = np.nan
              #####################################  

                west_WS_20[Date] = np.empty(shape=(1,))
                west_WS_20[Date][:] = np.nan    ## west > 20
                east_WS_20[Date] = np.empty(shape=(1,))
                east_WS_20[Date][:] = np.nan    ## east > 20

                precip_west_20[Date] = np.nan
                precip_east_20[Date] = np.nan

            else:
     #           print([Date])
                IDX2 = np.arange(np.array(obs_west_idx[Date]).shape[0])
                IDX3 = np.arange(np.array(obs_east_idx[Date]).shape[0])
                ######################
                ## 0 - 4 m/s
                if len(obs_west_idx[Date]) == 0.:
                        west_WS_0_4[Date] = np.empty(shape=(1,))
                        west_WS_0_4[Date][:] = np.nan    ## west 0-4

                        precip_west_0_4[Date] = np.nan
                else:
                        west_WS_0_4[Date] = IDX2[np.logical_and(    ## west 0-4
                                   obs_wind_speed[Date][obs_west_idx[Date],] >= 0, 
                                   obs_wind_speed[Date][obs_west_idx[Date],] < 4 )]
                        precip_west_0_4[Date] = np.nansum(obs_precip_diff[Date]
                                                                     [obs_west_idx[Date][west_WS_0_4[Date]]])

                if len(obs_east_idx[Date]) == 0.:
                        east_WS_0_4[Date] = np.empty(shape=(1,))
                        east_WS_0_4[Date][:] = np.nan    ## east 0-4

                        precip_east_0_4[Date] = np.nan
                else:
                        east_WS_0_4[Date] = IDX3[np.logical_and(
                                    obs_wind_speed[Date][obs_east_idx[Date],] >= 0,
                                    obs_wind_speed[Date][obs_east_idx[Date],] <4)]
                        precip_east_0_4[Date] = np.nansum(obs_precip_diff[Date]
                                                                      [obs_east_idx[Date][east_WS_0_4[Date]]])



                ######################
                ## 4 - 8 m/s
                if len(obs_west_idx[Date]) == 0.:
                        west_WS_4_8[Date] = np.empty(shape=(1,))
                        west_WS_4_8[Date][:] = np.nan    ## west 4-8

                        precip_west_4_8[Date] = np.nan
                else:
                        west_WS_4_8[Date] = IDX2[np.logical_and(    ## west 4-8
                                   obs_wind_speed[Date][obs_west_idx[Date],] >= 4, 
                                   obs_wind_speed[Date][obs_west_idx[Date],] < 8 )]
                        precip_west_4_8[Date] = np.nansum(obs_precip_diff[Date]
                                                                     [obs_west_idx[Date][west_WS_4_8[Date]]])
                if len(obs_east_idx[Date]) == 0.:
                        east_WS_4_8[Date] = np.empty(shape=(1,))
                        east_WS_4_8[Date][:] = np.nan    ## east 4-8

                        precip_east_4_8[Date] = np.nan
                else:
                        east_WS_4_8[Date] = IDX3[np.logical_and(
                                    obs_wind_speed[Date][obs_east_idx[Date],] >= 4,
                                    obs_wind_speed[Date][obs_east_idx[Date],] <8)]
                        precip_east_4_8[Date] = np.nansum(obs_precip_diff[Date]
                                                                      [obs_east_idx[Date][east_WS_4_8[Date]]])



                ######################
                ## 8 - 12 m/s
                if len(obs_west_idx[Date]) == 0.:
                        west_WS_8_12[Date] = np.empty(shape=(1,))
                        west_WS_8_12[Date][:] = np.nan    ## west 8-12

                        precip_west_8_12[Date] = np.nan
                else:
                        west_WS_8_12[Date] = IDX2[np.logical_and(    ## west 8-12
                                   obs_wind_speed[Date][obs_west_idx[Date],] >= 8, 
                                   obs_wind_speed[Date][obs_west_idx[Date],] < 12 )]
                        precip_west_8_12[Date] = np.nansum(obs_precip_diff[Date]
                                                                     [obs_west_idx[Date][west_WS_8_12[Date]]])
                if len(obs_east_idx[Date]) == 0.:
                        east_WS_8_12[Date] = np.empty(shape=(1,))
                        east_WS_8_12[Date][:] = np.nan    ## east 8-12

                        precip_east_8_12[Date] = np.nan
                else:
                        east_WS_8_12[Date] = IDX3[np.logical_and(
                                    obs_wind_speed[Date][obs_east_idx[Date],] >= 8,
                                    obs_wind_speed[Date][obs_east_idx[Date],] <12)]
                        precip_east_8_12[Date] = np.nansum(obs_precip_diff[Date]
                                                                      [obs_east_idx[Date][east_WS_8_12[Date]]])

                ######################
                ## 12 - 16 m/s
                if len(obs_west_idx[Date]) == 0.:
                        west_WS_12_16[Date] = np.empty(shape=(1,))
                        west_WS_12_16[Date][:] = np.nan    ## west 12-16

                        precip_west_12_16[Date] = np.nan
                else:
                        west_WS_12_16[Date] = IDX2[np.logical_and(    ## west 12-16
                                   obs_wind_speed[Date][obs_west_idx[Date],] >= 12, 
                                   obs_wind_speed[Date][obs_west_idx[Date],] < 16 )]
                        precip_west_12_16[Date] = np.nansum(obs_precip_diff[Date]
                                                                     [obs_west_idx[Date][west_WS_12_16[Date]]])
                if len(obs_east_idx[Date]) == 0.:
                        east_WS_12_16[Date] = np.empty(shape=(1,))
                        east_WS_12_16[Date][:] = np.nan    ## east 12-16

                        precip_east_12_16[Date] = np.nan
                else:
                        east_WS_12_16[Date] = IDX3[np.logical_and(
                                    obs_wind_speed[Date][obs_east_idx[Date],] >= 12,
                                    obs_wind_speed[Date][obs_east_idx[Date],] <16)]
                        precip_east_12_16[Date] = np.nansum(obs_precip_diff[Date]
                                                                      [obs_east_idx[Date][east_WS_12_16[Date]]])




               ######################
                ## 16 - 20 m/s
                if len(obs_west_idx[Date]) == 0.:
                        west_WS_16_20[Date] = np.empty(shape=(1,))
                        west_WS_16_20[Date][:] = np.nan    ## west 16-20

                        precip_west_16_20[Date] = np.nan
                else:
                        west_WS_16_20[Date] = IDX2[np.logical_and(    ## west 16-20
                                   obs_wind_speed[Date][obs_west_idx[Date],] >= 16, 
                                   obs_wind_speed[Date][obs_west_idx[Date],] < 20 )]
                        precip_west_16_20[Date] = np.nansum(obs_precip_diff[Date]
                                                                     [obs_west_idx[Date][west_WS_16_20[Date]]])
                if len(obs_east_idx[Date]) == 0.:
                        east_WS_16_20[Date] = np.empty(shape=(1,))
                        east_WS_16_20[Date][:] = np.nan    ## east 16-20

                        precip_east_16_20[Date] = np.nan
                else:
                        east_WS_16_20[Date] = IDX3[np.logical_and(
                                    obs_wind_speed[Date][obs_east_idx[Date],] >= 16,
                                    obs_wind_speed[Date][obs_east_idx[Date],] <20)]
                        precip_east_16_20[Date] = np.nansum(obs_precip_diff[Date]
                                                                      [obs_east_idx[Date][east_WS_16_20[Date]]])


                ######################
                ## 20 - inf m/s
                if len(obs_west_idx[Date]) == 0.:
                        west_WS_20[Date] = np.empty(shape=(1,))
                        west_WS_20[Date][:] = np.nan    ## west > 20

                        precip_west_20[Date] = np.nan
                else:
                        west_WS_20[Date] = IDX2[(       ## west >20
                                   obs_wind_speed[Date][obs_west_idx[Date],] >= 20)]
                        precip_west_20[Date] = np.nansum(obs_precip_diff[Date]
                                                                     [obs_west_idx[Date][west_WS_20[Date]]])
                if len(obs_east_idx[Date]) == 0.:
                        east_WS_20[Date] = np.empty(shape=(1,))
                        east_WS_20[Date][:] = np.nan    ## east > 20

                        precip_east_20[Date] = np.nan
                else:
                        east_WS_20[Date] = IDX3[(
                                    obs_wind_speed[Date][obs_east_idx[Date],] >= 20)]
                        precip_east_20[Date] = np.nansum(obs_precip_diff[Date]
                                                                      [obs_east_idx[Date][east_WS_20[Date]]])
                        
    vals_west_0_4 = get_array_of_dictionary(precip_west_0_4)
    vals_east_0_4 = get_array_of_dictionary(precip_east_0_4)
    
    vals_west_4_8 = get_array_of_dictionary(precip_west_4_8)
    vals_east_4_8 = get_array_of_dictionary(precip_east_4_8)
    
    vals_west_8_12 = get_array_of_dictionary(precip_west_8_12)
    vals_east_8_12 = get_array_of_dictionary(precip_east_8_12)
    
    vals_west_12_16 = get_array_of_dictionary(precip_west_12_16)
    vals_east_12_16 = get_array_of_dictionary(precip_east_12_16)
    
    vals_west_16_20 = get_array_of_dictionary(precip_west_16_20)
    vals_east_16_20 = get_array_of_dictionary(precip_east_16_20)
    
    vals_west_20 = get_array_of_dictionary(precip_west_20)
    vals_east_20 = get_array_of_dictionary(precip_east_20)
    

    return vals_west_0_4, vals_east_0_4, vals_west_4_8, vals_east_4_8 , vals_west_8_12, vals_east_8_12, vals_west_12_16, vals_east_12_16, vals_west_16_20, vals_east_16_20,vals_west_20, vals_east_20


# In[ ]:


def get_precip_amount_for_wind_speed_MEPS(meps_precip_diff, meps_wind_dir, meps_west_idx, meps_east_idx, meps_wind_speed):
    
    west_WS_0_4 = dict()
    east_WS_0_4 = dict()
    west_WS_4_8 = dict()
    east_WS_4_8 = dict()
    west_WS_8_12 = dict()
    east_WS_8_12 = dict()
    west_WS_12_16 = dict()
    east_WS_12_16 = dict()
    west_WS_16_20 = dict()
    east_WS_16_20 = dict()
    west_WS_20 = dict()
    east_WS_20 = dict()

    precip_west_0_4 = dict()
    precip_east_0_4 = dict()
    precip_west_4_8 = dict()
    precip_east_4_8 = dict()
    precip_west_8_12 = dict()
    precip_east_8_12 = dict()
    precip_west_12_16 = dict()
    precip_east_12_16 = dict()
    precip_west_16_20 = dict()
    precip_east_16_20 = dict()
    precip_west_20 = dict()
    precip_east_20 = dict()
    
    for month in m:
        if month == '11':
    #        t = np.arange(8,31)
            t = [9, 12, 13, 14, 15, 16, 17, 23, 24, 26,29,30]
        if month == '12': #or  or 
     #       t = np.arange(1,32)
            t = [7, 9, 10, 11, 12, 13, 15, 21, 22, 23, 24, 25, 26, 29, 30, 31]
        if month == '01':
            t = [2, 3, 6, 8, 9, 10, 11, 12, 13, 16, 17, 18, 19, 20, 23, 24, 25, 28, 29]
        if month == '02':
      #      t = np.arange(1,29)
            t = [2, 3, 4, 6]
        if month == '03':
            t = [9, 14, 16, 17, 19, 21, 24, 30, 31]
        if month == '11' or month == '12':
            year = '2016'
        if month == '01' or month == '02' or month == '03':
            year = '2017'
        
        for day in t:

            if day < 10:
                    day = '0%s' %day
            else:
                    day = '%s' %day

            Date = int(year+month+day)
            ############################   
            precip_west_0_4[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_west_0_4[Date][:] = np.nan
            precip_east_0_4[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_east_0_4[Date][:] = np.nan
            ############################   
            precip_west_4_8[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_west_4_8[Date][:] = np.nan
            precip_east_4_8[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_east_4_8[Date][:] = np.nan
            ############################   
            precip_west_8_12[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_west_8_12[Date][:] = np.nan
            precip_east_8_12[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_east_8_12[Date][:] = np.nan
            ############################   
            precip_west_12_16[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_west_12_16[Date][:] = np.nan
            precip_east_12_16[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_east_12_16[Date][:] = np.nan
            ############################   
            precip_west_16_20[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_west_16_20[Date][:] = np.nan
            precip_east_16_20[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_east_16_20[Date][:] = np.nan
            ############################   
            precip_west_20[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_west_20[Date][:] = np.nan
            precip_east_20[Date] = np.empty(shape=(meps_precip_diff[Date].shape[1],))
            precip_east_20[Date][:] = np.nan
        
            for ens_memb in range(10):
            

               

                
                IDX2 = np.array(np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])))

                IDX3 = np.array(np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])))

                ######################
                ## 0 - 4 m/s
                if len(meps_west_idx[Date][:,ens_memb]) == 0.:    ## west 0-4
                        west_WS_0_4[Date] = np.empty(shape=(1,))
                        west_WS_0_4[Date][:] = np.nan

                        precip_west_0_4[Date][ens_memb,] = np.nan
                else:
                        west_WS_0_4[Date] = IDX2[0, np.logical_and(   ## west 0-4
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] >=0,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] < 4)[0,:]]

                        precip_west_0_4[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][west_WS_0_4[Date],ens_memb])

                if len(meps_east_idx[Date][:,ens_memb]) == 0.:     ## east 0-4
                        east_WS_0_4[Date] = np.empty(shape=(1,))
                        east_WS_0_4[Date][:] = np.nan

                        precip_east_0_4[Date][ens_memb,] = np.nan

                else:
                        east_WS_0_4[Date] = IDX3[0, np.logical_and(   ## east 0-4
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] >=0,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] < 4)[0,:]]

                        precip_east_0_4[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][east_WS_0_4[Date],ens_memb])

                ######################
                ## 4 - 8 m/s
                if len(meps_west_idx[Date][:,ens_memb]) == 0.:   ## west 4-8
                        west_WS_4_8[Date] = np.empty(shape=(1,))
                        west_WS_4_8[Date][:] = np.nan

                        precip_west_4_8[Date][ens_memb,] = np.nan     
                else:
                        west_WS_4_8[Date] = IDX2[0, np.logical_and(  ## west 4-8
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] >=4,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] < 8)[0,:]]

                        precip_west_4_8[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][west_WS_4_8[Date],ens_memb])

                if len(meps_east_idx[Date][:,ens_memb]) == 0.:   ## east 4-8
                        east_WS_4_8[Date] = np.empty(shape=(1,))
                        east_WS_4_8[Date][:] = np.nan

                        precip_east_4_8[Date][ens_memb,] = np.nan     
                else:
                        east_WS_4_8[Date] = IDX3[0, np.logical_and(  ## east 4-8
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] >=4,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] < 8)[0,:]]

                        precip_east_4_8[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][east_WS_4_8[Date],ens_memb])

                ######################
                ## 8 - 12 m/s
                if len(meps_west_idx[Date][:,ens_memb]) == 0.:   ## west 8-12
                        west_WS_8_12[Date] = np.empty(shape=(1,))
                        west_WS_8_12[Date][:] = np.nan

                        precip_west_8_12[Date][ens_memb,] = np.nan     
                else:
                        west_WS_8_12[Date] = IDX2[0, np.logical_and(  ## west 8-12
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] >=8,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] < 12)[0,:]]

                        precip_west_8_12[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][west_WS_8_12[Date],ens_memb])

                if len(meps_east_idx[Date][:,ens_memb]) == 0.:   ## east 8-12
                        east_WS_8_12[Date] = np.empty(shape=(1,))
                        east_WS_8_12[Date][:] = np.nan

                        precip_east_8_12[Date][ens_memb,] = np.nan     
                else:
                        east_WS_8_12[Date] = IDX3[0, np.logical_and(  ## east 8-12
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] >=8,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] < 12)[0,:]]

                        precip_east_8_12[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][east_WS_8_12[Date],ens_memb]) 

                ######################
                ## 12 - 16 m/s
                if len(meps_west_idx[Date][:,ens_memb]) == 0.:   ## west 12-16
                        west_WS_12_16[Date] = np.empty(shape=(1,))
                        west_WS_12_16[Date][:] = np.nan

                        precip_west_12_16[Date][ens_memb,] = np.nan     
                else:
                        west_WS_12_16[Date] = IDX2[0, np.logical_and(  ## west 12-16
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] >=12,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] < 16)[0,:]]

                        precip_west_12_16[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][west_WS_12_16[Date],ens_memb])

                if len(meps_east_idx[Date][:,ens_memb]) == 0.:   ## east 12-16
                        east_WS_12_16[Date] = np.empty(shape=(1,))
                        east_WS_12_16[Date][:] = np.nan

                        precip_east_12_16[Date][ens_memb,] = np.nan     
                else:
                        east_WS_12_16[Date] = IDX3[0, np.logical_and(  ## east 12-16
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] >=12,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] < 16)[0,:]]

                        precip_east_12_16[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][east_WS_12_16[Date],ens_memb])

                ######################
                ## 16 - 20 m/s
                if len(meps_west_idx[Date][:,ens_memb]) == 0.:   ## west 16-20
                        west_WS_16_20[Date] = np.empty(shape=(1,))
                        west_WS_16_20[Date][:] = np.nan

                        precip_west_16_20[Date][ens_memb,] = np.nan     
                else:
                        west_WS_16_20[Date] = IDX2[0, np.logical_and(  ## west 16-20
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] >=16,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] < 20)[0,:]]

                        precip_west_16_20[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][west_WS_16_20[Date],ens_memb])

                if len(meps_east_idx[Date][:,ens_memb]) == 0.:   ## east 16-20
                        east_WS_16_20[Date] = np.empty(shape=(1,))
                        east_WS_16_20[Date][:] = np.nan

                        precip_east_16_20[Date][ens_memb,] = np.nan     
                else:
                        east_WS_16_20[Date] = IDX3[0, np.logical_and(  ## east 16-20
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] >=16,
                            meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] < 20)[0,:]]

                        precip_east_16_20[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][east_WS_16_20[Date],ens_memb])

                ######################
                ## 20 - inf m/s
                if len(meps_west_idx[Date][:,ens_memb]) == 0.:   ## west > 20
                        west_WS_20[Date] = np.empty(shape=(1,))
                        west_WS_20[Date][:] = np.nan

                        precip_west_20[Date][ens_memb,] = np.nan     
                else:
                        west_WS_20[Date] = IDX2[0,              ## west > 20
                                                (meps_wind_speed[Date][np.where(~np.isnan(meps_west_idx[Date][:,ens_memb])),ens_memb] >=20)[0,:]]

                        precip_west_20[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][west_WS_20[Date],ens_memb])

                if len(meps_east_idx[Date][:,ens_memb]) == 0.:   ## east > 20
                        east_WS_20[Date] = np.empty(shape=(1,))
                        east_WS_20[Date][:] = np.nan

                        precip_east_20[Date][ens_memb,] = np.nan     
                else:
                        east_WS_20[Date] = IDX3[0,              ## east > 20
                                                (meps_wind_speed[Date][np.where(~np.isnan(meps_east_idx[Date][:,ens_memb])),ens_memb] >=20)[0,:]]

                        precip_east_20[Date][ens_memb,]  = np.nansum(meps_precip_diff[Date][east_WS_20[Date],ens_memb])
    


    vals_west_0_4 = get_array_of_dictionary_MEPS(precip_west_0_4)
    vals_east_0_4 = get_array_of_dictionary_MEPS(precip_east_0_4)

    vals_west_4_8 = get_array_of_dictionary_MEPS(precip_west_4_8)
    vals_east_4_8 = get_array_of_dictionary_MEPS(precip_east_4_8)

    vals_west_8_12 = get_array_of_dictionary_MEPS(precip_west_8_12)
    vals_east_8_12 = get_array_of_dictionary_MEPS(precip_east_8_12)

    vals_west_12_16 = get_array_of_dictionary_MEPS(precip_west_12_16)
    vals_east_12_16 = get_array_of_dictionary_MEPS(precip_east_12_16)

    vals_west_16_20 = get_array_of_dictionary_MEPS(precip_west_16_20)
    vals_east_16_20 = get_array_of_dictionary_MEPS(precip_east_16_20)

    vals_west_20 = get_array_of_dictionary_MEPS(precip_west_20)
    vals_east_20 = get_array_of_dictionary_MEPS(precip_east_20)

    return vals_west_0_4, vals_east_0_4, vals_west_4_8, vals_east_4_8 , vals_west_8_12, vals_east_8_12, vals_west_12_16, vals_east_12_16, vals_west_16_20, vals_east_16_20,vals_west_20, vals_east_20


# In[ ]:




