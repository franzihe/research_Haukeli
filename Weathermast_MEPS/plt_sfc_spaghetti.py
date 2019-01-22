#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt

import matplotlib as mpl
mpl.style.use('ggplot')


# In[ ]:


### Define colorbar colors
champ = 255.
blue = np.array([1,74,159])/champ           # for the date
memb_col = np.array([99,99,99])/champ       # ensemble member color
vert_col = np.array([197,197,197])/champ    # vertical line for day marker
#dofe = np.array([64,180,233])/champ         # color for double fence measurement
dofe = np.array([125,98,179])/champ

fontsize = 30.
tick_fs = fontsize-2
label_fs = fontsize


# In[ ]:


def plt_variable(lead_time_sfc,wd_MEPS,WD,time_EM_mean, model_var_mean,var,xdays,title):
    
    fig = plt.figure(figsize=(20,6))
    ax = plt.axes()
# Vertical line to show end of day
    ax.axvline(0,color = vert_col, linewidth = 3)
    ax.axvline(24,color = vert_col, linewidth = 3)
    ax.axvline(48,color = vert_col, linewidth = 3)




### ensemble member
    for ens_memb in range(2,10):
        ax.plot(lead_time_sfc,wd_MEPS[:,:,ens_memb],color = memb_col,
           linestyle='-', label='_nolegend_')
    ax.plot(lead_time_sfc, wd_MEPS[:,:,1], color = memb_col,
           linestyle = '-', label = 'ensemble member')
    ax.plot(lead_time_sfc, wd_MEPS[:,:,0], 'k', linewidth = 4, label = 'deterministic')

### observation
    plt.plot(np.arange(0,np.asarray(WD).shape[0]), WD,  markersize=20,  linestyle='--', 
         label = 'observation', linewidth= 4)
### ensemble mean
    ax.plot(time_EM_mean[:], np.asarray(model_var_mean)[:], color='dodgerblue', linewidth = 3.5, 
            linestyle = '--', label = 'ensemble mean') 

### fine tuning
#    lgd = ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.37),
 #          fancybox=True, shadow=True, ncol=3, fontsize=label_fs)
  #  frame = lgd.get_frame()
   # frame.set_facecolor('white')

# xaxis
    a = lead_time_sfc[0:48]
    ax.set_xlim(-0.5,49-0.5)
    ax.set_xlabel('time', fontsize=label_fs)
    ax.set_xticks(np.arange(0,49,3))
  #  dates = pvert.dates_plt_00(hour, mm, dy, yr, ini_day)
    
    ax.set_xticklabels(xdays, rotation = 25, fontsize = tick_fs)
# title
    ax.set_title(title, fontsize=fontsize, color =blue )    
    
    if var == 'WD':
        # Horizontal line to show Wind direction
        ax.axhline(90,color=vert_col, linewidth= 3)
        ax.axhline(180,color=vert_col, linewidth= 3)
        ax.axhline(270,color=vert_col, linewidth= 3)
        ax.axhline(360,color=vert_col, linewidth= 3)
        # yaxis
        ax.set_ylabel('Wind direction', fontsize=label_fs)
        ax.set_ylim(-0.5,360)
        T = np.arange(0,361,45)
        ax.set_yticks(T)
        ax.set_yticklabels(['N', '', 'E', '','S', '', 'W', '' , 'N'],fontsize = tick_fs)
    elif var == 'WS':
        # yaxis
        ax.set_ylabel('Wind speed [m$\,$s$^{-1}$]', fontsize = label_fs)
        ax.set_ylim(0,30)
        ax.set_yticks(np.arange(0,32.5,2.5))
        ax.set_yticklabels([0, '', 5, '',10,'',15,'',20,'',25,'',30], fontsize=label_fs)
    elif var == 'T2':
        ax.axhline(0,color=vert_col, linewidth= 3)
        # yaxis
        ax.set_ylabel('Air Temperature [$^\circ$C]', fontsize=label_fs)
        ax.set_ylim(-9,6)
        T = np.arange(-9,7)
        ax.set_yticks(T)
        ax.set_yticklabels([-9, '' , '', -6, '' , '', -3, '' , '', 0, 
                            '' , '', 3, '' , '', 6], fontsize=tick_fs)
    elif var == 'SP':
        # yaxis
        ax.set_ylabel('Sea Level Pressure [hPa]', fontsize=label_fs)
        ax.set_ylim(975, 1040)
        ax.set_yticks(np.arange(975,1045,5))
        ax.set_yticklabels(['' , 980,'', '','', 1000, '','','', 1020, '','','', 1040], fontsize=tick_fs)
    elif var == 'PP':
        # yaxis
        ax.set_ylabel('Precipitation amount [mm]', fontsize=label_fs)
        ax.set_ylim(-0.5,80)
        ax.set_yticks(np.arange(0,90,5))
        ax.set_yticklabels([0, '',10,'',20,'',30,'',40,'',50,'',60,'',70,'',80,'',90],fontsize = tick_fs)


        

# tight layout
    plt.tight_layout()
    
#    return(lgd)


# In[ ]:




