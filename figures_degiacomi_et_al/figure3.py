# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 11:24:11 2021

@author: Tullio Degiacomi
"""

from netCDF4 import Dataset
import numpy as np

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import get_cmap,colors
import os
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, get_basemap, ALL_TIMES)


path_data = os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_500m/")

ncfile1 = Dataset(path_data + '/SIM_500m_domain1')
ncfile2 = Dataset(path_data + '/SIM_500m_domain2')
ncfile3 = Dataset(path_data + '/SIM_500m_domain3')

HGT = getvar(ncfile1,"HGT")
eta = getvar(ncfile1,"z")

import numpy as np

# setting X and Y for domain 1 with resolution of 3 km
cells1 = 400
X1 = np.arange(1.5,1201.5,3) 
Y1=X1


# definition of the domain for domain 2 and 3 in order to have KM on axis
i_start_2 = 163
i_start_3 = 41
j_start_3 = 51
start2 = (i_start_2-1)*3 
cells2 = 225
cells3 = 390
cells3_x = 975
cells3_y = 975

X2 = np.arange(start2,start2+cells2) 
Y2 = X2

#setting domain 3 for 500 m resolution
start3 = start2 + 16 
X3 = np.arange(start3,start3+cells3/2,0.5)
Y3 = X3

#setting domain 3 for 200 m of resolution
#start3_x = (start2) + 16
#start3_y = (start2) + 16
#X3 = np.arange(start3_x,start3_x + cells3_x/5,0.2) 
#Y3 = np.arange(start3_y,start3_y + cells3_y/5,0.2)

# code for cutting the colormp
import matplotlib.pyplot as plt
import matplotlib.colors as colors

################### Function to truncate color map ###################
def truncate_colormap(cmapIn='terrain', minval=0.0, maxval=1.0, n=100):
    '''truncate_colormap(cmapIn='terrain', minval=0.0, maxval=1.0, n=100)'''    
    cmapIn = plt.get_cmap(cmapIn)

    new_cmap = colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmapIn.name, a=minval, b=maxval),
        cmapIn(np.linspace(minval, maxval, n)))

    arr = np.linspace(0, 50, 100).reshape((10, 10))
    fig, ax = plt.subplots(ncols=2)
    ax[0].imshow(arr, interpolation='nearest', cmap=cmapIn)
    ax[1].imshow(arr, interpolation='nearest', cmap=new_cmap)

    return new_cmap

cmap_mod = truncate_colormap(minval=.4, maxval=.9)  # calls function to truncate colormap

#%%
fig = plt.figure(figsize=(18,18))
ax = plt.axes()

xlim_d01=ylim_d01 =np.zeros(2)
xlim_d01[0] = X1[0]
xlim_d01[1] = X1[-1]
ylim_d01[0] = Y1[0]
ylim_d01[1] = Y1[-1]
# defining domain 2 limits 
xlim_d02=ylim_d02 =np.zeros(2)
xlim_d02[0] = X2[0]
xlim_d02[1] = X2[-1]
ylim_d02[0] = Y2[0]
ylim_d02[1] = Y2[-1]
# defining domain 3 limits 
xlim_d03=np.zeros(2)
ylim_d03 =np.zeros(2)
xlim_d03[0] = X3[0]
xlim_d03[1] = X3[-1]
ylim_d03[0] = Y3[0]
ylim_d03[1] = Y3[-1]


 
# d01 box
ax.add_patch(mpl.patches.Rectangle((xlim_d01[0], ylim_d01[0]), xlim_d01[1]-xlim_d01[0], ylim_d01[1]-ylim_d01[0],
             fill=None, lw=4.5, edgecolor='black', zorder=10))
ax.text(xlim_d01[0]+(xlim_d01[1]-xlim_d01[0])*0.04, ylim_d01[0]+(ylim_d01[1]-ylim_d01[0])*0.94, 'D01, $\Delta$=3 km',
        size=44, color='black', zorder=10)

# d02 box
ax.add_patch(mpl.patches.Rectangle((xlim_d02[0], ylim_d02[0]), xlim_d02[1]-xlim_d02[0], ylim_d02[1]-ylim_d02[0],
             fill=None, lw=2.5, edgecolor='black', zorder=10))

ax.text(xlim_d02[0]-(xlim_d02[1]-xlim_d02[0])*0.2, ylim_d02[0]+(ylim_d02[1]-ylim_d02[0])*1.05, 'D02, $\Delta$=1 km',
        size=44, color='black', zorder=10)

# d03 box
ax.add_patch(mpl.patches.Rectangle((xlim_d03[0], ylim_d03[0]), xlim_d03[1]-xlim_d03[0], ylim_d03[1]-ylim_d03[0],
              fill=None, lw=2.5, edgecolor='black', zorder=10))
ax.text(xlim_d03[0]+(xlim_d03[1]-xlim_d03[0])*0.1, ylim_d03[0]+(ylim_d03[1]-ylim_d03[0])*0.78, 'D03 ',
        size=44, color='black', zorder=10)

plt.rcParams["figure.figsize"] = (18,14) 
levels = [100,300,500,1000,1250,1450]
levels_color =[0,100,300,500,1000,1250,1450,1500]
plt.contour(X1,Y1,HGT,levels=levels,colors=['black'],linewidths=0.7)
plt.contourf(X1,Y1,HGT,cmap= cmap_mod,levels=levels_color)
ax.set_xlabel('$\it{x}$ [km]',fontsize = 44)
ax.set_ylabel('$\it{y}$ [km]',fontsize = 44)
#ax.set_title('Nested domains',fontsize = 20)
plt.tick_params(length=16,width=2)
plt.yticks(fontsize=44)
plt.xticks(fontsize=44)

plt.savefig('figure3a',dpi=300)

#%%
fig = plt.figure(figsize=(14,10))
ax = plt.axes()
HGT = getvar(ncfile2,"HGT")

xlim_d01=ylim_d01 =np.zeros(2)
xlim_d01[0] = X1[0]
xlim_d01[1] = X1[-1]
ylim_d01[0] = Y1[0]
ylim_d01[1] = Y1[-1]
# defining domain 2 limits 
xlim_d02=ylim_d02 =np.zeros(2)
xlim_d02[0] = X2[0]
xlim_d02[1] = X2[-1]
ylim_d02[0] = Y2[0]
ylim_d02[1] = Y2[-1]
# defining domain 3 limits 
xlim_d03=np.zeros(2)
ylim_d03 =np.zeros(2)
xlim_d03[0] = X3[0]
xlim_d03[1] = X3[-1]
ylim_d03[0] = Y3[0]
ylim_d03[1] = Y3[-1]



# d02 box
ax.add_patch(mpl.patches.Rectangle((xlim_d02[0], ylim_d02[0]), xlim_d02[1]-xlim_d02[0], ylim_d02[1]-ylim_d02[0],
             fill=None, lw=4.5, edgecolor='black', zorder=10))

ax.text(xlim_d02[0]+(xlim_d02[1]-xlim_d02[0])*0.05, ylim_d02[0]+(ylim_d02[1]-ylim_d02[0])*0.952, 'D02, $\Delta$ = 1000 m',
        size=32, color='black', zorder=10)

#d03 box
ax.add_patch(mpl.patches.Rectangle((xlim_d03[0], ylim_d03[0]), xlim_d03[1]-xlim_d03[0], ylim_d03[1]-ylim_d03[0],
              fill=None, lw=2.5, edgecolor='black', zorder=10))
ax.text(xlim_d03[0]+(xlim_d03[1]-xlim_d03[0])*0.1, ylim_d03[0]+(ylim_d03[1]-ylim_d03[0])*0.94, 'D03, $\Delta$ = 500 m or 200 m',
        size=32, color='black', zorder=10)

plt.rcParams["figure.figsize"] = (14,10)
levels = [100,300,500,1000,1250,1450]
levels_color =np.array([0,100,300,500,1000,1250,1450,1500])
CS=plt.contour(X2,Y2,HGT,levels=levels,colors=['black'],linewidths=0.7)
plt.contourf(X2,Y2,HGT,cmap= cmap_mod,levels=levels_color)
ax.clabel(CS, inline=1, fontsize=24,fmt='%1d')

ax.set_xlabel('$\it{x}$ [km]',fontsize = 32)
ax.set_ylabel('$\it{y}$ [km]',fontsize = 32)
plt.tick_params(length=16,width=2)
plt.yticks(fontsize=32)
plt.xticks(fontsize=32)
cb = plt.colorbar(orientation = 'vertical')
cb.ax.tick_params(labelsize=30)
cb.set_label('m MSL', labelpad=-80, y=1.07, rotation=0,fontsize=32)
#plt.title('TERRAIN HEIGHT',y=1.02,fontsize=20)


plt.savefig('figure3b.png',dpi=300)
