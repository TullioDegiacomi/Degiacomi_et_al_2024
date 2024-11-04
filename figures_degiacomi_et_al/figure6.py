# -*- coding: utf-8 -*-
"""
Created on Tue Jan 25 10:14:56 2022

@author: Admin
"""


from netCDF4 import Dataset, MFDataset
import xarray as xr
from copy import copy
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.cm import get_cmap,colors
from matplotlib.colors import TwoSlopeNorm
import os
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, get_basemap, ALL_TIMES,CoordPair,vertcross, interpline)

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

#%% plot control simulation

path_data= os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_IDEALIZED_SOUNDING/CTRL_SIMULATION/")
SIM_CTRL = Dataset(path_data + '/CTRL_SIM_domain2')

# definition of the domain for domain 2 and 3 in order to have KM on axis
i_start_2 = 163
#i_start_3 = 35
#i_start_3 = 41
#j_start_3 = 51
#start2 = i_start_2*3   # multiplying for resolution of domain 1
start2 = (i_start_2-1)*3 #+0.5
cells2 = 225
cells3 = 390
#cells3_x = 975
#cells3_y = 975

X2 = np.arange(start2,start2+cells2) 
Y2 = X2

#setting domain 3 for 500 m resolution
start3 = (start2-0.5) + 16 + 0.25
X3 = np.arange(start3,start3+cells3/2,0.5)
Y3 = X3

#setting domain 3 for 200 m of resolution
#start3_x = (start2-0.5) + 16
#start3_y = (start2-0.5) + 16
#X3 = np.arange(start3_x,start3_x + cells3_x/5,0.2) 
#Y3 = np.arange(start3_y,start3_y + cells3_y/5,0.2)

#%%


# plot of QRain hour 6 and 9
z = getvar(SIM_CTRL,'z')
zagl = getvar(SIM_CTRL,'height_agl')
HGT = getvar(SIM_CTRL,'HGT')
QRAIN = getvar(SIM_CTRL,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000 = interplevel(QRAIN,z,2000)
RAINNC = getvar(SIM_CTRL,'RAINNC',timeidx = ALL_TIMES)

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_2000 = QRAIN_2000.assign_coords({'west_east' : QRAIN_2000.west_east.values +X2[0]+1})
QRAIN_2000 = QRAIN_2000.assign_coords({'south_north' : QRAIN_2000.south_north.values +X2[0]+1})
RAINNC = RAINNC.assign_coords({'west_east' : RAINNC.west_east.values +X2[0]+1})
RAINNC = RAINNC.assign_coords({'south_north' : RAINNC.south_north.values +X2[0]+1})


fig = plt.figure(figsize=(15,13))
plt.rcParams.update({'font.size': 22})

ax = plt.subplot(2, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, '(a)', fontsize=26)
plt.hlines(596,525,675, linestyle='dashed',color = 'black')
plt.vlines(606,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]',fontsize = 26)
plt.ylabel('$\it{y}$ [km]',fontsize = 26)

ax = plt.subplot(2, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000[9].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000[9].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, '(b)', fontsize=26)

#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]',fontsize = 26)
plt.ylabel('',fontsize = 26)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.55, 0.67, 0.015])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg kg$^{-1}$', labelpad=5, y=0.95, rotation=0,fontsize=26)

ax = plt.subplot(2, 2, 3)
plot2=(RAINNC[6]-RAINNC[5]).plot.contourf(ax=ax,levels=np.linspace(5,30,11),cmap='BuPu',add_colorbar=False)
#cb = plt.colorbar(plot)
#(RAINNC[6]-RAINNC[5]).plot.contour(ax=ax,levels=10,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, '(c)', fontsize=26)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]',fontsize = 26)
plt.ylabel('$\it{y}$ [km]',fontsize = 26)

ax = plt.subplot(2, 2, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
(RAINNC[9]-RAINNC[8]).plot.contourf(ax=ax,levels=np.linspace(5,30,11),cmap='BuPu',add_colorbar=False)
#cb = plt.colorbar(plot)
#(RAINNC[9]-RAINNC[8]).plot.contour(ax=ax,levels=10,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, '(d)', fontsize=26)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]',fontsize = 26)
plt.ylabel('',fontsize = 26)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.07, 0.67, 0.015])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
cb.set_label('mm h$^{-1}$', labelpad=0, y=0.95, rotation=0,fontsize=26)

plt.subplots_adjust(left=0.1,
                    bottom=0.15, 
                    right=0.87, 
                    top=0.95, 
                    wspace=0.25, 
                    hspace=0.55)


plt.savefig('figure6.png',dpi=300)

