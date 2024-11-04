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

# definition of the domain for domain 2 in order to have KM on axis
i_start_2 = 163
start2 = (i_start_2-1)*3
cells2 = 225


X2 = np.arange(start2,start2+cells2) 
Y2 = X2

z = getvar(SIM_CTRL,'z')
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
plt.rcParams.update({'font.size': 40})

ax = plt.subplot(1, 1, 1)
plot=(RAINNC[12]-RAINNC[6]).plot.contourf(ax=ax,levels=np.linspace(0,140,8),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]', fontsize = 44)
plt.ylabel('$\it{y}$ [km]', fontsize = 44)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.16, 0.08, 0.71, 0.015])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('mm', labelpad=0, y=0.95, rotation=0,fontsize=44)

plt.subplots_adjust(left=0.12,
                    bottom=0.2, 
                    right=0.91, 
                    top=0.93)

plt.savefig('figure7.png',dpi=300)

