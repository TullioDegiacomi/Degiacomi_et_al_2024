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


#%% plot RH profiles:
RH_CTRL = pd.read_csv("Parameters_CTRL_sounding.txt",delimiter= '\t')
CTRL = pd.read_csv("thetae_ctrl.txt",delimiter = '\t')
RH_INCREASED5 = pd.read_csv("Parameters_Sondaggio_RH_increased5.txt",delimiter = '\t')
RH_REDUCED5 = pd.read_csv("Parameters_Sondaggio_RH_reduced5.txt",delimiter = '\t')
RH_INCREASED5_under2600 = pd.read_csv("Parameters_Sondaggio_RH_increased5_under2600m.txt",delimiter = '\t')
RH_INCREASED5_over2300 = pd.read_csv("Parameters_Sondaggio_RH_increased5_over2300m.txt",delimiter = '\t')

RH_CTRL = RH_CTRL['RH']
z_CTRL = CTRL['Height']
RH_increased5 = RH_INCREASED5['RH']
z_increased5 = RH_INCREASED5['Height']
RH_reduced5 = RH_REDUCED5['RH']
z_reduced5 = RH_REDUCED5['Height']
RH_increased5_under2600 = RH_INCREASED5_under2600['RH']
RH_increased5_over2300 = RH_INCREASED5_over2300['RH']


fig = plt.figure(figsize=(18,8))
plt.rcParams.update({'font.size': 24})

ax = plt.subplot(1, 2, 1)
ax.plot(RH_CTRL,z_CTRL,linewidth=4, label = 'CTRL',color='green',alpha=0.7,linestyle='dashed')
ax.plot(RH_increased5,z_increased5,linewidth=4,label='RH$\_$INCR5',color='blue',alpha=0.5)
ax.plot(RH_reduced5,z_reduced5,linewidth=4,label = 'RH$\_$RED5',color='red',alpha=0.7)
plt.legend(fontsize=19,loc=3)
plt.xlabel('$\it{RH}$ [$\%$]',fontsize=28)
plt.ylabel('Height [m]', fontsize=28)
plt.ylim((50,6000))
plt.xlim((40,100))
plt.text(42., 5600, '(a)', fontsize=28)
plt.grid()

ax = plt.subplot(1, 2, 2)
ax.plot(RH_increased5,z_increased5,linewidth=5, label = 'RH$\_$INCR5',color='black',alpha=1.,linestyle='dashed')
ax.plot(RH_increased5_under2600,z_increased5,linewidth=5,label='RH$\_$INCR5$\_$LL',color='blue',alpha=0.5)
ax.plot(RH_increased5_over2300,z_reduced5,linewidth=5,label = 'RH$\_$INCR5$\_$UL',color='red',alpha=0.5)
plt.legend(fontsize=19,loc=3)
plt.xlabel('$\it{RH}$ [$\%$]',fontsize=28)
plt.ylabel('Height [m]',fontsize=28)
plt.ylim((50,6000))
plt.xlim((40,100))
plt.text(42., 5600, '(b)', fontsize=28)
plt.grid()

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.28, 
                    hspace=0.6)

plt.savefig('figure19.png', dpi=600)

