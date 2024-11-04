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

# profiles stability
N1_000001 = pd.read_csv("thetae_N1_000001.txt",delimiter = '\t')
N1_00004 = pd.read_csv("thetae_N1_00004.txt",delimiter = '\t')
N1_00015 = pd.read_csv("thetae_N1_00015.txt",delimiter = '\t')
CTRL = pd.read_csv("thetae_ctrl.txt",delimiter = '\t')
N3_00004 = pd.read_csv("thetae_N3_00004.txt",delimiter='\t')
N3_00012 = pd.read_csv("thetae_N3_00012.txt",delimiter='\t')
N3_00009 = pd.read_csv("thetae_N3_00009.txt",delimiter='\t')

z = N1_000001['Height']
thetae_N1_00004 = N1_00004['theta_e']
thetae_CTRL = CTRL['theta_e']
z_CTRL = CTRL['Height']
thetae_N1_00015 = N1_00015['theta_e']
thetae_N3_00004 = N3_00004['theta_e']
thetae_N1_000001 = N1_000001['theta_e']
thetae_N1_000001 = N1_000001['theta_e']
thetae_N3_00012 = N3_00012['theta_e']
thetae_N3_00009 = N3_00009['theta_e']

fig = plt.figure(figsize=(20,9))
plt.rcParams.update({'font.size': 28})

ax = plt.subplot(1, 2, 1)
ax.plot(thetae_CTRL,z_CTRL,linewidth=4, label = 'CTRL',color='green',alpha=0.7,linestyle='dashed')
ax.plot(thetae_N1_000001,z,linewidth=3,label='N1_000001',color='blue',alpha=0.6)
ax.plot(thetae_N1_00004,z,linewidth=3,label='N1_00004',alpha=0.6)
ax.plot(thetae_N1_00015,z,linewidth=3,label = 'N1_00015',color='red',alpha=0.6)

plt.legend(fontsize=20)
plt.xlabel(r'$\theta_e$ [K]',fontsize=34)
plt.ylabel('Height [m]',fontsize=34)
plt.ylim((50,6000))
plt.xlim((305,355))
plt.grid()
plt.text(306., 5600, '(a)', fontsize=34)

ax = plt.subplot(1,2,2)
ax.plot(thetae_CTRL,z_CTRL,linewidth=4, label = 'CTRL',color='green',alpha=0.7,linestyle='dashed')
ax.plot(thetae_N1_00004,z,linewidth=3,label='N1_00004',alpha=0.6)
ax.plot(thetae_N3_00004,z,linewidth=3,label = 'N1_00004_N3_00004',color='blue',alpha=0.6)
ax.plot(thetae_N3_00009,z,linewidth=3,label = 'N1_00004_N3_00009',color ='orange',alpha=0.6)
ax.plot(thetae_N3_00012,z,linewidth=3,label = 'N1_00004_N3_00012',color ='red',alpha=0.6)
plt.legend(fontsize=20,loc=1,bbox_to_anchor=(0.999,0.87))
plt.ylim((50,6000))
plt.xlim((305,345))
plt.xlabel(r'$\theta_e$ [K]',fontsize=34)
plt.ylabel('Height [m]',fontsize=34)
plt.grid()
plt.text(306., 5600, '(b)',  fontsize=34)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.28, 
                    hspace=0.6)

plt.savefig('figure15.png',dpi=600)
