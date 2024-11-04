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

# setting X and Y for domain 1 with resolution of 3 km
cells1 = 400
X1 = np.arange(1.5,1201.5,3) 
Y1=X1

# definition of the domain for domain 2 in order to have KM on axis
i_start_2 = 163
start2 = (i_start_2-1)*3
cells2 = 225

X2 = np.arange(start2,start2+cells2) 
Y2 = X2

path_data = os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_IDEALIZED_SOUNDING/SIM_IDEAL_WIND_DIRECTION/")
#SIM_CTRL = Dataset(path_data + '/CTRL_SIM_domain2')
SIM_210 = Dataset(path_data + '/SIM_wind_210degrees_dominio2')
SIM_UDINE_20 = Dataset(path_data + '/SIM_original_wind20degrees_dominio2')
SIM_SHEAR_TILTED = Dataset(path_data + '/SIM_shear_tilted_prova2_dominio2')

z = getvar(SIM_210,'z')      

ua = getvar(SIM_UDINE_20,'ua',timeidx=ALL_TIMES)
va = getvar(SIM_UDINE_20,'va',timeidx=ALL_TIMES)
v_2000 = interplevel(va, z, 2000)
u_2000= interplevel(ua, z, 2000)
v_4000 = interplevel(va, z, 4000)
u_4000= interplevel(ua, z, 4000)
v_5500 = interplevel(va, z, 5500)
u_5500= interplevel(ua, z, 5500)

HGT = getvar(SIM_210,'HGT')
QRAIN_210 = getvar(SIM_210,'QRAIN',timeidx=ALL_TIMES)
QRAIN_UDINE_20 = getvar(SIM_UDINE_20,'QRAIN',timeidx=ALL_TIMES)
QRAIN_SHEAR_TILTED = getvar(SIM_SHEAR_TILTED,'QRAIN',timeidx=ALL_TIMES)

QRAIN_210_2000 = interplevel(QRAIN_210,z,2000)
QRAIN_UDINE_20_2000 = interplevel(QRAIN_UDINE_20,z,2000)
QRAIN_SHEAR_TILTED_2000 = interplevel(QRAIN_SHEAR_TILTED,z,2000)

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
#QRAIN_CTRL_2000 = QRAIN_CTRL_2000.assign_coords({'west_east' : QRAIN_CTRL_2000.west_east.values +X2[0]+1})
#QRAIN_CTRL_2000 = QRAIN_CTRL_2000.assign_coords({'south_north' : QRAIN_CTRL_2000.south_north.values +X2[0]+1})
QRAIN_210_2000 = QRAIN_210_2000.assign_coords({'west_east' : QRAIN_210_2000.west_east.values +X2[0]+1})
QRAIN_210_2000 = QRAIN_210_2000.assign_coords({'south_north' : QRAIN_210_2000.south_north.values +X2[0]+1})
QRAIN_UDINE_20_2000 = QRAIN_UDINE_20_2000.assign_coords({'west_east' : QRAIN_UDINE_20_2000.west_east.values +X2[0]+1})
QRAIN_UDINE_20_2000 = QRAIN_UDINE_20_2000.assign_coords({'south_north' : QRAIN_UDINE_20_2000.south_north.values +X2[0]+1})
QRAIN_SHEAR_TILTED_2000 = QRAIN_SHEAR_TILTED_2000.assign_coords({'west_east' : QRAIN_SHEAR_TILTED_2000.west_east.values +X2[0]+1})
QRAIN_SHEAR_TILTED_2000 = QRAIN_SHEAR_TILTED_2000.assign_coords({'south_north' : QRAIN_SHEAR_TILTED_2000.south_north.values +X2[0]+1})
u_2000 = u_2000.assign_coords({'west_east' : u_2000.west_east.values +X2[0]+1})
u_2000 = u_2000.assign_coords({'south_north' : u_2000.south_north.values +X2[0]+1})
v_2000 = v_2000.assign_coords({'west_east' : v_2000.west_east.values +X2[0]+1})
v_2000 = v_2000.assign_coords({'south_north' : v_2000.south_north.values +X2[0]+1})
u_4000 = u_4000.assign_coords({'west_east' : u_4000.west_east.values +X2[0]+1})
u_4000 = u_4000.assign_coords({'south_north' : u_4000.south_north.values +X2[0]+1})
v_4000 = v_4000.assign_coords({'west_east' : v_4000.west_east.values +X2[0]+1})
v_4000 = v_4000.assign_coords({'south_north' : v_4000.south_north.values +X2[0]+1})
u_5500 = u_5500.assign_coords({'west_east' : u_5500.west_east.values +X2[0]+1})
u_5500 = u_5500.assign_coords({'south_north' : u_5500.south_north.values +X2[0]+1})
v_5500 = v_5500.assign_coords({'west_east' : v_5500.west_east.values +X2[0]+1})
v_5500 = v_5500.assign_coords({'south_north' : v_5500.south_north.values +X2[0]+1})

#%%
fig = plt.figure(figsize=(5,15))
plt.rcParams.update({'font.size': 14})

ax = plt.subplot(3, 1, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_210_2000[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_210_2000[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.text(529., 648, '(a)', fontsize=16)
props = dict(facecolor='gray', alpha=0.5)
plt.text(653, 652.5, '210$^\circ$', fontsize=16, verticalalignment='center',
          bbox=props)
#plt.vlines(585,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
plt.title('', y =1.02,fontsize=18)
plt.text(631, 543, 'avg. wind 0-5000 m', fontsize=8, verticalalignment='center', horizontalalignment='left')
plt.text(631, 551, 'shear 0-5000 m', fontsize=8, verticalalignment='center', horizontalalignment='left')
plt.arrow(614., 551,10, 0, width = 0.7,color='black')
plt.arrow(614., 543,10, 0, width = 0.7,color='green')
plt.arrow(545., 550,5, 10, width = 0.7,color='black')
plt.arrow(548., 550,5, 10, width = 0.7,color='g')

ax = plt.subplot(3, 1, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_UDINE_20_2000[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_UDINE_20_2000[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(b)',  fontsize=16)
#plt.vlines(589,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
props = dict(facecolor='gray', alpha=0.5)
plt.text(603.5, 652.5, 'UDINE_ROT20', fontsize=16, verticalalignment='center',
          bbox=props)
ax.quiver(X2[50],Y2[60],to_np(u_2000[9,50,60]),to_np(v_2000[9,50,60]),scale_units='xy', scale=1.2,color='orange',width=0.008,label = 'wind 2000m')
ax.quiver(X2[50],Y2[60],to_np(u_4000[9,50,60]),to_np(v_4000[9,50,60]),scale_units='xy', scale=1.2,color='red',width=0.008,label = 'wind 2000m')
ax.quiver(X2[50],Y2[60],to_np(u_5500[9,50,60]),to_np(v_5500[9,50,60]),scale_units='xy', scale=1.2,color='blue',width=0.008,label = 'wind 2000m')
plt.text(645, 555, 'wind 5500 m', fontsize=8, verticalalignment='center', horizontalalignment='left')
plt.text(645, 549, 'wind 4000 m', fontsize=8, verticalalignment='center', horizontalalignment='left')
plt.text(645, 543, 'wind 2000 m', fontsize=8, verticalalignment='center', horizontalalignment='left')
plt.arrow(628., 555,10, 0, width = 0.7,color='blue')
plt.arrow(628., 549,10, 0, width = 0.7,color='red')
plt.arrow(628., 543,10, 0, width = 0.7,color='orange')

ax = plt.subplot(3, 1, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_SHEAR_TILTED_2000[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_SHEAR_TILTED_2000[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(c)',  fontsize=16)
#plt.vlines(587,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
props = dict(facecolor='gray', alpha=0.5)
plt.text(596.5, 652.5, 'SHEAR_TILTED', fontsize=16, verticalalignment='center',
          bbox=props)
plt.text(631, 543, 'avg. wind 0-5000 m', fontsize=8, verticalalignment='center', horizontalalignment='left')
plt.text(631, 551, 'shear 0-5000 m', fontsize=8, verticalalignment='center', horizontalalignment='left')
plt.arrow(612., 551,10, 0, width = 0.7,color='black')
plt.arrow(612., 543,10, 0, width = 0.7,color='green')
plt.arrow(545., 550,10, 10, width = 0.7,color='black')
plt.arrow(545., 550,1, 10, width = 0.7,color='g')

from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.22, 0.15, 0.58, 0.005])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.ax.tick_params(labelsize=12)
cb.set_label('kg kg$^{-1}$',  fontsize=18, labelpad=0, y=0.92, rotation=0)

plt.subplots_adjust(left=0.2,
                    bottom=0.2,
                    right=0.84,
                    top=0.76,
                    wspace=0.2,
                    hspace=0.15)

plt.savefig('figure13.png', dpi=300)
