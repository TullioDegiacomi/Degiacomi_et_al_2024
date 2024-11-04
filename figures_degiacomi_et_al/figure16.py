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

#%% plot stability layer 1
path_data = os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_IDEALIZED_SOUNDING/SIM_IDEAL_N/")
#SIM_CTRL = Dataset(path_data + '/SIM_N1_000001_dominio2')
SIM_N1_000001 = Dataset(path_data + '/SIM_N1_000001_dominio2')
SIM_N1_00004 = Dataset(path_data + '/SIM_N1_00004_dominio2')
SIM_N1_00015 = Dataset(path_data + '/SIM_N1_00015_dominio2')

z = getvar(SIM_N1_000001,'z')      

HGT = getvar(SIM_N1_000001,'HGT')
QRAIN_N1_000001 = getvar(SIM_N1_000001,'QRAIN',timeidx=ALL_TIMES)
QRAIN_N1_00004 = getvar(SIM_N1_00004,'QRAIN',timeidx=ALL_TIMES)
QRAIN_N1_00015 = getvar(SIM_N1_00015,'QRAIN',timeidx=ALL_TIMES)

QRAIN_N1_000001_2000 = interplevel(QRAIN_N1_000001,z,2000)
QRAIN_N1_00004_2000 = interplevel(QRAIN_N1_00004,z,2000)
QRAIN_N1_00015_2000 = interplevel(QRAIN_N1_00015,z,2000)


HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
#QRAIN_CTRL_2000 = QRAIN_CTRL_2000.assign_coords({'west_east' : QRAIN_CTRL_2000.west_east.values +X2[0]+1})
#QRAIN_CTRL_2000 = QRAIN_CTRL_2000.assign_coords({'south_north' : QRAIN_CTRL_2000.south_north.values +X2[0]+1})
QRAIN_N1_000001_2000 = QRAIN_N1_000001_2000.assign_coords({'west_east' : QRAIN_N1_000001_2000.west_east.values +X2[0]+1})
QRAIN_N1_000001_2000 = QRAIN_N1_000001_2000.assign_coords({'south_north' : QRAIN_N1_000001_2000.south_north.values +X2[0]+1})
QRAIN_N1_00004_2000 = QRAIN_N1_00004_2000.assign_coords({'west_east' : QRAIN_N1_00004_2000.west_east.values +X2[0]+1})
QRAIN_N1_00004_2000 = QRAIN_N1_00004_2000.assign_coords({'south_north' : QRAIN_N1_00004_2000.south_north.values +X2[0]+1})
QRAIN_N1_00015_2000 = QRAIN_N1_00015_2000.assign_coords({'west_east' : QRAIN_N1_00015_2000.west_east.values +X2[0]+1})
QRAIN_N1_00015_2000 = QRAIN_N1_00015_2000.assign_coords({'south_north' : QRAIN_N1_00015_2000.south_north.values +X2[0]+1})


#%%
fig = plt.figure(figsize=(10,15))
plt.rcParams.update({'font.size': 14})

ax = plt.subplot(3, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N1_000001_2000[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_N1_000001_2000[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.text(529., 648, '(a)', fontsize=16)
props = dict(facecolor='gray', alpha=0.5)
plt.text(528, 547.5, 'N1_000001', fontsize=16, verticalalignment='center',
          bbox=props)
#plt.vlines(585,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
plt.title('6 h', y =1.02,fontsize=18)


ax = plt.subplot(3, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N1_000001_2000[9].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_N1_000001_2000[9].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, '(b)', fontsize=16)
#plt.vlines(589,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('', fontsize=18)
plt.title('9 h', y =1.02,fontsize=18)

ax = plt.subplot(3, 2, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N1_00004_2000[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_N1_00004_2000[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(c)',  fontsize=16)
#plt.vlines(589,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
props = dict(facecolor='gray', alpha=0.5)
plt.text(528, 547.5, 'N1_00004', fontsize=16, verticalalignment='center',
          bbox=props)


ax = plt.subplot(3, 2, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N1_00004_2000[9].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_N1_00004_2000[9].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(d)',  fontsize=16)
#plt.vlines(590,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('', fontsize=18)

ax = plt.subplot(3, 2, 5)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N1_00015_2000[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_N1_00015_2000[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(e)',  fontsize=16)
#plt.vlines(587,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
props = dict(facecolor='gray', alpha=0.5)
plt.text(528, 547.5, 'N1_00015', fontsize=16, verticalalignment='center',
          bbox=props)


ax = plt.subplot(3, 2, 6)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N1_00015_2000[9].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_N1_00015_2000[9].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(f)',  fontsize=16)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]', fontsize=18)
plt.ylabel('', fontsize=18)


from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.22, 0.14, 0.66, 0.01])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg kg$^{-1}$',  fontsize=18, labelpad=0, y=0.92, rotation=0)


plt.subplots_adjust(left=0.2,
                    bottom=0.2,
                    right=0.9,
                    top=0.76,
                    wspace=0.2,
                    hspace=0.15)


plt.savefig('figure16.png',dpi=300)
