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
cells3 = 390

X2 = np.arange(start2,start2+cells2) 
Y2 = X2

#%% PLOT SIM WITH IDEAL WINDS
path_data= os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_IDEALIZED_SOUNDING/SIM_IDEAL_WIND/")
SIM_U10 = Dataset(path_data + '/SIM_U10_dominio2')
SIM_U20 = Dataset(path_data + '/SIM_U20_dominio2')
SIM_U30 = Dataset(path_data + '/SIM_U30_dominio2')
SIM_shear = Dataset(path_data + '/SIM_SHEAR_FINALE_dominio2')

z = getvar(SIM_U10,'z')
HGT = getvar(SIM_U10,'HGT')
QRAIN = getvar(SIM_U10,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_U10 = interplevel(QRAIN,z,2000)

z = getvar(SIM_U20,'z')
HGT = getvar(SIM_U20,'HGT')
QRAIN = getvar(SIM_U20,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_U20 = interplevel(QRAIN,z,2000)

z = getvar(SIM_U30,'z')
HGT = getvar(SIM_U30,'HGT')
QRAIN = getvar(SIM_U30,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_U30 = interplevel(QRAIN,z,2000)

z = getvar(SIM_shear,'z')
HGT = getvar(SIM_shear,'HGT')
QRAIN = getvar(SIM_shear,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_shear = interplevel(QRAIN,z,2000)


# making the plot
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_2000_U10 = QRAIN_2000_U10.assign_coords({'west_east' : QRAIN_2000_U10.west_east.values +X2[0]+1})
QRAIN_2000_U10 = QRAIN_2000_U10.assign_coords({'south_north' : QRAIN_2000_U10.south_north.values +X2[0]+1})
QRAIN_2000_U20 = QRAIN_2000_U20.assign_coords({'west_east' : QRAIN_2000_U20.west_east.values +X2[0]+1})
QRAIN_2000_U20 = QRAIN_2000_U20.assign_coords({'south_north' : QRAIN_2000_U20.south_north.values +X2[0]+1})
QRAIN_2000_U30 = QRAIN_2000_U30.assign_coords({'west_east' : QRAIN_2000_U30.west_east.values +X2[0]+1})
QRAIN_2000_U30 = QRAIN_2000_U30.assign_coords({'south_north' : QRAIN_2000_U30.south_north.values +X2[0]+1})
QRAIN_2000_shear = QRAIN_2000_shear.assign_coords({'west_east' : QRAIN_2000_shear.west_east.values +X2[0]+1})
QRAIN_2000_shear = QRAIN_2000_shear.assign_coords({'south_north' : QRAIN_2000_shear.south_north.values +X2[0]+1})

#%%
fig = plt.figure(figsize=(10,15))
plt.rcParams.update({'font.size': 14})

ax = plt.subplot(4, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_U10[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_U10[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.text(529., 648, '(a)', fontsize=16)
props = dict(facecolor='gray', alpha=0.5)
plt.text(528, 548, 'V10', fontsize=16, verticalalignment='center',
          bbox=props)
plt.vlines(585,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
plt.title('6 h', y =1.02,fontsize=18)

ax = plt.subplot(4, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_U10[9].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_U10[9].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, '(b)', fontsize=16)
#plt.vlines(589,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('', fontsize=18)
plt.title('9 h', y =1.02,fontsize=18)

ax = plt.subplot(4, 2, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_U20[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_U20[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(c)',  fontsize=16)
plt.vlines(589,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
props = dict(facecolor='gray', alpha=0.5)
plt.text(528, 548, 'V20', fontsize=16, verticalalignment='center',
          bbox=props)


ax = plt.subplot(4, 2, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_U20[9].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_U20[9].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(d)',  fontsize=16)
#plt.vlines(590,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('', fontsize=18)

ax = plt.subplot(4, 2, 5)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_U30[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_U30[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(e)',  fontsize=16)
plt.vlines(587,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
props = dict(facecolor='gray', alpha=0.5)
plt.text(528, 548, 'V30', fontsize=16, verticalalignment='center',
          bbox=props)


ax = plt.subplot(4, 2, 6)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_U30[9].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_U30[9].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(f)',  fontsize=16)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('', fontsize=18)
plt.ylabel('', fontsize=18)

ax = plt.subplot(4, 2, 7)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_shear[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_shear[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(g)', fontsize=16)
plt.vlines(615,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
props = dict(facecolor='gray', alpha=0.5)
plt.text(528, 548, 'V10_SHEAR', fontsize=16, verticalalignment='center',
          bbox=props)



ax = plt.subplot(4, 2, 8)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_shear[9].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_shear[9].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('', y =1.03,fontsize=16)
plt.text(529., 648, '(h)', fontsize=16)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]', fontsize=18)
plt.ylabel('', fontsize=18)


from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.22, 0.14, 0.63, 0.01])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg kg$^{-1}$',  fontsize=18, labelpad=0, y=0.92, rotation=0)


plt.subplots_adjust(left=0.2,
                    bottom=0.2,
                    right=0.87,
                    top=0.92,
                    wspace=0.2,
                    hspace=0.15)


plt.savefig('figure9.png',dpi=300)
