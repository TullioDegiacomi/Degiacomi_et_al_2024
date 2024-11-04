# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 11:24:11 2021

@author: Tullio Degiacomi
DEFINIZIONE DEL DOMINIO:
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.cm import get_cmap,colors
import os
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, get_basemap, ALL_TIMES)

# making an automatic code for the precipitation plot normalized
import numpy as np

# setting X and Y for domain 1 with resolution of 3 km
cells1 = 400
X1 = np.arange(1.5,1201.5,3) 
Y1=X1


# definition of the domain for domain 2 and 3 in order to have KM on axis
i_start_2 = 163
#i_start_3 = 35
#i_start_3 = 41
#j_start_3 = 51
#start2 = i_start_2*3   # multiplying for resolution of domain 1
start2 = (i_start_2-1)*3 +0.5
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
#%% plot QRain sensitivity on resolution
path_data = os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_200m")
SIM_1KM = Dataset(path_data + '/SIM_200m_domain2')
SIM_200m = Dataset(path_data + '/SIM_200m_domain3')

HGT = getvar(SIM_1KM,'HGT')
z = getvar(SIM_1KM,'z')
z_200 = getvar(SIM_200m,'z')
QRAIN_1KM = getvar(SIM_1KM,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_1KM = interplevel(QRAIN_1KM,z,2000)
QRAIN_200m = getvar(SIM_200m,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_200m = interplevel(QRAIN_200m,z_200,2000)

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_2000_1KM = QRAIN_2000_1KM.assign_coords({'west_east' : QRAIN_2000_1KM.west_east.values +X2[0]+1})
QRAIN_2000_1KM = QRAIN_2000_1KM.assign_coords({'south_north' : QRAIN_2000_1KM.south_north.values +X2[0]+1})
QRAIN_2000_200m = QRAIN_2000_200m.assign_coords({'west_east' : (QRAIN_2000_200m.west_east.values/5 +X3[0]+1)})
QRAIN_2000_200m = QRAIN_2000_200m.assign_coords({'south_north' : (QRAIN_2000_200m.south_north.values/5 +X3[0]+1)})

#%%
# making the same for the simulation with 500 m of reslution
path_data = os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_500m")
SIM_1KM_500m = Dataset(path_data + '/SIM_500m_domain2')
SIM_500m = Dataset(path_data + '/SIM_500m_domain3')

HGT = getvar(SIM_1KM_500m,'HGT')
z = getvar(SIM_1KM_500m,'z')
z_500 = getvar(SIM_500m,'z')
QRAIN_1KM_500m = getvar(SIM_1KM_500m,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_1KM_500m = interplevel(QRAIN_1KM_500m,z,2000)
QRAIN_500m = getvar(SIM_500m,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_500m = interplevel(QRAIN_500m,z_500,2000)

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_2000_1KM_500m = QRAIN_2000_1KM_500m.assign_coords({'west_east' : QRAIN_2000_1KM_500m.west_east.values +X2[0]+1})
QRAIN_2000_1KM_500m = QRAIN_2000_1KM_500m.assign_coords({'south_north' : QRAIN_2000_1KM_500m.south_north.values +X2[0]+1})
QRAIN_2000_500m = QRAIN_2000_500m.assign_coords({'west_east' : (QRAIN_2000_500m.west_east.values/2 +X3[0]+1)})
QRAIN_2000_500m = QRAIN_2000_500m.assign_coords({'south_north' : (QRAIN_2000_500m.south_north.values/2 +X3[0]+1)})


#%%
fig = plt.figure(figsize=(15,10))
plt.rcParams.update({'font.size': 18})

ax = plt.subplot(2, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_1KM_500m[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_1KM_500m[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(560.5, 550, 'Resolution: 1000 m', fontsize=18,
          bbox={'facecolor': 'grey', 'alpha': 0.5})
#rect = patches.Rectangle((588.5, 560), 10, 10, linewidth=2, edgecolor='r', facecolor='none')
plt.text(529., 648, '(a)', fontsize=22)
#plt.arrow(593.5, 570,0.1, 10, width = 0.4,color='r')
#ax.add_patch(rect)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('',fontsize = 22)
plt.ylabel('$\it{y}$ [km]',fontsize = 22)

ax = plt.subplot(2, 2, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_500m[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_500m[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(560.5, 550, 'Resolution: 500 m', fontsize=18,
         bbox={'facecolor': 'grey', 'alpha': 0.5})
plt.text(529., 648, '(c)', fontsize=22)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]',fontsize = 22)
plt.ylabel('$\it{y}$ [km]',fontsize = 22)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.24, 0.09, 0.5, 0.015])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg kg$^{-1}$', labelpad=5, y=0.95, rotation=0,fontsize=22)

ax = plt.subplot(2, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_1KM[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_1KM[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(560.5, 550, 'Resolution: 1000 m', fontsize=18,
         bbox={'facecolor': 'grey', 'alpha': 0.5})
plt.text(529., 648, '(b)',  fontsize=22)
plt.hlines(592,525,675, linestyle='dashed',color = 'black')
#plt.hlines(608,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('',fontsize = 22)
plt.ylabel('',fontsize = 22)

ax = plt.subplot(2, 2, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_200m[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_200m[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.2)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(560.5, 550, 'Resolution: 200 m', fontsize=18,
         bbox={'facecolor': 'grey', 'alpha': 0.5})
plt.text(529., 648, '(d)', fontsize=22)
plt.hlines(592,525,675, linestyle='dashed',color = 'black')
#plt.hlines(608,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('$\it{x}$ [km]',fontsize = 22)
plt.ylabel('',fontsize = 22)

plt.subplots_adjust(left=0.2,
                    bottom=0.2, 
                    right=0.78, 
                    top=0.9, 
                    wspace=0.25, 
                    hspace=0.25)

plt.savefig('figure4.png',dpi=300)

