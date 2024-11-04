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

from matplotlib import ticker

# definition of the domain for domain 2 in order to have KM on axis
i_start_2 = 163
start2 = (i_start_2-1)*3
cells2 = 225

X2 = np.arange(start2,start2+cells2) 
Y2 = X2

#%% PLOT WIND
path_data= os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_IDEALIZED_SOUNDING/SIM_IDEAL_WIND_DIRECTION/")

SIM_210 = Dataset(path_data + '/SIM_wind_210degrees_dominio2')

HGT = getvar(SIM_210,'ter')
z = getvar(SIM_210,'z')
z_agl_SIM210 = getvar(SIM_210,'height_agl')

v_SIM210 = getvar(SIM_210, "va",timeidx=ALL_TIMES)
u_SIM210 = getvar(SIM_210, "ua",timeidx=ALL_TIMES)
wspd_SIM210 = getvar(SIM_210, "wspd",timeidx=ALL_TIMES)

# interpolating on levels
u_SIM210_100 = interplevel(u_SIM210,z_agl_SIM210,100)
v_SIM210_100 = interplevel(v_SIM210,z_agl_SIM210,100)
wspd_SIM210_100 = interplevel(wspd_SIM210,z_agl_SIM210,100)

u_SIM210_2000 = interplevel(u_SIM210,z,2000)
v_SIM210_2000 = interplevel(v_SIM210,z,2000)
wspd_SIM210_2000 = interplevel(wspd_SIM210,z,2000)

# changing coordinates to the wind
v_SIM210_100 = v_SIM210_100.assign_coords({'west_east' : v_SIM210_100.west_east.values +X2[0]+1})
v_SIM210_100 = v_SIM210_100.assign_coords({'south_north' : v_SIM210_100.south_north.values +X2[0]+1})
u_SIM210_100 = u_SIM210_100.assign_coords({'west_east' : u_SIM210_100.west_east.values +X2[0]+1})
u_SIM210_100 = u_SIM210_100.assign_coords({'south_north' : u_SIM210_100.south_north.values +X2[0]+1})
wspd_SIM210_100 = wspd_SIM210_100.assign_coords({'west_east' : wspd_SIM210_100.west_east.values +X2[0]+1})
wspd_SIM210_100 = wspd_SIM210_100.assign_coords({'south_north' : wspd_SIM210_100.south_north.values +X2[0]+1})

v_SIM210_2000 = v_SIM210_2000.assign_coords({'west_east' : v_SIM210_2000.west_east.values +X2[0]+1})
v_SIM210_2000 = v_SIM210_2000.assign_coords({'south_north' : v_SIM210_2000.south_north.values +X2[0]+1})
u_SIM210_2000 = u_SIM210_2000.assign_coords({'west_east' : u_SIM210_2000.west_east.values +X2[0]+1})
u_SIM210_2000 = u_SIM210_2000.assign_coords({'south_north' : u_SIM210_2000.south_north.values +X2[0]+1})
wspd_SIM210_2000 = wspd_SIM210_2000.assign_coords({'west_east' : wspd_SIM210_2000.west_east.values +X2[0]+1})
wspd_SIM210_2000 = wspd_SIM210_2000.assign_coords({'south_north' : wspd_SIM210_2000.south_north.values +X2[0]+1})

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})

#%% making the subplot
fig = plt.figure(figsize=(10,5))
plt.rcParams.update({'font.size': 14})
ax = plt.subplot(1, 2, 1)
#div[::2,::2].plot.contourf(levels=np.linspace(-3*10**-3,3*10**-3,41),cmap='Reds',x='west_east',y = 'south_north')
plot=wspd_SIM210_100[6].plot.contourf(levels=np.linspace(0,50,51),cmap='Reds',add_colorbar=False)
#cb=plt.colorbar(plot)
levels =[-200,-100,100,500,1000,1400,1499]
#CS=HGT[::2,::2].plot.contour(levels=levels,colors='k',alpha=0.5,linewidths=1.)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
qv = ax.quiver(u_SIM210_100[6,::5,::5].west_east, u_SIM210_100[6,::5,::5].south_north, u_SIM210_100[6,::5,::5], v_SIM210_100[6,::5,::5],color='black',scale =400.)
plt.title('', fontsize = 18)
ax.text(529, 648, '(a)', fontsize=18)
props = dict(facecolor='gray', alpha=0.8)
#plt.text(528, 547.5, '210$^\circ$', fontsize=18, verticalalignment='center',
#          bbox=props)
plt.xlabel('$\it{x}$ [km]', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
plt.xlim((525,675))     
plt.ylim((540,660))

ax = plt.subplot(1, 2, 2)             
plot=wspd_SIM210_2000[6].plot.contourf(levels=np.linspace(0,50,51),cmap='Reds',add_colorbar=False)
#cb=plt.colorbar(plot)
levels =[-200,-100,100,500,1000,1400,1499]
#CS=HGT[::2,::2].plot.contour(levels=levels,colors='k',alpha=0.5,linewidths=1.)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
qv = ax.quiver(u_SIM210_2000[6,::5,::5].west_east, u_SIM210_2000[6,::5,::5].south_north, u_SIM210_2000[6,::5,::5], v_SIM210_2000[6,::5,::5],color='black',scale = 400.)
plt.title('')
props = dict(facecolor='gray', alpha=0.8)
#plt.text(528, 547.5, 'V30', fontsize=18, verticalalignment='center',
#          bbox=props)
ax.text(529, 648, '(b)', fontsize=18)
plt.xlabel('$\it{x}$ [km]', fontsize=18)      
plt.ylabel('', fontsize=18)          
plt.xlim((525,675))        
plt.ylim((540,660))

cbar_ax = fig.add_axes([0.22, 0.20, 0.715, 0.01])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m s$^{-1}$',  fontsize=18, labelpad=0, y=0.95, rotation=0)

plt.subplots_adjust(left=0.2,
                    bottom=0.35,          
                    right=0.955,
                    top=0.9,
                    wspace=0.2,
                    hspace=0.15)

plt.savefig('figure14.png',dpi=300)
