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
path_data= os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_IDEALIZED_SOUNDING/SIM_IDEAL_WIND/")
SIM_U10 = Dataset(path_data + '/SIM_U10_dominio2')
SIM_U20 = Dataset(path_data + '/SIM_U20_dominio2')
SIM_U30 = Dataset(path_data + '/SIM_U30_dominio2')
SIM_shear = Dataset(path_data + '/SIM_SHEAR_FINALE_dominio2')

HGT = getvar(SIM_U10,'ter')
z = getvar(SIM_U10,'z')
z_agl_SIMU10 = getvar(SIM_U10,'height_agl')
z_agl_SIMU20 = getvar(SIM_U20,'height_agl')
z_agl_SIMU30 = getvar(SIM_U30,'height_agl')
z_agl_SIM_shear = getvar(SIM_shear,'height_agl')

v_SIMU10 = getvar(SIM_U10, "va",timeidx=6)
u_SIMU10 = getvar(SIM_U10, "ua",timeidx=6)
wspd_SIMU10 = getvar(SIM_U10, "wspd",timeidx=6)

v_SIMU20 = getvar(SIM_U20, "va",timeidx=6)       
u_SIMU20 = getvar(SIM_U20, "ua",timeidx=6)       
wspd_SIMU20 = getvar(SIM_U20, "wspd",timeidx=6)  

v_SIMU30 = getvar(SIM_U30, "va",timeidx=6)
u_SIMU30 = getvar(SIM_U30, "ua",timeidx=6)
wspd_SIMU30 = getvar(SIM_U30, "wspd",timeidx=6)

v_SIM_shear = getvar(SIM_shear, "va",timeidx=6)
u_SIM_shear = getvar(SIM_shear, "ua",timeidx=6)
wspd_SIM_shear = getvar(SIM_shear, "wspd",timeidx=6)


# interpolating on levels
u_SIMU10_100 = interplevel(u_SIMU10,z_agl_SIMU10,100)
v_SIMU10_100 = interplevel(v_SIMU10,z_agl_SIMU10,100)
wspd_SIMU10_100 = interplevel(wspd_SIMU10,z_agl_SIMU10,100)

u_SIMU20_100 = interplevel(u_SIMU20,z_agl_SIMU20,100)
v_SIMU20_100 = interplevel(v_SIMU20,z_agl_SIMU20,100)
wspd_SIMU20_100 = interplevel(wspd_SIMU20,z_agl_SIMU20,100)

u_SIMU30_100 = interplevel(u_SIMU30,z_agl_SIMU30,100)
v_SIMU30_100 = interplevel(v_SIMU30,z_agl_SIMU30,100)
wspd_SIMU30_100 = interplevel(wspd_SIMU30,z_agl_SIMU30,100)

u_SIM_shear_100 = interplevel(u_SIM_shear,z_agl_SIM_shear,100)
v_SIM_shear_100 = interplevel(v_SIM_shear,z_agl_SIM_shear,100)
wspd_SIM_shear_100 = interplevel(wspd_SIM_shear,z_agl_SIM_shear,100) 


# changing coordinates to the wind
v_SIMU10_100 = v_SIMU10_100.assign_coords({'west_east' : v_SIMU10_100.west_east.values +X2[0]+1})
v_SIMU10_100 = v_SIMU10_100.assign_coords({'south_north' : v_SIMU10_100.south_north.values +X2[0]+1})
u_SIMU10_100 = u_SIMU10_100.assign_coords({'west_east' : u_SIMU10_100.west_east.values +X2[0]+1})
u_SIMU10_100 = u_SIMU10_100.assign_coords({'south_north' : u_SIMU10_100.south_north.values +X2[0]+1})
wspd_SIMU10_100 = wspd_SIMU10_100.assign_coords({'west_east' : wspd_SIMU10_100.west_east.values +X2[0]+1})
wspd_SIMU10_100 = wspd_SIMU10_100.assign_coords({'south_north' : wspd_SIMU10_100.south_north.values +X2[0]+1})

v_SIMU20_100 = v_SIMU20_100.assign_coords({'west_east' : v_SIMU20_100.west_east.values +X2[0]+1})
v_SIMU20_100 = v_SIMU20_100.assign_coords({'south_north' : v_SIMU20_100.south_north.values +X2[0]+1})
u_SIMU20_100 = u_SIMU20_100.assign_coords({'west_east' : u_SIMU20_100.west_east.values +X2[0]+1})
u_SIMU20_100 = u_SIMU20_100.assign_coords({'south_north' : u_SIMU20_100.south_north.values +X2[0]+1})
wspd_SIMU20_100 = wspd_SIMU20_100.assign_coords({'west_east' : wspd_SIMU20_100.west_east.values +X2[0]+1})
wspd_SIMU20_100 = wspd_SIMU20_100.assign_coords({'south_north' : wspd_SIMU20_100.south_north.values +X2[0]+1})

v_SIMU30_100 = v_SIMU30_100.assign_coords({'west_east' : v_SIMU30_100.west_east.values +X2[0]+1})
v_SIMU30_100 = v_SIMU30_100.assign_coords({'south_north' : v_SIMU30_100.south_north.values +X2[0]+1})
u_SIMU30_100 = u_SIMU30_100.assign_coords({'west_east' : u_SIMU30_100.west_east.values +X2[0]+1})
u_SIMU30_100 = u_SIMU30_100.assign_coords({'south_north' : u_SIMU30_100.south_north.values +X2[0]+1})
wspd_SIMU30_100 = wspd_SIMU30_100.assign_coords({'west_east' : wspd_SIMU30_100.west_east.values +X2[0]+1})
wspd_SIMU30_100 = wspd_SIMU30_100.assign_coords({'south_north' : wspd_SIMU30_100.south_north.values +X2[0]+1})

v_SIM_shear_100 = v_SIM_shear_100.assign_coords({'west_east' : v_SIM_shear_100.west_east.values +X2[0]+1})
v_SIM_shear_100 = v_SIM_shear_100.assign_coords({'south_north' : v_SIM_shear_100.south_north.values +X2[0]+1})
u_SIM_shear_100 = u_SIM_shear_100.assign_coords({'west_east' : u_SIM_shear_100.west_east.values +X2[0]+1})
u_SIM_shear_100 = u_SIM_shear_100.assign_coords({'south_north' : u_SIM_shear_100.south_north.values +X2[0]+1})
wspd_SIM_shear_100 = wspd_SIM_shear_100.assign_coords({'west_east' : wspd_SIM_shear_100.west_east.values +X2[0]+1})
wspd_SIM_shear_100 = wspd_SIM_shear_100.assign_coords({'south_north' : wspd_SIM_shear_100.south_north.values +X2[0]+1})

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})

#%% making the subplot
#%%
fig = plt.figure(figsize=(10,10))
plt.rcParams.update({'font.size': 14})
ax = plt.subplot(2, 2, 1)
#div[::2,::2].plot.contourf(levels=np.linspace(-3*10**-3,3*10**-3,41),cmap='Reds',x='west_east',y = 'south_north')
plot=wspd_SIMU10_100[::2,::2].plot.contourf(levels=np.linspace(0,40,41),cmap='Reds',add_colorbar=False)
#cb=plt.colorbar(plot)
levels =[-200,-100,100,500,1000,1400,1499]
#CS=HGT[::2,::2].plot.contour(levels=levels,colors='k',alpha=0.5,linewidths=1.)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
qv = ax.quiver(u_SIMU10_100[::5,::5].west_east, u_SIMU10_100[::5,::5].south_north, u_SIMU10_100[::5,::5], v_SIMU10_100[::5,::5],color='black',scale =200.)
plt.title('')
ax.text(529, 648, '(a)', fontsize=18)
props = dict(facecolor='gray', alpha=0.8)
plt.text(528, 547.5, 'V10', fontsize=18, verticalalignment='center',
          bbox=props)
plt.xlabel('', fontsize=18)
plt.ylabel('$\it{y}$ [km]', fontsize=18)
plt.xlim((525,675))     
plt.ylim((540,660))

ax = plt.subplot(2, 2, 2)
plot=wspd_SIMU20_100[::2,::2].plot.contourf(levels=np.linspace(0,40,41),cmap='Reds',add_colorbar=False)
#cb=plt.colorbar(plot)
levels =[-200,-100,100,500,1000,1400,1499]
#CS=HGT[::2,::2].plot.contour(levels=levels,colors='k',alpha=0.5,linewidths=1.)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
qv = ax.quiver(u_SIMU20_100[::5,::5].west_east, u_SIMU20_100[::5,::5].south_north, u_SIMU20_100[::5,::5], v_SIMU20_100[::5,::5],color='black',scale = 200.)
plt.title('')
props = dict(facecolor='gray', alpha=0.8)
plt.text(528, 547.5, 'V20', fontsize=18, verticalalignment='center',
          bbox=props)
ax.text(529, 648, '(b)', fontsize=18)
plt.xlabel('', fontsize=18)
plt.ylabel('', fontsize=18)
plt.xlim((525,675))
plt.ylim((540,660))

ax = plt.subplot(2, 2, 3)             
plot=wspd_SIMU30_100[::2,::2].plot.contourf(levels=np.linspace(0,40,41),cmap='Reds',add_colorbar=False)
#cb=plt.colorbar(plot)
levels =[-200,-100,100,500,1000,1400,1499]
#CS=HGT[::2,::2].plot.contour(levels=levels,colors='k',alpha=0.5,linewidths=1.)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
qv = ax.quiver(u_SIMU30_100[::5,::5].west_east, u_SIMU30_100[::5,::5].south_north, u_SIMU30_100[::5,::5], v_SIMU30_100[::5,::5],color='black',scale = 200.)
plt.title('')
props = dict(facecolor='gray', alpha=0.8)
plt.text(528, 547.5, 'V30', fontsize=18, verticalalignment='center',
          bbox=props)
ax.text(529, 648, '(c)', fontsize=18)
plt.xlabel('$\it{x}$ [km]', fontsize=18)      
plt.ylabel('$\it{y}$ [km]', fontsize=18)          
plt.xlim((525,675))        
plt.ylim((540,660))

ax = plt.subplot(2, 2, 4)
plot=wspd_SIM_shear_100[::2,::2].plot.contourf(levels=np.linspace(0,40,41),cmap='Reds',add_colorbar=False)
#cb=plt.colorbar(plot)
levels =[-200,-100,100,500,1000,1400,1499]                  
#CS=HGT[::2,::2].plot.contour(levels=levels,colors='k',alpha=0.5,linewidths=1.)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
qv = ax.quiver(u_SIM_shear_100[::5,::5].west_east, u_SIM_shear_100[::5,::5].south_north, u_SIM_shear_100[::5,::5], v_SIM_shear_100[::5,::5],color='black',scale = 200.)
plt.title('')
props = dict(facecolor='gray', alpha=0.8)
plt.text(528, 548, 'V10_SHEAR', fontsize=18, verticalalignment='center',
          bbox=props)
ax.text(529, 648, '(d)', fontsize=18)
plt.xlabel('$\it{x}$ [km]', fontsize=18)
plt.ylabel('', fontsize=18)          
plt.xlim((525,675))
plt.ylim((540,660))     
cbar_ax = fig.add_axes([0.22, 0.12, 0.66, 0.01])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m s$^{-1}$',  fontsize=18, labelpad=0, y=0.95, rotation=0)

plt.subplots_adjust(left=0.2,
                    bottom=0.2,          
                    right=0.9,
                    top=0.75,
                    wspace=0.2,
                    hspace=0.15)

plt.savefig('figure11.png',dpi=300)
