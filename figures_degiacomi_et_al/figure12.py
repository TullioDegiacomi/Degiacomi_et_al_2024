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


import metpy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import SkewT,Hodograph
from metpy.units import units

SIM_TILTED = pd.read_csv("Parameters_Sondaggio_wind_shear_tilted_FINALE.txt",delimiter = '\t')
SIM_ROTATED = pd.read_csv("Parameters_Sondaggio_original_direction_sum20.txt",delimiter='\t')
SIM_SHEAR = pd.read_csv("Sondaggio_U10_double_shear_finale.txt",delimiter='\t')

u=SIM_ROTATED['u']
v=SIM_ROTATED['v']
hght=SIM_ROTATED['Height']
u2 = SIM_TILTED['u']
v2 = SIM_TILTED['v']
hght_2 =SIM_TILTED['Height']
u3 = SIM_SHEAR['u']
v3 = SIM_SHEAR['v']
hght_3 =SIM_SHEAR['Height']
# computing the wind speed profiles
speed = np.sqrt(u**2+v**2)
speed2 = np.sqrt(u2**2+v2**2)
speed3 = np.sqrt(u3**2+v3**2)


#%%
import matplotlib.gridspec as gridspec
u = u[0:985]
v = v[0:985]
u2 = u2[0:985]
v2 = v2[0:985]
u3 = u3[:985]
v3 = v3[:985]

fig = plt.figure(figsize=(14,12))
plt.rcParams.update({'font.size': 18})
gs = gridspec.GridSpec(1, 3, width_ratios=[1,1,1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])

#ax1 = plt.subplot(1, 2, 1)
h = Hodograph(ax=ax1,component_range=40.)
h.add_grid(increment=20)
cb = h.plot_colormapped(u, v, hght[0:985],linewidth=2.5)
#plt.colorbar(cb,shrink=0.7)
ax1.set_xlabel('$\it{u}$ [m s$^{-1}$]',labelpad=0, fontsize=22)
ax1.set_ylabel('$\it{v}$ [m s$^{-1}$]',labelpad=-5, fontsize=22)
ax1.set_title('UDINE_ROT20',y=1.03,fontsize=18)
ax1.text(-37,33,'(a)', fontsize=22)
plt.xticks()
plt.yticks()
h.wind_vectors(u[::100],v[::100],scale = 1.,width = 0.4)

#ax2 = plt.subplot(1, 2, 2)
h = Hodograph(ax=ax2,component_range=40.)
h.add_grid(increment=20)
cb = h.plot_colormapped(u2, v2, hght[0:985],linewidth=2.5)
#plt.colorbar(cb,shrink=0.7)
ax2.set_xlabel('$\it{u}$ [m s$^{-1}$]',labelpad=0, fontsize=22)
ax2.set_ylabel('$\it{v}$ [m s$^{-1}$]',labelpad=-5, fontsize=22)
ax2.set_title('SHEAR_TILTED',y=1.03,fontsize=18)
ax2.text(-37,33,'(b)', fontsize=22)
plt.xticks()
plt.yticks()
h.wind_vectors(u2[::100],v2[::100],scale = 1.,width = 0.4) 

#ax2 = plt.subplot(1, 2, 2)
ax3.plot(speed,hght,color='b',linewidth=2,alpha=0.8,label='CTRL')
ax3.plot(speed2,hght,color='g',alpha=0.7,linewidth=2,label='SHEAR_TILTED')
ax3.plot(speed3,hght_3,color='orange',linewidth=2,label='V10_SHEAR')
ax3.set_ylim((0,13000))
ax3.text(5,11800,'(c)', fontsize=22)
ax3.legend(bbox_to_anchor=(0.715,0.55),fontsize=14)
ax3.set_xlabel('Wind speed [m s$^{-1}$]', fontsize=22)
ax3.set_ylabel('$\it{z}$ [m]',labelpad=0, fontsize=22)
ax3.grid()

from matplotlib import ticker
#fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.17, 0.13, 0.4, 0.01])
cb = fig.colorbar(cb,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m MSL', labelpad=3, y=0.9, rotation=0,fontsize=22)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.5, 
                    wspace=0.4, 
                    hspace=0.3)

plt.savefig('figure12.png',dpi=600)
