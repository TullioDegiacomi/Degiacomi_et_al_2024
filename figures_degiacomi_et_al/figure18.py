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


path_data= os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_IDEALIZED_SOUNDING/SIM_IDEAL_N/")
SIM_N3_00004 = Dataset(path_data + '/SIM_N1_00004_N3_00004_dominio2')
SIM_N3_00009 = Dataset(path_data + '/SIM_N1_00004_N3_00009_dominio2')
SIM_N3_00012 = Dataset(path_data + '/SIM_N1_00004_N3_00012_dominio2')

z = getvar(SIM_N3_00004,'z')
HGT = getvar(SIM_N3_00004,'HGT')
QRAIN = getvar(SIM_N3_00004,'QRAIN',timeidx=ALL_TIMES)

z = getvar(SIM_N3_00009,'z')
HGT = getvar(SIM_N3_00009,'HGT')
QRAIN = getvar(SIM_N3_00009,'QRAIN',timeidx=ALL_TIMES)

z = getvar(SIM_N3_00012,'z')
HGT = getvar(SIM_N3_00012,'HGT')
QRAIN = getvar(SIM_N3_00012,'QRAIN',timeidx=ALL_TIMES)

w =getvar(SIM_N3_00009,"wa",timeidx=6)
v = getvar(SIM_N3_00009, "va",timeidx=6)  # default is in m/s
QCLOUD = getvar(SIM_N3_00009,"QCLOUD",timeidx=6)
theta_e = getvar(SIM_N3_00009,"eth",timeidx=6)

start_point = CoordPair(x=141, y=50)
end_point = CoordPair(x=141,y=150)

w_cross = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
v_cross = vertcross(v,z,start_point=start_point,end_point=end_point)
thetae_cross = vertcross(theta_e,z,start_point=start_point,end_point=end_point)
# make a copy of the w cross section
w_cross_filled = np.ma.copy(to_np(w_cross))
v_cross_filled = np.ma.copy(to_np(v_cross))
q_cross_filled = np.ma.copy(to_np(qcloud_cross))
ter_line = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)


for i in range(w_cross_filled.shape[-1]):
    column_vals = w_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    w_cross_filled[0:first_idx, i] = w_cross_filled[first_idx, i]

for i in range(v_cross_filled.shape[-1]):
    column_vals = v_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    v_cross_filled[0:first_idx, i] = v_cross_filled[first_idx, i]

for i in range(q_cross_filled.shape[-1]):
    column_vals = q_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    q_cross_filled[0:first_idx, i] = q_cross_filled[first_idx, i]
    
xs = np.arange(0, w_cross.shape[-1],1)
ys = to_np(w_cross.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

w_modified = np.zeros((100,101))
w_modified[:] = np.nan
w_modified[::2,::5]=w_cross_filled[::2,::5]

v_modified = np.zeros((100,101))
v_modified[:] = np.nan
v_modified[::2,::5]=v_cross_filled[::2,::5]

z = getvar(SIM_N3_00004, "z",units="m",timeidx=6)
w =getvar(SIM_N3_00004,"wa",timeidx=6)
v = getvar(SIM_N3_00004, "va",timeidx=6)  
thetae = getvar(SIM_N3_00004,"eth",timeidx=6)
QCLOUD = getvar(SIM_N3_00004,"QCLOUD",timeidx=6)
start_point = CoordPair(x=134,y=50)
end_point = CoordPair(x=134, y=150)
# Get the terrain heights along the cross section line
ter_line_2 = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)
# getting the cross sections of the variables
w_cross_2 = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross_2 = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
v_cross_2 = vertcross(v,z,start_point=start_point,end_point=end_point)
thetae_cross_2 = vertcross(thetae,z,start_point=start_point,end_point=end_point)

w_cross_filled_2 = np.ma.copy(to_np(w_cross_2))

for i in range(w_cross_filled_2.shape[-1]):
    column_vals = w_cross_filled_2[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    w_cross_filled_2[0:first_idx, i] = w_cross_filled_2[first_idx, i]


w_modified_2 = np.zeros((100,101))
w_modified_2[:] = np.nan
w_modified_2[::2,::5]=w_cross_2[::2,::5]
v_modified_2 = np.zeros((100,101))
v_modified_2[:] = np.nan
v_modified_2[::2,::5]=v_cross_2[::2,::5]

xs_2 = np.arange(0, w_cross_2.shape[-1], 1)
ys_2 = to_np(w_cross_2.coords["vertical"])
Xs_2, Ys_2 = np.meshgrid(xs_2,ys_2)


z = getvar(SIM_N3_00012, "z",units="m",timeidx=6)
w =getvar(SIM_N3_00012,"wa",timeidx=6)
v = getvar(SIM_N3_00012, "va",timeidx=6)
thetae = getvar(SIM_N3_00012,"eth",timeidx=6)
QCLOUD = getvar(SIM_N3_00012,"QCLOUD",timeidx=6)
start_point = CoordPair(x=118,y=50)
end_point = CoordPair(x=118, y=150)
# Get the terrain heights along the cross section line
ter_line_3 = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)
# getting the cross sections of the variables
w_cross_3 = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross_3 = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
v_cross_3 = vertcross(v,z,start_point=start_point,end_point=end_point)
thetae_cross_3 = vertcross(thetae,z,start_point=start_point,end_point=end_point)

w_cross_filled_3 = np.ma.copy(to_np(w_cross_3))

for i in range(w_cross_filled_3.shape[-1]):
    column_vals = w_cross_filled_3[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    w_cross_filled_3[0:first_idx, i] = w_cross_filled_3[first_idx, i]


w_modified_3 = np.zeros((100,101))
w_modified_3[:] = np.nan
w_modified_3[::2,::5]=w_cross_3[::2,::5]     
v_modified_3 = np.zeros((100,101))
v_modified_3[:] = np.nan
v_modified_3[::2,::5]=v_cross_3[::2,::5]

xs_3 = np.arange(0, w_cross_3.shape[-1], 1)
ys_3 = to_np(w_cross_3.coords["vertical"])          
Xs_3, Ys_3 = np.meshgrid(xs_3,ys_3)


# making the plot

#%%
fig = plt.figure(figsize=(10,15))
plt.rcParams.update({'font.size': 20})

ax = plt.subplot(3, 1, 1)
plot2 = ax.contourf(xs_2+ Y2[0]+51,ys_2,w_cross_filled_2,levels=np.linspace(-16,16,33),cmap='RdBu_r')
CS=ax.contour(xs_2+Y2[0]+51,ys_2,thetae_cross_2,levels=[316,318,320,322,324,326,328,330,332,334],colors=['black'])
props = dict(facecolor='gray', alpha=0.8)
plt.text(571, 8430, 'N1_00004_N3_00004', fontsize=24, verticalalignment='center',
          bbox=props)
ax.clabel(CS, inline=3, fontsize=18,fmt='%.0f')
ax.fill_between(xs_2+Y2[0]+51,-300,to_np(ter_line_2),facecolor="black")
ax.set_ylim((0,9000))
ax.set_xlabel('',fontsize = 24)
ax.set_ylabel('$\it{z}$ [m]',fontsize = 24)
ax.text(540,8000,'(a)', fontsize=24)

ax = plt.subplot(3, 1, 2)
plot2 = ax.contourf(xs+ Y2[0]+51,ys,w_cross_filled,levels=np.linspace(-16,16,33),cmap='RdBu_r')
CS=ax.contour(xs+Y2[0]+51,ys,thetae_cross,levels=[316,318,320,322,324,326,328,330,332,334],colors=['black'])
props = dict(facecolor='gray', alpha=0.8)
plt.text(571, 8430, 'N1_00004_N3_00009', fontsize=24, verticalalignment='center',
          bbox=props)
ax.clabel(CS, inline=2, fontsize=18,fmt='%.0f')
ax.fill_between(xs+Y2[0]+51,-300,to_np(ter_line),facecolor="black")
ax.set_ylim((0,9000))
ax.set_xlabel('',fontsize = 24)
ax.set_ylabel('$\it{z}$ [m]',fontsize = 24)
ax.text(540,8000,'(b)',fontsize=24)

ax = plt.subplot(3, 1, 3)
plot2 = ax.contourf(xs_3+ Y2[0]+51,ys_3,w_cross_filled_3,levels=np.linspace(-16,16,33),cmap='RdBu_r')
CS=ax.contour(xs_3+Y2[0]+51,ys_3,thetae_cross_3,levels=[316,318,320,322,324,326,328,330,332,334],colors=['black'])
props = dict(facecolor='gray', alpha=0.8)
plt.text(571, 8430, 'N1_00004_N3_00012', fontsize=24, verticalalignment='center',
          bbox=props)
ax.clabel(CS, inline=2, fontsize=18,fmt='%.0f')
ax.fill_between(xs_3+Y2[0]+51,-300,to_np(ter_line),facecolor="black")
ax.set_ylim((0,9000))
ax.set_xlabel('$\it{y}$ [km]',fontsize = 24)
ax.set_ylabel('$\it{z}$ [m]',fontsize = 24)
ax.text(540,8000,'(c)',fontsize=24)


from matplotlib import ticker

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.22, 0.13, 0.42, 0.01])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=6)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m s$^{-1}$',  fontsize=24, labelpad=0, y=0.95, rotation=0)

plt.subplots_adjust(left=0.2,
                    bottom=0.2,
                    right=0.66,             
                    top=0.95,
                    wspace=0.2,
                    hspace=0.2)

plt.savefig('figure18.png',dpi=300)
