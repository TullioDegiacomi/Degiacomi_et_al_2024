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


w =getvar(SIM_U20,"wa",timeidx=6)
v = getvar(SIM_U20, "va",timeidx=6)  # default is in m/s
QCLOUD = getvar(SIM_U20,"QCLOUD",timeidx=6)
theta_e = getvar(SIM_U20,"eth",timeidx=6)

start_point = CoordPair(x=102, y=50)
end_point = CoordPair(x=102,y=150)

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

z = getvar(SIM_U10, "z",units="m",timeidx=6)
w =getvar(SIM_U10,"wa",timeidx=6)
v = getvar(SIM_U10, "va",timeidx=6)  
thetae = getvar(SIM_U10,"eth",timeidx=6)
QCLOUD = getvar(SIM_U10,"QCLOUD",timeidx=6)
start_point = CoordPair(x=98,y=50)
end_point = CoordPair(x=98, y=150)
# Get the terrain heights along the cross section line
ter_line_2 = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)
# getting the cross sections of the variables
w_cross_2 = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross_2 = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
v_cross_2 = vertcross(v,z,start_point=start_point,end_point=end_point)
thetae_cross_2 = vertcross(thetae,z,start_point=start_point,end_point=end_point)

q_cross_filled_2 = np.ma.copy(to_np(qcloud_cross_2))

for i in range(q_cross_filled_2.shape[-1]):
    column_vals = q_cross_filled_2[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    q_cross_filled_2[0:first_idx, i] = q_cross_filled_2[first_idx, i]


w_modified_2 = np.zeros((100,101))
w_modified_2[:] = np.nan
w_modified_2[::2,::5]=w_cross_2[::2,::5]
v_modified_2 = np.zeros((100,101))
v_modified_2[:] = np.nan
v_modified_2[::2,::5]=v_cross_2[::2,::5]

xs_2 = np.arange(0, w_cross_2.shape[-1], 1)
ys_2 = to_np(w_cross_2.coords["vertical"])
Xs_2, Ys_2 = np.meshgrid(xs_2,ys_2)


z = getvar(SIM_U30, "z",units="m",timeidx=6)
w =getvar(SIM_U30,"wa",timeidx=6)
v = getvar(SIM_U30, "va",timeidx=6)
thetae = getvar(SIM_U30,"eth",timeidx=6)
QCLOUD = getvar(SIM_U30,"QCLOUD",timeidx=6)
start_point = CoordPair(x=100,y=50)
end_point = CoordPair(x=100, y=150)
# Get the terrain heights along the cross section line
ter_line_3 = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)
# getting the cross sections of the variables
w_cross_3 = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross_3 = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
v_cross_3 = vertcross(v,z,start_point=start_point,end_point=end_point)
thetae_cross_3 = vertcross(thetae,z,start_point=start_point,end_point=end_point)

q_cross_filled_3 = np.ma.copy(to_np(qcloud_cross_3))

for i in range(q_cross_filled_3.shape[-1]):
    column_vals = q_cross_filled_3[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    q_cross_filled_3[0:first_idx, i] = q_cross_filled_3[first_idx, i]


w_modified_3 = np.zeros((100,101))
w_modified_3[:] = np.nan
w_modified_3[::2,::5]=w_cross_3[::2,::5]     
v_modified_3 = np.zeros((100,101))
v_modified_3[:] = np.nan
v_modified_3[::2,::5]=v_cross_3[::2,::5]

xs_3 = np.arange(0, w_cross_3.shape[-1], 1)
ys_3 = to_np(w_cross_3.coords["vertical"])          
Xs_3, Ys_3 = np.meshgrid(xs_3,ys_3)


z = getvar(SIM_shear, "z",units="m",timeidx=6)
w =getvar(SIM_shear,"wa",timeidx=6)
v = getvar(SIM_shear, "va",timeidx=6)
thetae = getvar(SIM_shear,"eth",timeidx=6)
QCLOUD = getvar(SIM_shear,"QCLOUD",timeidx=6)
start_point = CoordPair(x=128,y=50)
end_point = CoordPair(x=128, y=150)
# Get the terrain heights along the cross section line
ter_line_4 = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)
# getting the cross sections of the variables
w_cross_4 = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross_4 = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
v_cross_4 = vertcross(v,z,start_point=start_point,end_point=end_point)
thetae_cross_4 = vertcross(thetae,z,start_point=start_point,end_point=end_point)

q_cross_filled_4 = np.ma.copy(to_np(qcloud_cross_4))

for i in range(q_cross_filled_4.shape[-1]):
    column_vals = q_cross_filled_4[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    q_cross_filled_4[0:first_idx, i] = q_cross_filled_4[first_idx, i]


w_modified_4 = np.zeros((100,101))
w_modified_4[:] = np.nan
w_modified_4[::2,::5]=w_cross_4[::2,::5]
v_modified_4 = np.zeros((100,101))
v_modified_4[:] = np.nan
v_modified_4[::2,::5]=v_cross_4[::2,::5]

xs_4 = np.arange(0, w_cross_4.shape[-1], 1)
ys_4 = to_np(w_cross_4.coords["vertical"])
Xs_4, Ys_4 = np.meshgrid(xs_4,ys_4)


# making the plot


#%%
fig = plt.figure(figsize=(15,10))
plt.rcParams.update({'font.size': 20})

ax = plt.subplot(2, 2, 1)
plot2 = ax.contourf(xs_2+ Y2[0]+51,ys_2,q_cross_filled_2,levels=np.linspace(0,0.0024,25),cmap='Blues')
ax.quiver(Xs_2+Y2[0]+51,Ys_2,v_modified_2/8,w_modified_2,scale=100.,alpha=1.)
CS=ax.contour(xs_2+Y2[0]+51,ys_2,thetae_cross_2,levels=[322,324,326,328],colors=['black'])
props = dict(facecolor='gray', alpha=0.8)
plt.text(623.5, 8430, 'V10', fontsize=24, verticalalignment='center',
          bbox=props)
ax.clabel(CS, inline=3, fontsize=18,fmt='%.0f')
ax.fill_between(xs_2+Y2[0]+51,-300,to_np(ter_line_2),facecolor="black")
ax.set_ylim((0,9000))
ax.set_xlabel('',fontsize = 24)
ax.set_ylabel('$\it{z}$ [m]',fontsize = 24)
ax.text(540,8000,'(a)', fontsize=24)

ax = plt.subplot(2, 2, 2)
plot2 = ax.contourf(xs+ Y2[0]+51,ys,q_cross_filled,levels=np.linspace(0,0.0024,25),cmap='Blues')
ax.quiver(Xs+Y2[0]+51,Ys,v_modified/8,w_modified,scale=100.,alpha=1.)
CS=ax.contour(xs+Y2[0]+51,ys,thetae_cross,levels=[322,324,326,328],colors=['black'])
props = dict(facecolor='gray', alpha=0.8)
plt.text(623.5, 8430, 'V20', fontsize=24, verticalalignment='center',
          bbox=props)
ax.clabel(CS, inline=2, fontsize=18,fmt='%.0f')
ax.fill_between(xs+Y2[0]+51,-300,to_np(ter_line),facecolor="black")
ax.set_ylim((0,9000))
ax.set_xlabel('',fontsize = 24)
ax.set_ylabel('',fontsize = 24)
ax.text(540,8000,'(b)',fontsize=24)

ax = plt.subplot(2, 2, 3)
plot2 = ax.contourf(xs_3+ Y2[0]+51,ys_3,q_cross_filled_3,levels=np.linspace(0,0.0024,25),cmap='Blues')
ax.quiver(Xs_3+Y2[0]+51,Ys_3,v_modified_3/8,w_modified_3,scale=100.,alpha=1.)
CS=ax.contour(xs_3+Y2[0]+51,ys_3,thetae_cross_3,levels=[322,324,326,328],colors=['black'])
#CS=ax.contour(xs_3+Y2[0]+51,ys_3,w_cross_3,levels=[0.5,1,4,8,14],colors=['black'])
props = dict(facecolor='gray', alpha=0.8)
plt.text(623.5, 8430, 'V30', fontsize=24, verticalalignment='center',
          bbox=props)
ax.clabel(CS, inline=2, fontsize=18,fmt='%.0f')
ax.fill_between(xs_3+Y2[0]+51,-300,to_np(ter_line),facecolor="black")
ax.set_ylim((0,9000))
ax.set_xlabel('$\it{y}$ [km]',fontsize = 24)
ax.set_ylabel('$\it{z}$ [m]',fontsize = 24)
ax.text(540,8000,'(c)',fontsize=24)

ax = plt.subplot(2, 2, 4)
plot2 = ax.contourf(xs_4+ Y2[0]+51,ys_4,q_cross_filled_4,levels=np.linspace(0,0.0024,25),cmap='Blues')
ax.quiver(Xs_4+Y2[0]+51,Ys_4,v_modified_4/8,w_modified_4,scale=100.,alpha=1.)
CS=ax.contour(xs_4+Y2[0]+51,ys_4,thetae_cross_4,levels=[322,324,326,328],colors=['black'])
#CS=ax.contour(xs_4+Y2[0]+51,ys_4,w_cross_4,levels=[0.5,1,4,8,14],colors=['black'])
props = dict(facecolor='gray', alpha=0.8)
plt.text(597.5, 8430, 'V10_SHEAR', fontsize=24, verticalalignment='center',
          bbox=props)
ax.clabel(CS, inline=2, fontsize=18,fmt='%.0f')
ax.fill_between(xs_4+Y2[0]+51,-300,to_np(ter_line),facecolor="black")
ax.set_ylim((0,9000))
ax.set_xlabel('$\it{y}$ [km]',fontsize = 24)
ax.set_ylabel('',fontsize = 24)
ax.text(540,8000,'(d)',fontsize=24)


from matplotlib import ticker

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.1, 0.6, 0.01])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=6)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg kg$^{-1}$',  fontsize=24, labelpad=0, y=0.95, rotation=0)

plt.subplots_adjust(left=0.2,
                    bottom=0.2,
                    right=0.9,             
                    top=0.95,
                    wspace=0.2,
                    hspace=0.2)

plt.savefig('figure10.png',dpi=300)
