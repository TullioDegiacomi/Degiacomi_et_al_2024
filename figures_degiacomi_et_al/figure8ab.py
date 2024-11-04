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
import matplotlib.font_manager
from matplotlib.cm import get_cmap,colors
from matplotlib.colors import TwoSlopeNorm
import os
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, get_basemap, ALL_TIMES)

from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, get_basemap, ALL_TIMES,CoordPair,vertcross, interpline)


font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 16,
        }

matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


#%% plot struttura bande tramite cross section orizzontale-->Roll Vortices

path_data= os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_IDEALIZED_SOUNDING/CTRL_SIMULATION/")
SIM_CTRL = Dataset(path_data + '/CTRL_SIM_domain2')


# definition of the domain for domain 2 in order to have KM on axis
i_start_2 = 163         
start2 = (i_start_2-1)*3 
cells2 = 225            
cells3 = 390

X2 = np.arange(start2,start2+cells2)
Y2 = X2

w =getvar(SIM_CTRL,"wa",timeidx=6)
u = getvar(SIM_CTRL, "ua",timeidx=6)  # default is in m/s
QCLOUD = getvar(SIM_CTRL,"QCLOUD",timeidx=6)
z = getvar(SIM_CTRL,'z')
HGT = getvar(SIM_CTRL,'HGT')

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})

start_point = CoordPair(x=90, y=110)
end_point = CoordPair(x=135,y=110)

w_cross = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
u_cross = vertcross(u,z,start_point=start_point,end_point=end_point)
# make a copy of the w cross section
w_cross_filled = np.ma.copy(to_np(w_cross))
u_cross_filled = np.ma.copy(to_np(u_cross))
q_cross_filled = np.ma.copy(to_np(qcloud_cross))
ter_line = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)

fig = plt.figure(figsize=(15,10))


for i in range(w_cross_filled.shape[-1]):
    column_vals = w_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    w_cross_filled[0:first_idx, i] = w_cross_filled[first_idx, i]

for i in range(u_cross_filled.shape[-1]):
    column_vals = u_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    u_cross_filled[0:first_idx, i] = u_cross_filled[first_idx, i]

for i in range(q_cross_filled.shape[-1]):
    column_vals = q_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    q_cross_filled[0:first_idx, i] = q_cross_filled[first_idx, i]
    
xs = np.arange(0, w_cross.shape[-1],1)+X2[0]
ys = to_np(w_cross.coords["vertical"])

w_modified = np.zeros((100,46))
w_modified[:] = np.nan
w_modified[::2,::]=w_cross_filled[::2,::]

u_modified = np.zeros((100,46))
u_modified[:] = np.nan
u_modified[::2,::]=u_cross_filled[::2,::]

#%% CROSS SECTION LUNGO LA BANDA
z = getvar(SIM_CTRL, "z",units="m",timeidx=6)
w =getvar(SIM_CTRL,"wa",timeidx=6)
v = getvar(SIM_CTRL, "va",timeidx=6)  # default is in m/s
#theta = getvar(ncfile,"theta",timeidx=6)
#thetae = getvar(ncfile,"theta_e",timeidx=6)
QCLOUD = getvar(SIM_CTRL,"QCLOUD",timeidx=6)
start_point = CoordPair(x=120,y=50)
end_point = CoordPair(x=120, y=130)
# Get the terrain heights along the cross section line
ter_line_2 = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)
# getting the cross sections of the variables
w_cross_2 = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross_2 = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
v_cross_2 = vertcross(v,z,start_point=start_point,end_point=end_point)

q_cross_filled_2 = np.ma.copy(to_np(qcloud_cross_2))

for i in range(q_cross_filled_2.shape[-1]):
    column_vals = q_cross_filled_2[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    q_cross_filled_2[0:first_idx, i] = q_cross_filled_2[first_idx, i]


w_modified_2 = np.zeros((100,81))
w_modified_2[:] = np.nan
w_modified_2[::2,::5]=w_cross_2[::2,::5]
v_modified_2 = np.zeros((100,81))
v_modified_2[:] = np.nan
v_modified_2[::2,::5]=v_cross_2[::2,::5]

xs_2 = np.arange(0, w_cross_2.shape[-1], 1)
ys_2 = to_np(w_cross_2.coords["vertical"])
Xs_2, Ys_2 = np.meshgrid(xs_2,ys_2)

#making the subplot
plt.rcParams.update({'font.size': 20})
fig, axs = plt.subplots(1,2,figsize=(18, 6))

norm = TwoSlopeNorm(vmin=w_cross_filled.min(),vcenter=0, vmax=w_cross_filled.max())
levs = [10**(-6)]
plot = axs[0].contourf(xs+90,ys,w_cross_filled,norm=norm,levels=30, cmap='RdBu_r')
#cb = fig.colorbar(plot,orientation='horizontal')
#cb.ax.tick_params(labelsize=18)
#cb.set_label('m/s', labelpad=-40, y=1.05, rotation=0,fontsize=18)
axs[0].quiver(xs+90,ys,u_modified/2,w_modified,scale=100.,color='black',alpha=1.)
axs[0].set_ylim((0,9000))
axs[0].set_xlabel('$\it{x}$ [km]',fontsize = 24)
axs[0].set_ylabel('$\it{z}$ [m]',fontsize = 24)

axs[0].contour(xs+90,ys,q_cross_filled,levels=levs,linestyles='dashed',linewidths=1.8,colors='black',alpha=0.8)
axs[0].fill_between(xs+90,0,to_np(ter_line),facecolor="black",alpha=0.93)
axs[0].text(578,8200,'(a)', fontsize=24)

#
levs_w = [0.3,0.5,1,2,4,8,12,15,20,25]
plot2 = axs[1].contourf(xs_2+ Y2[0]+51,ys_2,q_cross_filled_2,levels=np.linspace(0,0.0024,25),cmap='Blues')
axs[1].quiver(Xs_2+Y2[0]+51,Ys_2,v_modified_2/8,w_modified_2,scale=120.,alpha=1)
axs[1].fill_between(xs_2+Y2[0]+51,-300,to_np(ter_line_2),facecolor="black")
CS=axs[1].contour(xs_2+Y2[0]+51,ys_2,w_cross_2,levels=levs_w,colors='black',linewidths=0.8)
axs[1].clabel(CS, inline=1, fontsize=14,fmt='%.1f')
axs[1].set_ylim((0,9000))
axs[1].set_xlabel('$\it{y}$ [km]',fontsize = 24)
axs[1].set_ylabel('',fontsize = 24)
axs[1].text(540,8200,'(b)', fontsize=24)

from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.21, 0.12, 0.3, 0.03])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m s$^{-1}$', rotation=0,fontsize=24)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.59, 0.12, 0.3, 0.03])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=4)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg kg$^{-1}$', labelpad=-1, y=0.95, rotation=0,fontsize=24)

plt.subplots_adjust(left=0.2,
                    bottom=0.3, 
                    right=0.9, 
                    top=0.95, 
                    wspace=0.2, 
                    hspace=0.5)

plt.savefig('figure8ab.png', dpi=300)

