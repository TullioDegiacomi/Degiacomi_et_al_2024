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


path_data = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/Simulazioni_Coriolis/SIM_500m_005K")
path_output = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini")
#path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIMULAZIONI_FINALI")

ncfile1 = Dataset(path_data + '/SIM_500m_005K_dominio1')
ncfile2 = Dataset(path_data + '/SIM_500m_005K_dominio2')
ncfile3 = Dataset(path_data + '/SIM_500m_005K_dominio3')

# Get the latitude and longitude
TIMES = getvar(ncfile1,"times",timeidx=0)

# Estraggo orografia
HGT = getvar(ncfile1,"HGT")
eta = getvar(ncfile1,"z")
#%% 
import numpy as np

# setting X and Y for domain 1 with resolution of 3 km
cells1 = 400
X1 = np.arange(1.5,1201.5,3) 
Y1=X1


# definition of the domain for domain 2 and 3 in order to have KM on axis
i_start_2 = 163
#i_start_3 = 35
i_start_3 = 41
j_start_3 = 51
#start2 = i_start_2*3   # multiplying for resolution of domain 1
start2 = (i_start_2-1)*3 +0.5
cells2 = 225
cells3 = 390
cells3_x = 975
cells3_y = 975

X2 = np.arange(start2,start2+cells2) 
Y2 = X2

#setting domain 3 for 500 m resolution
# start3 = (start2-0.5) + 16 + 0.25
# X3 = np.arange(start3,start3+cells3/2,0.5)
# Y3 = X3

#setting domain 3 for 200 m of resolution
start3_x = (start2-0.5) + 16
start3_y = (start2-0.5) + 16
X3 = np.arange(start3_x,start3_x + cells3_x/5,0.2) 
Y3 = np.arange(start3_y,start3_y + cells3_y/5,0.2)

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
from matplotlib.colors import DivergingNorm, TwoSlopeNorm
import os
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, get_basemap, ALL_TIMES,CoordPair,vertcross, interpline)

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
#%% plot QRain sensitivity on resolution
path_data= os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/Simulazioni_Coriolis/SIM_200m_005K")
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIMULAZIONI_DOMINIO_CENTRATO/SIM_200m_noCoriolis_sounding_sud")
SIM_1KM = Dataset(path_data + '/SIM_200m_noCoriolis_sounding_sud_dominio2')
SIM_200m = Dataset(path_data + '/SIM_200m_noCoriolis_sounding_sud_dominio3')

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
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIMULAZIONI_DOMINIO_CENTRATO/SIM_500m_noCoriolis_sounding_sud")
SIM_1KM_500m = Dataset(path_data + '/SIM_500m_noCoriolis_sounding_sud_dominio2')
SIM_500m = Dataset(path_data + '/SIM_500m_noCoriolis_sounding_sud_dominio3')

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
plt.rcParams.update({'font.size': 20})

ax = plt.subplot(2, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_1KM_500m[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_1KM_500m[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(560.5, 550, 'Resolution: 1000m', style='italic', fontsize=22,
          bbox={'facecolor': 'grey', 'alpha': 0.5})
#rect = patches.Rectangle((588.5, 560), 10, 10, linewidth=2, edgecolor='r', facecolor='none')
plt.text(529., 648, 'a)', style='italic', fontsize=22)
#plt.arrow(593.5, 570,0.1, 10, width = 0.4,color='r')
#ax.add_patch(rect)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(2, 2, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_500m[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_500m[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(560.5, 550, 'Resolution: 500 m', style='italic', fontsize=22,
         bbox={'facecolor': 'grey', 'alpha': 0.5})
plt.text(529., 648, 'b)', style='italic', fontsize=22)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0., 0.72, 0.015])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=5, y=0.95, rotation=0,fontsize=24)

ax = plt.subplot(2, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_1KM[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_1KM[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(560.5, 550, 'Resolution: 1000 m', style='italic', fontsize=22,
         bbox={'facecolor': 'grey', 'alpha': 0.5})
plt.text(529., 648, 'c)', style='italic', fontsize=22)
plt.hlines(592,525,675, linestyle='dashed',color = 'black')
#plt.hlines(608,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(2, 2, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_200m[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_200m[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(560.5, 550, 'Resolution: 200 m', style='italic', fontsize=22,
         bbox={'facecolor': 'grey', 'alpha': 0.5})
plt.text(529., 648, 'd)', style='italic', fontsize=22)
plt.hlines(592,525,675, linestyle='dashed',color = 'black')
#plt.hlines(608,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.25, 
                    hspace=0.25)

plt.savefig('')

#%% plot control simulation

path_data= os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/CTRL_SIMULATION/")
SIM_CTRL = Dataset(path_data + '/SIM_wrfinput_prova_dominio2')

# plot of QRain hour 5 and 6
z = getvar(SIM_CTRL,'z')
HGT = getvar(SIM_CTRL,'HGT')
QRAIN = getvar(SIM_CTRL,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000 = interplevel(QRAIN,z,2000)
RAINNC = getvar(SIM_CTRL,'RAINNC',timeidx = ALL_TIMES)

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_2000 = QRAIN_2000.assign_coords({'west_east' : QRAIN_2000.west_east.values +X2[0]+1})
QRAIN_2000 = QRAIN_2000.assign_coords({'south_north' : QRAIN_2000.south_north.values +X2[0]+1})
RAINNC = RAINNC.assign_coords({'west_east' : RAINNC.west_east.values +X2[0]+1})
RAINNC = RAINNC.assign_coords({'south_north' : RAINNC.south_north.values +X2[0]+1})


fig = plt.figure(figsize=(15,13))
plt.rcParams.update({'font.size': 22})

ax = plt.subplot(2, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000[5].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000[5].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'a)', style='italic', fontsize=22)
plt.text(613.,555, 'B', style='italic',fontsize=22)
plt.arrow(616, 563,0.1, 11, width = 0.5,color='black')
plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(2, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'b)', style='italic', fontsize=22)

#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.5, 0.72, 0.015])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=5, y=0.95, rotation=0,fontsize=22)

ax = plt.subplot(2, 2, 3)
plot2=(RAINNC[5]-RAINNC[4]).plot.contourf(ax=ax,levels=np.linspace(5,30,11),cmap='BuPu',add_colorbar=False)
#cb = plt.colorbar(plot)
#(RAINNC[5]-RAINNC[4]).plot.contour(ax=ax,levels=10,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'c)', style='italic', fontsize=22)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(2, 2, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
(RAINNC[6]-RAINNC[5]).plot.contourf(ax=ax,levels=np.linspace(5,30,11),cmap='BuPu',add_colorbar=False)
#cb = plt.colorbar(plot)
#QRAIN_2000_200m[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'd)', style='italic', fontsize=22)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.02, 0.72, 0.015])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
cb.set_label('mm/h', labelpad=0, y=0.95, rotation=0,fontsize=22)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.25, 
                    hspace=0.5)

#%%
# prova subplot

plt.rcParams.update({'font.size': 20})
fig, axs = plt.subplots(2,2,figsize=(15, 13))

levels=[0.0004,0.0008,0.0012,0.002,0.003]
plot=axs[0,0].contourf(X2,Y2,QRAIN_2000[5],levels=levels,cmap='Blues')
axs[0,0].contour(X2,Y2,QRAIN_2000[5],levels=levels,colors='black',alpha=0.5)
axs[0,0].contour(X2,Y2,HGT,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
axs[0,0].set_xlim((525,675))
axs[0,0].set_ylim((540,660))
axs[0,0].set_xlabel('x (km)')
axs[0,0].set_ylabel('y (km)')




plot2=axs[0,1].contourf(X2,Y2,QRAIN_2000[6],levels=levels,cmap='Blues')
axs[0,1].contour(X2,Y2,QRAIN_2000[6],levels=levels,colors='black',alpha=0.5)
axs[0,1].contour(X2,Y2,HGT,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
axs[0,1].set_xlim((525,675))
axs[0,1].set_ylim((540,660))
axs[0,1].set_xlabel('x (km)')
axs[0,1].set_ylabel('y (km)')

plot3 = axs[1,0].contourf(X2,Y2,RAINNC[5]-RAINNC[4],levels=np.linspace(5,30,11),cmap = 'BuPu')
axs[1,0].contour(X2,Y2,HGT,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
axs[1,0].set_xlim((525,675))
axs[1,0].set_ylim((540,660))
axs[1,0].set_xlabel('x (km)')
axs[1,0].set_ylabel('y (km)')

axs[1,1].contourf(X2,Y2,RAINNC[6]-RAINNC[5],levels=np.linspace(5,30,11),cmap = 'BuPu')
axs[1,1].contour(X2,Y2,HGT,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
axs[1,1].set_xlim((525,675))
axs[1,1].set_ylim((540,660))
axs[1,1].set_xlabel('x (km)')
axs[1,1].set_ylabel('y (km)')



fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.5, 0.72, 0.015])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=5, y=0.95, rotation=0,fontsize=18)

cbar_ax2 = fig.add_axes([0.13, 0.02, 0.72, 0.015])
cb = fig.colorbar(plot3,cax=cbar_ax2,orientation='horizontal')
cb.set_label('mm/h', labelpad=5, y=0.95, rotation=0,fontsize=18)


plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.25, 
                    hspace=0.5)


#%% plot struttura bande tramite cross section orizzontale-->Roll Vortices

w =getvar(SIM_CTRL,"wa",timeidx=5)
u = getvar(SIM_CTRL, "ua",timeidx=5)  # default is in m/s
QCLOUD = getvar(SIM_CTRL,"QCLOUD",timeidx=5)

start_point = CoordPair(x=90, y=105)
end_point = CoordPair(x=135,y=105)

w_cross = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
u_cross = vertcross(u,z,start_point=start_point,end_point=end_point)
# make a copy of the w cross section
w_cross_filled = np.ma.copy(to_np(w_cross))
u_cross_filled = np.ma.copy(to_np(u_cross))
q_cross_filled = np.ma.copy(to_np(qcloud_cross))
ter_line = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)


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
z = getvar(SIM_CTRL, "z",units="m",timeidx=5)
w =getvar(SIM_CTRL,"wa",timeidx=5)
v = getvar(SIM_CTRL, "va",timeidx=5)  # default is in m/s
#theta = getvar(ncfile,"theta",timeidx=5)
#thetae = getvar(ncfile,"theta_e",timeidx=5)
QCLOUD = getvar(SIM_CTRL,"QCLOUD",timeidx=5)
start_point = CoordPair(x=128,y=50)
end_point = CoordPair(x=128, y=130)
# Get the terrain heights along the cross section line
ter_line_2 = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)
# getting the cross sections of the variables
w_cross_2 = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross_2 = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)*1000
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
plt.rcParams.update({'font.size': 22})
fig, axs = plt.subplots(1,2,figsize=(18, 6))

norm = DivergingNorm(vmin=w_cross_filled.min(),vcenter=0, vmax=w_cross_filled.max())
levs = [10**(-6)]
plot = axs[0].contourf(xs+90,ys,w_cross_filled,norm=norm,levels=30, cmap='RdBu_r')
#cb = fig.colorbar(plot,orientation='horizontal')
#cb.ax.tick_params(labelsize=18)
#cb.set_label('m/s', labelpad=-40, y=1.05, rotation=0,fontsize=18)
axs[0].quiver(xs+90,ys,u_modified/2,w_modified,scale=100.,color='black',alpha=0.5)
axs[0].set_ylim((0,9000))
axs[0].set_xlabel('x [km]',fontsize = 22)
axs[0].set_ylabel('z [m]',fontsize = 22)
axs[0].contour(xs+90,ys,q_cross_filled,levels=levs,linestyles='dashed',linewidths=1.8,colors='black',alpha=0.8)
axs[0].fill_between(xs+90,0,to_np(ter_line),facecolor="black",alpha=0.93)
axs[0].text(578,8200,'e)',style='italic', fontsize=26)


#
levs_w = [0.3,0.5,1,2,4,8,12,15,20,25]
plot2 = axs[1].contourf(xs_2+ Y2[0]+51,ys_2,q_cross_filled_2,levels=np.linspace(0,0.0025*1000,26),cmap='Blues')
axs[1].quiver(Xs_2+Y2[0]+51,Ys_2,v_modified_2/8,w_modified_2,scale=120.,alpha=0.6)
axs[1].fill_between(xs_2+Y2[0]+51,-300,to_np(ter_line_2),facecolor="black")
CS=axs[1].contour(xs_2+Y2[0]+51,ys_2,w_cross_2,levels=levs_w,colors='black',linewidths=0.8)
axs[1].clabel(CS, inline=1, fontsize=18,fmt='%.1f')
axs[1].set_ylim((0,9000))
axs[1].set_xlabel('y [km]',fontsize = 22)
axs[1].set_ylabel('z [m]',fontsize = 22)
axs[1].text(540,8200,'f)',style='italic', fontsize=26)

from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.095, -0.06, 0.365, 0.03])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m/s', labelpad=5, y=0.95, rotation=0,fontsize=22)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.54, -0.06, 0.365, 0.03])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('$10^{-3}$ kg/kg', labelpad=5, y=0.95, rotation=0,fontsize=22)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.25, 
                    hspace=0.5)


#%% PLOT SIM WITH IDEAL WINDS
path_data= os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_WIND/")
path_data2 = os.path.abspath('C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_SHEAR/')
SIM_U10 = Dataset(path_data + '/SIM_U10_dominio2')
SIM_U20 = Dataset(path_data + '/SIM_U20_dominio2')
SIM_U27 = Dataset(path_data + '/SIM_U30_dominio2')
SIM_shear = Dataset(path_data2 + '/SIM_SHEAR_FINALE_dominio2')

z = getvar(SIM_U10,'z')
HGT = getvar(SIM_U10,'HGT')
QRAIN = getvar(SIM_U10,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_U10 = interplevel(QRAIN,z,2000)

z = getvar(SIM_U20,'z')
HGT = getvar(SIM_U20,'HGT')
QRAIN = getvar(SIM_U20,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_U20 = interplevel(QRAIN,z,2000)

z = getvar(SIM_U27,'z')
HGT = getvar(SIM_U27,'HGT')
QRAIN = getvar(SIM_U27,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_U27 = interplevel(QRAIN,z,2000)

z = getvar(SIM_shear,'z')
HGT = getvar(SIM_shear,'HGT')
QRAIN = getvar(SIM_shear,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_shear = interplevel(QRAIN,z,2000)

# calcolo le cross section per SIM_U20
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

# stessa cosa per SIM_U10
z = getvar(SIM_U10, "z",units="m",timeidx=6)
w =getvar(SIM_U10,"wa",timeidx=6)
v = getvar(SIM_U10, "va",timeidx=6)  
thetae = getvar(SIM_U10,"eth",timeidx=6)
QCLOUD = getvar(SIM_U10,"QCLOUD",timeidx=6)
start_point = CoordPair(x=102,y=50)
end_point = CoordPair(x=102, y=150)
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

# making the plot
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_2000_U10 = QRAIN_2000_U10.assign_coords({'west_east' : QRAIN_2000_U10.west_east.values +X2[0]+1})
QRAIN_2000_U10 = QRAIN_2000_U10.assign_coords({'south_north' : QRAIN_2000_U10.south_north.values +X2[0]+1})
QRAIN_2000_U20 = QRAIN_2000_U20.assign_coords({'west_east' : QRAIN_2000_U20.west_east.values +X2[0]+1})
QRAIN_2000_U20 = QRAIN_2000_U20.assign_coords({'south_north' : QRAIN_2000_U20.south_north.values +X2[0]+1})
QRAIN_2000_U27 = QRAIN_2000_U27.assign_coords({'west_east' : QRAIN_2000_U27.west_east.values +X2[0]+1})
QRAIN_2000_U27 = QRAIN_2000_U27.assign_coords({'south_north' : QRAIN_2000_U27.south_north.values +X2[0]+1})
QRAIN_2000_shear = QRAIN_2000_shear.assign_coords({'west_east' : QRAIN_2000_shear.west_east.values +X2[0]+1})
QRAIN_2000_shear = QRAIN_2000_shear.assign_coords({'south_north' : QRAIN_2000_shear.south_north.values +X2[0]+1})

#%%
fig = plt.figure(figsize=(20,24))
plt.rcParams.update({'font.size': 28})

ax = plt.subplot(3, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003]
plot=QRAIN_2000_U10[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_U10[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'a)', style='italic', fontsize=34)
plt.vlines(591,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_V10', y =1.05,fontsize=30)

ax = plt.subplot(3, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003]
plot=QRAIN_2000_U20[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_U20[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_V20', y =1.05,fontsize=30)
plt.text(529., 648, 'b)', style='italic', fontsize=34)
plt.vlines(591,540,660, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(3, 2, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003]
plot=QRAIN_2000_U27[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_U27[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_V30', y =1.05,fontsize=30)
plt.text(529., 648, 'c)', style='italic', fontsize=34)

#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(3, 2, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003]
plot=QRAIN_2000_shear[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_shear[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_V10_SHEAR', y =1.05,fontsize=30)
plt.text(529., 648, 'd)', style='italic', fontsize=34)


#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')


ax = plt.subplot(3, 2, 5)
plot2 = ax.contourf(xs_2+ Y2[0]+51,ys_2,q_cross_filled_2,levels=np.linspace(0,0.0025,26),cmap='Greys')
ax.quiver(Xs_2+Y2[0]+51,Ys_2,v_modified_2/8,w_modified_2,scale=60.,alpha=0.6)
CS=ax.contour(xs_2+Y2[0]+51,ys_2,thetae_cross_2,levels=[322,325,327.8,328.3],colors=['black'])
ax.clabel(CS, inline=3, fontsize=18,fmt='%.1f')
ax.fill_between(xs_2+Y2[0]+51,-300,to_np(ter_line_2),facecolor="black")
ax.set_ylim((0,6000))
ax.set_xlabel('y [km]',fontsize = 24)
ax.set_ylabel('z [m]',fontsize = 24)

ax.text(540,5500,'e)',style='italic', fontsize=34)

ax = plt.subplot(3, 2, 6)
plot2 = ax.contourf(xs+ Y2[0]+51,ys,q_cross_filled,levels=np.linspace(0,0.002,21),cmap='Greys')
ax.quiver(Xs+Y2[0]+51,Ys,v_modified/8,w_modified,scale=110.,alpha=0.6)
CS=ax.contour(xs+Y2[0]+51,ys,thetae_cross,levels=[322,325,327.8,328.3],colors=['black'])
ax.clabel(CS, inline=2, fontsize=18,fmt='%.1f')
ax.fill_between(xs+Y2[0]+51,-300,to_np(ter_line),facecolor="black")
ax.set_ylim((0,6000))
ax.set_xlabel('y [km]',fontsize = 24)
ax.set_ylabel('z [m]',fontsize = 24)
ax.text(540,5500,'f)',style='italic', fontsize=34)



from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.342, 0.72, 0.009])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=0, y=0.92, rotation=0)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.03, 0.72, 0.009])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg/kg', labelpad=0, y=0.95, rotation=0)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.3, 
                    hspace=0.53)

#%% PLOT WAVE PATTERN SIM_V27
SIM_U25 = Dataset(path_data + '/SIM_U30_dominio2')
SIM_U30 = Dataset(path_data + '/SIM_U30_dominio2')


# calcolo le cross section per SIM_U27
w =getvar(SIM_U25,"wa",timeidx=6)
v = getvar(SIM_U25, "va",timeidx=6)  # default is in m/s
QCLOUD = getvar(SIM_U25,"QCLOUD",timeidx=6)
theta_e = getvar(SIM_U25,"eth",timeidx=6)

start_point = CoordPair(x=104, y=50)
end_point = CoordPair(x=104,y=150)

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
w_modified[::3,::3]=w_cross_filled[::3,::3]

v_modified = np.zeros((100,101))
v_modified[:] = np.nan
v_modified[::3,::3]=v_cross_filled[::3,::3]

fig = plt.figure(figsize=(13,19))
plt.rcParams.update({'font.size': 28})
ax = plt.subplot(2, 1, 1)
plot2 = ax.contourf(xs+ Y2[0]+51,ys,w_cross_filled,levels=np.linspace(-10,10,21),cmap='RdBu_r')
cb = plt.colorbar(plot2)
tick_locator = ticker.MaxNLocator(nbins=10)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m/s', labelpad=-60, y=1.08, rotation=0)
ax.quiver(Xs+Y2[0]+51,Ys,v_modified/5,w_modified,scale=190.,alpha=0.6)
CS=ax.contour(xs+Y2[0]+51,ys,thetae_cross,levels=[320,323,324,325,326,327,329,331,335,345,370,400,450],colors=['black'],linewidths=1.1)
ax.clabel(CS, inline=1, fontsize=22,fmt='%.1f')
ax.fill_between(xs+Y2[0]+51,-300,to_np(ter_line),facecolor="black")
ax.set_ylim((0,15000))
ax.set_xlabel('y [km]')
ax.set_ylabel('z [m]')
ax.text(540,14000,'b)',style='italic', fontsize=34)
ax.set_title('SIM_V30',y=1.03,fontsize=28)




#%% PLOT SONDAGGI CON WIND ROTATION
import metpy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import SkewT,Hodograph
from metpy.units import units

SIM_TILTED = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_Wind_direction/Sondaggio_wind_shear_tilted/Parameters_Sondaggio_wind_shear_tilted_FINALE.txt",delimiter = '\t')
SIM_ROTATED = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_Wind_direction/Original_direction_sum20/Parameters_Sondaggio_original_direction_sum20.txt",delimiter='\t')
SIM_SHEAR = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_wind_profile/Sondaggio_U10_double_shear_finale.txt",delimiter='\t')

u=SIM_ROTATED['u']
v=SIM_ROTATED['v']
hght=SIM_ROTATED['Height']
u2 = SIM_TILTED['u']
v2 = SIM_TILTED['v']
hght_2 =SIM_TILTED['Height']
u3 = SIM_SHEAR['u']
v3 = SIM_SHEAR['v']
hght_3 =SIM_SHEAR['Height']
# computing the wnd speed profiles
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
ax1.set_xlabel('u [m/s]',labelpad=0)
ax1.set_ylabel('v [m/s]',labelpad=-5)
ax1.set_title('SIM_UDINE_ROT20',y=1.03,fontsize=18)
ax1.text(-37,33,'a)',style='italic', fontsize=20)
plt.xticks()
plt.yticks()
h.wind_vectors(u[::100],v[::100],scale = 1.,width = 0.4)

#ax2 = plt.subplot(1, 2, 2)
h = Hodograph(ax=ax2,component_range=40.)
h.add_grid(increment=20)
cb = h.plot_colormapped(u2, v2, hght[0:985],linewidth=2.5)
#plt.colorbar(cb,shrink=0.7)
ax2.set_xlabel('u [m/s]',labelpad=0)
ax2.set_ylabel('v [m/s]',labelpad=-5)
ax2.set_title('SIM_SHEAR_TILTED',y=1.03,fontsize=18)
ax2.text(-37,33,'b)',style='italic', fontsize=20)
plt.xticks()
plt.yticks()
h.wind_vectors(u2[::100],v2[::100],scale = 1.,width = 0.4) 

#ax2 = plt.subplot(1, 2, 2)
ax3.plot(speed,hght,color='b',linewidth=2,alpha=0.8,label='SIM_CTRL')
ax3.plot(speed2,hght,color='g',alpha=0.7,linewidth=2,label='SIM_SHEAR_TILTED')
ax3.plot(speed3,hght_3,color='orange',linewidth=2,label='SIM_V10_SHEAR')
ax3.set_ylim((0,13000))
ax3.text(5,11800,'c)',style='italic', fontsize=20)
ax3.legend(bbox_to_anchor=(0.815,0.55),fontsize=14)
ax3.set_xlabel('Wind speed [m/s]')
ax3.set_ylabel('z [m]',labelpad=0)
ax3.grid()
from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.32, -0.22, 0.4, 0.01])
cb = fig.colorbar(cb,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m a.s.l', labelpad=3, y=0.9, rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=-0.15, 
                    right=0.95, 
                    top=0.1, 
                    wspace=0.3, 
                    hspace=0.5)


#%% SUBPLOT QRAIN IDEAL WIND DIRECTION
path_data= os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_WIND_DIRECTION/")
SIM_210 = Dataset(path_data + '/SIM_wind_210degrees_dominio2')
SIM_UDINE = Dataset(path_data + '/SIM_original_wind20degrees_dominio2')
SIM_shear = Dataset(path_data2 + '/SIM_SHEAR_ROT_FINALE_dominio2')

z = getvar(SIM_210,'z')
HGT = getvar(SIM_210,'HGT')
QRAIN = getvar(SIM_210,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_210 = interplevel(QRAIN,z,2000)

z = getvar(SIM_UDINE,'z')
HGT_UDINE = getvar(SIM_UDINE,'HGT')
QRAIN = getvar(SIM_UDINE,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_UDINE = interplevel(QRAIN,z,2000)
ua = getvar(SIM_UDINE,'ua',timeidx=ALL_TIMES)
va = getvar(SIM_UDINE,'va',timeidx=ALL_TIMES)
v_2000 = interplevel(va, z, 2000)
u_2000= interplevel(ua, z, 2000)
v_4000 = interplevel(va, z, 4000)
u_4000= interplevel(ua, z, 4000)
v_5500 = interplevel(va, z, 5500)
u_5500= interplevel(ua, z, 5500)

z = getvar(SIM_shear,'z')
HGT = getvar(SIM_shear,'HGT')
QCLOUD = getvar(SIM_shear,'QCLOUD',timeidx=ALL_TIMES)
qrain =  getvar(SIM_shear,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_shear = interplevel(qrain,z,2000)
QRAIN_2000_shear = interplevel(QCLOUD,z,4000)
#QRAIN_2000_shear = QCLOUD_2000_shear + QRAIN_2000_shear

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_2000_210 = QRAIN_2000_210.assign_coords({'west_east' : QRAIN_2000_210.west_east.values +X2[0]+1})
QRAIN_2000_210 = QRAIN_2000_210.assign_coords({'south_north' : QRAIN_2000_210.south_north.values +X2[0]+1})
QRAIN_2000_UDINE = QRAIN_2000_UDINE.assign_coords({'west_east' : QRAIN_2000_UDINE.west_east.values +X2[0]+1})
QRAIN_2000_UDINE = QRAIN_2000_UDINE.assign_coords({'south_north' : QRAIN_2000_UDINE.south_north.values +X2[0]+1})
QRAIN_2000_shear = QRAIN_2000_shear.assign_coords({'west_east' : QRAIN_2000_shear.west_east.values +X2[0]+1})
QRAIN_2000_shear = QRAIN_2000_shear.assign_coords({'south_north' : QRAIN_2000_shear.south_north.values +X2[0]+1})
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
# making the plot
fig = plt.figure(figsize=(14,18))
plt.rcParams.update({'font.size': 24})

ax = plt.subplot(3, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003]
plot=QRAIN_2000_210[3].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_210[3].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'a)', style='italic', fontsize=26)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_210°', y =1.05,fontsize=24)
plt.arrow(545., 550,5, 10, width = 0.7,color='black')
plt.arrow(548., 550,5, 10, width = 0.7,color='g')

ax = plt.subplot(3, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003]
plot=QRAIN_2000_210[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_210[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_210°', y =1.05,fontsize=24)
plt.text(529., 648, 'b)', style='italic', fontsize=26)
plt.arrow(545., 550,5, 10, width = 0.7,color='black')
plt.arrow(548., 550,5, 10, width = 0.7,color='g')

#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(3, 2, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003]
plot=QRAIN_2000_UDINE[3].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_UDINE[3].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_UDINE_ROT20', y =1.05,fontsize=24)
plt.text(529., 648, 'c)', style='italic', fontsize=26)
ax.quiver(X2[50],Y2[60],to_np(u_2000[3,50,60]),to_np(v_2000[3,50,60]),scale_units='xy', scale=1.2,color='orange',width=0.008,label = 'wind 2000m')
ax.quiver(X2[50],Y2[60],to_np(u_4000[3,50,60]),to_np(v_4000[3,50,60]),scale_units='xy', scale=1.2,color='red',width=0.008,label = 'wind 2000m')
ax.quiver(X2[50],Y2[60],to_np(u_5500[3,50,60]),to_np(v_5500[3,50,60]),scale_units='xy', scale=1.2,color='blue',width=0.008,label = 'wind 2000m')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(3, 2, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003]
plot=QRAIN_2000_UDINE[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_UDINE[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.quiver(X2[50],Y2[60],to_np(u_2000[6,50,60]),to_np(v_2000[6,50,60]),scale_units='xy', scale=1.2,color='orange',width=0.008,label = 'wind 2000m')
ax.quiver(X2[50],Y2[60],to_np(u_4000[6,50,60]),to_np(v_4000[6,50,60]),scale_units='xy', scale=1.2,color='red',width=0.008,label = 'wind 2000m')
ax.quiver(X2[50],Y2[60],to_np(u_5500[6,50,60]),to_np(v_5500[6,50,60]),scale_units='xy', scale=1.2,color='blue',width=0.008,label = 'wind 2000m')
plt.title('SIM_UDINE_ROT20', y =1.05,fontsize=24)
plt.text(529., 648, 'd)', style='italic', fontsize=26)


#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(3, 2, 5)
levels=[0.0006,0.001,0.0015,0.002,0.003]
plot2=QRAIN_2000_shear[4].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_shear[4].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_SHEAR_TILTED', y =1.05,fontsize=24)
plt.text(529., 648, 'e)', style='italic', fontsize=26)
plt.arrow(545., 550,10, 10, width = 0.7,color='black')
plt.arrow(545., 550,1, 10, width = 0.7,color='g')

#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(3, 2, 6)
levels=[0.0006,0.001,0.0015,0.002,0.003]
plot2=QRAIN_2000_shear[5].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_2000_shear[5].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_SHEAR_TILTED', y =1.05,fontsize=24)
plt.text(529., 648, 'f)', style='italic', fontsize=26)
plt.arrow(545., 550,10, 10, width = 0.7,color='black')
plt.arrow(545., 550,1, 10, width = 0.7,color='g')


#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')



from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.14, 0.344, 0.72, 0.0075])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=0, y=0.92, rotation=0,fontsize=24)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.14, 0.03, 0.72, 0.0075])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg/kg', labelpad=0, y=0.95, rotation=0,fontsize=24)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.37, 
                    hspace=0.63)

#%% PLOT PROFILI VERTICALI TEMPERATURA POTENZIALE EQUIVALENTE

# profile on stability
N1_000001 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_N_profile/Stable_layer1/Sondaggio_N1_00001/Parameters_Sondaggio_N1_00001.txt",delimiter = '\t')
N1_00015 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_N_profile/Stable_layer1/Sondaggio_N1_00015/Parameters_Sondaggio_N1_00015.txt",delimiter = '\t')
CTRL = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_N_profile/CTRL_SIMULATION/Parameters_CTRL_sounding.txt",delimiter = '\t')
N3_00004 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_N_profile/Modified_layer3/Sondaggio_N1_00004_N3_00004/Parameters_Sondaggio_N3_00004.txt",delimiter = '\t')
N1_00004 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_N_profile/Stable_layer1/Sondaggio_N1_00004/Parameters_Sondaggio_N1_00004.txt",delimiter = '\t')
T300 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_N_profile/SIM_T300/Parameters_sounding_T300.txt", delimiter='\t')
T290 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_N_profile/SIM_T290/Parameters_sounding_T290.txt", delimiter='\t')
N3_00012 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_N_profile/Modified_layer3/Sondaggio_N1_00004_N3_00012/Parameters_Sondaggio_N3_00012.txt",delimiter='\t')
N3_0001 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_N_profile/Modified_layer3/Parameters_sondaggio_N3_0001.txt",delimiter='\t')

thetae_N1_00004 = N1_00004['theta_e']
z = N1_000001['Height']
thetae_CTRL = CTRL['theta_e']
z_CTRL = CTRL['Height']
thetae_N1_00015 = N1_00015['theta_e']
thetae_N3_00004 = N3_00004['theta_e']
thetae_T290 = T290['theta_e']
thetae_T300 = T300['theta_e']
thetae_N1_000001 = N1_000001['theta_e']
thetae_N3_00012 = N3_00012['theta_e']
thetae_N3_0001 = N3_0001['theta_e']

fig = plt.figure(figsize=(20,9))
plt.rcParams.update({'font.size': 28})

ax = plt.subplot(1, 2, 1)
# ax.plot(thetae_CTRL,z_CTRL,linewidth=4, label = 'CTRL',color='green',alpha=0.7,linestyle='dashed')
# ax.plot(thetae_T290,z,linewidth=4,label='$T_s=290 K$',color='blue',alpha=0.3)
# ax.plot(thetae_T300,z,linewidth=4,label = '$T_s=300K$',color='red',alpha=0.7)
# ax.plot(thetae_N1_00015,z,linewidth=3,label = '$N_{1}^{2} = 0.00015 \ s^{-2}$',color='red',alpha=0.5)
# ax.plot(thetae_N1_000001,z,linewidth=3,label='$N_{1}^{2} = 0.000001 \ s^{-2}$',color='blue',alpha=0.6)
# ax.plot(thetae_N1_00004,z,linewidth=3,label='$N_{1}^{2} = 0.00004 \ s^{-2}$',alpha=0.9)
ax.plot(thetae_CTRL,z_CTRL,linewidth=4, label = 'SIM_CTRL',color='green',alpha=0.7,linestyle='dashed')
ax.plot(thetae_T290,z,linewidth=4,label='SIM_T290',color='blue',alpha=0.3)
ax.plot(thetae_T300,z,linewidth=4,label = 'SIM_T300',color='red',alpha=0.7)
ax.plot(thetae_N1_00015,z,linewidth=3,label = 'SIM_N1_00015',color='red',alpha=0.5)
ax.plot(thetae_N1_000001,z,linewidth=3,label='SIM_N1_000001',color='blue',alpha=0.6)
ax.plot(thetae_N1_00004,z,linewidth=3,label='SIM_N1_00004',alpha=0.9)



plt.legend(fontsize=20)
plt.xlabel(r'$\theta_e$ [K]')
plt.ylabel('Heigth [m]')
plt.ylim((50,6000))
plt.xlim((305,355))
plt.grid()
plt.text(306., 5600, 'a)', style='italic', fontsize=34)

ax = plt.subplot(1,2,2)
# ax.plot(thetae_CTRL,z_CTRL,linewidth=3, label = 'CTRL',color='green',alpha=0.7,linestyle='dashed')
# ax.plot(thetae_N3_00004,z,linewidth=4,label = '$N_{3}^{2} = 0.00004 \ s^{-2}$',color='blue',alpha=0.5)
# ax.plot(thetae_N3_00012,z,linewidth=3,label = '$N_{3}^{2} = 0.00012 \ s^{-2}$',color ='red',alpha=0.5)
# ax.plot(thetae_N1_00004,z,linewidth=3,label='$N_{3}^{2} = 0.00008 \ s^{-2}$',alpha=0.9)
ax.plot(thetae_CTRL,z_CTRL,linewidth=3, label = 'SIM_CTRL',color='green',alpha=0.7,linestyle='dashed')
ax.plot(thetae_N3_00004,z,linewidth=4,label = 'SIM_N1_00004_N3_00004',color='blue',alpha=0.5)
ax.plot(thetae_N3_00012,z,linewidth=3,label = 'SIM_N1_00004_N3_00012',color ='red',alpha=0.5)
ax.plot(thetae_N1_00004,z,linewidth=3,label='SIM_N1_00004',alpha=0.9)
#ax.plot(thetae_N3_0001,z,linewidth=4,label = '$N_{3}^{2} = 0.0001 \ s^{-2}$',color ='blue',alpha=0.5)
plt.legend(fontsize=20,loc=1,bbox_to_anchor=(0.999,0.87))
plt.ylim((50,6000))
plt.xlim((305,345))
plt.xlabel(r'$\theta_e$ [K]')
plt.ylabel('Heigth [m]')
plt.grid()
plt.text(306., 5600, 'b)', style='italic', fontsize=34)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.28, 
                    hspace=0.6)

#%% PLOT QRAIN SENSITIVITà ALLA TEMPERATURA SUPERFICIALE
path_data= os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/CTRL_SIMULATION/")
SIM_CTRL = Dataset(path_data + '/SIM_wrfinput_prova_dominio2')
path_data = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_T_QSNOW/")
SIM_T290 = Dataset(path_data+'/SIM_T290_dominio2')
SIM_T300 = Dataset(path_data+'/SIM_T300_dominio2')

QRAIN_CTRL = getvar(SIM_CTRL,'QRAIN',timeidx=ALL_TIMES)
HGT = getvar(SIM_CTRL,'HGT')
z = getvar(SIM_CTRL,'z')
QRAIN_T300 = getvar(SIM_T300,'QRAIN',timeidx=ALL_TIMES)
QSNOW_T290 = getvar(SIM_T290,'QSNOW',timeidx=ALL_TIMES)

QRAIN_CTRL_2000 = interplevel(QRAIN_CTRL,z,2000)
QRAIN_T300_2000 = interplevel(QRAIN_T300,z,2000)
QSNOW_T290_3000 = interplevel(QSNOW_T290,z,3000)

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_CTRL_2000 = QRAIN_CTRL_2000.assign_coords({'west_east' : QRAIN_CTRL_2000.west_east.values +X2[0]+1})
QRAIN_CTRL_2000 = QRAIN_CTRL_2000.assign_coords({'south_north' : QRAIN_CTRL_2000.south_north.values +X2[0]+1})
QRAIN_T300_2000 = QRAIN_T300_2000.assign_coords({'west_east' : QRAIN_T300_2000.west_east.values +X2[0]+1})
QRAIN_T300_2000 = QRAIN_T300_2000.assign_coords({'south_north' : QRAIN_T300_2000.south_north.values +X2[0]+1})
QSNOW_T290_3000 = QSNOW_T290_3000.assign_coords({'west_east' : QSNOW_T290_3000.west_east.values +X2[0]+1})
QSNOW_T290_3000 = QSNOW_T290_3000.assign_coords({'south_north' : QSNOW_T290_3000.south_north.values +X2[0]+1})

#%%
# making the plot
fig = plt.figure(figsize=(20,4))
plt.rcParams.update({'font.size': 24})

ax = plt.subplot(1, 3, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QSNOW_T290_3000[8].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QSNOW_T290_3000[8].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.3)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'a)', style='italic', fontsize=26)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_290', y =1.05,fontsize=24)
#plt.hlines(597,525,675, linestyle='dashed',color = 'black')
plt.vlines(600,525,675, linestyle='dashed',color='black')

ax = plt.subplot(1, 3, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_CTRL_2000[5].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_CTRL_2000[5].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.3)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_CTRL', y =1.05,fontsize=24)
plt.text(529., 648, 'b)', style='italic', fontsize=26)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
#plt.hlines(597,525,675, linestyle='dashed',color = 'black')
plt.vlines(600,525,675, linestyle='dashed',color='black')

ax = plt.subplot(1, 3, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_T300_2000[5].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_T300_2000[5].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.3)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_T300', y =1.05,fontsize=24)
plt.text(529., 648, 'c)', style='italic', fontsize=26)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
#plt.hlines(597,525,675, linestyle='dashed',color = 'black')
plt.vlines(600,525,675, linestyle='dashed',color='black')


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.14, -0.05, 0.72, 0.032])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg/kg', labelpad=0, y=0.95, rotation=0,fontsize=24)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.32, 
                    hspace=0.6)

#%% provo a fare il plot diviso in 2 righe
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(12,10))
plt.rcParams.update({'font.size': 20})
gs = gridspec.GridSpec(4, 4)

ax1 = plt.subplot(gs[:2, :2])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QSNOW_T290_3000[8].plot.contourf(ax=ax1,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QSNOW_T290_3000[8].plot.contour(ax=ax1,levels=levels,colors='black',alpha=0.3)
HGT.plot.contour(ax=ax1,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'a)', style='italic', fontsize=22)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_T290', y =1.05,fontsize=20)
plt.hlines(597,525,675, linestyle='dashed',color = 'black')

ax2 = plt.subplot(gs[:2, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_CTRL_2000[5].plot.contourf(ax=ax2,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_CTRL_2000[5].plot.contour(ax=ax2,levels=levels,colors='black',alpha=0.3)
HGT.plot.contour(ax=ax2,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_CTRL', y =1.05,fontsize=20)
plt.text(529., 648, 'b)', style='italic', fontsize=22)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.hlines(597,525,675, linestyle='dashed',color = 'black')

ax3 = plt.subplot(gs[2:4, 1:3])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_T300_2000[5].plot.contourf(ax=ax3,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_T300_2000[5].plot.contour(ax=ax3,levels=levels,colors='black',alpha=0.3)
HGT.plot.contour(ax=ax3,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_T300', y =1.05,fontsize=20)
plt.text(529., 648, 'c)', style='italic', fontsize=22)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.hlines(597,525,675, linestyle='dashed',color = 'black')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.14, 0.095, 0.72, 0.017])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg/kg', labelpad=0, y=0.95, rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.82, 
                    hspace=1.5)

#%% PLOT OF MOIST BRUNT VAISALA FREQUENCY
# EXTRACTING VARIABLES NEEDED ALONG CROSS SECTION FOR CTRL_SIM
Theta = getvar(SIM_CTRL,"theta",timeidx=1)
T = getvar(SIM_CTRL,"tk",timeidx=1)
z = getvar(SIM_CTRL,"z", timeidx=1)
RH = getvar(SIM_CTRL,"rh",timeidx=1)
QVAPOR = getvar(SIM_CTRL,"QVAPOR",timeidx=1)
QRAIN = getvar(SIM_CTRL,"QRAIN",timeidx=1)
QCLOUD = getvar(SIM_CTRL,"QCLOUD",timeidx=1)
QGRAUP = getvar(SIM_CTRL,"QGRAUP",timeidx=1)
QICE = getvar(SIM_CTRL,"QICE",timeidx=1)
#QSNOW = getvar(SIM_CTRL,"QSNOW",timeidx=1)
Theta_e = getvar(SIM_CTRL,"eth",timeidx=1)
# EXTRACTING VARIABLES NEEDED ALONG CROSS SECTION FOR SIM_T290
Theta_T290 = getvar(SIM_T290,"theta",timeidx=1)
T_T290 = getvar(SIM_T290,"tk",timeidx=1)
z_T290 = getvar(SIM_T290,"z", timeidx=1)
RH_T290 = getvar(SIM_T290,"rh",timeidx=1)
QVAPOR_T290 = getvar(SIM_T290,"QVAPOR",timeidx=1)
QRAIN_T290 = getvar(SIM_T290,"QRAIN",timeidx=1)
QCLOUD_T290 = getvar(SIM_T290,"QCLOUD",timeidx=1)
QGRAUP_T290 = getvar(SIM_T290,"QGRAUP",timeidx=1)
QICE_T290 = getvar(SIM_T290,"QICE",timeidx=1)
QSNOW_T290 = getvar(SIM_T290,"QSNOW",timeidx=1)
Theta_e_T290 = getvar(SIM_T290,"eth",timeidx=1)
# EXTRACTING VARIABLES NEEDED ALONG CROSS SECTION FOR SIM_T300
Theta_T300 = getvar(SIM_T300,"theta",timeidx=1)
T_T300 = getvar(SIM_T300,"tk",timeidx=1)
z_T300 = getvar(SIM_T300,"z", timeidx=1)
RH_T300 = getvar(SIM_T300,"rh",timeidx=1)
QVAPOR_T300 = getvar(SIM_T300,"QVAPOR",timeidx=1)
QRAIN_T300 = getvar(SIM_T300,"QRAIN",timeidx=1)
QCLOUD_T300 = getvar(SIM_T300,"QCLOUD",timeidx=1)
QGRAUP_T300 = getvar(SIM_T300,"QGRAUP",timeidx=1)
QICE_T300 = getvar(SIM_T300,"QICE",timeidx=1)
QSNOW_T300 = getvar(SIM_T300,"QSNOW",timeidx=1)
Theta_e_T300 = getvar(SIM_T300,"eth",timeidx=1)

g = 9.81
L= 2501000 # J/kg latent heat of vaporization
R = 287 # dry air gas constant
Rw = 461.5
eps = 287/461.5 # ratio between gas constant of dry and moist air
cp = 1004
q_s = 100*QVAPOR/RH
q_s_290 = 100*QVAPOR_T290/RH_T290
q_s_300 = 100*QVAPOR_T300/RH_T300

# applichiamola sulla cross section
start_point = CoordPair(x=113,y=50)
end_point = CoordPair(x=113, y=125)

ter_line = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables fro SIM_CTRL
qcloud_cross = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
thetae_cross = vertcross(Theta_e,z,start_point=start_point,end_point=end_point)
theta_cross = vertcross(Theta,z,start_point=start_point,end_point=end_point)
qvapor_cross = vertcross(QVAPOR,z,start_point=start_point,end_point=end_point)
RH_cross = vertcross(RH,z,start_point=start_point,end_point=end_point)
T_cross = vertcross(T,z,start_point=start_point,end_point=end_point)
qrain_cross = vertcross(QRAIN,z,start_point=start_point,end_point=end_point)
qgraup_cross = vertcross(QGRAUP,z,start_point=start_point,end_point=end_point)
qice_cross = vertcross(QICE,z,start_point=start_point,end_point=end_point)
#qsnow_cross = vertcross(QSNOW,z,start_point=start_point,end_point=end_point)

q_s_cross = 100*qvapor_cross/RH_cross
q_w_cross = qrain_cross + q_s_cross + qcloud_cross + qgraup_cross + qice_cross #+ qsnow_cross

# getting the cross sections of the variables for SIM_T290
qcloud_cross_T290 = vertcross(QCLOUD_T290,z,start_point=start_point,end_point=end_point)
thetae_cross_T290 = vertcross(Theta_e_T290,z,start_point=start_point,end_point=end_point)
theta_cross_T290 = vertcross(Theta_T290,z,start_point=start_point,end_point=end_point)
qvapor_cross_T290 = vertcross(QVAPOR_T290,z,start_point=start_point,end_point=end_point)
RH_cross_T290 = vertcross(RH_T290,z,start_point=start_point,end_point=end_point)
T_cross_T290 = vertcross(T_T290,z,start_point=start_point,end_point=end_point)
qrain_cross_T290 = vertcross(QRAIN_T290,z,start_point=start_point,end_point=end_point)
qgraup_cross_T290 = vertcross(QGRAUP_T290,z,start_point=start_point,end_point=end_point)
qice_cross_T290 = vertcross(QICE_T290,z,start_point=start_point,end_point=end_point)
qsnow_cross_T290 = vertcross(QSNOW_T290,z,start_point=start_point,end_point=end_point)

q_s_cross_T290 = 100*qvapor_cross_T290/RH_cross_T290
q_w_cross_T290 = qrain_cross_T290 + q_s_cross_T290 + qcloud_cross_T290 + qgraup_cross_T290 + qice_cross_T290 + qsnow_cross_T290

# getting the cross sections of the variables for SIM_T300
qcloud_cross_T300 = vertcross(QCLOUD_T300,z,start_point=start_point,end_point=end_point)
thetae_cross_T300 = vertcross(Theta_e_T300,z,start_point=start_point,end_point=end_point)
theta_cross_T300 = vertcross(Theta_T300,z,start_point=start_point,end_point=end_point)
qvapor_cross_T300 = vertcross(QVAPOR_T300,z,start_point=start_point,end_point=end_point)
RH_cross_T300 = vertcross(RH_T300,z,start_point=start_point,end_point=end_point)
T_cross_T300 = vertcross(T_T300,z,start_point=start_point,end_point=end_point)
qrain_cross_T300 = vertcross(QRAIN_T300,z,start_point=start_point,end_point=end_point)
qgraup_cross_T300 = vertcross(QGRAUP_T300,z,start_point=start_point,end_point=end_point)
qice_cross_T300 = vertcross(QICE_T300,z,start_point=start_point,end_point=end_point)
qsnow_cross_T300 = vertcross(QSNOW_T300,z,start_point=start_point,end_point=end_point)

q_s_cross_T300 = 100*qvapor_cross_T300/RH_cross_T300
q_w_cross_T300 = qrain_cross_T300 + q_s_cross_T300 + qcloud_cross_T300 + qgraup_cross_T300 + qice_cross_T300 + qsnow_cross_T300

xs = np.arange(0, thetae_cross.shape[-1], 1)
ys = to_np(thetae_cross.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

dz = np.diff(Ys,axis=0)
dq_s_dz = q_s_cross.differentiate('vertical')
dq_s_dz_T290 = q_s_cross_T290.differentiate('vertical')
dq_s_dz_T300 = q_s_cross_T300.differentiate('vertical')

log_Theta = np.log(theta_cross)
log_Theta_T290 = np.log(theta_cross_T290)
log_Theta_T300 = np.log(theta_cross_T300)
dTheta_dz = log_Theta.differentiate('vertical')
dTheta_dz_T290 = log_Theta_T290.differentiate('vertical')
dTheta_dz_T300 = log_Theta_T290.differentiate('vertical')

dq_w_dz = q_w_cross.differentiate('vertical')
dq_w_dz_T290 = q_w_cross_T290.differentiate('vertical')
dq_w_dz_T300 = q_w_cross_T300.differentiate('vertical')

Nm_squared = g*((1+(L*q_s_cross/(R*T_cross)))/(1+((eps*q_s_cross*L**2)/(cp*R*T_cross**2))))*(dTheta_dz+L/(cp*T_cross)*(dq_s_dz)) - g*dq_w_dz
Nm_squared_T290 = g*((1+(L*q_s_cross_T290/(R*T_cross_T290)))/(1+((eps*q_s_cross_T290*L**2)/(cp*R*T_cross_T290**2))))*(dTheta_dz_T290+L/(cp*T_cross_T290)*(dq_s_dz_T290)) - g*dq_w_dz_T290
Nm_squared_T300 = g*((1+(L*q_s_cross_T300/(R*T_cross_T300)))/(1+((eps*q_s_cross_T300*L**2)/(cp*R*T_cross_T300**2))))*(dTheta_dz_T300+L/(cp*T_cross_T300)*(dq_s_dz_T300)) - g*dq_w_dz_T300

# interpolo i dati su una griglia motlo fitta per vedere meglio il plot
vertical = Nm_squared.vertical.values
horizontal = Nm_squared.cross_line_idx.values
NM_squared = xr.DataArray(Nm_squared, coords=[vertical, horizontal])
NM_squared_T290 = xr.DataArray(Nm_squared_T290, coords=[vertical, horizontal])
NM_squared_T300 = xr.DataArray(Nm_squared_T300, coords=[vertical, horizontal])
TER_LINE = xr.DataArray(ter_line, coords=[horizontal])

z_coord = np.linspace(0,NM_squared.vertical.max(),1000)
y_coord= np.linspace(NM_squared.cross_line_idx.min(),NM_squared.cross_line_idx.max(),1000)
NM_squared_interp = NM_squared.interp(vertical = z_coord,cross_line_idx = y_coord)
NM_squared_interp_T290 = NM_squared_T290.interp(vertical = z_coord,cross_line_idx = y_coord)
NM_squared_interp_T300 = NM_squared_T300.interp(vertical = z_coord,cross_line_idx = y_coord)
ter_line_interp = TER_LINE.interp(line_idx = y_coord)

QCROSS = xr.DataArray(qcloud_cross, coords=[vertical, horizontal])
QCROSS_T290 = xr.DataArray(qcloud_cross_T290, coords=[vertical, horizontal])
QCROSS_T300 = xr.DataArray(qcloud_cross_T300, coords=[vertical, horizontal])

QCROSS = QCROSS.interp(vertical = z_coord,cross_line_idx = y_coord)
QCROSS_T290 = QCROSS_T290.interp(vertical = z_coord,cross_line_idx = y_coord)
QCROSS_T300 = QCROSS_T300.interp(vertical = z_coord,cross_line_idx = y_coord)


xs = to_np(NM_squared_interp.cross_line_idx)
ys = to_np(NM_squared_interp.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

# metto nan fuori dalla cloud per plottare solo dove ho saturazione            
NM_squared_interp = NM_squared_interp.where(QCROSS>10**(-5))
NM_squared_interp_T290 = NM_squared_interp_T290.where(QCROSS_T290>10**(-5))
NM_squared_interp_T300 = NM_squared_interp_T300.where(QCROSS_T300>10**(-5))

NM_squared_interp_T300 = NM_squared_interp_T300.assign_coords({'cross_line_idx' : NM_squared_interp_T300.cross_line_idx.values+Y2[0]+50})
QCROSS_T300 = QCROSS_T300.assign_coords({'cross_line_idx' : QCROSS_T300.cross_line_idx.values +Y2[0]+50})
NM_squared_interp = NM_squared_interp.assign_coords({'cross_line_idx' : NM_squared_interp.cross_line_idx.values+Y2[0]+50})
QCROSS = QCROSS.assign_coords({'cross_line_idx' : QCROSS.cross_line_idx.values +Y2[0]+50})
NM_squared_interp_T290 = NM_squared_interp_T290.assign_coords({'cross_line_idx' : NM_squared_interp_T290.cross_line_idx.values+Y2[0]+50})
QCROSS_T290 = QCROSS_T290.assign_coords({'cross_line_idx' : QCROSS_T290.cross_line_idx.values +Y2[0]+50})


#NM_squared_interp_T300 = NM_squared_interp_T300.assign_coords({'south_north' : NM_squared_interp_T300.south_north.values +X2[0]+1})
#%% making the subplot
fig = plt.figure(figsize=(20,5))
plt.rcParams.update({'font.size': 24})

ax = plt.subplot(1, 3, 1)
plot=(NM_squared_interp_T290*10000).plot(levels=np.linspace(-2,2,41),cmap='RdBu_r',add_colorbar=False)
ax.fill_between(xs+Y2[0]+50,0,to_np(ter_line_interp),facecolor="black")
ax.contour(Xs+Y2[0]+50,Ys,QCROSS_T290,levels=[10**(-5)],colors='Gray')
plt.ylabel('z [m]')
plt.xlabel('y [km]')
plt.title('T$_s$ = 290 K',fontsize=26,y=1.03)
plt.ylim((0,7000))
plt.text(538., 6300, 'a)', style='italic', fontsize=26)

ax = plt.subplot(1, 3, 2)
plot=(NM_squared_interp*10000).plot(levels=np.linspace(-2,2,41),cmap='RdBu_r',add_colorbar=False)
ax.fill_between(xs+Y2[0]+50,0,to_np(ter_line_interp),facecolor="black")
ax.contour(Xs+Y2[0]+50,Ys,QCROSS,levels=[10**(-5)],colors='Gray')
plt.title('Control Simulation',fontsize=26,y=1.03)
plt.ylabel('z [m]')
plt.xlabel('y [km]')
plt.ylim((0,7000))
plt.text(538., 6300, 'b)', style='italic', fontsize=26)

ax = plt.subplot(1, 3, 3)
plot=(NM_squared_interp_T300*10000).plot(levels=np.linspace(-2,2,41),cmap='RdBu_r',add_colorbar=False)
ax.fill_between(xs+Y2[0]+50,0,to_np(ter_line_interp),facecolor="black")
ax.contour(Xs+Y2[0]+50,Ys,QCROSS_T300,levels=[10**(-5)],colors='Gray')
plt.ylabel('z [m]')
plt.xlabel('y [km]')
plt.title('T$_s$ = 300 K',fontsize=26,y=1.03)
plt.ylim((0,7000))
plt.text(538., 6300, 'c)', style='italic', fontsize=26)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.14, -0.01, 0.72, 0.032])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
cb.update_ticks()
cb.set_label('10$^{-4}\  s^{-2}$', labelpad=0, y=0.95, x=0.53, rotation=0,fontsize=24)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.32, 
                    hspace=0.3)

#%% plot stability layer 1
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_N/SIM/")
SIM_CTRL = Dataset(path_data + '/SIM_wrfinput_prova_dominio2')
SIM_N1_00004 = Dataset(path_data + '/SIM_N1_00004_dominio2')
SIM_N1_00015 = Dataset(path_data + '/SIM_N1_00015_dominio2')
SIM_N1_000001 = Dataset(path_data + '/SIM_N1_000001_dominio2')

z = getvar(SIM_CTRL,'z')
HGT = getvar(SIM_CTRL,'HGT')
QRAIN_CTRL = getvar(SIM_CTRL,'QRAIN',timeidx=5)
QRAIN_N1_00004 = getvar(SIM_N1_00004,'QRAIN',timeidx=5)
QRAIN_N1_00015 = getvar(SIM_N1_00015,'QRAIN',timeidx=5)
QRAIN_N1_000001 = getvar(SIM_N1_000001,'QRAIN',timeidx=2)

QRAIN_CTRL_2000 = interplevel(QRAIN_CTRL,z,2000)
QRAIN_N1_00004_2000 = interplevel(QRAIN_N1_00004,z,2000)
QRAIN_N1_00015_2000 = interplevel(QRAIN_N1_00015,z,2000)
QRAIN_N1_000001_2000 = interplevel(QRAIN_N1_000001,z,2000)

# changing the coordinates
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_CTRL_2000 = QRAIN_CTRL_2000.assign_coords({'west_east' : QRAIN_CTRL_2000.west_east.values +X2[0]+1})
QRAIN_CTRL_2000 = QRAIN_CTRL_2000.assign_coords({'south_north' : QRAIN_CTRL_2000.south_north.values +X2[0]+1})
QRAIN_N1_00004_2000 = QRAIN_N1_00004_2000.assign_coords({'west_east' : QRAIN_N1_00004_2000.west_east.values +X2[0]+1})
QRAIN_N1_00004_2000 = QRAIN_N1_00004_2000.assign_coords({'south_north' : QRAIN_N1_00004_2000.south_north.values +X2[0]+1})
QRAIN_N1_00015_2000 = QRAIN_N1_00015_2000.assign_coords({'west_east' : QRAIN_N1_00015_2000.west_east.values +X2[0]+1})
QRAIN_N1_00015_2000 = QRAIN_N1_00015_2000.assign_coords({'south_north' : QRAIN_N1_00015_2000.south_north.values +X2[0]+1})
QRAIN_N1_000001_2000 = QRAIN_N1_000001_2000.assign_coords({'west_east' : QRAIN_N1_000001_2000.west_east.values +X2[0]+1})
QRAIN_N1_000001_2000 = QRAIN_N1_000001_2000.assign_coords({'south_north' : QRAIN_N1_000001_2000.south_north.values +X2[0]+1})



# making the subplots
fig = plt.figure(figsize=(15,13))
plt.rcParams.update({'font.size': 26})
ax = plt.subplot(2, 2, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_CTRL_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_CTRL_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(535, 655, 'c)', style='italic',fontsize=22)
plt.title('SIM_CTRL',fontsize=28,y=1.02)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(2, 2, 2)
plot=QRAIN_N1_00004_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_N1_00004_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(535, 655, 'b)', style='italic',fontsize=22)
plt.title('SIM_N1_00004',fontsize=28,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))

ax = plt.subplot(2, 2, 4)
plot=QRAIN_N1_00015_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_N1_00015_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(535, 655, 'd)', style='italic',fontsize=22)
plt.title('SIM_N1_00015',fontsize=28,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))

ax = plt.subplot(2, 2, 1)
plot=QRAIN_N1_000001_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_N1_000001_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(535, 655, 'a)', style='italic',fontsize=22)
plt.title('SIM_N1_000001',fontsize=28,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.09, 0.5, 0.015])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=5, y=0.95, rotation=0,fontsize=26)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.27, 
                    hspace=0.4)

#%% provo a fare il plot diviso in 2 righe
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(15,10))
plt.rcParams.update({'font.size': 20})
gs = gridspec.GridSpec(4, 4)

ax1 = plt.subplot(gs[:2, :2])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N1_000001_2000.plot.contourf(ax=ax1, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_N1_000001_2000.plot.contour(ax=ax1, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax1,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax1.text(530, 645, 'a)', style='italic',fontsize=22)
plt.title('SIM_N1_000001',fontsize=24,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))

'''
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QSNOW_T290_3000[8].plot.contourf(ax=ax1,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QSNOW_T290_3000[8].plot.contour(ax=ax1,levels=levels,colors='black',alpha=0.3)
HGT.plot.contour(ax=ax1,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'a)', style='italic', fontsize=22)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_T290', y =1.05,fontsize=20)
plt.hlines(597,525,675, linestyle='dashed',color = 'black')
'''

ax2 = plt.subplot(gs[:2, 2:])
plot=QRAIN_N1_00004_2000.plot.contourf(ax=ax2, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_N1_00004_2000.plot.contour(ax=ax2, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax2,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax2.text(530, 645, 'b)', style='italic',fontsize=22)
plt.title('SIM_N1_00004',fontsize=24,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))


ax3 = plt.subplot(gs[2:4, 1:3])
plot=QRAIN_N1_00015_2000.plot.contourf(ax=ax3, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_N1_00015_2000.plot.contour(ax=ax3, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax3,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax3.text(530, 645, 'c)', style='italic',fontsize=22)
plt.title('SIM_N1_00015',fontsize=24,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))


#fig.subplots_adjust(right=0.8)
#cbar_ax = fig.add_axes([0.14, 0.095, 0.72, 0.017])
#cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
#cb.update_ticks()
#cb.set_label('kg/kg', labelpad=0, y=0.95, rotation=0,fontsize=20)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0., 0.72, 0.015])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=5, y=0.95, rotation=0,fontsize=24)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.65, 
                    hspace=1.15)



#%% faccio subplot con roll vortices

HGT = getvar(SIM_N1_000001,"HGT")
QCLOUD_N1_000001 = getvar(SIM_N1_000001,"QCLOUD",timeidx=2)
z = getvar(SIM_N1_000001, "z",units="m")
w_N1_000001 =getvar(SIM_N1_000001,"wa",timeidx=2)
u_N1_000001 = getvar(SIM_N1_000001, "ua",timeidx=2) 
v_N1_000001 = getvar(SIM_N1_000001, "va",timeidx=2) 
 

HGT = getvar(SIM_CTRL,"HGT")
QCLOUD_CTRL = getvar(SIM_CTRL,"QCLOUD",timeidx=5)
z = getvar(SIM_N1_000001, "z",units="m")
w_CTRL =getvar(SIM_CTRL,"wa",timeidx=5)
u_CTRL = getvar(SIM_CTRL, "ua",timeidx=5)
v_CTRL = getvar(SIM_CTRL, "va",timeidx=5)

start_point = CoordPair(x=90, y=90)
end_point = CoordPair(x=126,y=90)

w_cross_N1_000001 = vertcross(w_N1_000001,z,start_point=start_point,end_point=end_point)
qcloud_cross_N1_000001 = vertcross(QCLOUD_N1_000001,z,start_point=start_point,end_point=end_point)
u_cross_N1_000001 = vertcross(u_N1_000001,z,start_point=start_point,end_point=end_point)  

w_cross_CTRL = vertcross(w_CTRL,z,start_point=start_point,end_point=end_point)
qcloud_cross_CTRL = vertcross(QCLOUD_CTRL,z,start_point=start_point,end_point=end_point)
u_cross_CTRL = vertcross(u_CTRL,z,start_point=start_point,end_point=end_point)
ter_line = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)

xs = np.arange(0, w_cross_CTRL.shape[-1],1)+X2[0]
ys = to_np(w_cross_CTRL.coords["vertical"]) 

#making SN cross sections for SIM_N100001
start_point = CoordPair(x=134, y=50)
end_point = CoordPair(x=134,y=150)

ter_line_SN = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)
# getting the cross sections of the variables
w_cross_SN_N1_000001 = vertcross(w_N1_000001,z,start_point=start_point,end_point=end_point)
qcloud_cross_SN_N1_000001 = vertcross(QCLOUD_N1_000001,z,start_point=start_point,end_point=end_point)
v_cross_SN_N1_000001 = vertcross(v_N1_000001,z,start_point=start_point,end_point=end_point)
#thetae_cross_SN = vertcross(thetae,z,start_point=start_point,end_point=end_point)
start_point = CoordPair(x=128, y=50)
end_point = CoordPair(x=128,y=150)
w_cross_SN_CTRL = vertcross(w_CTRL,z,start_point=start_point,end_point=end_point)
qcloud_cross_SN_CTRL = vertcross(QCLOUD_CTRL,z,start_point=start_point,end_point=end_point)
v_cross_SN_CTRL = vertcross(v_CTRL,z,start_point=start_point,end_point=end_point)

xs_SN = np.arange(0, w_cross_SN_CTRL.shape[-1],1)+X2[0]
ys_SN = to_np(w_cross_SN_CTRL.coords["vertical"]) 
Xs_SN, Ys_SN = np.meshgrid(xs_SN,ys_SN)

#%%
plt.rcParams.update({'font.size': 20})
fig, axs = plt.subplots(2,2,figsize=(14, 12))

wspd_N1_000001 = np.sqrt(w_cross_N1_000001**2+u_cross_N1_000001**2)
norm = DivergingNorm(vmin=-5,vcenter=0, vmax=10)
levs = [10**(-6),10**(-4)]
plot = axs[0,0].contourf(xs+90,ys,w_cross_N1_000001,norm=norm,levels=30, cmap='RdBu_r')
axs[0,0].quiver(xs[::]+90,ys[::],u_cross_N1_000001/(3*wspd_N1_000001),w_cross_N1_000001/wspd_N1_000001,scale=25.,color='black',alpha=0.6)
axs[0,0].set_ylim((0,6500))
axs[0,0].set_xlabel('x [km]',fontsize = 22)
axs[0,0].set_ylabel('z [m]',fontsize = 22)
axs[0,0].contour(xs+90,ys,qcloud_cross_N1_000001,levels=levs,linestyles='dashed',linewidths=1.8,colors='black',alpha=0.8)
axs[0,0].fill_between(xs+90,0,to_np(ter_line),facecolor="black",alpha=0.93)
axs[0,0].text(577,6000,'a)',style='italic', fontsize=22)


wspd_CTRL = np.sqrt(w_cross_CTRL**2+u_cross_CTRL**2)
plot = axs[0,1].contourf(xs+90,ys,w_cross_CTRL,norm=norm,levels=30, cmap='RdBu_r')
axs[0,1].quiver(xs+90,ys,u_cross_CTRL/(2*wspd_CTRL),w_cross_CTRL/wspd_CTRL,scale=30.,color='black',alpha=0.6)
axs[0,1].set_ylim((0,6500))
axs[0,1].set_xlabel('x [km]',fontsize = 22)
axs[0,1].set_ylabel('z [m]',fontsize = 22)
axs[0,1].contour(xs+90,ys,qcloud_cross_CTRL,levels=levs,linestyles='dashed',linewidths=1.8,colors='black',alpha=0.8)
axs[0,1].fill_between(xs+90,0,to_np(ter_line),facecolor="black",alpha=0.93)
axs[0,1].text(577,6000,'b)',style='italic', fontsize=22)

axs[1,1].fill_between(xs_SN+50,0,to_np(ter_line_SN),facecolor="black")
axs[1,1].set_ylim((0,10000))
plot2= axs[1,1].contourf(xs_SN+50,ys_SN,qcloud_cross_SN_CTRL,levels=10,cmap="Greys")
#plt.quiver(Xs+20,Ys,v_modified/8,w_modified,scale=120.,alpha=0.5)
axs[1,1].streamplot(Xs_SN+50,Ys_SN,v_cross_SN_CTRL/1000,w_cross_SN_CTRL,density=2,arrowsize=1)
axs[1,1].set_xlabel('y [km]',fontsize = 22)
axs[1,1].contour(xs_SN+50,ys_SN,w_cross_SN_CTRL,levels =[1,2,4,8,10,15],colors='black',alpha=0.7,linewidths=1.)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
axs[1,1].set_ylabel('z [m]',fontsize = 22)
plt.title('',fontsize=20)
axs[1,1].text(540,9000,'d)',style='italic', fontsize=22)

axs[1,0].fill_between(xs_SN+50,0,to_np(ter_line_SN),facecolor="black")
axs[1,0].set_ylim((0,10000))
plot2= axs[1,0].contourf(xs_SN+50,ys_SN,qcloud_cross_SN_N1_000001,levels=10,cmap="Greys")
#plt.quiver(Xs+20,Ys,v_modified/8,w_modified,scale=120.,alpha=0.5)
axs[1,0].streamplot(Xs_SN+50,Ys_SN,v_cross_SN_N1_000001/1000,w_cross_SN_N1_000001,density=2,arrowsize=1)
axs[1,0].contour(xs_SN+50,ys_SN,w_cross_SN_N1_000001,levels =[1,2,4,8,10,15],colors='black',alpha=0.7,linewidths=1.)
plt.yticks(fontsize=20)
plt.xticks(fontsize=20)
axs[1,0].set_ylabel('z [m]',fontsize = 22)
axs[1,0].set_xlabel('y [km]',fontsize = 22)
plt.title('',fontsize=20)
axs[1,0].text(540,9000,'c)',style='italic', fontsize=22)

from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.2, 0.53, 0.6, 0.013])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m/s', labelpad=-60, y=0.92, rotation=0,fontsize=18)

cbar_ax = fig.add_axes([0.2, 0.11, 0.6, 0.013])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg/kg', labelpad=-60, y=0.95, rotation=0,fontsize=18)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.35, 
                    hspace=0.6)

#%% plot sensitività layer 3:
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_N/SIM/")
SIM_N1_00004 = Dataset(path_data + '/SIM_N1_00004_dominio2')
SIM_N3_00004 = Dataset(path_data + '/SIM_N1_00004_N3_00004_dominio2')
SIM_N3_00012 = Dataset(path_data + '/SIM_N1_00004_N3_00012_dominio2')

z = getvar(SIM_N1_00004,'z')
HGT = getvar(SIM_N1_00004,'HGT')
QRAIN_N1_00004 = getvar(SIM_N1_00004,'QRAIN',timeidx=5)
QRAIN_N3_00004 = getvar(SIM_N3_00004,'QRAIN',timeidx=5)
QRAIN_N3_00012 = getvar(SIM_N3_00012,'QRAIN',timeidx=9)

QRAIN_N1_00004_2000 = interplevel(QRAIN_N1_00004,z,2000)
QRAIN_N3_00004_2000 = interplevel(QRAIN_N3_00004,z,2000)
QRAIN_N3_00012_2000 = interplevel(QRAIN_N3_00012,z,2000)

# changing the coordinates
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_N1_00004_2000 = QRAIN_N1_00004_2000.assign_coords({'west_east' : QRAIN_N1_00004_2000.west_east.values +X2[0]+1})
QRAIN_N1_00004_2000 = QRAIN_N1_00004_2000.assign_coords({'south_north' : QRAIN_N1_00004_2000.south_north.values +X2[0]+1})
QRAIN_N3_00004_2000 = QRAIN_N3_00004_2000.assign_coords({'west_east' : QRAIN_N3_00004_2000.west_east.values +X2[0]+1})
QRAIN_N3_00004_2000 = QRAIN_N3_00004_2000.assign_coords({'south_north' : QRAIN_N3_00004_2000.south_north.values +X2[0]+1})
QRAIN_N3_00012_2000 = QRAIN_N3_00012_2000.assign_coords({'west_east' : QRAIN_N3_00012_2000.west_east.values +X2[0]+1})
QRAIN_N3_00012_2000 = QRAIN_N3_00012_2000.assign_coords({'south_north' : QRAIN_N3_00012_2000.south_north.values +X2[0]+1})

# faccio le cross sections
#HGT_SN = getvar(SIM_N3_00004,"HGT")
QCLOUD_N3_00004 = getvar(SIM_N3_00004,"QCLOUD",timeidx=5)
z = getvar(SIM_N3_00004, "z",units="m")
w_N3_00004 =getvar(SIM_N3_00004,"wa",timeidx=5)
u_N3_00004 = getvar(SIM_N3_00004, "ua",timeidx=5) 
v_N3_00004 = getvar(SIM_N3_00004, "va",timeidx=5) 
 
#HGT = getvar(SIM_N3_00012,"HGT")
QCLOUD_N3_00012 = getvar(SIM_N3_00012,"QCLOUD",timeidx=9)
z = getvar(SIM_N3_00004, "z",units="m")
w_N3_00012 =getvar(SIM_N3_00012,"wa",timeidx=9)
u_N3_00012 = getvar(SIM_N3_00012, "ua",timeidx=9) 
v_N3_00012 = getvar(SIM_N3_00012, "va",timeidx=9)

#making SN cross sections
start_point = CoordPair(x=86, y=50)
end_point = CoordPair(x=86,y=150)

ter_line_SN = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)
# getting the cross sections of the variables
w_cross_SN_N3_00004 = vertcross(w_N3_00004,z,start_point=start_point,end_point=end_point)
qcloud_cross_SN_N3_00004 = vertcross(QCLOUD_N3_00004,z,start_point=start_point,end_point=end_point)
v_cross_SN_N3_00004 = vertcross(v_N3_00004,z,start_point=start_point,end_point=end_point)
#thetae_cross_SN = vertcross(thetae,z,start_point=start_point,end_point=end_point)
start_point = CoordPair(x=119, y=50)
end_point = CoordPair(x=119,y=150)
w_cross_SN_N3_00012 = vertcross(w_N3_00012,z,start_point=start_point,end_point=end_point)
qcloud_cross_SN_N3_00012 = vertcross(QCLOUD_N3_00012,z,start_point=start_point,end_point=end_point)
v_cross_SN_N3_00012 = vertcross(v_N3_00012,z,start_point=start_point,end_point=end_point)

xs_SN = np.arange(0, w_cross_SN_N3_00012.shape[-1],1)+X2[0]
ys_SN = to_np(w_cross_SN_N3_00012.coords["vertical"]) 
Xs_SN, Ys_SN = np.meshgrid(xs_SN,ys_SN)


import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'


# making the subplots
fig = plt.figure(figsize=(15,12))
plt.rcParams.update({'font.size': 22})
ax = plt.subplot(2, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N3_00004_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_N3_00004_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(530, 645, 'a)', style='italic',fontsize=22)
plt.title('SIM_N1_00004_N3_00004',fontsize=22,y=1.05)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.vlines(573,533,665,linestyle='dashed',color='black')
plt.xlim((525,675))
plt.ylim((540,660))


ax = plt.subplot(2, 2, 2)
plot=QRAIN_N3_00012_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_N3_00012_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(530, 645, 'b)', style='italic',fontsize=22)
plt.title('SIM_N1_00004_N3_00012',fontsize=22,y=1.05)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.vlines(606,533,665,linestyle='dashed',color='black')
plt.xlim((525,675))
plt.ylim((540,660))


ax = plt.subplot(2, 2, 3)
norm = DivergingNorm(vmin=w_cross_SN_N3_00004.min(),vcenter=0, vmax=w_cross_SN_N3_00004.max())
ax.fill_between(xs_SN+50,0,to_np(ter_line_SN),facecolor="black")
ax.set_ylim((0,10000))
plot2= ax.contourf(xs_SN+50,ys_SN,w_cross_SN_N3_00004,levels=20,cmap="RdBu_r",norm=norm)
#plt.quiver(Xs+20,Ys,v_modified/8,w_modified,scale=120.,alpha=0.5)
ax.streamplot(Xs_SN+50,Ys_SN,v_cross_SN_N3_00004/1000,w_cross_SN_N3_00004,density=2,arrowsize=1)
ax.contour(xs_SN+50,ys_SN,qcloud_cross_SN_N3_00004,levels =[10**(-5)],colors='black',alpha=0.7,linewidths=2.)
plt.xlabel('y [km]')
plt.ylabel('z [m]')
plt.title('',fontsize=20)
ax.text(540,9000,'c)',style='italic', fontsize=22)


ax = plt.subplot(2, 2, 4)
ax.fill_between(xs_SN+50,0,to_np(ter_line_SN),facecolor="black")
ax.set_ylim((0,10000))
plot3= ax.contourf(xs_SN+50,ys_SN,w_cross_SN_N3_00012,levels=20,cmap="RdBu_r",norm=norm)
#plt.quiver(Xs+20,Ys,v_modified/8,w_modified,scale=120.,alpha=0.5)
ax.streamplot(Xs_SN+50,Ys_SN,v_cross_SN_N3_00012/1000,w_cross_SN_N3_00012,density=2,arrowsize=1)
ax.contour(xs_SN+50,ys_SN,qcloud_cross_SN_N3_00012,levels =[10**(-5)],colors='black',alpha=0.7,linewidths=2.)
plt.xlabel('y [km]')
plt.ylabel('z [m]')
ax.text(540,9000,'d)',style='italic', fontsize=22)


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.54, 0.5, 0.012])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=5, y=0.95, rotation=0,fontsize=20)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.12, 0.5, 0.012])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
cb.set_label('m/s', labelpad=5, y=0.95, rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.35, 
                    hspace=0.47)

#%% plot RH profiles:
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_N/SIM/")
SIM_CTRL = Dataset(path_data+'/SIM_wrfinput_prova_dominio2')
RH_INCREASED5 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_RH_profile/Sondaggio_RH_increased5/Parameters_Sondaggio_RH_increased5.txt",delimiter = '\t')
RH_REDUCED5 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_RH_profile/Sondaggio_RH_reduced5/Parameters_Sondaggio_RH_reduced5.txt",delimiter = '\t')
RH_INCREASED5_under2600 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_RH_profile/Sondaggio_RH_increased5_under2600m/Parameters_Sondaggio_RH_increased5_under2600m.txt",delimiter = '\t')
RH_INCREASED5_over2300 = pd.read_csv("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_RH_profile/Sondaggio_RH_increased5_over2300m/Parameters_Sondaggio_RH_increased5_over2300m.txt",delimiter = '\t')

RH_CTRL = getvar(SIM_CTRL,'rh',timeidx=0)
z = getvar(SIM_CTRL,'z')
RH_increased5 = RH_INCREASED5['RH']
z_increased5 = RH_INCREASED5['Height']
RH_reduced5 = RH_REDUCED5['RH']
z_reduced5 = RH_REDUCED5['Height']
RH_increased5_under2600 = RH_INCREASED5_under2600['RH']
RH_increased5_over2300 = RH_INCREASED5_over2300['RH']


fig = plt.figure(figsize=(18,8))
plt.rcParams.update({'font.size': 24})

ax = plt.subplot(1, 2, 1)
ax.plot(RH_CTRL[:,50,50],z[:,50,50],linewidth=4, label = 'SIM$\_$CTRL',color='green',alpha=0.7,linestyle='dashed')
ax.plot(RH_increased5,z_increased5,linewidth=4,label='SIM$\_$RH$\_$INCR5',color='blue',alpha=0.5)
ax.plot(RH_reduced5,z_reduced5,linewidth=4,label = 'SIM$\_$RH$\_$RED5',color='red',alpha=0.7)
plt.legend(fontsize=19,loc=3)
plt.xlabel('RH $\%$')
plt.ylabel('Heigth [m]')
plt.ylim((50,6000))
plt.xlim((40,100))
plt.text(42., 5600, 'a)', style='italic', fontsize=26)
plt.grid()

ax = plt.subplot(1, 2, 2)
ax.plot(RH_increased5,z_increased5,linewidth=5, label = 'SIM$\_$RH$\_$INCR5',color='black',alpha=1.,linestyle='dashed')
ax.plot(RH_increased5_under2600,z_increased5,linewidth=5,label='SIM$\_$RH$\_$INCR5$\_$LL',color='blue',alpha=0.5)
ax.plot(RH_increased5_over2300,z_reduced5,linewidth=5,label = 'SIM$\_$RH$\_$INCR5$\_$UL',color='red',alpha=0.5)
plt.legend(fontsize=19,loc=3)
plt.xlabel('RH $\%$')
plt.ylabel('Heigth [m]')
plt.ylim((50,6000))
plt.xlim((40,100))
plt.text(42., 5600, 'b)', style='italic', fontsize=26)
plt.grid()

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.28, 
                    hspace=0.6)

#%% plot QRain sensitivity to relative humidity

path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_RH/")
SIM_RH_reduced5 = Dataset(path_data + '/SIM_RH_reduced5_dominio2')
SIM_RH_increased5 = Dataset(path_data + '/SIM_RH_increased5_dominio2')
SIM_RH_increased5_under2600m = Dataset(path_data + '/SIM_RH_increased5_under2600m_dominio2')
SIM_RH_increased5_over2300m = Dataset(path_data + '/SIM_RH_increased5_over2300m_dominio2')


z = getvar(SIM_RH_reduced5,'z')
HGT = getvar(SIM_RH_reduced5,'HGT')
QRAIN_RH_reduced5 = getvar(SIM_RH_reduced5,'QRAIN',timeidx=5)
QRAIN_RH_increased5 = getvar(SIM_RH_increased5,'QRAIN',timeidx=3)
QRAIN_RH_increased5_under2600m = getvar(SIM_RH_increased5_under2600m,'QRAIN',timeidx=3)
QRAIN_RH_increased5_over2300m = getvar(SIM_RH_increased5_over2300m,'QRAIN',timeidx=3)

QRAIN_RH_reduced5_2000 = interplevel(QRAIN_RH_reduced5,z,2000)
QRAIN_RH_increased5_2000 = interplevel(QRAIN_RH_increased5,z,2000)
QRAIN_RH_increased5_under2600m_2000 = interplevel(QRAIN_RH_increased5_under2600m,z,2000)
QRAIN_RH_increased5_over2300m_2000 = interplevel(QRAIN_RH_increased5_over2300m,z,2000)

#%%
# changing the coordinates
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_RH_reduced5_2000 = QRAIN_RH_reduced5_2000.assign_coords({'west_east' : QRAIN_RH_reduced5_2000.west_east.values +X2[0]+1})
QRAIN_RH_reduced5_2000 = QRAIN_RH_reduced5_2000.assign_coords({'south_north' : QRAIN_RH_reduced5_2000.south_north.values +X2[0]+1})
QRAIN_RH_increased5_2000 = QRAIN_RH_increased5_2000.assign_coords({'west_east' : QRAIN_RH_increased5_2000.west_east.values +X2[0]+1})
QRAIN_RH_increased5_2000 = QRAIN_RH_increased5_2000.assign_coords({'south_north' : QRAIN_RH_increased5_2000.south_north.values +X2[0]+1})
QRAIN_RH_increased5_under2600m_2000 = QRAIN_RH_increased5_under2600m_2000.assign_coords({'west_east' : QRAIN_RH_increased5_under2600m_2000.west_east.values +X2[0]+1})
QRAIN_RH_increased5_under2600m_2000 = QRAIN_RH_increased5_under2600m_2000.assign_coords({'south_north' : QRAIN_RH_increased5_under2600m_2000.south_north.values +X2[0]+1})
QRAIN_RH_increased5_over2300m_2000 = QRAIN_RH_increased5_over2300m_2000.assign_coords({'west_east' : QRAIN_RH_increased5_over2300m_2000.west_east.values +X2[0]+1})
QRAIN_RH_increased5_over2300m_2000 = QRAIN_RH_increased5_over2300m_2000.assign_coords({'south_north' : QRAIN_RH_increased5_over2300m_2000.south_north.values +X2[0]+1})

#%%
# making the subplots
fig = plt.figure(figsize=(15,12))
plt.rcParams.update({'font.size': 20})
ax = plt.subplot(2, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_RH_reduced5_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_RH_reduced5_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(530, 645, 'a)', style='italic',fontsize=22)
plt.title('SIM_RH_RED5',fontsize=22,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))

ax = plt.subplot(2, 2, 2)
plot=QRAIN_RH_increased5_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_RH_increased5_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(530, 645, 'b)', style='italic',fontsize=22)
plt.title('SIM_RH_INCR5',fontsize=22,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))

ax = plt.subplot(2, 2, 3)
plot=QRAIN_RH_increased5_under2600m_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_RH_increased5_under2600m_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(530, 645, 'c)', style='italic',fontsize=22)
plt.title('SIM_RH_INCR5_LL',fontsize=22,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))

ax = plt.subplot(2, 2, 4)
plot=QRAIN_RH_increased5_over2300m_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_RH_increased5_over2300m_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(530, 645, 'd)', style='italic',fontsize=22)
plt.title('SIM_RH_INCR5_UL',fontsize=22,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((525,675))
plt.ylim((540,660))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.12, 0.5, 0.015])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=5, y=0.95, rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.27, 
                    hspace=0.33)


#%% OROGRAFIA SPORCAAAAAAAAA
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
path_data2 = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_HGT_0_50 = Dataset(path_data + '/SIM_HGT_pert_0_50_each_cell_dominio2')
SIM_HGT_100_150_dominio2 = Dataset(path_data2 + '/SIM_HGT_100_150_dominio2')
SIM_HGT_100_150_dominio3 = Dataset(path_data2 + '/SIM_HGT_100_150_dominio3')
SIM_HGT_50_100 = Dataset(path_data + '/SIM_HGT_pert_50_100_each_cell_dominio2')
SIM_HGT_100_150_first = Dataset(path_data + '/SIM_HGT_100_150_firstpart_dominio2')
SIM_U20_HGT_100_150 = Dataset(path_data + '/SIM_U20_HGT_100_150_each_cell_dominio2')


QRAIN_HGT_0_50 = getvar(SIM_HGT_0_50,'QRAIN',timeidx=2)
z = getvar(SIM_HGT_0_50,'z')
z_dominio2 = getvar(SIM_HGT_100_150_dominio2,'z')
z_dominio3 = getvar(SIM_HGT_100_150_dominio3,'z')
HGT_0_50 = getvar(SIM_HGT_0_50,'HGT')
QRAIN_HGT_50_100 = getvar(SIM_HGT_50_100,'QRAIN',timeidx=2)
HGT_50_100 = getvar(SIM_HGT_50_100,'HGT')
QRAIN_HGT_100_150_dominio2 = getvar(SIM_HGT_100_150_dominio2,'QRAIN',timeidx=2)
HGT_100_150_dominio2 = getvar(SIM_HGT_100_150_dominio2,'HGT')
QRAIN_HGT_100_150_dominio3 = getvar(SIM_HGT_100_150_dominio3,'QRAIN',timeidx=2)
HGT_100_150_dominio3 = getvar(SIM_HGT_100_150_dominio3,'HGT')
QRAIN_HGT_100_150_first = getvar(SIM_HGT_100_150_first,'QRAIN',timeidx=2)
HGT_100_150_first = getvar(SIM_HGT_100_150_first,'HGT')
QRAIN_U20_HGT_100_150 = getvar(SIM_U20_HGT_100_150,'QRAIN',timeidx=5)
HGT_U20_100_150 = getvar(SIM_U20_HGT_100_150,'HGT')

QRAIN_HGT_0_50_2000 = interplevel(QRAIN_HGT_0_50,z,2000)
QRAIN_HGT_50_100_2000 = interplevel(QRAIN_HGT_50_100,z,2000)
QRAIN_HGT_100_150_dominio2_2000 = interplevel(QRAIN_HGT_100_150_dominio2,z_dominio2,2000)
QRAIN_HGT_100_150_first_2000 = interplevel(QRAIN_HGT_100_150_first,z,2000)
QRAIN_U20_HGT_100_150_2000 = interplevel(QRAIN_U20_HGT_100_150,z,2000)
QRAIN_HGT_100_150_dominio3_2000 = interplevel(QRAIN_HGT_100_150_dominio3,z_dominio3,2000)

# changing coordinates
HGT_0_50 = HGT_0_50.assign_coords({'west_east' : HGT_0_50.west_east.values +X2[0]+1})
HGT_0_50 = HGT_0_50.assign_coords({'south_north' : HGT_0_50.south_north.values +X2[0]+1})
QRAIN_HGT_0_50_2000 = QRAIN_HGT_0_50_2000.assign_coords({'west_east' : QRAIN_HGT_0_50_2000.west_east.values +X2[0]+1})
QRAIN_HGT_0_50_2000 = QRAIN_HGT_0_50_2000.assign_coords({'south_north' : QRAIN_HGT_0_50_2000.south_north.values +X2[0]+1})
HGT_50_100 = HGT_50_100.assign_coords({'west_east' : HGT_50_100.west_east.values +X2[0]+1})
HGT_50_100 = HGT_50_100.assign_coords({'south_north' : HGT_50_100.south_north.values +X2[0]+1})
QRAIN_HGT_50_100_2000 = QRAIN_HGT_50_100_2000.assign_coords({'west_east' : QRAIN_HGT_50_100_2000.west_east.values +X2[0]+1})
QRAIN_HGT_50_100_2000 = QRAIN_HGT_50_100_2000.assign_coords({'south_north' : QRAIN_HGT_50_100_2000.south_north.values +X2[0]+1})
HGT_100_150_dominio2 = HGT_100_150_dominio2.assign_coords({'west_east' : HGT_100_150_dominio2.west_east.values +X2[0]+1})
HGT_100_150_dominio2 = HGT_100_150_dominio2.assign_coords({'south_north' : HGT_100_150_dominio2.south_north.values +X2[0]+1})
QRAIN_HGT_100_150_dominio2_2000 = QRAIN_HGT_100_150_dominio2_2000.assign_coords({'west_east' : QRAIN_HGT_100_150_dominio2_2000.west_east.values +X2[0]+1})
QRAIN_HGT_100_150_dominio2_2000 = QRAIN_HGT_100_150_dominio2_2000.assign_coords({'south_north' : QRAIN_HGT_100_150_dominio2_2000.south_north.values +X2[0]+1})
HGT_100_150_dominio3 = HGT_100_150_dominio3.assign_coords({'west_east' : HGT_100_150_dominio3.west_east.values/2 +X3[0]+1})
HGT_100_150_dominio3 = HGT_100_150_dominio3.assign_coords({'south_north' : HGT_100_150_dominio3.south_north.values/2 +X3[0]+1})
QRAIN_HGT_100_150_dominio3_2000 = QRAIN_HGT_100_150_dominio3_2000.assign_coords({'west_east' : QRAIN_HGT_100_150_dominio3_2000.west_east.values/2 +X3[0]+1})
QRAIN_HGT_100_150_dominio3_2000 = QRAIN_HGT_100_150_dominio3_2000.assign_coords({'south_north' : QRAIN_HGT_100_150_dominio3_2000.south_north.values/2 +X3[0]+1})


HGT_U20_100_150 = HGT_U20_100_150.assign_coords({'west_east' : HGT_U20_100_150.west_east.values +X2[0]+1})
HGT_U20_100_150 = HGT_U20_100_150.assign_coords({'south_north' : HGT_U20_100_150.south_north.values +X2[0]+1})
QRAIN_U20_HGT_100_150_2000 = QRAIN_U20_HGT_100_150_2000.assign_coords({'west_east' : QRAIN_U20_HGT_100_150_2000.west_east.values +X2[0]+1})
QRAIN_U20_HGT_100_150_2000 = QRAIN_U20_HGT_100_150_2000.assign_coords({'south_north' : QRAIN_U20_HGT_100_150_2000.south_north.values +X2[0]+1})
HGT_100_150_first = HGT_100_150_first.assign_coords({'west_east' : HGT_100_150_first.west_east.values +X2[0]+1})
HGT_100_150_first = HGT_100_150_first.assign_coords({'south_north' : HGT_100_150_first.south_north.values +X2[0]+1})
QRAIN_HGT_100_150_first_2000 = QRAIN_HGT_100_150_first_2000.assign_coords({'west_east' : QRAIN_HGT_100_150_first_2000.west_east.values +X2[0]+1})
QRAIN_HGT_100_150_first_2000 = QRAIN_HGT_100_150_first_2000.assign_coords({'south_north' : QRAIN_HGT_100_150_first_2000.south_north.values +X2[0]+1})

#%% making the plot
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(13,15))
plt.rcParams.update({'font.size': 16})
gs = gridspec.GridSpec(6, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
ax1 = plt.subplot(gs[:2, :2])
plot=QRAIN_HGT_0_50_2000[45:180,45:180].plot.contourf(ax=ax1, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_HGT_0_50_2000[45:180,45:180].plot.contour(ax=ax1, levels=levels,colors = 'black',alpha = 0.5)  
HGT_0_50[45:180,45:180].plot.contour(ax=ax1,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax1.text(537, 645, 'a)', style='italic',fontsize=22)
plt.title('SIM_HGT_0_50',fontsize=18,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((545,655))

ax2 = plt.subplot(gs[:2, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_HGT_50_100_2000[45:180,45:180].plot.contourf(ax=ax2, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_HGT_50_100_2000[45:180,45:180].plot.contour(ax=ax2, levels=levels,colors = 'black',alpha = 0.5)  
HGT_50_100[45:180,45:180].plot.contour(ax=ax2,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax2.text(537, 645, 'b)', style='italic',fontsize=22)
plt.title('SIM_HGT_50_100',fontsize=18,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((545,655))

ax3 = plt.subplot(gs[2:4, :2])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_HGT_100_150_dominio2_2000[45:180,45:180].plot.contourf(ax=ax3, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_HGT_100_150_dominio2_2000[45:180,45:180].plot.contour(ax=ax3, levels=levels,colors = 'black',alpha = 0.5)  
HGT_100_150_dominio2[45:180,45:180].plot.contour(ax=ax3,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax3.text(537, 645, 'c)', style='italic',fontsize=22)
plt.title('SIM_HGT_100_150 (Res=1km)',fontsize=18,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((545,655))

ax4 = plt.subplot(gs[2:4, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_HGT_100_150_dominio3_2000.plot.contourf(ax=ax4, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_HGT_100_150_dominio3_2000.plot.contour(ax=ax4, levels=levels,colors = 'black',alpha = 0.5)  
HGT_100_150_dominio3.plot.contour(ax=ax4,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax4.text(537, 645, 'd)', style='italic',fontsize=22)
plt.title('SIM_HGT_100_150 (Res=500m)',fontsize=18,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((545,655))

ax5 = plt.subplot(gs[4:6, 1:3])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_HGT_100_150_first_2000[45:180,45:180].plot.contourf(ax=ax5, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_HGT_100_150_first_2000[45:180,45:180].plot.contour(ax=ax5, levels=levels,colors = 'black',alpha = 0.5)  
HGT_100_150_first[45:180,45:180].plot.contour(ax=ax5,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax5.text(537, 645, 'e)', style='italic',fontsize=22)
plt.title('SIM_HGT_100_150_FP',fontsize=18,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((545,655))

# ax6 = plt.subplot(gs[4:6, 2:4])
# levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
# plot=QRAIN_U20_HGT_100_150_2000[45:180,45:180].plot.contourf(ax=ax6, levels=levels,cmap='Blues',add_colorbar=False)  
# #cb = plt.colorbar(plot,shrink=0.9)
# QRAIN_U20_HGT_100_150_2000[45:180,45:180].plot.contour(ax=ax6, levels=levels,colors = 'black',alpha = 0.5)  
# HGT_U20_100_150[45:180,45:180].plot.contour(ax=ax6,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
# ax6.text(537, 645, 'f)', style='italic',fontsize=22)
# plt.title('SIM_U20_HGT_100_150',fontsize=18,y=1.01)
# plt.xlabel('x [km]')
# plt.ylabel('y [km]')
# plt.xlim((535,665))
# plt.ylim((545,655))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.04, 0.5, 0.008])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=5, y=0.95,x=0.55,rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.6, 
                    hspace=1.)

#%% plto sezione cumulata più distribuzione precipitazione lungo il ridge:
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_HGT_0_50 = Dataset(path_data + '/SIM_HGT_pert_0_50_each_cell_dominio2')
SIM_HGT_100_150 = Dataset(path_data + '/SIM_HGT_100_150_each_cell_dominio2')

RAINNC_HGT_0_50 = getvar(SIM_HGT_0_50,'RAINNC',timeidx=ALL_TIMES)
RAINNC_HGT_100_150 = getvar(SIM_HGT_100_150,'RAINNC',timeidx=ALL_TIMES)

# computing the west-east averaged rainfall over the ridge
RAIN_FINALE_HGT_100_150 = RAINNC_HGT_100_150[5]-RAINNC_HGT_100_150[0]
RAIN_FINALE_HGT_0_50 = RAINNC_HGT_0_50[5]-RAINNC_HGT_0_50[0]

RAIN_mean_HGT_100_150 = RAIN_FINALE_HGT_100_150[55:160,50:175].mean(dim='west_east')
RAIN_mean_HGT_0_50 = RAIN_FINALE_HGT_0_50[55:160,50:175].mean(dim='west_east')

# changing coordinates
RAINNC_HGT_100_150 = RAINNC_HGT_100_150.assign_coords({'west_east' : RAINNC_HGT_100_150.west_east.values +X2[0]+1})
RAINNC_HGT_100_150 = RAINNC_HGT_100_150.assign_coords({'south_north' : RAINNC_HGT_100_150.south_north.values +X2[0]+1})
RAINNC_HGT_0_50 = RAINNC_HGT_0_50.assign_coords({'west_east' : RAINNC_HGT_0_50.west_east.values +X2[0]+1})
RAINNC_HGT_0_50 = RAINNC_HGT_0_50.assign_coords({'south_north' : RAINNC_HGT_0_50.south_north.values +X2[0]+1})

RAIN_mean_HGT_100_150 = RAIN_mean_HGT_100_150.assign_coords({'south_north' : RAIN_mean_HGT_100_150.south_north.values +X2[0]+1+55})
RAIN_mean_HGT_0_50 = RAIN_mean_HGT_0_50.assign_coords({'south_north' : RAIN_mean_HGT_0_50.south_north.values +X2[0]+1+55})


fig = plt.figure(figsize=(18,18))
plt.rcParams.update({'font.size': 32})
ax = plt.subplot(2, 1, 1)
(RAINNC_HGT_0_50[5,110,50:175]-RAINNC_HGT_0_50[0,110,50:175]).plot(linewidth=5,alpha=0.7,color='red',label='Range 0-50 m')
(RAINNC_HGT_100_150[5,110,50:175]-RAINNC_HGT_100_150[0,110,50:175]).plot(linewidth=5.,alpha=0.4,color='blue',label='Range 100-150 m')
plt.grid()
plt.title('')
plt.xlabel('x [km]')
plt.ylabel('RAIN [mm]')
plt.legend(loc=2)
#ax.text(663,165,'a)',style='italic',fontsize=40)

# ax = plt.subplot(2, 1, 2)
# RAIN_mean_HGT_0_50.plot(linewidth=5,color='red',alpha=0.7,label='HGT_0_50')
# RAIN_mean_HGT_100_150.plot(linewidth=5,color='blue',alpha=0.4,label='HGT_100_150')
# plt.grid()
# plt.xlabel('y [km]')
# plt.ylabel('RAIN [mm]')
# plt.legend(loc='best')
# ax.text(538,52,'b)',style='italic',fontsize=40)

# plt.subplots_adjust(left=0.1,
#                     bottom=0.1, 
#                     right=0.95, 
#                     top=0.9, 
#                     wspace=0.8, 
#                     hspace=0.3)

#%% plot TRIGGER MECHANISM
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_HGT_0_50 = Dataset(path_data + '/SIM_HGT_pert_0_50_each_cell_dominio2')
SIM_HGT_50_100 = Dataset(path_data + '/SIM_HGT_100_150_each_cell_dominio2')

z_HGT_0_50 = getvar(SIM_HGT_0_50, "z",units="m",timeidx=2)
w_HGT_0_50 =getvar(SIM_HGT_0_50,"wa",timeidx=2)
v_HGT_0_50 = getvar(SIM_HGT_0_50, "va",timeidx=2)  # default is in m/s
QCLOUD_HGT_0_50 = getvar(SIM_HGT_0_50,"QCLOUD",timeidx=2)
HGT_0_50 = getvar(SIM_HGT_0_50,"HGT")

start_point = CoordPair(x=101,y=50)
end_point = CoordPair(x=101, y=130)

# Get the terrain heights along the cross section line
ter_line_HGT_0_50 = interpline(HGT_0_50,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables
w_cross_HGT_0_50 = vertcross(w_HGT_0_50,z_HGT_0_50,start_point=start_point,end_point=end_point)
qcloud_cross_HGT_0_50 = vertcross(QCLOUD_HGT_0_50,z_HGT_0_50,start_point=start_point,end_point=end_point)
v_cross_HGT_0_50 = vertcross(v_HGT_0_50,z_HGT_0_50,start_point=start_point,end_point=end_point)

w_modified_HGT_0_50 = np.zeros((100,81))
w_modified_HGT_0_50[:] = np.nan
w_modified_HGT_0_50[::2,::5]=w_cross_HGT_0_50[::2,::5]

v_modified_HGT_0_50 = np.zeros((100,81))
v_modified_HGT_0_50[:] = np.nan
v_modified_HGT_0_50[::2,::5]=v_cross_HGT_0_50[::2,::5]

z_HGT_50_100 = getvar(SIM_HGT_50_100, "z",units="m")
w_HGT_50_100 =getvar(SIM_HGT_50_100,"wa",timeidx=2)
v_HGT_50_100 = getvar(SIM_HGT_50_100, "va",timeidx=2)  # default is in m/s
QCLOUD_HGT_50_100 = getvar(SIM_HGT_50_100,"QCLOUD",timeidx=2)
HGT_50_100 = getvar(SIM_HGT_50_100,"HGT")

start_point = CoordPair(x=101,y=50)
end_point = CoordPair(x=101, y=130)

# Get the terrain heights along the cross section line
ter_line_HGT_50_100 = interpline(HGT_50_100,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables
w_cross_HGT_50_100 = vertcross(w_HGT_50_100,z_HGT_50_100,start_point=start_point,end_point=end_point)
qcloud_cross_HGT_50_100 = vertcross(QCLOUD_HGT_50_100,z_HGT_50_100,start_point=start_point,end_point=end_point)
v_cross_HGT_50_100 = vertcross(v_HGT_50_100,z_HGT_50_100,start_point=start_point,end_point=end_point)

w_modified_HGT_50_100 = np.zeros((100,81))
w_modified_HGT_50_100[:] = np.nan
w_modified_HGT_50_100[::2,::5]=w_cross_HGT_50_100[::2,::5]

v_modified_HGT_50_100 = np.zeros((100,81))
v_modified_HGT_50_100[:] = np.nan
v_modified_HGT_50_100[::2,::5]=v_cross_HGT_50_100[::2,::5]


xs = np.arange(0, w_cross_HGT_50_100.shape[-1], 1)
ys = to_np(w_cross_HGT_50_100.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)


#%% making the plot

fig = plt.figure(figsize=(14,12))
plt.rcParams.update({'font.size': 22})
ax = plt.subplot(2, 2, 1)
levs = [0.25,0.5,1,2,4,8,12]
levs_w = [-0.03,0,0.25,0.5,1,2,4,8,12]
plot=plt.contourf(xs+Y2[0]+50,ys,w_cross_HGT_0_50,levels=levs,cmap="Blues")
CS=plt.contour(xs+Y2[0]+50,ys,w_cross_HGT_0_50,levels=levs_w,colors='black',linewidths=0.8)
plt.clabel(CS, inline=1, fontsize=13,fmt='%.1f')
plt.quiver(Xs+Y2[0]+50,Ys,v_modified_HGT_0_50/8,w_modified_HGT_0_50,scale=120.,alpha=0.5)
plt.contour(xs+ Y2[0]+50,ys,qcloud_cross_HGT_0_50,levels=[10**(-5)],colors = 'Black',linewidths=2.5)
plt.fill_between(xs+Y2[0]+50,0,to_np(ter_line_HGT_0_50),facecolor="black")
plt.ylim((0,7000))
plt.title('Perturbations in range 0-50 m', y =1.02,fontsize=24)
plt.xlabel('y [km]')
plt.ylabel('z [m]')
plt.yticks()
plt.xticks()
plt.text(538,6500,'a)',style='italic')

ax = plt.subplot(2, 2, 2)
levs_w = [-0.03,0,0.25,0.5,1,2,4,8,12]
plt.contourf(xs+Y2[0]+50,ys,w_cross_HGT_50_100,levels=levs,cmap="Blues")
CS=plt.contour(xs+Y2[0]+50,ys,w_cross_HGT_50_100,levels=levs_w,colors='black',linewidths=0.8)
plt.clabel(CS, inline=True, inline_spacing=-5,levels=[0.5,1,2,4,8],fontsize=13,fmt='%.1f')
plt.quiver(Xs+Y2[0]+50,Ys,v_modified_HGT_50_100/8,w_modified_HGT_50_100,scale=120.,alpha=0.5)
plt.contour(xs+ Y2[0]+50,ys,qcloud_cross_HGT_50_100,levels=[10**(-5)],colors = 'Black',linewidths=2.5)
plt.fill_between(xs+Y2[0]+50,0,to_np(ter_line_HGT_50_100),facecolor="black")
plt.ylim((0,7000))
plt.xlabel('y [km]')
plt.ylabel('z [m]')
plt.title('Perturbations in range 100-150 m', y =1.02,fontsize=24)
plt.yticks()
plt.xticks()
plt.text(538,6500,'b)',style='italic')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.44, 0.5, 0.008])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('m/s', labelpad=5, y=0.95,x=0.55,rotation=0,fontsize=22)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.4, 
                    hspace=0.1)

#%% plot QRAIN simulazioni a vento costante
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_U20_HGT_100_150 = Dataset(path_data + '/SIM_U20_HGT_100_150_each_cell_dominio2')
SIM_U20_HGT_0_50 = Dataset(path_data + '/SIM_U20_HGT_0_50_each_cell_dominio2')
#extracting variables
z=getvar(SIM_U20_HGT_100_150,'z')
QRAIN_U20_HGT_100_150 = getvar(SIM_U20_HGT_100_150,'QRAIN',timeidx=5)
HGT_U20_100_150 = getvar(SIM_U20_HGT_100_150,'HGT')
QRAIN_U20_HGT_100_150_2000 = interplevel(QRAIN_U20_HGT_100_150,z,2000)

QRAIN_U20_HGT_0_50 = getvar(SIM_U20_HGT_0_50,'QRAIN',timeidx=5)
HGT_U20_0_50 = getvar(SIM_U20_HGT_0_50,'HGT')
QRAIN_U20_HGT_0_50_2000 = interplevel(QRAIN_U20_HGT_0_50,z,2000)


# changing coordinates
HGT_U20_100_150 = HGT_U20_100_150.assign_coords({'west_east' : HGT_U20_100_150.west_east.values +X2[0]+1})
HGT_U20_100_150 = HGT_U20_100_150.assign_coords({'south_north' : HGT_U20_100_150.south_north.values +X2[0]+1})
QRAIN_U20_HGT_100_150_2000 = QRAIN_U20_HGT_100_150_2000.assign_coords({'west_east' : QRAIN_U20_HGT_100_150_2000.west_east.values +X2[0]+1})
QRAIN_U20_HGT_100_150_2000 = QRAIN_U20_HGT_100_150_2000.assign_coords({'south_north' : QRAIN_U20_HGT_100_150_2000.south_north.values +X2[0]+1})

HGT_U20_0_50 = HGT_U20_0_50.assign_coords({'west_east' : HGT_U20_0_50.west_east.values +X2[0]+1})
HGT_U20_0_50 = HGT_U20_0_50.assign_coords({'south_north' : HGT_U20_0_50.south_north.values +X2[0]+1})
QRAIN_U20_HGT_0_50_2000 = QRAIN_U20_HGT_0_50_2000.assign_coords({'west_east' : QRAIN_U20_HGT_0_50_2000.west_east.values +X2[0]+1})
QRAIN_U20_HGT_0_50_2000 = QRAIN_U20_HGT_0_50_2000.assign_coords({'south_north' : QRAIN_U20_HGT_0_50_2000.south_north.values +X2[0]+1})

#%% making the plot
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(14,12))
plt.rcParams.update({'font.size': 20})
gs = gridspec.GridSpec(4, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
ax1 = plt.subplot(gs[:2, :2])
plot=QRAIN_U20_HGT_0_50_2000[45:180,45:180].plot.contourf(ax=ax1, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_U20_HGT_0_50_2000[45:180,45:180].plot.contour(ax=ax1, levels=levels,colors = 'black',alpha = 0.5)  
HGT_U20_0_50[45:180,45:180].plot.contour(ax=ax1,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax1.text(537, 645, 'a)', style='italic',fontsize=22)
plt.title('SIM_V20_HGT_0_50',fontsize=20,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((545,655))

ax2 = plt.subplot(gs[:2, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_U20_HGT_100_150_2000[45:180,45:180].plot.contourf(ax=ax2, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_U20_HGT_100_150_2000[45:180,45:180].plot.contour(ax=ax2, levels=levels,colors = 'black',alpha = 0.5)  
HGT_U20_100_150[45:180,45:180].plot.contour(ax=ax2,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax2.text(537, 645, 'b)', style='italic',fontsize=22)
plt.title('SIM_V20_HGT_100_150',fontsize=20,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((545,655))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.48, 0.5, 0.008])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=5, y=0.95,x=0.55,rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.6, 
                    hspace=1.)

#%% PLOT QRAIN AND MECHANISM SINGLE BUMP
path_data2 = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_BUMP_A = Dataset(path_data2+'/SIM_BUMP_Y180_NOPERT_dominio2')
SIM_BUMP_B = Dataset(path_data2+'/SIM_BUMP_Y187_NOPERT_dominio2')

HGT_BUMP_A = getvar(SIM_BUMP_A,'HGT')
z_BUMP_A = getvar(SIM_BUMP_A,'z')
z_agl_A = getvar(SIM_BUMP_A,'height_agl')
HGT_BUMP_B = getvar(SIM_BUMP_B,'HGT')
z_BUMP_B = getvar(SIM_BUMP_B,'z')
z_agl_B = getvar(SIM_BUMP_B,'height_agl')
QRAIN_BUMP_A = getvar(SIM_BUMP_A,'QRAIN',timeidx=3)
QRAIN_BUMP_B = getvar(SIM_BUMP_B,'QRAIN',timeidx=2)
QRAIN_BUMP_A_2000 = interplevel(QRAIN_BUMP_A,z_BUMP_A,2000)
QRAIN_BUMP_B_2000 = interplevel(QRAIN_BUMP_B,z_BUMP_B,2000)
wspd = getvar(SIM_BUMP_A,'wspd',timeidx=3)
#%%
# changing the coordinates
HGT_BUMP_A = HGT_BUMP_A.assign_coords({'west_east' : HGT_BUMP_A.west_east.values +X2[0]+1})
HGT_BUMP_A = HGT_BUMP_A.assign_coords({'south_north' : HGT_BUMP_A.south_north.values +X2[0]+1})
HGT_BUMP_B = HGT_BUMP_B.assign_coords({'west_east' : HGT_BUMP_B.west_east.values +X2[0]+1})
HGT_BUMP_B = HGT_BUMP_B.assign_coords({'south_north' : HGT_BUMP_B.south_north.values +X2[0]+1})
QRAIN_BUMP_A_2000 = QRAIN_BUMP_A_2000.assign_coords({'west_east' : QRAIN_BUMP_A_2000.west_east.values +X2[0]+1})
QRAIN_BUMP_A_2000 = QRAIN_BUMP_A_2000.assign_coords({'south_north' : QRAIN_BUMP_A_2000.south_north.values +X2[0]+1})
QRAIN_BUMP_B_2000 = QRAIN_BUMP_B_2000.assign_coords({'west_east' : QRAIN_BUMP_B_2000.west_east.values +X2[0]+1})
QRAIN_BUMP_B_2000 = QRAIN_BUMP_B_2000.assign_coords({'south_north' : QRAIN_BUMP_B_2000.south_north.values +X2[0]+1})

# gravity wave mechanism (as Fuhrer and Schar)
QCLOUD_BUMP_A = getvar(SIM_BUMP_A, "QCLOUD",timeidx=3)
w_BUMP_A = getvar(SIM_BUMP_A, "wa",timeidx=3)
QCLOUD_BUMP_B = getvar(SIM_BUMP_B, "QCLOUD",timeidx=2)
w_BUMP_B = getvar(SIM_BUMP_B, "wa",timeidx=2)
QCLOUD_BUMP_A_1500 = interplevel(QCLOUD_BUMP_A,z_agl_A,2000)
QCLOUD_BUMP_B_1500 = interplevel(QCLOUD_BUMP_B,z_agl_B,1500)
w_BUMP_A_1500 = interplevel(w_BUMP_A,z_agl_A,2000)
w_BUMP_B_1500 = interplevel(w_BUMP_B,z_agl_B,1500)

# changing coordinates
QCLOUD_BUMP_A_1500 = QCLOUD_BUMP_A_1500.assign_coords({'west_east' : QCLOUD_BUMP_A_1500.west_east.values +X2[0]+1})
QCLOUD_BUMP_A_1500 = QCLOUD_BUMP_A_1500.assign_coords({'south_north' : QCLOUD_BUMP_A_1500.south_north.values +X2[0]+1})
QCLOUD_BUMP_B_1500 = QCLOUD_BUMP_B_1500.assign_coords({'west_east' : QCLOUD_BUMP_B_1500.west_east.values +X2[0]+1})
QCLOUD_BUMP_B_1500 = QCLOUD_BUMP_B_1500.assign_coords({'south_north' : QCLOUD_BUMP_B_1500.south_north.values +X2[0]+1})
w_BUMP_A_1500 = w_BUMP_A_1500.assign_coords({'west_east' : w_BUMP_A_1500.west_east.values +X2[0]+1})
w_BUMP_A_1500 = w_BUMP_A_1500.assign_coords({'south_north' : w_BUMP_A_1500.south_north.values +X2[0]+1})
w_BUMP_B_1500 = w_BUMP_B_1500.assign_coords({'west_east' : w_BUMP_B_1500.west_east.values +X2[0]+1})
w_BUMP_B_1500 = w_BUMP_B_1500.assign_coords({'south_north' : w_BUMP_B_1500.south_north.values +X2[0]+1})


#%% making the plot

fig = plt.figure(figsize=(18,15))
plt.rcParams.update({'font.size': 28})
gs = gridspec.GridSpec(4, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
ax1 = plt.subplot(gs[:2, :2])
plot=QRAIN_BUMP_A_2000.plot.contourf(ax=ax1, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_BUMP_A_2000.plot.contour(ax=ax1, levels=levels,colors = 'black',alpha = 0.5)  
CS=HGT_BUMP_A.plot.contour(ax=ax1,levels = [100,450,1000,1400,1499],colors='black',alpha=0.5)
#plt.clabel(CS, inline=1, inline_spacing=-20,fontsize=15,fmt='%.1f')
ax1.text(538, 530, 'a)', style='italic',fontsize=32)
plt.title('SIM_BUMP_A',fontsize=30,y=1.02)
ax1.text(645, 530, 't=3 h',fontsize=32)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
#plt.vlines(599.5,525,665, linestyle='dashed',color = 'black')
plt.xlim((535,665))
plt.ylim((525,655))

ax2 = plt.subplot(gs[:2, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_BUMP_B_2000.plot.contourf(ax=ax2, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_BUMP_B_2000.plot.contour(ax=ax2, levels=levels,colors = 'black',alpha = 0.5)  
HGT_BUMP_B[45:180,45:180].plot.contour(ax=ax2,levels = [100,450,1000,1400,1499],colors='black',alpha=0.5)
ax2.text(538, 530, 'b)', style='italic',fontsize=32)
ax2.text(645, 530, 't=2 h',fontsize=32)
plt.title('SIM_BUMP_B',fontsize=30,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((525,655))

ax3 = plt.subplot(gs[2:4, :2])
#levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=QCLOUD_BUMP_A_1500.plot.contourf(ax=ax3, levels=np.linspace(0.,0.0015,10),cmap='Greys',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
CS=w_BUMP_A_1500.plot.contour(ax=ax3,levels=np.linspace(-1,1,8),colors='black',linewidths=1.3)
plt.clabel(CS, inline=1, inline_spacing=-10,fontsize=22,fmt='%.1f')
HGT_BUMP_A.plot.contour(ax=ax3,levels = [100,450,1000,1400],colors='black',alpha=0.5,linestyles='dashed',linewidths=5.)
ax3.text(558, 528, 'c)', style='italic',fontsize=32)
plt.title('SIM_BUMP_A',fontsize=30,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((555,645))
plt.ylim((525,600))

ax4 = plt.subplot(gs[2:4, 2:])
#levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=QCLOUD_BUMP_B_1500.plot.contourf(ax=ax4, levels=np.linspace(0.0001,0.0015,20),cmap='Greys',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
CS=w_BUMP_B_1500.plot.contour(ax=ax4,levels=np.linspace(-1,1,9),colors='black',linewidths=1.3)
plt.clabel(CS, inline=1, inline_spacing=-10,fontsize=22,fmt='%.1f')
HGT_BUMP_B.plot.contour(ax=ax4,levels = [100,450,1000,1400],colors='black',alpha=0.5,linestyles='dashed',linewidths=5.)
ax4.text(558, 528, 'd)', style='italic',fontsize=32)
plt.title('SIM_BUMP_B',fontsize=30,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((555,645))
plt.ylim((525,600))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.5, 0.5, 0.01])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=-75, y=0.95,x=0.55,rotation=0,fontsize=26)

from matplotlib import ticker
cbar_ax = fig.add_axes([0.25, 0.02, 0.5, 0.01])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg/kg', labelpad=5, y=0.95,x=0.55,rotation=0,fontsize=26)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.8, 
                    hspace=1.6)


#%% PLOT QRAIN AND MECHANISM SINUSOIDAL COLUMN
path_data2 = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_COLUMN_A = Dataset(path_data2+'/SIM_COLUMN_Y180_NO_BUMPETTINO')
SIM_COLUMN_B = Dataset(path_data2+'/SIM_COLUMN_Y185_NO_BUMPETTINO')

HGT_COLUMN_A = getvar(SIM_COLUMN_A,'HGT')
z_COLUMN_A = getvar(SIM_COLUMN_A,'z')
z_agl_COLUMN_A = getvar(SIM_COLUMN_A,'height_agl')
HGT_COLUMN_B = getvar(SIM_COLUMN_B,'HGT')
z_COLUMN_B = getvar(SIM_COLUMN_B,'z')
z_agl_COLUMN_B = getvar(SIM_COLUMN_B,'height_agl')
QRAIN_COLUMN_A = getvar(SIM_COLUMN_A,'QRAIN',timeidx=2)
QRAIN_COLUMN_B = getvar(SIM_COLUMN_B,'QRAIN',timeidx=3)
QRAIN_COLUMN_A_2000 = interplevel(QRAIN_COLUMN_A,z_COLUMN_A,2000)
QRAIN_COLUMN_B_2000 = interplevel(QRAIN_COLUMN_B,z_COLUMN_B,2000)

# changing the coordinates
HGT_COLUMN_A = HGT_COLUMN_A.assign_coords({'west_east' : HGT_COLUMN_A.west_east.values +X2[0]+1})
HGT_COLUMN_A = HGT_COLUMN_A.assign_coords({'south_north' : HGT_COLUMN_A.south_north.values +X2[0]+1})
HGT_COLUMN_B = HGT_COLUMN_B.assign_coords({'west_east' : HGT_COLUMN_B.west_east.values +X2[0]+1})
HGT_COLUMN_B = HGT_COLUMN_B.assign_coords({'south_north' : HGT_COLUMN_B.south_north.values +X2[0]+1})
QRAIN_COLUMN_A_2000 = QRAIN_COLUMN_A_2000.assign_coords({'west_east' : QRAIN_COLUMN_A_2000.west_east.values +X2[0]+1})
QRAIN_COLUMN_A_2000 = QRAIN_COLUMN_A_2000.assign_coords({'south_north' : QRAIN_COLUMN_A_2000.south_north.values +X2[0]+1})
QRAIN_COLUMN_B_2000 = QRAIN_COLUMN_B_2000.assign_coords({'west_east' : QRAIN_COLUMN_B_2000.west_east.values +X2[0]+1})
QRAIN_COLUMN_B_2000 = QRAIN_COLUMN_B_2000.assign_coords({'south_north' : QRAIN_COLUMN_B_2000.south_north.values +X2[0]+1})

# gravity wave mechanism (as Fuhrer and Schar)
QCLOUD_COLUMN_A = getvar(SIM_COLUMN_A, "QCLOUD",timeidx=2)
w_COLUMN_A = getvar(SIM_COLUMN_A, "wa",timeidx=2)
QCLOUD_COLUMN_B = getvar(SIM_COLUMN_B, "QCLOUD",timeidx=3)
w_COLUMN_B = getvar(SIM_COLUMN_B, "wa",timeidx=3)
QCLOUD_COLUMN_A_1500 = interplevel(QCLOUD_COLUMN_A,z_agl_COLUMN_A,1500)
QCLOUD_COLUMN_B_1500 = interplevel(QCLOUD_COLUMN_B,z_agl_COLUMN_B,1500)
w_COLUMN_A_1500 = interplevel(w_COLUMN_A,z_agl_COLUMN_A,1500)
w_COLUMN_B_1500 = interplevel(w_COLUMN_B,z_agl_COLUMN_B,1500)

# changing coordinates
QCLOUD_COLUMN_A_1500 = QCLOUD_COLUMN_A_1500.assign_coords({'west_east' : QCLOUD_COLUMN_A_1500.west_east.values +X2[0]+1})
QCLOUD_COLUMN_A_1500 = QCLOUD_COLUMN_A_1500.assign_coords({'south_north' : QCLOUD_COLUMN_A_1500.south_north.values +X2[0]+1})
QCLOUD_COLUMN_B_1500 = QCLOUD_COLUMN_B_1500.assign_coords({'west_east' : QCLOUD_COLUMN_B_1500.west_east.values +X2[0]+1})
QCLOUD_COLUMN_B_1500 = QCLOUD_COLUMN_B_1500.assign_coords({'south_north' : QCLOUD_COLUMN_B_1500.south_north.values +X2[0]+1})
w_COLUMN_A_1500 = w_COLUMN_A_1500.assign_coords({'west_east' : w_COLUMN_A_1500.west_east.values +X2[0]+1})
w_COLUMN_A_1500 = w_COLUMN_A_1500.assign_coords({'south_north' : w_COLUMN_A_1500.south_north.values +X2[0]+1})
w_COLUMN_B_1500 = w_COLUMN_B_1500.assign_coords({'west_east' : w_COLUMN_B_1500.west_east.values +X2[0]+1})
w_COLUMN_B_1500 = w_COLUMN_B_1500.assign_coords({'south_north' : w_COLUMN_B_1500.south_north.values +X2[0]+1})

#%% making the plot
import matplotlib.gridspec as gridspec
fig = plt.figure(figsize=(20,18))
plt.rcParams.update({'font.size': 34})
gs = gridspec.GridSpec(4, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
ax1 = plt.subplot(gs[:2, :2])
plot=QRAIN_COLUMN_A_2000.plot.contourf(ax=ax1, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_COLUMN_A_2000.plot.contour(ax=ax1, levels=levels,colors = 'black',alpha = 0.5)  
HGT_COLUMN_A.plot.contour(ax=ax1,levels = [-100,100,500,1000,1400,1499],colors='black',alpha=0.5)
ax1.text(538, 642, 'a)', style='italic',fontsize=32)
plt.title('SIM_COLUMN_A',fontsize=30,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((525,655))

ax2 = plt.subplot(gs[:2, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_COLUMN_B_2000.plot.contourf(ax=ax2, levels=levels,cmap='Blues',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
QRAIN_COLUMN_B_2000.plot.contour(ax=ax2, levels=levels,colors = 'black',alpha = 0.5)  
HGT_COLUMN_B[45:180,45:180].plot.contour(ax=ax2,levels = [-100,100,500,1000,1400,1499],colors='black',alpha=0.5)
ax2.text(538, 642, 'b)', style='italic',fontsize=32)
plt.title('SIM_COLUMN_B',fontsize=30,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((525,655))

ax3 = plt.subplot(gs[2:4, :2])
#levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=QCLOUD_COLUMN_A_1500.plot.contourf(ax=ax3, levels=np.linspace(0.,0.0015,11),cmap='Greys',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
CS=w_COLUMN_A_1500.plot.contour(ax=ax3,levels=np.linspace(-1,1,8),colors='black',linewidths=1.3)
plt.clabel(CS, inline=1,inline_spacing=-10, fontsize=24,fmt='%.1f')
HGT_COLUMN_A.plot.contour(ax=ax3,levels = [-100,100,500,1000,1400],colors='black',alpha=0.5,linewidths=5.)
ax3.text(558, 528, 'c)', style='italic',fontsize=32)
plt.title('SIM_COLUMN_A',fontsize=30,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((555,645))
plt.ylim((525,600))

ax4 = plt.subplot(gs[2:4, 2:])
#levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=QCLOUD_COLUMN_B_1500.plot.contourf(ax=ax4, levels=np.linspace(0.,0.0015,11),cmap='Greys',add_colorbar=False)  
#cb = plt.colorbar(plot,shrink=0.9)
CS=w_COLUMN_B_1500.plot.contour(ax=ax4,levels=np.linspace(-1,1,9),colors='black',linewidths=1.3)
plt.clabel(CS, inline=1, inline_spacing=-10,fontsize=24,fmt='%.1f')
HGT_COLUMN_B.plot.contour(ax=ax4,levels = [-100,100,500,1000,1400],colors='black',alpha=0.5,linewidths=5.)
ax4.text(600, 528, 'd)', style='italic',fontsize=32)
plt.title('SIM_COLUMN_B',fontsize=30,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((540,605))
plt.ylim((525,600))

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.25, 0.5, 0.5, 0.01])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('kg/kg', labelpad=-89, y=0.98,x=0.55,rotation=0,fontsize=34)

from matplotlib import ticker
cbar_ax = fig.add_axes([0.25, 0.02, 0.5, 0.01])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg/kg', labelpad=5, y=0.95,x=0.55,rotation=0,fontsize=34)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.8, 
                    hspace=1.6)

#%% PLOT FLOW DEFLECTION VALLEY LATERALI:
path_data2 = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_COLUMN_B = Dataset(path_data2+'/SIM_COLUMN_Y185_NO_BUMPETTINO')

HGT_COLUMN_B = getvar(SIM_COLUMN_B,'ter')
z_COLUMN_B = getvar(SIM_COLUMN_B,'z')
z_agl_COLUMN_B = getvar(SIM_COLUMN_B,'height_agl')
QCLOUD_COLUMN_B = getvar(SIM_COLUMN_B,'QCLOUD',timeidx=3)
w_COLUMN_B = getvar(SIM_COLUMN_B, "wa",timeidx=3)
v_COLUMN_B = getvar(SIM_COLUMN_B, "va",timeidx=3)
u_COLUMN_B = getvar(SIM_COLUMN_B, "ua",timeidx=3)
wspd_COLUMN_B = getvar(SIM_COLUMN_B, "wspd",timeidx=3)
# interpolating on levels
QCLOUD_COLUMN_B_1500 = interplevel(QCLOUD_COLUMN_B,z_agl_COLUMN_B,1500)
w_COLUMN_B_1500 = interplevel(w_COLUMN_B,z_agl_COLUMN_B,1500)
u_COLUMN_B_50 = interplevel(u_COLUMN_B,z_agl_COLUMN_B,100)
v_COLUMN_B_50 = interplevel(v_COLUMN_B,z_agl_COLUMN_B,100)
wspd_COLUMN_B_50 = interplevel(wspd_COLUMN_B,z_agl_COLUMN_B,100)
# computing horizontal convergence
from numpy import diff
du_dx = np.array(np.gradient(u_COLUMN_B_50,axis=1))
dv_dy = np.array(np.gradient(v_COLUMN_B_50,axis=0))
div = (du_dx + dv_dy)

# changing coordinates to the wind
v_COLUMN_B_50 = v_COLUMN_B_50.assign_coords({'west_east' : v_COLUMN_B_50.west_east.values +X2[0]+1})
v_COLUMN_B_50 = v_COLUMN_B_50.assign_coords({'south_north' : v_COLUMN_B_50.south_north.values +X2[0]+1})
u_COLUMN_B_50 = u_COLUMN_B_50.assign_coords({'west_east' : u_COLUMN_B_50.west_east.values +X2[0]+1})
u_COLUMN_B_50 = u_COLUMN_B_50.assign_coords({'south_north' : u_COLUMN_B_50.south_north.values +X2[0]+1})
wspd_COLUMN_B_50 = wspd_COLUMN_B_50.assign_coords({'west_east' : wspd_COLUMN_B_50.west_east.values +X2[0]+1})
wspd_COLUMN_B_50 = wspd_COLUMN_B_50.assign_coords({'south_north' : wspd_COLUMN_B_50.south_north.values +X2[0]+1})
# changing coordinates to qcloud, w and hgt
QCLOUD_COLUMN_B_1500 = QCLOUD_COLUMN_B_1500.assign_coords({'west_east' : QCLOUD_COLUMN_B_1500.west_east.values +X2[0]+1})
QCLOUD_COLUMN_B_1500 = QCLOUD_COLUMN_B_1500.assign_coords({'south_north' : QCLOUD_COLUMN_B_1500.south_north.values +X2[0]+1})
w_COLUMN_B_1500 = w_COLUMN_B_1500.assign_coords({'west_east' : w_COLUMN_B_1500.west_east.values +X2[0]+1})
w_COLUMN_B_1500 = w_COLUMN_B_1500.assign_coords({'south_north' : w_COLUMN_B_1500.south_north.values +X2[0]+1})
HGT_COLUMN_B = HGT_COLUMN_B.assign_coords({'west_east' : HGT_COLUMN_B.west_east.values +X2[0]+1})
HGT_COLUMN_B = HGT_COLUMN_B.assign_coords({'south_north' : HGT_COLUMN_B.south_north.values +X2[0]+1})

# making divergence a data array
import xarray as xr
Xs, Ys = np.meshgrid(HGT_COLUMN_B.west_east.values,HGT_COLUMN_B.south_north.values)

div = xr.DataArray(data=-div,dims=["x", "y"],
    coords=dict(
        west_east=(["x", "y"], Xs),
        south_north=(["x", "y"], Ys),),
    attrs=dict(
        description="Divergence.",
        units="s$^-1$",
    ),)
#%% making the subplot
fig = plt.figure(figsize=(18,14))
plt.rcParams.update({'font.size': 22})
gs = gridspec.GridSpec(4, 4)
ax1 = plt.subplot(gs[:2, :2])
plot=div[::2,::2].plot.contourf(levels=np.linspace(-3,3,41),cmap='seismic',x='west_east',y = 'south_north',add_colorbar=False)
#wspd_50[3,::2,::2].plot.contourf(levels=np.linspace(0,20,21),cmap='Reds')
cb=plt.colorbar(plot)
cb.set_label('$10^{-3} \ s^{-1}$',rotation=0,y=1.08,labelpad=-50)
levels =[-200,-100,100,500,1000,1400,1499]
CS=HGT_COLUMN_B[::2,::2].plot.contour(levels=levels,colors='k',alpha=0.5,linewidths=4.)
qv = ax1.quiver(u_COLUMN_B_50[::3,::3].west_east, u_COLUMN_B_50[::3,::3].south_north, u_COLUMN_B_50[::3,::3], v_COLUMN_B_50[::3,::3],color='black',scale = 150.)
plt.title('')
ax1.text(522, 597, 'a)', style='italic',fontsize=32)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((530,600))
plt.ylim((519,599))
#%%
fig = plt.figure(figsize=(18,14))
plt.rcParams.update({'font.size': 22})
gs = gridspec.GridSpec(4, 4)
ax1 = plt.subplot(gs[:2, :2])
#div[::2,::2].plot.contourf(levels=np.linspace(-3*10**-3,3*10**-3,41),cmap='seismic',x='west_east',y = 'south_north')
plot=wspd_COLUMN_B_50[::2,::2].plot.contourf(levels=np.linspace(0,20,21),cmap='Reds',add_colorbar=False)
cb=plt.colorbar(plot)
cb.set_label('m/s',rotation=0,y=1.08,labelpad=-40)
levels =[-200,-100,100,500,1000,1400,1499]
CS=HGT_COLUMN_B[::2,::2].plot.contour(levels=levels,colors='k',alpha=0.5,linewidths=4.)
qv = ax1.quiver(u_COLUMN_B_50[::3,::3].west_east, u_COLUMN_B_50[::3,::3].south_north, u_COLUMN_B_50[::3,::3], v_COLUMN_B_50[::3,::3],color='black',scale = 150.)
plt.title('')
ax1.text(522, 597, 'b)', style='italic',fontsize=32)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((530,600))
plt.ylim((519,599))
#%%
fig = plt.figure(figsize=(18,14))
plt.rcParams.update({'font.size': 22})
gs = gridspec.GridSpec(4, 4)
ax1 = plt.subplot(gs[:2, :2])
#levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=QCLOUD_COLUMN_B_1500.plot.contourf(ax=ax1, levels=np.linspace(0.,0.0015,11),cmap='Greys',add_colorbar=False)  
cb=plt.colorbar(plot2)
cb.set_label('kg/kg',rotation=0,y=1.11,labelpad=-40)
CS=w_COLUMN_B_1500.plot.contour(ax=ax1,levels=np.linspace(-1,1,8),colors='black',linewidths=1.2)
plt.clabel(CS, inline=1, inline_spacing=-5,fontsize=18,fmt='%.1f')
HGT_COLUMN_B.plot.contour(ax=ax1,levels = [-100,100,500,1000,1400],colors='black',alpha=0.5,linewidths=5.)
ax1.text(522, 597, 'c)', style='italic',fontsize=32)
plt.title('',fontsize=30,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((530,600))
plt.ylim((519,599))

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=1.2, 
                    hspace=1.6)

#%% PLOT PRECIPITAZIONI SIM_REAL
import matplotlib.colors as mcolors
import cartopy.crs as crs
import cartopy
import wrf
# Hourly Precipitation
path_data= os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_REAL_FINALE/")
ncfile = Dataset(path_data + '/SIM_ORIGINAL_SLP_ORIGINAL_WIND_ore18')
# Get the latitude and longitude
lat = getvar(ncfile, "lat")
lon = getvar(ncfile, "lon")
z= getvar(ncfile,'z')
TIMES = getvar(ncfile,"times",timeidx=ALL_TIMES)
HGT = getvar(ncfile,"HGT")
RAINNC = getvar(ncfile,"RAINNC",timeidx = ALL_TIMES)

# defining the colormap
nws_precip_colors = [
    (1.0, 1.0, 1.0),
    "#04e9e7",  # 0.01 - 0.10 inches
    "#019ff4",  # 0.10 - 0.25 inches
    "#0300f4",  # 0.25 - 0.50 inches
    "#02fd02",  # 0.50 - 0.75 inches
    "#01c501",  # 0.75 - 1.00 inches
    "#008e00",  # 1.00 - 1.50 inches
    "#fdf802",  # 1.50 - 2.00 inches
    "#e5bc00",  # 2.00 - 2.50 inches
    "#fd9500",  # 2.50 - 3.00 inches
    "#fd0000",  # 3.00 - 4.00 inches
    "#d40000",  # 4.00 - 5.00 inches
    "#bc0000",  # 5.00 - 6.00 inches
    "#f800fd",  # 6.00 - 8.00 inches
    "#9854c6",  # 8.00 - 10.00 inches
    #"#fdfdfd"   # 10.00+
]
precip_colormap = mcolors.ListedColormap(nws_precip_colors)

#%% modifico orografia per mettere coastlines
LANDMASK = ncfile['LANDMASK']
# modifying orography
HGT = HGT.where(HGT>0,0.1)
HGT.plot.contourf()
#%% cumulata
from matplotlib.patches import Rectangle
RAIN_hour10 = RAINNC[10]-RAINNC[9]
RAIN_hour9 = RAINNC[9]-RAINNC[8]
cart_proj=wrf.get_cartopy(HGT)
plt.rcParams.update({'font.size': 24})
plt.figure(figsize=(22,12))
ax = plt.axes(projection=cart_proj)
transform = crs.PlateCarree()._as_mpl_transform(ax)
ax.set_aspect('auto')
#levels = [10,15,20,22,25,30,35,40,45,50,75,80,100,150,200,300]
levels=[5,7,9,10,12,14,16,18,21,24,27,30,40,50]
plot = RAIN_hour10.plot.contourf(ax=ax,cmap = precip_colormap,levels=levels,x='XLONG',y='XLAT',transform = crs.PlateCarree(),add_colorbar=False)
cb = plt.colorbar(plot,shrink=0.8)
cb.set_label('mm/h', labelpad=-45, y=1.05, rotation=0)
HGT.plot.contour(ax=ax,levels=[200,600,1000,1500,2000,2500,3000],colors='black',alpha=0.7,linewidths=0.5,x='XLONG',y='XLAT',transform = crs.PlateCarree())
#plt.title('Hourly_Precipitation hour 9',y=1.03)
plt.title('')
lon_Venezia = 12.22
lat_Venezia = 45.46
#plt.plot(lon_Venezia,lat_Venezia,'bo',label='Venice',transform = crs.PlateCarree())
#ax.annotate("Venice", (lon_Venezia-0.2, lat_Venezia+0.08),xycoords=transform,fontsize=30,fontweight='bold')
grd = ax.gridlines(draw_labels=True,x_inline=False, y_inline=False)
grd.top_labels = grd.right_labels = False
ax.add_feature(cartopy.feature.COASTLINE)
#ax.annotate("b)", (10.05,47.1),xycoords=transform,fontsize=30)
import matplotlib.patches as mpatches
ax.add_patch(mpatches.Rectangle(xy=[10.65, 45.6], width=0.28, height=1.1,
                                     facecolor='blue',
                                     alpha=0.2,
                                     angle=-30,
                                     transform=crs.PlateCarree()))
ax.annotate('Adige valley and Sarca valley effect',(10.1,45.1),xycoords=transform,fontsize=30,fontweight='bold')
plt.arrow(10.8, 45.2,0, 0.1, width = 0.02,color='black',transform=crs.PlateCarree())

#%%
cart_proj=wrf.get_cartopy(HGT)
plt.rcParams.update({'font.size': 24})
plt.figure(figsize=(22,12))
ax = plt.axes(projection=cart_proj)
transform = crs.PlateCarree()._as_mpl_transform(ax)
ax.set_aspect('auto')
#levels = [10,15,20,22,25,30,35,40,45,50,75,80,100,150,200,300]
levels=[5,7,9,10,12,14,16,18,21,24,27,30,40,50]
plot = RAIN_hour9.plot.contourf(ax=ax,cmap = precip_colormap,levels=levels,x='XLONG',y='XLAT',transform = crs.PlateCarree(),add_colorbar=False)
cb = plt.colorbar(plot,shrink=0.8)
cb.set_label('mm/h', labelpad=-45, y=1.05, rotation=0)
HGT.plot.contour(ax=ax,levels=[200,600,1000,1500,2000,2500,3000],colors = 'black',alpha=0.7,linewidths=0.5,x='XLONG',y='XLAT',transform = crs.PlateCarree())
plt.title('')
grd = ax.gridlines(draw_labels=True,x_inline=False, y_inline=False)
grd.top_labels = grd.right_labels = False
ax.add_feature(cartopy.feature.COASTLINE)
#ax.annotate("a)", (10.05,47.1),xycoords=transform,fontsize=30)

# # disegno le due cross section
# #line 4
# point1 = [12.184, 45.522]
# point2 = [12.019, 46.632]
# #start_point = CoordPair(lat=45.268, lon=13.283)
# #end_point = CoordPair(lat=46.893, lon=12.051)
# x_values = [point1[0], point2[0]]
# y_values = [point1[1], point2[1]]
# ax.plot(x_values, y_values,label='Cross 1',color='black',linewidth=3.,alpha=0.7,transform=crs.PlateCarree())
# #plt.legend(loc = 'lower right')
# #line5
# point1 = [12.1, 45.41]
# point2 = [11.95, 46.693]
# x_values = [point1[0], point2[0]]
# y_values = [point1[1], point2[1]]
# plt.plot(x_values, y_values,label='Cross 2',color='black',linewidth=3.,alpha=0.7,transform=crs.PlateCarree())
# #plt.legend(loc = 'lower right')
# ax.annotate("A", (12., 45.31),xycoords=transform,fontsize=30)
# ax.annotate("B", (11.85, 46.693),xycoords=transform,fontsize=30)
# ax.annotate("C", (12.184, 45.462),xycoords=transform,fontsize=30)
# ax.annotate("D", (12.019, 46.648),xycoords=transform,fontsize=30)

#%% making cross section AB
plt.rcParams.update({'font.size': 34})
start_point = CoordPair(lat=45.522, lon=12.184)
end_point = CoordPair(lat=46.632, lon=12.019)

z= getvar(ncfile,'z')  # non usare mai height_agl!!!!
u = getvar(ncfile,'ua',timeidx = 9)
v = getvar(ncfile,'va',timeidx = 9)
w = getvar(ncfile,'wa',timeidx = 9)
q_rain = getvar(ncfile,'QRAIN',timeidx=9)
q_graup = getvar(ncfile,'QGRAUP',timeidx=9)
q_snow = getvar(ncfile,'QSNOW',timeidx=9)
q_cloud = getvar(ncfile,'QCLOUD',timeidx=9)
theta_e = getvar(ncfile,'eth',timeidx=9)

# Compute the vertical cross-section interpolation.  Also, include the
# lat/lon points along the cross-section.
ter_line = interpline(HGT,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
w_cross = vertcross(w,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
v_cross = vertcross(v,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
u_cross = vertcross(u,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
#q_cross = vertcross(q_cloud+q_rain+q_graup+q_snow,z,start_point=start_point,end_point=end_point)
q_cross = vertcross(q_cloud,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
theta_e_cross = vertcross(theta_e,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)

xs = np.arange(0, w_cross.shape[-1], 1)
ys = to_np(w_cross.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)
ys = to_np(w_cross.coords["vertical"])

fig = plt.figure(figsize=(22,11))
ax = plt.axes()
plt.fill_between(xs,0,to_np(ter_line),facecolor="black")
#plot = plt.contourf(xs,ys,q_cross,levels=np.linspace(0.00005,0.003,21),cmap=new_cmap)
plot = plt.contourf(xs,ys,w_cross,levels=np.linspace(-5, 5, 21), extend='both',cmap = 'bwr')
#plot = plt.contourf(xs,ys,theta_e_cross,cmap='bwr',levels = np.linspace(310,330,21))
cb = plt.colorbar(plot)
cb.set_label('m/s', labelpad=-80, y=1.08, rotation=0)
ax.streamplot(Xs,Ys,v_cross/1000,w_cross,density=2,arrowsize=2.3,linewidth=2.7)
# Set the x-ticks to use latitude and longitude labels.
coord_pairs = to_np(ter_line.coords["xy_loc"])
x_ticks = np.arange(coord_pairs.shape[0])
x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}")
            for pair in to_np(coord_pairs)]
ax.set_xticks(x_ticks[::20])
ax.set_xticklabels(x_labels[::20], rotation=45)
ax.set_xlabel('lat,lon [°N,°E]')
ax.set_ylabel('z [m]')
ax.set_ylim((0,8000))
ax.annotate("b)", (4,7500.20),fontsize=35)
ax.grid()

#%%
plt.rcParams.update({'font.size': 34})
start_point = CoordPair(lat=45.41, lon=12.1)
end_point = CoordPair(lat=46.693, lon=11.95)

# Compute the vertical cross-section interpolation.  Also, include the
# lat/lon points along the cross-section.
ter_line = interpline(HGT,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
w_cross = vertcross(w,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
v_cross = vertcross(v,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
u_cross = vertcross(u,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
#q_cross = vertcross(q_cloud+q_rain+q_graup+q_snow,z,start_point=start_point,end_point=end_point)
q_cross = vertcross(q_cloud,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)
theta_e_cross = vertcross(theta_e,z,wrfin=ncfile,start_point=start_point,end_point=end_point,latlon=True, meta=True)

xs = np.arange(0, w_cross.shape[-1], 1)
ys = to_np(w_cross.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)
ys = to_np(w_cross.coords["vertical"])

fig = plt.figure(figsize=(22,11))
plt.rcParams.update({'font.size': 34})
ax = plt.axes()
plt.fill_between(xs,0,to_np(ter_line),facecolor="black")
#plot = plt.contourf(xs,ys,q_cross,levels=np.linspace(0.00005,0.003,21),cmap=new_cmap)
plot = plt.contourf(xs,ys,w_cross,levels=np.linspace(-5, 5, 21), extend='both',cmap = 'bwr')
#plot = plt.contourf(xs,ys,theta_e_cross,cmap='bwr',levels = np.linspace(310,330,21))
cb = plt.colorbar(plot)
cb.set_label('m/s', labelpad=-80, y=1.08, rotation=0)
ax.streamplot(Xs,Ys,v_cross/1000,w_cross,density=2,arrowsize=2.3,linewidth=2.7)
# Set the x-ticks to use latitude and longitude labels.
coord_pairs = to_np(ter_line.coords["xy_loc"])
x_ticks = np.arange(coord_pairs.shape[0])
x_labels = [pair.latlon_str(fmt="{:.2f}, {:.2f}")
            for pair in to_np(coord_pairs)]
ax.set_xticks(x_ticks[::20])
ax.set_xticklabels(x_labels[::20], rotation=45)
ax.set_xlabel('lat,lon [°N,°E]')
ax.set_ylabel('z [m]')
ax.set_ylim((0,8000))
ax.annotate("a)", (4,7500.20),fontsize=35)
ax.grid()

#%% zoom sulla parte centrale del dominio per fare confronto con le simulazioni di Giovannini:

longitude = RAIN_hour9.XLONG[0,:]
latitude = RAIN_hour9.XLAT[:,0]

latitude = latitude.where((latitude>45.3) & (latitude<47))
longitude = longitude.where(longitude<13.5)
longitude = longitude.where(longitude>10.5)

#levels=[3,7,9,12,14,16,18,20,24,27,30,35,40,50,70]
levels=[5,7,9,10,12,14,16,18,21,24,27,30,40,50]
fig = plt.figure(figsize=(22,12))
plt.rcParams.update({'font.size': 24})
ax = plt.axes(projection=cart_proj)
#ax.set_extent([latitude.min(),latitude.max(),longitude.min(),longitude.max()])
ax.set_extent([longitude.min(),longitude.max(),latitude.min(),latitude.max()])
HGT.plot.contour(ax=ax,levels=10,colors='black',alpha=0.5,x='XLONG',y='XLAT',transform = crs.PlateCarree())
plot=RAIN_hour9.plot.contourf(ax=ax,cmap = precip_colormap,levels=levels,x='XLONG',y='XLAT',transform = crs.PlateCarree(),add_colorbar=False)
cb = plt.colorbar(plot,shrink=0.8)
cb.set_label('mm/h', labelpad=-45, y=1.05, rotation=0)
#pu, pv = u_200[::5,::5], v_200[::5,::5]
#qv = ax.quiver(to_np(pu.XLONG), to_np(pu.XLAT), to_np(pu), to_np(pv), transform=crs.PlateCarree(),color='blue',scale=650.)
ax.add_feature(cartopy.feature.COASTLINE)
grd = ax.gridlines(draw_labels=True,x_inline=False, y_inline=False)
grd.top_labels = grd.right_labels = False

  
    


