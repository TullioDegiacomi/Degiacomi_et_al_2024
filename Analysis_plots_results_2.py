# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 11:24:11 2021

@author: Tullio Degiacomi
DEFINIZIONE DEL DOMINIO
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
# making an automatic code for the precipitation plot normalized
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
Created on Wed Feb  2 16:01:28 2022

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

import matplotlib.gridspec as gridspec
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'
#%% plot appendice rumore numerico e hourly_Prec ore 6 SIM_200m e SIM_500m su cui applico Fourier
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIMULAZIONI_DOMINIO_CENTRATO/SIM_200m_noCoriolis_sounding_sud")
SIM_1KM = Dataset(path_data + '/SIM_200m_noCoriolis_sounding_sud_dominio2')
SIM_200m = Dataset(path_data + '/SIM_200m_noCoriolis_sounding_sud_dominio3')

HGT = getvar(SIM_1KM,'HGT')
z = getvar(SIM_1KM,'z')
z_200 = getvar(SIM_200m,'z')
RAINNC_1KM = getvar(SIM_1KM,'RAINNC',timeidx=ALL_TIMES)
RAINNC_200m = getvar(SIM_200m,'RAINNC',timeidx=ALL_TIMES)
QRAIN_200m = getvar(SIM_200m,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_200m = interplevel(QRAIN_200m,z_200,2000)

HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
RAINNC_1KM = RAINNC_1KM.assign_coords({'west_east' : RAINNC_1KM.west_east.values +X2[0]+1})
RAINNC_1KM = RAINNC_1KM.assign_coords({'south_north' : RAINNC_1KM.south_north.values +X2[0]+1})
QRAIN_2000_200m = QRAIN_2000_200m.assign_coords({'west_east' : (QRAIN_2000_200m.west_east.values/5 +X3[0]+1)})
QRAIN_2000_200m = QRAIN_2000_200m.assign_coords({'south_north' : (QRAIN_2000_200m.south_north.values/5 +X3[0]+1)})
RAINNC_200m = RAINNC_200m.assign_coords({'west_east' : (RAINNC_200m.west_east.values/5 +X3[0]+1)})
RAINNC_200m = RAINNC_200m.assign_coords({'south_north' : (RAINNC_200m.south_north.values/5 +X3[0]+1)})

#%% extracting data for sim500m
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIMULAZIONI_DOMINIO_CENTRATO/SIM_500m_noCoriolis_sounding_sud")
SIM_1KM_500m = Dataset(path_data + '/SIM_500m_noCoriolis_sounding_sud_dominio2')
SIM_500m = Dataset(path_data + '/SIM_500m_noCoriolis_sounding_sud_dominio3')

z_500 = getvar(SIM_500m,'z')
RAINNC_1KM_500m = getvar(SIM_1KM_500m,'RAINNC',timeidx=ALL_TIMES)
RAINNC_500m = getvar(SIM_500m,'RAINNC',timeidx=ALL_TIMES)
QRAIN_500m = getvar(SIM_500m,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000_500m = interplevel(QRAIN_500m,z_500,2000)

RAINNC_1KM_500m = RAINNC_1KM_500m.assign_coords({'west_east' : RAINNC_1KM_500m.west_east.values +X2[0]+1})
RAINNC_1KM_500m = RAINNC_1KM_500m.assign_coords({'south_north' : RAINNC_1KM_500m.south_north.values +X2[0]+1})
QRAIN_2000_500m = QRAIN_2000_500m.assign_coords({'west_east' : (QRAIN_2000_500m.west_east.values/2 +X3[0]+1)})
QRAIN_2000_500m = QRAIN_2000_500m.assign_coords({'south_north' : (QRAIN_2000_500m.south_north.values/2 +X3[0]+1)})
RAINNC_500m = RAINNC_500m.assign_coords({'west_east' : (RAINNC_500m.west_east.values/2 +X3[0]+1)})
RAINNC_500m = RAINNC_500m.assign_coords({'south_north' : (RAINNC_500m.south_north.values/2 +X3[0]+1)})

#%%
(RAINNC_500m[7]-RAINNC_500m[6]).plot.contourf(levels=np.linspace(5,30,11),cmap='BuPu',add_colorbar=False)
#%% making the plot
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(12,18))
plt.rcParams.update({'font.size': 20})
gs = gridspec.GridSpec(6, 4)

ax1 = plt.subplot(gs[:2, :2])
plot=(RAINNC_1KM[6]-RAINNC_1KM[5]).plot.contourf(ax=ax1,levels=np.linspace(5,35,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(ax=ax1,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'a)', style='italic', fontsize=22)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_200m: Res=1km', y =1.05,fontsize=20)
plt.hlines(592,525,675, linestyle='dashed',color = 'black')
plt.hlines(608,525,675, linestyle='dashed',color = 'black')

ax2 = plt.subplot(gs[:2, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=(RAINNC_200m[6]-RAINNC_200m[5]).plot.contourf(ax=ax2,levels=np.linspace(5,35,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(ax=ax2,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_200m: Res=200m', y =1.05,fontsize=20)
plt.text(529., 648, 'b)', style='italic', fontsize=22)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.hlines(592,525,675, linestyle='dashed',color = 'black')
plt.hlines(608,525,675, linestyle='dashed',color = 'black')

ax3 = plt.subplot(gs[2:4, :2])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=(RAINNC_1KM_500m[7]-RAINNC_1KM_500m[6]).plot.contourf(ax=ax3,levels=np.linspace(5,35,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(ax=ax3,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'c)', style='italic', fontsize=22)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_500m: Res=1km', y =1.05,fontsize=20)
plt.hlines(595,525,675, linestyle='dashed',color = 'black')

ax4 = plt.subplot(gs[2:4, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=(RAINNC_500m[7]-RAINNC_500m[6]).plot.contourf(ax=ax4,levels=np.linspace(5,35,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(ax=ax4,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_500m: Res=500m', y =1.05,fontsize=20)
plt.text(529., 648, 'd)', style='italic', fontsize=22)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.hlines(595,525,675, linestyle='dashed',color = 'black')

ax5 = plt.subplot(gs[4:6, 1:3])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000_200m[2].plot.contourf(ax=ax5,levels=levels,cmap='Blues',add_colorbar=False)
QRAIN_2000_200m[2].plot.contour(ax=ax5,levels=levels,colors='black',alpha=0.3)
HGT.plot.contour(ax=ax5,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_200m: Res=200m', y =1.05,fontsize=20)
plt.text(529., 648, 'e)', style='italic', fontsize=22)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.14, 0.15, 0.72, 0.008])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg/kg', labelpad=0, x=0.52, rotation=0,fontsize=20)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.14, 0.42, 0.72, 0.008])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
cb.update_ticks()
cb.set_label('mm/h', labelpad=-60, x=0.52, rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.82, 
                    hspace=3.8)

#%% confronto updrafts SIM_CTRL per fare vedere che sono roll vortices:
path_data= os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/CTRL_SIMULATION/")
SIM_CTRL = Dataset(path_data + '/SIM_wrfinput_prova_dominio2')

z = getvar(SIM_CTRL,'z')
HGT = getvar(SIM_CTRL,'HGT')

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
    
# making the second cross section to make the comparison
start_point = CoordPair(x=90, y=95)
end_point = CoordPair(x=135,y=95)

w_cross_2 = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross_2 = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
u_cross_2 = vertcross(u,z,start_point=start_point,end_point=end_point)
# make a copy of the w cross section
w_cross_filled_2 = np.ma.copy(to_np(w_cross_2))
u_cross_filled_2 = np.ma.copy(to_np(u_cross_2))
q_cross_filled_2 = np.ma.copy(to_np(qcloud_cross_2))
ter_line_2 = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)


for i in range(w_cross_filled_2.shape[-1]):
    column_vals = w_cross_filled_2[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    w_cross_filled_2[0:first_idx, i] = w_cross_filled_2[first_idx, i]

for i in range(u_cross_filled_2.shape[-1]):
    column_vals = u_cross_filled_2[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    u_cross_filled_2[0:first_idx, i] = u_cross_filled_2[first_idx, i]

for i in range(q_cross_filled_2.shape[-1]):
    column_vals = q_cross_filled_2[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    q_cross_filled_2[0:first_idx, i] = q_cross_filled_2[first_idx, i]
    
xs = np.arange(0, w_cross.shape[-1],1)+X2[0]
ys = to_np(w_cross.coords["vertical"])

w_modified_2 = np.zeros((100,46))
w_modified_2[:] = np.nan
w_modified_2[::2,::]=w_cross_filled_2[::2,::]

u_modified_2 = np.zeros((100,46))
u_modified_2[:] = np.nan
u_modified_2[::2,::]=u_cross_filled_2[::2,::]

#%% making the plot    
plt.rcParams.update({'font.size': 22})
fig, axs = plt.subplots(1,2,figsize=(14, 6))

wspd_1 = np.sqrt(w_cross**2+u_cross**2)
norm = DivergingNorm(vmin=w_cross_filled.min(),vcenter=0, vmax=w_cross_filled.max())
levs = [10**(-6),10**(-4)]
plot = axs[0].contourf(xs+90,ys,w_cross_filled,norm=norm,levels=30, cmap='RdBu_r')
axs[0].quiver(xs+90,ys,u_modified/2,w_modified,scale=95.,color='black',alpha=0.7,width=0.0025)
axs[0].set_ylim((0,9000))
axs[0].set_xlabel('x [km]',fontsize = 22)
axs[0].set_ylabel('z [m]',fontsize = 22)
axs[0].contour(xs+90,ys,q_cross_filled,levels=levs,linestyles='dashed',linewidths=1.8,colors='black',alpha=0.8)
axs[0].fill_between(xs+90,0,to_np(ter_line),facecolor="black",alpha=0.93)
axs[0].text(578,8200,'a)',style='italic', fontsize=26)

plot = axs[1].contourf(xs+90,ys,w_cross_filled_2,norm=norm,levels=30, cmap='RdBu_r')
axs[1].quiver(xs+90,ys,u_modified_2/2,w_modified_2,scale=95.,color='black',alpha=0.7,width=0.0025)
axs[1].set_ylim((0,9000))
axs[1].set_xlabel('x [km]',fontsize = 22)
axs[1].set_ylabel('z [m]',fontsize = 22)
axs[1].contour(xs+90,ys,q_cross_filled_2,levels=levs,linestyles='dashed',linewidths=1.8,colors='black',alpha=0.8)
axs[1].fill_between(xs+90,0,to_np(ter_line_2),facecolor="black",alpha=0.93)
axs[1].text(578,8200,'b)',style='italic', fontsize=26)

from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.2, 0.03, 0.6, 0.025])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=5)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('m/s', labelpad=5, y=0.12, rotation=0,fontsize=22)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.35, 
                    hspace=0.6)


#%% PLOT CUMULATA SIM VELOCITà IDEALIZZATE E WIND PROFILE SIM_SHEAR
path_data= os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_WIND/")
path_data2 = os.path.abspath('C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_SHEAR/')
SIM_U10 = Dataset(path_data + '/SIM_U10_dominio2')
SIM_U20 = Dataset(path_data + '/SIM_U20_dominio2')
SIM_U27 = Dataset(path_data + '/SIM_U25_dominio2')
SIM_shear = Dataset(path_data2 + '/SIM_SHEAR_FINALE_dominio2')

z = getvar(SIM_U10,'z')
HGT = getvar(SIM_U10,'HGT')
RAINNC_SIM_U10 = getvar(SIM_U10,'RAINNC',timeidx=ALL_TIMES)
RAINNC_SIM_U20 = getvar(SIM_U20,'RAINNC',timeidx=ALL_TIMES)
RAINNC_SIM_U27 = getvar(SIM_U27,'RAINNC',timeidx=ALL_TIMES)
RAINNC_SIM_shear = getvar(SIM_shear,'RAINNC',timeidx=ALL_TIMES)

# modifying coordinates
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
RAINNC_SIM_U10 = RAINNC_SIM_U10.assign_coords({'west_east' : RAINNC_SIM_U10.west_east.values +X2[0]+1})
RAINNC_SIM_U10 = RAINNC_SIM_U10.assign_coords({'south_north' : RAINNC_SIM_U10.south_north.values +X2[0]+1})
RAINNC_SIM_U20 = RAINNC_SIM_U20.assign_coords({'west_east' : RAINNC_SIM_U20.west_east.values +X2[0]+1})
RAINNC_SIM_U20 = RAINNC_SIM_U20.assign_coords({'south_north' : RAINNC_SIM_U20.south_north.values +X2[0]+1})
RAINNC_SIM_U27 = RAINNC_SIM_U27.assign_coords({'west_east' : RAINNC_SIM_U27.west_east.values +X2[0]+1})
RAINNC_SIM_U27 = RAINNC_SIM_U27.assign_coords({'south_north' : RAINNC_SIM_U27.south_north.values +X2[0]+1})
RAINNC_SIM_shear = RAINNC_SIM_shear.assign_coords({'west_east' : RAINNC_SIM_shear.west_east.values +X2[0]+1})
RAINNC_SIM_shear = RAINNC_SIM_shear.assign_coords({'south_north' : RAINNC_SIM_shear.south_north.values +X2[0]+1})
#%%
# making the plot of the wind profile:
data = pd.read_csv('C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_Wind_Profile/CTRL_SIMULATION/Parameters_CTRL_sounding.txt',delimiter = '\t')
#v_shear = data['v']
v = data['v']
fig = plt.figure(figsize=(12,10))
plt.rcParams.update({'font.size': 40})
plt.plot(v,data.Height,linewidth=4)
wind_shear1 = (20-10)/300
wind_shear2 = (30-20)/1200
wind_shear3= (45-30)/2500
v_shear = np.zeros(len(data['v']))
v_shear[0]=10

for i in range(1,len(v)):
    if(data['Height'][i]<400):
        v_shear[i]=v_shear[i-1]+wind_shear1*(data['Height'][i]-data['Height'][i-1])
    elif(data['Height'][i]<1200):
        v_shear[i]=v_shear[i-1]+wind_shear2*(data['Height'][i]-data['Height'][i-1])
    elif(data['Height'][i]<2500):
        v_shear[i]=v_shear[i-1]
    elif(data['Height'][i]<5000):
        v_shear[i] = v_shear[i-1]+wind_shear3*0.87*(data['Height'][i]-data['Height'][i-1])
    elif(data['Height'][i]<12000):
        v_shear[i]=v_shear[i-1]
    else:    
        v_shear[i] = 20
        
# plt.plot(v,data.Height,linewidth=2) 
plt.plot(v_shear,data.Height,linewidth=4)  
plt.ylabel('Height',fontsize=40)
plt.xlabel('V (m/s)',fontsize=40)
plt.ylim((0,15000))
plt.xlim((5,45))
plt.grid() 

#%% making the plot
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(12,18))
plt.rcParams.update({'font.size': 20})
gs = gridspec.GridSpec(6, 4)

ax1 = plt.subplot(gs[:2, :2])
plot=(RAINNC_SIM_U10[6]-RAINNC_SIM_U10[5]).plot.contourf(ax=ax1,levels=np.linspace(0,30,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(ax=ax1,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'a)', style='italic', fontsize=22)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_V10', y =1.05,fontsize=20)

ax2 = plt.subplot(gs[:2, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=(RAINNC_SIM_U20[6]-RAINNC_SIM_U20[5]).plot.contourf(ax=ax2,levels=np.linspace(0,30,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(ax=ax2,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_V20', y =1.05,fontsize=20)
plt.text(529., 648, 'b)', style='italic', fontsize=22)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax3 = plt.subplot(gs[2:4, :2])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=(RAINNC_SIM_U27[6]-RAINNC_SIM_U27[5]).plot.contourf(ax=ax3,levels=np.linspace(0,30,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(ax=ax3,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'c)', style='italic', fontsize=22)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_V25', y =1.05,fontsize=20)

ax4 = plt.subplot(gs[2:4, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=(RAINNC_SIM_shear[7]-RAINNC_SIM_shear[6]).plot.contourf(ax=ax4,levels=np.linspace(0,30,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(ax=ax4,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('SIM_V10_SHEAR', y =1.05,fontsize=20)
plt.text(529., 648, 'd)', style='italic', fontsize=22)
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax5 = plt.subplot(gs[4:6, 1:3])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plt.text(6, 13000, 'e)', style='italic', fontsize=22)
plt.plot(v,data.Height,linewidth=2)
plt.plot(v_shear,data.Height,linewidth=2)  
plt.ylabel('z [m]')
plt.xlabel('Wind [m/s]')
plt.ylim((0,15000))
plt.xlim((5,45))
plt.grid() 


fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.14, 0.41, 0.72, 0.008])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
cb.update_ticks()
cb.set_label('mm/h', labelpad=0, x=0.52, rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.82, 
                    hspace=3.8)

#%% PLOT IDEAL WIND DIRECTION SECTION
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_WIND_DIRECTION/")
SIM_210 = Dataset(path_data + '/SIM_wind_210degrees_dominio2')
SIM_UDINE = Dataset(path_data + '/SIM_original_wind20degrees_dominio2')
SIM_UDINE_LONG_RIDGE = Dataset(path_data + '/SIM_original_wind_sum20degrees_dominio2')

# extracting variables for SIM_210
z = getvar(SIM_210, 'height_agl',units="m")
HGT = getvar(SIM_210,'HGT')
u_SIM_210 = getvar(SIM_210,'ua',timeidx=ALL_TIMES)
v_SIM_210 = getvar(SIM_210,'va',timeidx=ALL_TIMES)
wspd_SIM_210 = getvar(SIM_210,'wspd',timeidx=ALL_TIMES)
u_SIM_210_200 = interplevel(u_SIM_210, z, 200)
v_SIM_210_200 = interplevel(v_SIM_210, z,200)
wspd_SIM_210_200 = interplevel(wspd_SIM_210, z, 200)
RAINNC_SIM_210 = getvar(SIM_210,'RAINNC',timeidx=ALL_TIMES)
# extracting variables for SIM_UDINE
u_SIM_UDINE = getvar(SIM_UDINE,'ua',timeidx=ALL_TIMES)
v_SIM_UDINE = getvar(SIM_UDINE,'va',timeidx=ALL_TIMES)
wspd_SIM_UDINE = getvar(SIM_UDINE,'wspd',timeidx=ALL_TIMES)
u_SIM_UDINE_200 = interplevel(u_SIM_UDINE, z, 200)
v_SIM_UDINE_200 = interplevel(v_SIM_UDINE, z,200)
wspd_SIM_UDINE_200 = interplevel(wspd_SIM_UDINE, z, 200)
RAINNC_SIM_UDINE = getvar(SIM_UDINE,'RAINNC',timeidx=ALL_TIMES)

# extracting variables for long ridge
QRAIN_SIM_UDINE_LONG_RIDGE = getvar(SIM_UDINE_LONG_RIDGE,'QRAIN',timeidx=ALL_TIMES)
HGT_LONG = getvar(SIM_UDINE_LONG_RIDGE,'HGT')
z_LONG = getvar(SIM_UDINE_LONG_RIDGE,'z')
QRAIN_SIM_UDINE_LONG_RIDGE_2000 = interplevel(QRAIN_SIM_UDINE_LONG_RIDGE,z_LONG,2000)

#%%
# changing coordinates fro SIM_210
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
u_SIM_210_200 = u_SIM_210_200.assign_coords({'west_east' : u_SIM_210_200.west_east.values +X2[0]+1})
u_SIM_210_200 = u_SIM_210_200.assign_coords({'south_north' : u_SIM_210_200.south_north.values +X2[0]+1})
v_SIM_210_200 = v_SIM_210_200.assign_coords({'west_east' : v_SIM_210_200.west_east.values +X2[0]+1})
v_SIM_210_200 = v_SIM_210_200.assign_coords({'south_north' : v_SIM_210_200.south_north.values +X2[0]+1})
wspd_SIM_210_200 = wspd_SIM_210_200.assign_coords({'west_east' : wspd_SIM_210_200.west_east.values +X2[0]+1})
wspd_SIM_210_200 = wspd_SIM_210_200.assign_coords({'south_north' : wspd_SIM_210_200.south_north.values +X2[0]+1})
RAINNC_SIM_210 = RAINNC_SIM_210.assign_coords({'west_east' : RAINNC_SIM_210.west_east.values +X2[0]+1})
RAINNC_SIM_210 = RAINNC_SIM_210.assign_coords({'south_north' : RAINNC_SIM_210.south_north.values +X2[0]+1})

#changing the coordinate for SIM_UDINE
u_SIM_UDINE_200 = u_SIM_UDINE_200.assign_coords({'west_east' : u_SIM_UDINE_200.west_east.values +X2[0]+1})
u_SIM_UDINE_200 = u_SIM_UDINE_200.assign_coords({'south_north' : u_SIM_UDINE_200.south_north.values +X2[0]+1})
v_SIM_UDINE_200 = v_SIM_UDINE_200.assign_coords({'west_east' : v_SIM_UDINE_200.west_east.values +X2[0]+1})
v_SIM_UDINE_200 = v_SIM_UDINE_200.assign_coords({'south_north' : v_SIM_UDINE_200.south_north.values +X2[0]+1})
wspd_SIM_UDINE_200 = wspd_SIM_UDINE_200.assign_coords({'west_east' : wspd_SIM_UDINE_200.west_east.values +X2[0]+1})
wspd_SIM_UDINE_200 = wspd_SIM_UDINE_200.assign_coords({'south_north' : wspd_SIM_UDINE_200.south_north.values +X2[0]+1})
RAINNC_SIM_UDINE = RAINNC_SIM_UDINE.assign_coords({'west_east' : RAINNC_SIM_UDINE.west_east.values +X2[0]+1})
RAINNC_SIM_UDINE = RAINNC_SIM_UDINE.assign_coords({'south_north' : RAINNC_SIM_UDINE.south_north.values +X2[0]+1})

# changing coordinates for sim_long_ridge
HGT_LONG = HGT_LONG.assign_coords({'west_east' : HGT_LONG.west_east.values +X2_long[0]})
HGT_LONG = HGT_LONG.assign_coords({'south_north' : HGT_LONG.south_north.values +X2_long[0]})
QRAIN_SIM_UDINE_LONG_RIDGE_2000 = QRAIN_SIM_UDINE_LONG_RIDGE_2000.assign_coords({'west_east' : QRAIN_SIM_UDINE_LONG_RIDGE_2000.west_east.values +X2_long[0]})
QRAIN_SIM_UDINE_LONG_RIDGE_2000 = QRAIN_SIM_UDINE_LONG_RIDGE_2000.assign_coords({'south_north' : QRAIN_SIM_UDINE_LONG_RIDGE_2000.south_north.values +X2_long[0]})

#%% making the plot

fig = plt.figure(figsize=(14,10))
plt.rcParams.update({'font.size': 34})

ax = plt.subplot(1, 1, 1)
levels =[100,500,1000,1400,1499]
HGT.plot.contour(levels=levels,colors=['black'],alpha=0.5,x = 'west_east',y = 'south_north')
plot=wspd_SIM_210_200[3].plot.contourf(levels=np.linspace(0,20,41),cmap='bwr',add_colorbar=False)
qv = ax.quiver(u_SIM_210_200[3,::12,::12].west_east, u_SIM_210_200[3,::12,::12].south_north, u_SIM_210_200[3,::12,::12], v_SIM_210_200[3,::12,::12],color='black',scale=250.)
plt.title('South-north speed 20 m above ground hour 5 SIM_V5',fontsize =34)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_210°', y =1.05,fontsize=34)
plt.text(493., 692, 'a)', style='italic', fontsize=34)

'''
ax = plt.subplot(3, 2, 2)
levels =[100,500,1000,1400,1499]
HGT.plot.contour(levels=levels,colors=['black'],alpha=0.5,x = 'west_east',y = 'south_north')
plot=wspd_SIM_UDINE_200[3].plot.contourf(levels=np.linspace(0,20,41),cmap='bwr',add_colorbar=False)
qv = ax.quiver(u_SIM_UDINE_200[3,::10,::10].west_east, u_SIM_UDINE_200[3,::10,::10].south_north, u_SIM_UDINE_200[3,::10,::10], v_SIM_UDINE_200[3,::10,::10],color='black',scale=250.)
plt.title('SIM_UDINE_ROT20',y =1.05,fontsize=26)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.text(493., 692, 'b)', style='italic', fontsize=26)

ax = plt.subplot(3, 2, 3)
plot2=(RAINNC_SIM_210[3]-RAINNC_SIM_210[2]).plot.contourf(levels=np.linspace(0,20,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'c)', style='italic', fontsize=26)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_210°', y =1.05,fontsize=26)

ax = plt.subplot(3, 2, 4)
plot2=(RAINNC_SIM_UDINE[3]-RAINNC_SIM_UDINE[2]).plot.contourf(levels=np.linspace(0,20,11),cmap='BuPu',add_colorbar=False)
HGT.plot.contour(levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'd)', style='italic', fontsize=26)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_UDINE_ROT20',y =1.05,fontsize=26)

ax = plt.subplot(3, 2, 5)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot3=QRAIN_SIM_UDINE_LONG_RIDGE_2000[3].plot.contourf(levels=levels,cmap='Blues',add_colorbar=False)
QRAIN_SIM_UDINE_LONG_RIDGE_2000[3].plot.contour(levels=levels,colors='black',alpha=0.5)
HGT_LONG.plot.contour(levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(510., 648, 'e)', style='italic', fontsize=26)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((500,700))
plt.ylim((535,665))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_UDINE_ROT20_LONG',y =1.05,fontsize=26)

ax = plt.subplot(3, 2, 6)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot3=QRAIN_SIM_UDINE_LONG_RIDGE_2000[6].plot.contourf(levels=levels,cmap='Blues',add_colorbar=False)
QRAIN_SIM_UDINE_LONG_RIDGE_2000[6].plot.contour(levels=levels,colors='black',alpha=0.5)
HGT_LONG.plot.contour(levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(510., 648, 'f)', style='italic', fontsize=26)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((500,700))
plt.ylim((535,665))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_UDINE_ROT20_LONG',y =1.05,fontsize=26)


'''
from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.0000017, 0.72, 0.025])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
cb.set_label('m/s', labelpad=0, y=0.92,x=0.53, rotation=0,fontsize=34)

'''
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.342, 0.72, 0.012])
cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
cb.set_label('mm/h', labelpad=0, y=0.92,x=0.53, rotation=0,fontsize=24)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.13, 0.03, 0.72, 0.012])
cb = fig.colorbar(plot3,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
#cb.update_ticks()
cb.set_label('kg/kg', labelpad=0, y=0.95, x=0.53, rotation=0,fontsize=24)
'''

plt.subplots_adjust(left=0.1,
                    bottom=0.15, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.37, 
                    hspace=0.75)



#%% PLOT E CALCOLO VELOCITà ROLL VORTICES SIM IDEAL T_SURF
path_data= os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/CTRL_SIMULATION/")
SIM_CTRL = Dataset(path_data + '/SIM_wrfinput_prova_dominio2')
path_data = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_T_QSNOW/")
SIM_T290 = Dataset(path_data+'/SIM_T290_dominio2')
SIM_T300 = Dataset(path_data+'/SIM_T300_dominio2')

# Extracting orography
HGT = getvar(SIM_CTRL,"HGT")
QCLOUD = getvar(SIM_CTRL,"QCLOUD",timeidx=5)
QCLOUD_T300 = getvar(SIM_T300,"QCLOUD",timeidx=5)
QCLOUD_T290 = getvar(SIM_T290,"QCLOUD",timeidx=8)
z = getvar(SIM_CTRL, "z",units="m")
w =getvar(SIM_CTRL,"wa",timeidx=5)
w_T300 =getvar(SIM_T300,"wa",timeidx=5)
w_T290 =getvar(SIM_T290,"wa",timeidx=8)
u = getvar(SIM_CTRL, "ua",timeidx=5)
u_T300 = getvar(SIM_T300, "ua",timeidx=5)
u_T290 = getvar(SIM_T290, "ua",timeidx=8)

start_point = CoordPair(x=75, y=108)
end_point = CoordPair(x=130,y=108)

#getting cross section values
w_cross = vertcross(w,z,start_point=start_point,end_point=end_point)
w_cross_T290 = vertcross(w_T290,z,start_point=start_point,end_point=end_point)
w_cross_T300 = vertcross(w_T300,z,start_point=start_point,end_point=end_point)
qcloud_cross = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
qcloud_cross_T290 = vertcross(QCLOUD_T290,z,start_point=start_point,end_point=end_point)
qcloud_cross_T300 = vertcross(QCLOUD_T300,z,start_point=start_point,end_point=end_point)
u_cross = vertcross(u,z,start_point=start_point,end_point=end_point)
u_cross_T290 = vertcross(u_T290,z,start_point=start_point,end_point=end_point)
u_cross_T300 = vertcross(u_T300,z,start_point=start_point,end_point=end_point)

xs = np.arange(0, w_cross.shape[-1],1)+X2[0]
ys = to_np(w_cross.coords["vertical"])

# Get the terrain heights along the cross section line
ter_line = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)

#%%
fig = plt.figure(figsize=(14,9))
plt.rcParams.update({'font.size': 18})

wspd = np.sqrt(w_cross_T290**2+u_cross_T290**2)
ax = plt.subplot(2, 2, 1)
norm = DivergingNorm(vmin=w_cross_T300.min(),vcenter=0, vmax=w_cross_T300.max())
levs = [10**(-6),10**(-4)]
plot = ax.contourf(xs+75,ys,w_cross_T290,norm=norm,levels=30, cmap='RdBu_r')
ax.quiver(xs[::2]+75,ys[::3],u_cross_T290[::3,::2]/(4*wspd[::3,::2]),w_cross_T290[::3,::2]/wspd[::3,::2],scale=25.,color='black',alpha=0.8,width=0.003)
ax.set_ylim((0,10000))
ax.set_xlabel('x [km]',fontsize = 22)
ax.set_ylabel('z [m]',fontsize = 22)
ax.contour(xs+75,ys,qcloud_cross_T290,levels=levs,linestyles='dashed',linewidths=1.4,colors='black',alpha=0.8)
ax.fill_between(xs+75,0,to_np(ter_line),facecolor="black",alpha=0.93)
ax.text(563,9200,'a)',style='italic', fontsize=22)
plt.title('SIM_T290',y=1.03)

wspd = np.sqrt(w_cross**2+u_cross**2)
ax = plt.subplot(2, 2, 2)
levs = [10**(-6),10**(-4)]
plot = ax.contourf(xs+75,ys,w_cross,norm=norm,levels=30, cmap='RdBu_r')
#cb = fig.colorbar(plot,orientation='horizontal')
#cb.ax.tick_params(labelsize=18)
#cb.set_label('m/s', labelpad=-40, y=1.05, rotation=0,fontsize=18)
ax.quiver(xs[::2]+75,ys[::3],u_cross[::3,::2]/(4*wspd[::3,::2]),w_cross[::3,::2]/wspd[::3,::2],scale=25.,color='black',alpha=0.8,width=0.003)
ax.set_ylim((0,10000))
ax.set_xlabel('x [km]',fontsize = 22)
ax.set_ylabel('z [m]',fontsize = 22)
ax.contour(xs+75,ys,qcloud_cross,levels=levs,linestyles='dashed',linewidths=1.4,colors='black',alpha=0.8)
ax.fill_between(xs+75,0,to_np(ter_line),facecolor="black",alpha=0.93)
ax.text(563,9200,'b)',style='italic', fontsize=22)
plt.title('SIM_CTRL',y=1.03)

wspd = np.sqrt(w_cross_T300**2+u_cross_T300**2)
ax = plt.subplot(2, 2, 3)
levs = [10**(-6),10**(-4)]
plot = ax.contourf(xs+75,ys,w_cross_T300,norm=norm,levels=30, cmap='RdBu_r')
#cb = fig.colorbar(plot,orientation='horizontal')
#cb.ax.tick_params(labelsize=18)
#cb.set_label('m/s', labelpad=-40, y=1.05, rotation=0,fontsize=18)
ax.quiver(xs[::2]+75,ys[::3],u_cross_T300[::3,::2]/(4*wspd[::3,::2]),w_cross_T300[::3,::2]/(wspd[::3,::2]),scale=25.,color='black',alpha=0.8,width=0.003)
ax.set_ylim((0,12000))
ax.set_xlabel('x [km]',fontsize = 22)
ax.set_ylabel('z [m]',fontsize = 22)
ax.contour(xs+75,ys,qcloud_cross_T300,levels=levs,linestyles='dashed',linewidths=1.4,colors='black',alpha=0.8)
ax.fill_between(xs+75,0,to_np(ter_line),facecolor="black",alpha=0.93)
ax.text(563,11200,'c)',style='italic', fontsize=24)
plt.title('SIM_T300',y=1.03)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.5, 0., 0.02, 0.4])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='vertical')
cb.set_label('m/s', labelpad=-50, y=1.12, rotation=0,fontsize=22)

plt.subplots_adjust(left=0.1,
                    bottom=0., 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.45, 
                    hspace=0.4)

#%% PLOT QRAIN ORA 5 PER SIM_N1_000001 + Moist-Brunt_N1_00004
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_N/SIM/")
SIM_N1_000001 = Dataset(path_data + '/SIM_N1_000001_dominio2')
SIM_N1_00004 = Dataset(path_data + '/SIM_N1_00004_dominio2')

z = getvar(SIM_N1_000001,'z')
HGT = getvar(SIM_N1_000001,'HGT')
QRAIN_N1_000001 = getvar(SIM_N1_000001,'QRAIN',timeidx=4)
QRAIN_N1_000001_2000 = interplevel(QRAIN_N1_000001,z,2000)
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_N1_000001_2000 = QRAIN_N1_000001_2000.assign_coords({'west_east' : QRAIN_N1_000001_2000.west_east.values +X2[0]+1})
QRAIN_N1_000001_2000 = QRAIN_N1_000001_2000.assign_coords({'south_north' : QRAIN_N1_000001_2000.south_north.values +X2[0]+1})

# EXTRACTING VARIABLES NEEDED ALONG CROSS SECTION
Theta = getvar(SIM_N1_00004,"theta",timeidx=1)
T = getvar(SIM_N1_00004,"tk",timeidx=1)
z = getvar(SIM_N1_00004,"z", timeidx=1)
RH = getvar(SIM_N1_00004,"rh",timeidx=1)
QVAPOR = getvar(SIM_N1_00004,"QVAPOR",timeidx=1)
QRAIN = getvar(SIM_N1_00004,"QRAIN",timeidx=1)
QCLOUD = getvar(SIM_N1_00004,"QCLOUD",timeidx=1)
QGRAUP = getvar(SIM_N1_00004,"QGRAUP",timeidx=1)
QICE = getvar(SIM_N1_00004,"QICE",timeidx=1)
#QSNOW = getvar(ncfile,"QSNOW",timeidx=1)
Theta_e = getvar(SIM_N1_00004,"eth",timeidx=1)
HGT_N1_00004 = getvar(SIM_N1_00004,"HGT")

# computation of the approximated moist Brunt-Vaisala frequency
g = 9.81
L= 2501000 # J/kg latent heat of vaporization
R = 287 # dry air gas constant
Rw = 461.5
#R = 8.31
eps = 287/461.5 # ratio between gas constant of dry and moist air
cp = 1004
#per usare la fromula di DURRAN che usa anche MIglietta ti devi ricavare il saturated mixing ratio e il total mixing ratio dalla cross section
q_s = 100*QVAPOR/RH

# applichiamola sulla cross section
start_point = CoordPair(x=113,y=50)
end_point = CoordPair(x=113, y=125)

ter_line = interpline(HGT_N1_00004,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables
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

xs = np.arange(0, thetae_cross.shape[-1], 1)
ys = to_np(thetae_cross.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

dz = np.diff(Ys,axis=0)
dq_s_dz = q_s_cross.differentiate('vertical')
log_Theta = np.log(theta_cross)
dTheta_dz = log_Theta.differentiate('vertical')
dq_w_dz = q_w_cross.differentiate('vertical')
Nm_squared = g*((1+(L*q_s_cross/(R*T_cross)))/(1+((eps*q_s_cross*L**2)/(cp*R*T_cross**2))))*(dTheta_dz+L/(cp*T_cross)*(dq_s_dz)) - g*dq_w_dz

vertical = Nm_squared.vertical.values
horizontal = Nm_squared.cross_line_idx.values
NM_squared = xr.DataArray(Nm_squared, coords=[vertical, horizontal])
TER_LINE = xr.DataArray(ter_line, coords=[horizontal])

z_coord = np.linspace(0,NM_squared.vertical.max(),1000)
y_coord= np.linspace(NM_squared.cross_line_idx.min(),NM_squared.cross_line_idx.max(),1000)
NM_squared_interp = NM_squared.interp(vertical = z_coord,cross_line_idx = y_coord)
ter_line_interp = TER_LINE.interp(line_idx = y_coord)

QCROSS = xr.DataArray(qcloud_cross, coords=[vertical, horizontal])
QCROSS = QCROSS.interp(vertical = z_coord,cross_line_idx = y_coord)

xs = to_np(NM_squared_interp.cross_line_idx)
ys = to_np(NM_squared_interp.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

NM_squared_interp = NM_squared_interp.where(QCROSS>10**(-5))


#%% making the plot
fig = plt.figure(figsize=(28,10))
plt.rcParams.update({'font.size': 32})
ax = plt.subplot(2, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N1_000001_2000.plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
cb = plt.colorbar(plot,shrink=0.9,orientation='vertical')
cb.set_label('kg/kg',labelpad=-50, y=1.12, rotation=0)
QRAIN_N1_000001_2000.plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT.plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax.text(495, 695, 'a)', style='italic',fontsize=40)
plt.title('SIM_N1_000001, t=4h',fontsize=38,y=1.04)
plt.xlabel('x [km]')
plt.ylabel('y [km]')


ax = plt.subplot(2, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=ax.contourf(Xs+Y2[0]+50,Ys,NM_squared_interp*10000,levels=np.linspace(-2,2,41),cmap='RdBu_r')
cb = plt.colorbar(plot,shrink=0.9,orientation='vertical')
ax.fill_between(xs+Y2[0]+50,0,to_np(ter_line_interp),facecolor="black")
ax.contour(Xs+Y2[0]+50,Ys,QCROSS,levels=[10**(-5)],colors='Gray')
plt.title('SIM_N1_00004',fontsize=38,y=1.04)
plt.ylabel('z [m]')
plt.xlabel('y [km]')
plt.ylim((0,7000))
plt.text(538., 6300, 'b)', style='italic', fontsize=40)
cb.set_label('$10^{-4} \ s^{-2}$',labelpad=-50, y=1.12, rotation=0)

plt.subplots_adjust(left=0.1,
                    bottom=-0.6, 
                    right=0.9, 
                    top=0.9, 
                    wspace=0.25, 
                    hspace=0.15)

#%% plot moist Brunt e cloud SIM_N1_00015
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_N/SIM/")
SIM_N1_00015 = Dataset(path_data + '/SIM_N1_00015_dominio2')

#Moist Brunt cross section:
Theta = getvar(SIM_N1_00015,"theta",timeidx=1)
T = getvar(SIM_N1_00015,"tk",timeidx=1)
z = getvar(SIM_N1_00015,"z", timeidx=1)
RH = getvar(SIM_N1_00015,"rh",timeidx=1)
QVAPOR = getvar(SIM_N1_00015,"QVAPOR",timeidx=1)
QRAIN = getvar(SIM_N1_00015,"QRAIN",timeidx=1)
QCLOUD = getvar(SIM_N1_00015,"QCLOUD",timeidx=1)
QGRAUP = getvar(SIM_N1_00015,"QGRAUP",timeidx=1)
QICE = getvar(SIM_N1_00015,"QICE",timeidx=1)
#QSNOW = getvar(ncfile,"QSNOW",timeidx=1)
Theta_e = getvar(SIM_N1_00015,"eth",timeidx=1)
HGT_N1_00004 = getvar(SIM_N1_00015,"HGT")

# computation of the approximated moist Brunt-Vaisala frequency
g = 9.81
L= 2501000 # J/kg latent heat of vaporization
R = 287 # dry air gas constant
Rw = 461.5
#R = 8.31
eps = 287/461.5 # ratio between gas constant of dry and moist air
cp = 1004
#per usare la fromula di DURRAN che usa anche MIglietta ti devi ricavare il saturated mixing ratio e il total mixing ratio dalla cross section
q_s = 100*QVAPOR/RH

# applichiamola sulla cross section
start_point = CoordPair(x=113,y=50)
end_point = CoordPair(x=113, y=125)

ter_line = interpline(HGT_N1_00004,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables
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

xs = np.arange(0, thetae_cross.shape[-1], 1)
ys = to_np(thetae_cross.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

dz = np.diff(Ys,axis=0)
dq_s_dz = q_s_cross.differentiate('vertical')
log_Theta = np.log(theta_cross)
dTheta_dz = log_Theta.differentiate('vertical')
dq_w_dz = q_w_cross.differentiate('vertical')
Nm_squared = g*((1+(L*q_s_cross/(R*T_cross)))/(1+((eps*q_s_cross*L**2)/(cp*R*T_cross**2))))*(dTheta_dz+L/(cp*T_cross)*(dq_s_dz)) - g*dq_w_dz

vertical = Nm_squared.vertical.values
horizontal = Nm_squared.cross_line_idx.values
NM_squared = xr.DataArray(Nm_squared, coords=[vertical, horizontal])
TER_LINE = xr.DataArray(ter_line, coords=[horizontal])

z_coord = np.linspace(0,NM_squared.vertical.max(),1000)
y_coord= np.linspace(NM_squared.cross_line_idx.min(),NM_squared.cross_line_idx.max(),1000)
NM_squared_interp = NM_squared.interp(vertical = z_coord,cross_line_idx = y_coord)
ter_line_interp = TER_LINE.interp(line_idx = y_coord)

QCROSS = xr.DataArray(qcloud_cross, coords=[vertical, horizontal])
QCROSS = QCROSS.interp(vertical = z_coord,cross_line_idx = y_coord)

xs = to_np(NM_squared_interp.cross_line_idx)
ys = to_np(NM_squared_interp.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

NM_squared_interp_N1_00015 = NM_squared_interp.where(QCROSS>10**(-5))
NM_squared_interp_N1_00015 = NM_squared_interp_N1_00015.assign_coords({'cross_line_idx' : NM_squared_interp_N1_00015.cross_line_idx.values+Y2[0]+50})
QCROSS = QCROSS.assign_coords({'cross_line_idx' : QCROSS.cross_line_idx.values +Y2[0]+50})

# Moist Brunt for N1_00004_N3_00012
SIM_N1_00004_N3_00012 = Dataset(path_data + '/SIM_N1_00004_N3_00012_dominio2')

#Moist Brunt cross section:
Theta = getvar(SIM_N1_00004_N3_00012,"theta",timeidx=1)
T = getvar(SIM_N1_00004_N3_00012,"tk",timeidx=1)
z = getvar(SIM_N1_00004_N3_00012,"z", timeidx=1)
RH = getvar(SIM_N1_00004_N3_00012,"rh",timeidx=1)
QVAPOR = getvar(SIM_N1_00004_N3_00012,"QVAPOR",timeidx=1)
QRAIN = getvar(SIM_N1_00004_N3_00012,"QRAIN",timeidx=1)
QCLOUD = getvar(SIM_N1_00004_N3_00012,"QCLOUD",timeidx=1)
QGRAUP = getvar(SIM_N1_00004_N3_00012,"QGRAUP",timeidx=1)
QICE = getvar(SIM_N1_00004_N3_00012,"QICE",timeidx=1)
#QSNOW = getvar(ncfile,"QSNOW",timeidx=1)
Theta_e = getvar(SIM_N1_00004_N3_00012,"eth",timeidx=1)
HGT_N1_00004 = getvar(SIM_N1_00004_N3_00012,"HGT")

# computation of the approximated moist Brunt-Vaisala frequency
g = 9.81
L= 2501000 # J/kg latent heat of vaporization
R = 287 # dry air gas constant
Rw = 461.5
#R = 8.31
eps = 287/461.5 # ratio between gas constant of dry and moist air
cp = 1004
#per usare la fromula di DURRAN che usa anche MIglietta ti devi ricavare il saturated mixing ratio e il total mixing ratio dalla cross section
q_s = 100*QVAPOR/RH

# applichiamola sulla cross section
start_point = CoordPair(x=113,y=50)
end_point = CoordPair(x=113, y=125)

ter_line = interpline(HGT_N1_00004,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables
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

xs = np.arange(0, thetae_cross.shape[-1], 1)
ys = to_np(thetae_cross.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

dz = np.diff(Ys,axis=0)
dq_s_dz = q_s_cross.differentiate('vertical')
log_Theta = np.log(theta_cross)
dTheta_dz = log_Theta.differentiate('vertical')
dq_w_dz = q_w_cross.differentiate('vertical')
Nm_squared = g*((1+(L*q_s_cross/(R*T_cross)))/(1+((eps*q_s_cross*L**2)/(cp*R*T_cross**2))))*(dTheta_dz+L/(cp*T_cross)*(dq_s_dz)) - g*dq_w_dz

vertical = Nm_squared.vertical.values
horizontal = Nm_squared.cross_line_idx.values
NM_squared = xr.DataArray(Nm_squared, coords=[vertical, horizontal])
TER_LINE = xr.DataArray(ter_line, coords=[horizontal])

z_coord = np.linspace(0,NM_squared.vertical.max(),1000)
y_coord= np.linspace(NM_squared.cross_line_idx.min(),NM_squared.cross_line_idx.max(),1000)
NM_squared_interp = NM_squared.interp(vertical = z_coord,cross_line_idx = y_coord)
ter_line_interp = TER_LINE.interp(line_idx = y_coord)

QCROSS_N1_00004_N3_00012 = xr.DataArray(qcloud_cross, coords=[vertical, horizontal])
QCROSS_N1_00004_N3_00012 = QCROSS_N1_00004_N3_00012.interp(vertical = z_coord,cross_line_idx = y_coord)

xs = to_np(NM_squared_interp.cross_line_idx)
ys = to_np(NM_squared_interp.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

NM_squared_interp_N1_00004_N3_00012 = NM_squared_interp.where(QCROSS_N1_00004_N3_00012>10**(-5))


# cross section cloud
z = getvar(SIM_N1_00015, "z",units="m",timeidx=5)
w =getvar(SIM_N1_00015,"wa",timeidx=5)
v = getvar(SIM_N1_00015, "va",timeidx=5)  # default is in m/s
theta = getvar(SIM_N1_00015,"theta",timeidx=5)
thetae = getvar(SIM_N1_00015,"theta_e",timeidx=5)
QCLOUD = getvar(SIM_N1_00015,"QCLOUD",timeidx=5)

start_point = CoordPair(x=95,y=50)
end_point = CoordPair(x=95, y=150)

# Get the terrain heights along the cross section line
ter_line = interpline(HGT_N1_00004,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables
w_cross = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
qice = vertcross(QICE,z,start_point=start_point,end_point=end_point)
v_cross = vertcross(v,z,start_point=start_point,end_point=end_point)
theta_cross = vertcross(theta,z,start_point=start_point,end_point=end_point)
thetae_cross = vertcross(thetae,z,start_point=start_point,end_point=end_point)


q_cross_filled = np.ma.copy(to_np(qcloud_cross))

for i in range(q_cross_filled.shape[-1]):
    column_vals = q_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    q_cross_filled[0:first_idx, i] = q_cross_filled[first_idx, i]


w_modified = np.zeros((100,101))
w_modified[:] = np.nan
w_modified[::2,::5]=w_cross[::2,::5]

v_modified = np.zeros((100,101))
v_modified[:] = np.nan
v_modified[::2,::5]=v_cross[::2,::5]

xs_2 = np.arange(0, w_cross.shape[-1], 1)
ys_2 = to_np(w_cross.coords["vertical"])
Xs_2, Ys_2 = np.meshgrid(xs_2,ys_2)

#%% making the plot
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(24,20))
plt.rcParams.update({'font.size': 32})
gs = gridspec.GridSpec(4, 4)

ax1 = plt.subplot(gs[:2, :2])
plt.ylim((0,8000))
plot=ax1.contourf(xs_2+Y2[0]+50,ys_2,q_cross_filled,levels=10,cmap="Greys")
cb = plt.colorbar(plot,shrink=0.9,orientation='vertical')
cb.set_label('kg/kg', labelpad=-40, y=1.12, rotation=0)
ax1.fill_between(xs_2+Y2[0]+50,0,to_np(ter_line),facecolor="black")
plt.xlabel('y [km]')
ax1.streamplot(Xs_2+Y2[0]+50,Ys_2,v_cross/1000,w_cross,density=2,arrowsize=1.9,linewidth=2.)
plt.ylabel('z [m]')
plt.title('SIM_N1_00015',y=1.03)
plt.text(538., 7300, 'a)', style='italic', fontsize=40)

ax2 = plt.subplot(gs[:2, 2:])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=(NM_squared_interp_N1_00015*10000).plot.contourf(levels=np.linspace(-2,2,41),cmap='RdBu_r',add_colorbar=False)
cb = plt.colorbar(plot2,shrink=0.9,orientation='vertical')
cb.set_label('kg/kg', labelpad=-40, y=1.08, rotation=0)
ax2.fill_between(xs+Y2[0]+50,0,to_np(ter_line_interp),facecolor="black")
ax2.contour(Xs+Y2[0]+50,Ys,QCROSS,levels=[10**(-5)],colors='Gray')
plt.title('SIM_N1_00015',y=1.03)
plt.ylabel('z [m]')
plt.xlabel('y [km]')
plt.ylim((0,7000))
plt.text(538., 6300, 'b)', style='italic', fontsize=40)
cb.set_label('$10^{-4} \ s^{-2}$',labelpad=-60, y=1.12, rotation=0)

ax3 = plt.subplot(gs[2:4, 1:3])
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=ax3.contourf(Xs+Y2[0]+50,Ys,NM_squared_interp_N1_00004_N3_00012*10000,levels=np.linspace(-2,2,41),cmap='RdBu_r')
cb = plt.colorbar(plot2,shrink=0.9,orientation='vertical')
cb.set_label('kg/kg', labelpad=-40, y=1.08, rotation=0)
ax3.fill_between(xs+Y2[0]+50,0,to_np(ter_line_interp),facecolor="black")
ax3.contour(Xs+Y2[0]+50,Ys,QCROSS_N1_00004_N3_00012,levels=[10**(-5)],colors='Gray')
plt.title('SIM_N1_00004_N3_00012',y=1.03)
plt.ylabel('z [m]')
plt.xlabel('y [km]')
plt.ylim((0,7000))
plt.text(538., 6300, 'c)', style='italic', fontsize=40)
cb.set_label('$10^{-4} \ s^{-2}$',labelpad=-60, y=1.12, rotation=0)

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=1., 
                    hspace=1.)



#%% cross secition SIM_RH_INCR5_over2300m
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_RH/")
SIM_RH_increased5_above2300m = Dataset(path_data + '/SIM_RH_increased5_over2300m_dominio2')
SIM_RH_increased5 = Dataset(path_data + '/SIM_RH_increased5_dominio2')
SIM_RH_increased5_under2600m = Dataset(path_data + '/SIM_RH_increased5_under2600m_dominio2')


z = getvar(SIM_RH_increased5,'z')
HGT = getvar(SIM_RH_increased5,'HGT')
QRAIN = getvar(SIM_RH_increased5,'QRAIN',timeidx=ALL_TIMES)
QRAIN_2000 = interplevel(QRAIN,z,2000)
QRAIN_under2600m = getvar(SIM_RH_increased5_under2600m,'QRAIN',timeidx=ALL_TIMES)
QRAIN_under2600m_2000 = interplevel(QRAIN_under2600m,z,2000)
QRAIN_above2300m = getvar(SIM_RH_increased5_above2300m,'QRAIN',timeidx=ALL_TIMES)
QRAIN_above2300m_2000 = interplevel(QRAIN_above2300m,z,2000)

# varying coordinates
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_2000 = QRAIN_2000.assign_coords({'west_east' : QRAIN_2000.west_east.values +X2[0]+1})
QRAIN_2000 = QRAIN_2000.assign_coords({'south_north' : QRAIN_2000.south_north.values +X2[0]+1})
QRAIN_under2600m_2000 = QRAIN_under2600m_2000.assign_coords({'west_east' : QRAIN_under2600m_2000.west_east.values +X2[0]+1})
QRAIN_under2600m_2000 = QRAIN_under2600m_2000.assign_coords({'south_north' : QRAIN_under2600m_2000.south_north.values +X2[0]+1})
QRAIN_above2300m_2000 = QRAIN_above2300m_2000.assign_coords({'west_east' : QRAIN_above2300m_2000.west_east.values +X2[0]+1})
QRAIN_above2300m_2000 = QRAIN_above2300m_2000.assign_coords({'south_north' : QRAIN_above2300m_2000.south_north.values +X2[0]+1})

# ora faccio cross section SIM_RH_INCR5_above2300m
z = getvar(SIM_RH_increased5_above2300m, "z",units="m",timeidx=5)
w =getvar(SIM_RH_increased5_above2300m,"wa",timeidx=5)
v = getvar(SIM_RH_increased5_above2300m, "va",timeidx=5)  # default is in m/s
theta = getvar(SIM_RH_increased5_above2300m,"theta",timeidx=5)
thetae = getvar(SIM_RH_increased5_above2300m,"theta_e",timeidx=5)
QCLOUD = getvar(SIM_RH_increased5_under2600m,"QCLOUD",timeidx=5)

start_point = CoordPair(x=120,y=20)
end_point = CoordPair(x=120, y=150)

# Get the terrain heights along the cross section line
ter_line = interpline(HGT,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables
w_cross = vertcross(w,z,start_point=start_point,end_point=end_point)
qcloud_cross = vertcross(QCLOUD,z,start_point=start_point,end_point=end_point)
v_cross = vertcross(v,z,start_point=start_point,end_point=end_point)
theta_cross = vertcross(theta,z,start_point=start_point,end_point=end_point)
thetae_cross = vertcross(thetae,z,start_point=start_point,end_point=end_point)


q_cross_filled = np.ma.copy(to_np(qcloud_cross))

for i in range(q_cross_filled.shape[-1]):
    column_vals = q_cross_filled[:,i]
    first_idx = int(np.transpose((column_vals > -200).nonzero())[0])
    q_cross_filled[0:first_idx, i] = q_cross_filled[first_idx, i]


w_modified = np.zeros((100,131))
w_modified[:] = np.nan
w_modified[::2,::5]=w_cross[::2,::5]

v_modified = np.zeros((100,131))
v_modified[:] = np.nan
v_modified[::2,::5]=v_cross[::2,::5]

xs = np.arange(0, w_cross.shape[-1], 1)
ys = to_np(w_cross.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

#%% making the plot
fig = plt.figure(figsize=(18,12))
plt.rcParams.update({'font.size': 28})

ax = plt.subplot(2, 2, 1)
levels =[100,500,1000,1400,1499]
HGT.plot.contour(levels=levels,colors=['black'],alpha=0.5,x = 'west_east',y = 'south_north')
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_2000[4].plot.contourf(levels=levels,cmap='Blues',add_colorbar=False)
cb = plt.colorbar(plot,shrink=0.9,orientation='vertical')
cb.set_label('Kg/kg', labelpad=-40, y=1.15, rotation=0,fontsize=26)
QRAIN_2000[4].plot.contour(levels=levels,colors='black',alpha=0.5)
plt.title('South-north speed 20 m above ground hour 5 SIM_V5',fontsize =15)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_RH_INCR5', y =1.05,fontsize=28)
plt.text(493., 695, 'a)', style='italic', fontsize=32)

ax = plt.subplot(2, 2, 2)
levels =[100,500,1000,1400,1499]
HGT.plot.contour(levels=levels,colors=['black'],alpha=0.5,x = 'west_east',y = 'south_north')
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_under2600m_2000[6].plot.contourf(levels=levels,cmap='Blues',add_colorbar=False)
cb = plt.colorbar(plot,shrink=0.9,orientation='vertical')
cb.set_label('Kg/kg', labelpad=-40, y=1.15, rotation=0,fontsize=26)
QRAIN_under2600m_2000[6].plot.contour(levels=levels,colors='black',alpha=0.5)
plt.title('SIM_RH_INCR5_LL',y =1.05,fontsize=28)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.text(530., 650, 'b)', style='italic', fontsize=32)
plt.xlim((525,675))
plt.ylim((540,660))

ax = plt.subplot(2, 2, 3)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_above2300m_2000[5].plot.contourf(levels=levels,cmap='Blues',add_colorbar=False)
cb = plt.colorbar(plot,shrink=0.9,orientation='vertical')
cb.set_label('Kg/kg', labelpad=-40, y=1.15, rotation=0,fontsize=26)
QRAIN_above2300m_2000[5].plot.contour(levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'c)', style='italic', fontsize=32)
plt.vlines(619,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_RH_INCR5_UL', y =1.05,fontsize=28)

ax = plt.subplot(2, 2, 4)
plt.ylim((0,8000))
plot=ax.contourf(xs+Y2[0]+50,ys,w_cross,levels=np.linspace(-8,8,9),cmap="RdBu_r")
cb = plt.colorbar(plot,shrink=0.9,orientation='vertical')
cb.set_label('m/s', labelpad=-40, y=1.12, rotation=0,fontsize=26)
ax.fill_between(xs+Y2[0]+50,0,to_np(ter_line),facecolor="black")
plt.xlabel('y [km]',fontsize = 16)
ax.streamplot(Xs+Y2[0]+50,Ys,v_cross/1000,w_cross,density=2,arrowsize=1.7)
levs_w = [10**(-5)]
CS=plt.contour(xs+Y2[0]+50,ys,q_cross_filled,levels=levs_w,colors='black',linewidths=1.1)
#plt.clabel(CS, inline=1, fontsize=16,fmt='%.1f')
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.ylabel('z [m]',fontsize = 16)
plt.title('SIM_RH_INCR5_UL',fontsize=28,y=1.05)
plt.text(538., 7300, 'd)', style='italic', fontsize=32)


# from matplotlib import ticker
# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.13, 0.647, 0.72, 0.012])
# cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
# cb.set_label('m/s', labelpad=0, y=0.92,x=0.53, rotation=0,fontsize=16)


plt.subplots_adjust(left=0.1,
                    bottom=0.05, 
                    right=0.95, 
                    top=0.95, 
                    wspace=0.3, 
                    hspace=0.5)

#%% subplot con Fourier RAINNC e cumulative Precipitation
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_HGT_0_50 = Dataset(path_data + '/SIM_HGT_pert_0_50_each_cell_dominio2')
SIM_HGT_100_150 = Dataset(path_data + '/SIM_HGT_100_150_each_cell_dominio2')
SIM_HGT_100_150_first = Dataset(path_data + '/SIM_HGT_100_150_firstpart_dominio2')

RAINNC = getvar(SIM_HGT_100_150,"RAINNC",timeidx = ALL_TIMES) 
# Fourier terza ora
from scipy.fft import fft, fftfreq
from scipy.signal import blackman
# number of sample points
N = 125
#sample spacing
dx = 1
y_dominio3 = to_np(RAINNC[3,110,50:175]-RAINNC[2,110,50:175])
xf_dominio3 = fftfreq(N,dx)[:N//2]
y1 = y_dominio3-y_dominio3.mean()
yf_dominio3 = fft(y1)

# precipitazioni cumulate
RAINNC_HGT_0_50 = getvar(SIM_HGT_0_50,'RAINNC',timeidx=ALL_TIMES)
RAINNC_HGT_100_150 = getvar(SIM_HGT_100_150,'RAINNC',timeidx=ALL_TIMES)
HGT_0_50 = getvar(SIM_HGT_0_50,'HGT')
HGT_100_150 = getvar(SIM_HGT_100_150,'HGT')
QRAIN_HGT_100_150 = getvar(SIM_HGT_100_150,'QRAIN',timeidx=5)
QRAIN_HGT_100_150_2000 = interplevel(QRAIN_HGT_100_150,z,2000)
RAINNC_HGT_100_150_first = getvar(SIM_HGT_100_150_first,'RAINNC',timeidx=ALL_TIMES)
HGT_100_150_first = getvar(SIM_HGT_100_150_first,'HGT')

# changing coordinates
HGT_0_50 = HGT_0_50.assign_coords({'west_east' : HGT_0_50.west_east.values +X2[0]+1})
HGT_0_50 = HGT_0_50.assign_coords({'south_north' : HGT_0_50.south_north.values +X2[0]+1})
HGT_100_150 = HGT_100_150.assign_coords({'west_east' : HGT_100_150.west_east.values +X2[0]+1})
HGT_100_150 = HGT_100_150.assign_coords({'south_north' : HGT_100_150.south_north.values +X2[0]+1})
RAINNC_HGT_100_150 = RAINNC_HGT_100_150.assign_coords({'west_east' : RAINNC_HGT_100_150.west_east.values +X2[0]+1})
RAINNC_HGT_100_150 = RAINNC_HGT_100_150.assign_coords({'south_north' : RAINNC_HGT_100_150.south_north.values +X2[0]+1})
RAINNC_HGT_0_50 = RAINNC_HGT_0_50.assign_coords({'west_east' : RAINNC_HGT_0_50.west_east.values +X2[0]+1})
RAINNC_HGT_0_50 = RAINNC_HGT_0_50.assign_coords({'south_north' : RAINNC_HGT_0_50.south_north.values +X2[0]+1})
QRAIN_HGT_100_150_2000 = QRAIN_HGT_100_150_2000.assign_coords({'west_east' : QRAIN_HGT_100_150_2000.west_east.values +X2[0]+1})
QRAIN_HGT_100_150_2000 = QRAIN_HGT_100_150_2000.assign_coords({'south_north' : QRAIN_HGT_100_150_2000.south_north.values +X2[0]+1})
HGT_100_150_first = HGT_100_150_first.assign_coords({'west_east' : HGT_100_150_first.west_east.values +X2[0]+1})
HGT_100_150_first = HGT_100_150_first.assign_coords({'south_north' : HGT_100_150_first.south_north.values +X2[0]+1})
RAINNC_HGT_100_150_first = RAINNC_HGT_100_150_first.assign_coords({'west_east' : RAINNC_HGT_100_150_first.west_east.values +X2[0]+1})
RAINNC_HGT_100_150_first = RAINNC_HGT_100_150_first.assign_coords({'south_north' : RAINNC_HGT_100_150_first.south_north.values +X2[0]+1})


#%%
# plotting Fourier
fig = plt.figure(figsize=(15,18))
plt.rcParams.update({'font.size': 20})
ax = plt.subplot(3, 2, 3)
plt.rcParams["figure.figsize"] = (22,12)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.plot(xf_dominio3,2.0/N*np.abs(yf_dominio3[0:N//2]),linewidth=3.,color='blue',alpha=0.6)
plt.xlim((0,0.5))
#plt.ylim((0.1,10))
plt.xlabel('frequency (km$^{-1}$)', fontsize=20)
plt.ylabel('Amplitude',fontsize = 20)
plt.text(0.01, 8.6, 'c)', style='italic', fontsize=28)

ax = plt.subplot(3, 2, 1)
levels =[100,500,1000,1400,1499]
HGT_0_50.plot.contour(levels=levels,colors=['black'],alpha=0.3,x = 'west_east',y = 'south_north')
#levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=(RAINNC_HGT_0_50[5]-RAINNC_HGT_0_50[0]).plot.contourf(levels=np.linspace(5,120,116),cmap='BuPu',add_colorbar=False)
#cb = plt.colorbar(plot,shrink=0.9,orientation='vertical')
#cb.set_label('Kg/kg', labelpad=-40, y=1.15, rotation=0,fontsize=18)
#QRAIN_2000[4].plot.contour(levels=levels,colors='black',alpha=0.5)
plt.title('South-north speed 20 m above ground hour 5 SIM_V5',fontsize =15)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_HGT_0_50', y =1.02,fontsize=20)
plt.text(538., 642, 'a)', style='italic', fontsize=28)
plt.xlim((535,665))
plt.ylim((545,655))

ax = plt.subplot(3, 2, 2)
levels =[100,500,1000,1400,1499]
HGT_0_50.plot.contour(levels=levels,colors=['black'],alpha=0.3,x = 'west_east',y = 'south_north')
#levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=(RAINNC_HGT_100_150[5]-RAINNC_HGT_100_150[0]).plot.contourf(levels=np.linspace(0,150,116),cmap='BuPu',add_colorbar=False)
#cb = plt.colorbar(plot,shrink=0.9,orientation='vertical')
#cb.set_label('Kg/kg', labelpad=-40, y=1.15, rotation=0,fontsize=18)
#QRAIN_2000[4].plot.contour(levels=levels,colors='black',alpha=0.5)
plt.title('South-north speed 20 m above ground hour 5 SIM_V5',fontsize =15)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_HGT_100_150', y =1.02,fontsize=20)
plt.text(538., 642, 'b)', style='italic', fontsize=28)
plt.xlim((535,665))
plt.ylim((545,655))
plt.hlines(596,525,675, linestyle='dashed',color = 'black',alpha=0.85)

ax = plt.subplot(3, 2, 4)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot2=QRAIN_HGT_100_150_2000[45:180,45:180].plot.contourf(ax=ax, levels=levels,cmap='Blues',add_colorbar=False)  
cb = plt.colorbar(plot2,shrink=0.9)
cb.set_label('kg/kg', labelpad=-40, y=1.15, rotation=0,fontsize=18)
QRAIN_HGT_100_150_2000[45:180,45:180].plot.contour(ax=ax, levels=levels,colors = 'black',alpha = 0.5)  
HGT_100_150[45:180,45:180].plot.contour(ax=ax,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax3.text(537, 645, 'c)', style='italic',fontsize=22)
plt.title('SIM_HGT_100_150',fontsize=20,y=1.02)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((545,655))
plt.text(538., 642, 'd)', style='italic', fontsize=28)

ax = plt.subplot(3, 2, 5)
levels =[100,500,1000,1400,1499]
HGT_100_150_first.plot.contour(levels=levels,colors=['black'],alpha=0.3,x = 'west_east',y = 'south_north')
#levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot3=(RAINNC_HGT_100_150_first[3]-RAINNC_HGT_100_150_first[2]).plot.contourf(levels=np.linspace(4,36,20),cmap='BuPu',add_colorbar=False)
#QRAIN_2000[4].plot.contour(levels=levels,colors='black',alpha=0.5)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_HGT_100_150_FP', y =1.05,fontsize=20)
plt.text(567., 615, 'e)', style='italic', fontsize=28)
plt.xlim((565,625))
plt.ylim((580,620))

ax = plt.subplot(3, 2, 6)
levels =[100,500,1000,1400,1499]
HGT_100_150_first.plot.contour(levels=levels,colors=['black'],alpha=0.3,x = 'west_east',y = 'south_north')
#levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot4=(RAINNC_HGT_100_150_first[4]-RAINNC_HGT_100_150_first[3]).plot.contourf(levels=np.linspace(4,36,20),cmap='BuPu',add_colorbar=False)
#QRAIN_2000[4].plot.contour(levels=levels,colors='black',alpha=0.5)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.title('SIM_HGT_100_150_FP', y =1.03,fontsize=20)
plt.text(567., 615, 'f)', style='italic', fontsize=28)
plt.xlim((565,625))
plt.ylim((580,620))



from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.2, 0.645, 0.6, 0.013])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=10)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('mm', labelpad=-5, y=0.82,x=0.53, rotation=0,fontsize=20)


from matplotlib import ticker
fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.2, 0.045, 0.6, 0.013])
cb = fig.colorbar(plot3,cax=cbar_ax,orientation='horizontal')
tick_locator = ticker.MaxNLocator(nbins=10)
cb.locator = tick_locator
cb.update_ticks()
cb.set_label('mm/h', labelpad=0, y=0.82,x=0.53, rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.2, 
                    hspace=0.5)

#%% PLOT SINGLE BUMP: 
path_data2 = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_BUMP_A = Dataset(path_data2+'/SIM_BUMP_Y180_NOPERT_dominio2')
SIM_BUMP_B = Dataset(path_data2+'/SIM_BUMP_Y187_NOPERT_dominio2') 

HGT_BUMP_A = getvar(SIM_BUMP_A,'HGT')
z_BUMP_A = getvar(SIM_BUMP_A,'z')
HGT_BUMP_B = getvar(SIM_BUMP_B,'HGT')
z_BUMP_B = getvar(SIM_BUMP_B,'z')   
RAINNC_BUMP_A = getvar(SIM_BUMP_A,'RAINNC',timeidx=ALL_TIMES)
RAINNC_BUMP_B = getvar(SIM_BUMP_B,'RAINNC',timeidx=ALL_TIMES)

# changing coordinates   
HGT_BUMP_A = HGT_BUMP_A.assign_coords({'west_east' : HGT_BUMP_A.west_east.values +X2[0]+1})
HGT_BUMP_A = HGT_BUMP_A.assign_coords({'south_north' : HGT_BUMP_A.south_north.values +X2[0]+1})
HGT_BUMP_B = HGT_BUMP_B.assign_coords({'west_east' : HGT_BUMP_B.west_east.values +X2[0]+1})
HGT_BUMP_B = HGT_BUMP_B.assign_coords({'south_north' : HGT_BUMP_B.south_north.values +X2[0]+1})    
RAINNC_BUMP_A = RAINNC_BUMP_A.assign_coords({'west_east' : RAINNC_BUMP_A.west_east.values +X2[0]+1})
RAINNC_BUMP_A = RAINNC_BUMP_A.assign_coords({'south_north' : RAINNC_BUMP_A.south_north.values +X2[0]+1})
RAINNC_BUMP_B = RAINNC_BUMP_B.assign_coords({'west_east' : RAINNC_BUMP_B.west_east.values +X2[0]+1})
RAINNC_BUMP_B = RAINNC_BUMP_B.assign_coords({'south_north' : RAINNC_BUMP_B.south_north.values +X2[0]+1})

# south-north cross section to underline the mechanism:
w_BUMP_A = getvar(SIM_BUMP_A,'wa',timeidx=3)
QCLOUD_BUMP_A = getvar(SIM_BUMP_A,'QCLOUD',timeidx=3)
v_BUMP_A = getvar(SIM_BUMP_A,'va',timeidx=3)  
w_BUMP_B = getvar(SIM_BUMP_B,'wa',timeidx=2)
QCLOUD_BUMP_B = getvar(SIM_BUMP_B,'QCLOUD',timeidx=2)  
v_BUMP_B = getvar(SIM_BUMP_B,'va',timeidx=2)

start_point = CoordPair(x=112,y=30)
end_point = CoordPair(x=112, y=140)

# Get the terrain heights along the cross section line
ter_line_BUMP_A = interpline(HGT_BUMP_A,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables
w_cross_BUMP_A = vertcross(w_BUMP_A,z_BUMP_A,start_point=start_point,end_point=end_point)
qcloud_cross_BUMP_A = vertcross(QCLOUD_BUMP_A,z_BUMP_A,start_point=start_point,end_point=end_point)
v_cross_BUMP_A = vertcross(v_BUMP_A,z_BUMP_A,start_point=start_point,end_point=end_point)

w_modified_BUMP_A = np.zeros((100,111))
w_modified_BUMP_A[:] = np.nan
w_modified_BUMP_A[::2,::5]=w_cross_BUMP_A[::2,::5]

v_modified_BUMP_A = np.zeros((100,111))
v_modified_BUMP_A[:] = np.nan
v_modified_BUMP_A[::2,::5]=v_cross_BUMP_A[::2,::5]

xs = np.arange(0, w_cross_BUMP_A.shape[-1], 1)
ys = to_np(w_cross_BUMP_A.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

# stessa cross section per SIM_BUMP_B
start_point = CoordPair(x=112,y=30)
end_point = CoordPair(x=112, y=140)

# Get the terrain heights along the cross section line
ter_line_BUMP_B = interpline(HGT_BUMP_B,start_point=start_point,end_point=end_point,latlon=True)

# getting the cross sections of the variables
w_cross_BUMP_B = vertcross(w_BUMP_B,z_BUMP_B,start_point=start_point,end_point=end_point)
qcloud_cross_BUMP_B = vertcross(QCLOUD_BUMP_B,z_BUMP_B,start_point=start_point,end_point=end_point)
v_cross_BUMP_B = vertcross(v_BUMP_B,z_BUMP_B,start_point=start_point,end_point=end_point)

w_modified_BUMP_B = np.zeros((100,111))
w_modified_BUMP_B[:] = np.nan
w_modified_BUMP_B[::2,::5]=w_cross_BUMP_B[::2,::5]

v_modified_BUMP_B = np.zeros((100,111))
v_modified_BUMP_B[:] = np.nan
v_modified_BUMP_B[::2,::5]=v_cross_BUMP_B[::2,::5]

xs = np.arange(0, w_cross_BUMP_B.shape[-1], 1)
ys = to_np(w_cross_BUMP_B.coords["vertical"])
Xs, Ys = np.meshgrid(xs,ys)

# calcolo della convergenza per SIM_BUMP_A
u_BUMP_A = getvar(SIM_BUMP_A,'ua',timeidx=3)
z_agl_BUMP_A = getvar(SIM_BUMP_A,'height_agl')
u_BUMP_A_50 = interplevel(u_BUMP_A, z_agl_BUMP_A,50)
v_BUMP_A_50 = interplevel(v_BUMP_A, z_agl_BUMP_A,50)
z_agl_BUMP_A = getvar(SIM_BUMP_A,'height_agl')
from numpy import diff
du_dx = np.array(np.gradient(u_BUMP_A_50,axis=1))
dv_dy = np.array(np.gradient(v_BUMP_A_50,axis=0))
div = (du_dx + dv_dy)/1000 

# changing wind coordinates for quiver plot
u_BUMP_A_50 = u_BUMP_A_50.assign_coords({'west_east' : u_BUMP_A_50.west_east.values +X2[0]+1})
u_BUMP_A_50 = u_BUMP_A_50.assign_coords({'south_north' : u_BUMP_A_50.south_north.values +X2[0]+1})
v_BUMP_A_50 = v_BUMP_A_50.assign_coords({'west_east' : v_BUMP_A_50.west_east.values +X2[0]+1})
v_BUMP_A_50 = v_BUMP_A_50.assign_coords({'south_north' : v_BUMP_A_50.south_north.values +X2[0]+1})

#%% making the plot
fig = plt.figure(figsize=(18,14))
plt.rcParams.update({'font.size': 24})
gs = gridspec.GridSpec(4, 4)

ax1 = plt.subplot(gs[:2, :2])
plt.contourf(X2,Y2,(RAINNC_BUMP_A[3]-RAINNC_BUMP_A[2]),levels=np.linspace(5,90,18),cmap='BuPu')  
cb=plt.colorbar()
cb.set_label('mm/h',labelpad=-48, y=1.1,rotation=0)
HGT_BUMP_A.plot.contour(ax=ax1,levels = [100,400,1000,1400,1499],colors='black',alpha=0.5)
ax1.text(538, 644, 'a)', style='italic',fontsize=32)
plt.title('SIM_BUMP_A',fontsize=28,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((525,655))

ax2 = plt.subplot(gs[:2, 2:])
plt.contourf(X2,Y2,(RAINNC_BUMP_B[3]-RAINNC_BUMP_B[2]),levels=np.linspace(5,30,11),cmap='BuPu')  
cb=plt.colorbar()
cb.set_label('mm/h',labelpad=-48, y=1.1,rotation=0)
HGT_BUMP_B[45:180,45:180].plot.contour(ax=ax2,levels = [100,500,1000,1400,1499],colors='black',alpha=0.5)
ax2.text(538, 644, 'b)', style='italic',fontsize=32)
plt.title('SIM_BUMP_B',fontsize=28,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((525,655))

ax5 = plt.subplot(gs[2:4, 2:])
#plot=(RAINNC_BUMP_A[5]-RAINNC_BUMP_A[0]).plot.contourf(ax=ax5, levels=120,cmap='BuPu',add_colorbar=False)  
plt.contourf(X2,Y2,(RAINNC_BUMP_A[5]-RAINNC_BUMP_A[0]),levels=np.linspace(5,355,36),cmap='BuPu')
cb=plt.colorbar()
cb.set_label('mm/h',labelpad=-48, y=1.1,rotation=0)
HGT_BUMP_A.plot.contour(ax=ax5,levels = [100,400,1000,1400,1499],colors='black',alpha=0.5)
ax5.text(538, 644, 'd)', style='italic',fontsize=32)
plt.title('SIM_BUMP_A',fontsize=28,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((525,655))

ax5 = plt.subplot(gs[2:4, :2])
plot5=plt.contourf(X2,Y2,-div,levels=np.linspace(-10**-2,10**-2,41),cmap='seismic')
cb=plt.colorbar(plot5)
ax5.quiver(u_BUMP_A_50[::3,::2].west_east, u_BUMP_A_50[::3,::2].south_north, u_BUMP_A_50[::3,::2], v_BUMP_A_50[::3,::2]/2,color='black',scale=80.,width=0.003)
cb.set_label('$s^{-1}$',labelpad=-48, y=1.13,rotation=0)
HGT_BUMP_A.plot.contour(ax=ax5,levels = [100,400,1000,1400,1499],colors='black',alpha=0.5)
ax5.text(581, 575, 'c)', style='italic',fontsize=32)
plt.title('SIM_BUMP_A',fontsize=28,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((580,620))
plt.ylim((520,580))

# from matplotlib import ticker
# fig.subplots_adjust(right=0.8)
# cbar_ax = fig.add_axes([0.2, 0.347, 0.6, 0.011])
# cb = fig.colorbar(plot2,cax=cbar_ax,orientation='horizontal')
# #tick_locator = ticker.MaxNLocator(nbins=10)
# #cb.locator = tick_locator
# #cb.update_ticks()
# cb.set_label('kg/kg', labelpad=-75, y=0.82,x=0.53, rotation=0,fontsize=24)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.8, 
                    hspace=0.8)

#%% faccio i plot singoli della wave
# SIM_BUMP_A
fig = plt.figure(figsize=(18,12))
plt.rcParams.update({'font.size': 30})
plt.ylim((0,9000))
plt.contourf(xs+Y2[0]+30,ys,qcloud_cross_BUMP_A*1000,levels=np.linspace(0.,2.2,23),cmap="Greys")
cb = plt.colorbar(shrink=1.,orientation='vertical')
cb.set_label('$10^{-3}$ kg/kg', labelpad=-40, y=1.07, rotation=0,fontsize=28)
plt.fill_between(xs+Y2[0]+30,0,to_np(ter_line_BUMP_A),facecolor="black")
plt.xlabel('y [km]')
#plt.contour(xs+Y2[0]+30,ys,qcloud_cross_BUMP_A,levels=[0.00001],colors='royalblue',linewidths=3.)
plt.streamplot(Xs+Y2[0]+30,Ys,v_cross_BUMP_A/1000,w_cross_BUMP_A,density=2,arrowsize=2.3,linewidth=2.)
levs_w = [-0.5,-0.2,0.2,0.5,1,4,8,12,15,20]
CS=plt.contour(xs+Y2[0]+30,ys,w_cross_BUMP_A,levels=levs_w,colors='black',linewidths=1.5)
plt.clabel(CS, inline=1, fontsize=24,fmt='%.1f')
plt.ylabel('z [m]')
plt.title('SIM_BUMP_A',fontsize=30,y=1.02)

#%% SIM_BUMP_B
fig = plt.figure(figsize=(18,12))
plt.rcParams.update({'font.size': 30})
plt.ylim((0,9000))
plt.contourf(xs+Y2[0]+30,ys,qcloud_cross_BUMP_B*1000,levels=np.linspace(0.,2.2,23),cmap="Greys")
cb = plt.colorbar(shrink=1.,orientation='vertical')
cb.set_label('$10^{-3}$ kg/kg', labelpad=-40, y=1.07, rotation=0,fontsize=28)
plt.fill_between(xs+Y2[0]+30,0,to_np(ter_line_BUMP_B),facecolor="black")
plt.contour(xs+Y2[0]+30,ys,qcloud_cross_BUMP_B,levels=[0.00001],colors='royalblue',linewidths=3.)
plt.xlabel('y [km]')
plt.streamplot(Xs+Y2[0]+30,Ys,v_cross_BUMP_B/1200,w_cross_BUMP_B,density=2,arrowsize=2.3,linewidth=2.3)
levs_w = [-0.5,-0.2,0.2,0.5,1,4,8,12,15,20]
CS=plt.contour(xs+Y2[0]+30,ys,w_cross_BUMP_B,levels=levs_w,colors='black',linewidths=1.5)
plt.clabel(CS, inline=1, fontsize=24,fmt='%.1f')
plt.ylabel('z [m]')
plt.title('SIM_BUMP_B',fontsize=30,y=1.02)
plt.text(520., 8200, 'd)', style='italic', fontsize=30)

#%% PLOT APPENDICE MULTIPLE BUMP-->cumulata erza ora
path_data2 = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_HGT/")
SIM_COLUMN_A = Dataset(path_data2+'/SIM_COLUMN_Y180_NOPERT_dominio2')
SIM_COLUMN_B = Dataset(path_data2+'/SIM_COLUMN_Y185_NOPERT_1KM')

HGT_COLUMN_A = getvar(SIM_COLUMN_A,'HGT')
HGT_COLUMN_B = getvar(SIM_COLUMN_B,'HGT') 
RAINNC_COLUMN_A = getvar(SIM_COLUMN_A,'RAINNC',timeidx=ALL_TIMES)
RAINNC_COLUMN_B = getvar(SIM_COLUMN_B,'RAINNC',timeidx=ALL_TIMES)

# changing coordinates   
HGT_COLUMN_A = HGT_COLUMN_A.assign_coords({'west_east' : HGT_COLUMN_A.west_east.values +X2[0]+1})
HGT_COLUMN_A = HGT_COLUMN_A.assign_coords({'south_north' : HGT_COLUMN_A.south_north.values +X2[0]+1})
HGT_COLUMN_B = HGT_COLUMN_B.assign_coords({'west_east' : HGT_COLUMN_B.west_east.values +X2[0]+1})
HGT_COLUMN_B = HGT_COLUMN_B.assign_coords({'south_north' : HGT_COLUMN_B.south_north.values +X2[0]+1})    
RAINNC_COLUMN_A = RAINNC_COLUMN_A.assign_coords({'west_east' : RAINNC_COLUMN_A.west_east.values +X2[0]+1})
RAINNC_COLUMN_A = RAINNC_COLUMN_A.assign_coords({'south_north' : RAINNC_COLUMN_A.south_north.values +X2[0]+1})
RAINNC_COLUMN_B = RAINNC_COLUMN_B.assign_coords({'west_east' : RAINNC_COLUMN_B.west_east.values +X2[0]+1})
RAINNC_COLUMN_B = RAINNC_COLUMN_B.assign_coords({'south_north' : RAINNC_COLUMN_B.south_north.values +X2[0]+1})

fig = plt.figure(figsize=(18,14))
plt.rcParams.update({'font.size': 24})
gs = gridspec.GridSpec(4, 4)

ax1 = plt.subplot(gs[:2, :2])
plt.contourf(X2,Y2,(RAINNC_COLUMN_A[3]-RAINNC_COLUMN_A[2]),levels=np.linspace(5,60,12),cmap='BuPu')  
cb=plt.colorbar()
cb.set_label('mm/h',labelpad=-48, y=1.1,rotation=0)
HGT_COLUMN_A.plot.contour(ax=ax1,levels = [-100,100,250,400,1000,1400,1499],colors='black',alpha=0.5)
ax1.text(538, 644, 'a)', style='italic',fontsize=32)
plt.title('SIM_COLUMN_A',fontsize=28,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((525,655))

ax2 = plt.subplot(gs[:2, 2:])
plt.contourf(X2,Y2,(RAINNC_COLUMN_B[3]-RAINNC_COLUMN_B[2]),levels=np.linspace(5,60,12),cmap='BuPu')  
cb=plt.colorbar()
cb.set_label('mm/h',labelpad=-48, y=1.1,rotation=0)
HGT_COLUMN_B[45:180,45:180].plot.contour(ax=ax2,levels = [-100,100,250,500,1000,1400,1499],colors='black',alpha=0.5)
ax2.text(538, 644, 'b)', style='italic',fontsize=32)
plt.title('SIM_COLUMN_B',fontsize=28,y=1.01)
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.xlim((535,665))
plt.ylim((525,655))

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.6, 
                    hspace=0.8)

#%% plot profili verticali presi a Udine e Padova
path_data= os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_REAL_FINALE/")
SIM_REAL = Dataset(path_data + '/SIM_ORIGINAL_SLP_ORIGINAL_WIND_ore18')

lat = getvar(SIM_REAL, "lat")
lon = getvar(SIM_REAL, "lon")
HGT = getvar(SIM_REAL,'HGT')
z= getvar(SIM_REAL,'z')
TIMES = getvar(SIM_REAL,"times",timeidx=ALL_TIMES)
theta_e = getvar(SIM_REAL,'eth',timeidx=ALL_TIMES)
RH = getvar(SIM_REAL,'rh',timeidx=ALL_TIMES)
wspd = getvar(SIM_REAL,'wspd',timeidx=ALL_TIMES)

#%% plotto profili verticali variabili prime sei ore Udine:
# prendo le coordinate di Venezia che è lontano dalle montagne e al centro del dominio:
import matplotlib.gridspec as gridspec

fig = plt.figure(figsize=(24,20))
plt.rcParams.update({'font.size': 35})
gs = gridspec.GridSpec(4, 4)

ax1 = plt.subplot(gs[:2, :2])
plt.plot(theta_e[0,:,91,168],z[:,91,168],label='hour0',linewidth=2.)
plt.plot(theta_e[1,:,91,168],z[:,91,168],label='hour2',linewidth=2.)
plt.plot(theta_e[2,:,91,168],z[:,91,168],label='hour3',linewidth=2.)
plt.plot(theta_e[3,:,91,168],z[:,91,168],label='hour4',linewidth=2.)
plt.plot(theta_e[5,:,91,168],z[:,91,168],label='hour5',linewidth=2.)
plt.plot(theta_e[10,:,91,168],z[:,91,168],label='hour10',linewidth=2.)
plt.xlim((310,335))
plt.ylim((0,6000))
ax1.text(310.5, 5400, 'a)', style='italic',fontsize=40)
plt.xlabel('$\Theta_e$ [K]')
plt.ylabel('z [m]')

ax2 = plt.subplot(gs[:2, 2:])
plt.plot(RH[0,:,91,168],z[:,91,168],label='t = 0 h',linewidth=2.)
plt.plot(RH[2,:,91,168],z[:,91,168],label='t = 2 h',linewidth=2.)
plt.plot(RH[3,:,91,168],z[:,91,168],label='t = 3 h',linewidth=2.)
plt.plot(RH[5,:,91,168],z[:,91,168],label='t = 5 h',linewidth=2.)
plt.plot(RH[10,:,91,168],z[:,91,168],label='t = 10 h',linewidth=2.)
plt.legend(bbox_to_anchor =(0.4, -0.35))
plt.xlim((50,100))
plt.ylim((0,6000))
ax2.text(51, 5400, 'b)', style='italic',fontsize=40)
plt.xlabel('RH %')
plt.ylabel('z [m]')

ax3 = plt.subplot(gs[2:4, :2])
plt.plot(wspd[0,:,91,168],z[:,91,168],label='hour0',linewidth=2.)
plt.plot(wspd[2,:,91,168],z[:,91,168],label='hour2',linewidth=2.)
plt.plot(wspd[3,:,91,168],z[:,91,168],label='hour3',linewidth=2.)
plt.plot(wspd[5,:,91,168],z[:,91,168],label='hour5',linewidth=2.)
plt.plot(wspd[10,:,91,168],z[:,91,168],label='hour10',linewidth=2.)
ax3.text(4, 5400, 'c)', style='italic',fontsize=40)
plt.ylim((0,6000))
plt.xlabel('Wind speed [m/s]')
plt.ylabel('z [m]')

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=1.1, 
                    hspace=1.3)


#%% PLOT DELLA SIM_N3_00001

path_data = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/SIM_IDEAL_N/")
SIM_N3_0001 = Dataset(path_data + '/SIM_N1_00004_N3_0001_dominio2_FINALE')


z = getvar(SIM_N3_0001,'z')
HGT = getvar(SIM_N3_0001,'HGT')
QRAIN_N3_0001 = getvar(SIM_N3_0001,'QRAIN',timeidx=ALL_TIMES)
QRAIN_N3_0001_2000 = interplevel(QRAIN_N3_0001,z,2000)
HGT = HGT.assign_coords({'west_east' : HGT.west_east.values +X2[0]+1})
HGT = HGT.assign_coords({'south_north' : HGT.south_north.values +X2[0]+1})
QRAIN_N3_0001_2000 = QRAIN_N3_0001_2000.assign_coords({'west_east' : QRAIN_N3_0001_2000.west_east.values +X2[0]+1})
QRAIN_N3_0001_2000 = QRAIN_N3_0001_2000.assign_coords({'south_north' : QRAIN_N3_0001_2000.south_north.values +X2[0]+1})

# making the plot

fig = plt.figure(figsize=(14,16))
plt.rcParams.update({'font.size': 20})

ax = plt.subplot(2, 2, 1)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N3_0001_2000[6].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_N3_0001_2000[6].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'a)', style='italic', fontsize=26)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')

ax = plt.subplot(2, 2, 2)
levels=[0.0004,0.0008,0.0012,0.002,0.003,0.004]
plot=QRAIN_N3_0001_2000[7].plot.contourf(ax=ax,levels=levels,cmap='Blues',add_colorbar=False)
#cb = plt.colorbar(plot)
QRAIN_N3_0001_2000[7].plot.contour(ax=ax,levels=levels,colors='black',alpha=0.5)
HGT.plot.contour(ax=ax,levels=[100,500,1000,1400,1499],colors='black',alpha=0.5)
plt.title('')
plt.text(529., 648, 'b)', style='italic', fontsize=26)
#plt.hlines(593,525,675, linestyle='dashed',color = 'black')
plt.xlim((525,675))
plt.ylim((540,660))
plt.xlabel('x [km]')
plt.ylabel('y [km]')
plt.vlines(619,525,675,linestyle='dashed',color = 'black')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.14, 0.6, 0.72, 0.008])
cb = fig.colorbar(plot,cax=cbar_ax,orientation='horizontal')
#tick_locator = ticker.MaxNLocator(nbins=5)
#cb.locator = tick_locator
cb.update_ticks()
cb.set_label('kg/kg', labelpad=0, x=0.52, rotation=0,fontsize=20)

plt.subplots_adjust(left=0.1,
                    bottom=0.1, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.3, 
                    hspace=1.3)