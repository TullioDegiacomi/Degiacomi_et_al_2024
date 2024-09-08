# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 11:24:11 2021

@author: Tullio Degiacomi
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap,colors
import os
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, get_basemap, ALL_TIMES)


#path_data = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/Simulazioni_Coriolis/SIM_200m_005K")
#path_output = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Immagini_dominio_centrato/Immagini_200m_noCoriolis_sounding_sud/Immagini_200m_dominio2")
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIMULAZIONI_DOMINIO_CENTRATO/SIM_200m_noCoriolis_sounding_sud")
#path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/CTRL_SIMULATION/")
ncfile = Dataset(path_data + '/SIM_200m_noCoriolis_sounding_sud_dominio3')

# Get the latitude and longitude
lat = getvar(ncfile, "lat")
lon = getvar(ncfile, "lon")
TIMES = getvar(ncfile,"times",timeidx=ALL_TIMES)


# Estraggo orografia
HGT = getvar(ncfile,"HGT")

#%%
# Estraggo e calcolo precipitazioni 
RAINNC = getvar(ncfile,"RAINNC",timeidx = ALL_TIMES) 
QRAIN = getvar(ncfile,"QRAIN",timeidx = ALL_TIMES)

z = getvar(ncfile,"z",units="m")
q = interplevel(QRAIN,z,2000)
#q = q.values
# set dimensions of domain
imax = 975
jmax = 975
RAIN_H = np.arange((len(TIMES)-1)*imax*jmax).reshape((len(TIMES)-1, imax,jmax))

#%%
#RAIN_H = [None]* (len(TIMES)-2)
for i in range(len(TIMES)-2):
    RAIN_H[i]=RAINNC[i+1]-RAINNC[i]

# compute last hour precipitation
i=len(RAINNC)-1
#RAIN = RAINNC[12]-RAINNC[11]
#RAIN=RAIN.values

#%%
# plotting q_rain
levels =[100,500,1000,1400,1499]
plt.contour(X3,Y3,HGT,levels=levels,colors=['black'],alpha=0.5)

# setting plot parameters
plt.rcParams["figure.figsize"] = (18,14)
levels_color=np.linspace(RAIN_H[5].min()+5,RAIN_H[5].max(),50)
levels=np.linspace(RAIN_H[5].min()+5,RAIN_H[5].max(),6)
#levels_color = [0.0006,0.0008,0.0012,0.0020,0.0030,0.0040]
#q[5].plot.contourf(levels=levels_color,cmap='Blues')
#HGT.plot.contour(levels=levels,colors='black',linewidths=0.8)
plt.contourf(X3,Y3,RAIN_H[5],levels=levels_color,cmap='Blues')
cb = plt.colorbar(orientation='vertical',format='%d')
cb.ax.tick_params(labelsize=16)
cb.set_label('mm/h', labelpad=-40, y=1.05, rotation=0,fontsize=18)
#plt.contour(X3,Y3,RAIN_H[i],levels=levels,colors=['black'])
plt.xlabel('x (Km)',fontsize = 18)
plt.ylabel('y (Km)',fontsize = 18)
#plt.title('Hourly precipitation hour ' + str(i+1),fontsize = 20)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
    #plt.xlim((150,250))
    #plt.ylim((150,250))

#%% definition of the domain mesh
# # definition of the domain for domain 2 and 3 in order to have KM on axis

cells1 = 400
X1 = np.arange(1.5,1201.5,3) 
Y1=X1


# definition of the domain for domain 2 and 3 in order to have KM on axis
i_start_2 = 163
#i_start_3 = 35
i_start_3 = 41
j_start_3 = 51
#start2 = i_start_2*3   # multiplying for resolution of domain 1
start2 = (i_start_2-1)*3 + 0.5
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
start3_x = start2 + 15 +0.2
start3_y = start2 + 15 +0.2
X3 = np.arange(start3_x,start3_x + cells3_x/5 -0.1,0.2) 
Y3 = np.arange(start3_y,start3_y + cells3_y/5-0.1,0.2)

#%% computing the fourier transform
from scipy.fft import fft, fftfreq
from scipy.signal import blackman

# number of sample points
N = 625
#sample spacing
dx = 0.2

# per dominio 3 ho tagliato d 45 a 670
#y = RAIN_H[6,116,50:175]          # per dominio 500 m ho tagliato a 146 per avere corrispondenza con il taglio a 108 del dominio 2
#y_dominio2 = RAIN_H[6,110,50:175]
# taglio per dominio a 200m
y_dominio3 = RAIN_H[5,448,179:804]
#y_dominio3 = RAIN_H[6][565][175:800]
#y_dominio3 = RAIN_H[5,185,69:319]
#yf = fft(y)
xf_dominio3 = fftfreq(N,dx)[:N//2]
#xf = fftfreq(N,dx)[:N//2]


#%%
#let's clean a bit the rain signal y
# for i in range(len(y)):
#     if(y[i]<30):
#         y[i]=0
y1 = y_dominio3-y_dominio3.mean()

# for i in range(len(y1)):
#       if(y1[i]<10):
#           y1[i]=y1.min()
plt.plot(X3[0:N]+69/2,y1)
#plt.ylim((y1.min(),60))

# #%%
# import matplotlib.pyplot as plt

# sp = np.fft.fft(y)
# freq = np.fft.fftfreq(y.shape[-1])

# plt.plot(sp.real)
# #plt.xlim((0,0.5))
# #plt.clf()

#%%
w = blackman(N)
yf_dominio3 = fft(y1)
# plotting Fourier
plt.rcParams["figure.figsize"] = (22,12)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.subplot(1,2,1)
plt.plot(xf_dominio3,2.0/N*np.abs(yf_dominio3[0:N//2]))
plt.xlim((0,0.5))
#plt.ylim((0.1,10))
plt.xlabel('frequency (km$^{-1}$)', fontsize=18)
plt.ylabel('Amplitude',fontsize = 18)
plt.title("Fourier spectrum rainbands hour 7",fontsize=18)

#plotting rainbands
plt.subplot(1,2,2)
plt.plot(X3[0:N]+69/2,y1)
plt.ylim((y1.min(),45))
#plt.xlim((550,590))
#plt.ylim((0,0.001))
plt.title('Hourly prec signal hour  7',fontsize=18)
#plt.ylabel('q$_r$ [Kg/Kg]',fontsize = 18)
plt.ylabel('Hourly Prec [mm/hr]',fontsize = 18)
plt.xlabel('x [Km]', fontsize = 18)
plt.tick_params(axis='both', which='major', labelsize=16)
#plt.savefig(path_output + '/Immagini_Nesting_3_finali/Immagini_500m_kmopt4_dominio3/Fourier_migliori/Fourier_qrain_hour'+str(j)+'_tagliatoY3[130]eqrain0.0006.png')
# #plt.clf()     


#%% same procedure for the domain 2
# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 11:24:11 2021

@author: Tullio Degiacomi
"""

from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap,colors
import os
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, get_basemap, ALL_TIMES)


#path_data = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/SIMULAZIONI/")
#path_output = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Immagini_dominio_centrato/Immagini_200m_noCoriolis/Immagini_200m_dominio2")
path_data = os.path.abspath("E:/TESI/TESI/WRF_DATA/SIMULAZIONI/SIMULAZIONI_DOMINIO_CENTRATO/SIM_200m_noCoriolis_sounding_sud")

ncfile = Dataset(path_data + '/SIM_200m_noCoriolis_sounding_sud_dominio2')

# Get the latitude and longitude
lat = getvar(ncfile, "lat")
lon = getvar(ncfile, "lon")
TIMES = getvar(ncfile,"times",timeidx=ALL_TIMES)


# Estraggo orografia
HGT = getvar(ncfile,"HGT")

#%%
# Estraggo e calcolo precipitazioni 
RAINNC = getvar(ncfile,"RAINNC",timeidx = ALL_TIMES) 
QRAIN = getvar(ncfile,"QRAIN",timeidx = ALL_TIMES)

z = getvar(ncfile,"z",units="m")
q = interplevel(QRAIN,z,2000)
q = q.values
# set dimensions of domain
imax = 225
jmax = 225
RAIN_H = np.arange((len(TIMES)-1)*imax*jmax).reshape((len(TIMES)-1, imax,jmax))

#%%
for i in range(len(TIMES)-2):
    RAIN_H[i]=RAINNC[i+1]-RAINNC[i]

# compute last hour precipitation
i=len(RAINNC)-1
#RAIN = RAINNC[12]-RAINNC[11]
#RAIN=RAIN.values

#%%
# plotting q_rain
levels =[100,500,1000,1400,1499]
plt.contour(X2,Y2,HGT,levels=levels,colors=['black'],alpha=0.5)

    # setting plot parameters
plt.rcParams["figure.figsize"] = (16,14)
levels_color=np.linspace(RAIN_H[5].min()+5,RAIN_H[5].max(),50)
levels=np.linspace(RAIN_H[5].min()+5,RAIN_H[5].max(),6)
plt.contourf(X2,Y2,RAIN_H[5],levels=levels_color,cmap='BuPu')
cb = plt.colorbar(orientation='vertical',format='%d')
cb.ax.tick_params(labelsize=16)
cb.set_label('mm/h', labelpad=-40, y=1.05, rotation=0,fontsize=18)
#plt.contour(X3,Y3,RAIN_H[i],levels=levels,colors=['black'])
plt.xlabel('x (Km)',fontsize = 18)
plt.ylabel('y (Km)',fontsize = 18)
plt.title('Hourly precipitation hour ' + str(i+1),fontsize = 20)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
    #plt.xlim((150,250))
    #plt.ylim((150,250))

#%% computing the fourier transform
from scipy.fft import fft, fftfreq
from scipy.signal import blackman

# number of sample points
N = 125
#sample spacing
dx = 1

# per dominio 2 ho tagliato d 45 a 670
y_dominio2 = RAIN_H[5,105,50:175]         
#y = RAIN_H[6,174,68:318]
#yf = fft(y)
xf_dominio2 = fftfreq(N,dx)[:N//2]
#xf = fftfreq(N,dx)[:N//2]

#%%
#let's clean a bit the rain signal y
# for i in range(len(y)):
#     if(y[i]<0.000):
#         y[i]=0
y = y_dominio2-y_dominio2.mean()

plt.plot(X2[0:N]+50,y)

# #%%
# import matplotlib.pyplot as plt

# sp = np.fft.fft(y)
# freq = np.fft.fftfreq(y.shape[-1])

# plt.plot(sp.real)
# #plt.xlim((0,0.5))
# #plt.clf()

#%%
yf_dominio2 = fft(y)
# plotting Fourier
plt.rcParams["figure.figsize"] = (14,10)
plt.tick_params(axis='both', which='major', labelsize=16)
plt.subplot(2,1,1)
plt.plot(xf_dominio2,0.4/N*np.abs(yf_dominio2[0:N//2]))
plt.xlim((0,0.5))
#plt.ylim((0.1,10))
plt.xlabel('frequency (km$^{-1}$)', fontsize=18)
plt.ylabel('Amplitude',fontsize = 18)
plt.title("Fourier spectrum rainbands hour 6",fontsize=18)

#plotting rainbands
plt.subplot(2,1,2)
plt.plot(X2[0:N]+50,y)
#plt.xlim((550,590))
#plt.ylim((0,0.001))
plt.title('Hourly prec signal hour  6',fontsize=18)
#plt.ylabel('q$_r$ [Kg/Kg]',fontsize = 18)
plt.ylabel('Hourly Prec [mm/hr]',fontsize = 18)
plt.xlabel('x [Km]', fontsize = 18)
plt.tick_params(axis='both', which='major', labelsize=16)
#plt.savefig(path_output + '/Immagini_Nesting_3_finali/Immagini_500m_kmopt4_dominio3/Fourier_migliori/Fourier_qrain_hour'+str(j)+'_tagliatoY3[130]eqrain0.0006.png')
# #plt.clf()     


#%% plot finale con sovrapposizione
yf_dominio2 = fft(y)
# plotting Fourier
plt.rcParams["figure.figsize"] = (16,16)
plt.rcParams.update({'font.size': 38})
plt.tick_params(axis='both', which='major')
plt.subplot(2,1,1)
plt.plot(xf_dominio2,2.0/N*np.abs(yf_dominio2[0:N//2]),alpha=0.8,linewidth=4.5,label='Spectrum 1km')
plt.plot(xf_dominio3[:62],0.4/N*np.abs(yf_dominio3[:62]),alpha=0.8,linewidth=4.5,label = 'Spectrum 200m')
plt.xlim((0,0.5))
plt.xlabel('Frequency (km$^{-1}$)')
plt.ylabel('Amplitude')
#plt.title("Fourier spectrum rainbands hour 6",fontsize=18)
plt.legend(fontsize=34)
plt.xticks()
plt.yticks()
plt.ylim((0,7))
plt.text(0.007, 6.3, 'a)', style='italic', fontsize=38)


#plotting rainbands
x1 = np.arange(539, 664, 1)
x2 = np.arange(539, 664, 0.2)
plt.subplot(2,1,2)
plt.plot(x1,y_dominio2,linewidth=4.5,label='Hourly Prec. 1km')
plt.plot(x2,y_dominio3,linewidth=4.5,label='Hourly Prec. 200m',alpha=0.7)
plt.legend(fontsize=34,loc=1)
#plt.title('Hourly prec signal hour  6',fontsize=18)
#plt.ylabel('q$_r$ [Kg/Kg]',fontsize = 18)
plt.ylabel('Hourly Prec [mm/h]')
plt.xlabel('x [Km]')
plt.tick_params(axis='both', which='major')
plt.ylim((0,65))
plt.text(535, 57, 'b)', style='italic', fontsize=38)
#plt.savefig(path_output + '/Immagini_Nesting_3_finali/Immagini_500m_kmopt4_dominio3/Fourier_migliori/Fourier_qrain_hour'+str(j)+'_tagliatoY3[130]eqrain0.0006.png')
# #plt.clf()  

plt.subplots_adjust(left=0.1,
                    bottom=0.2, 
                    right=0.95, 
                    top=0.9, 
                    wspace=0.82, 
                    hspace=0.4)

