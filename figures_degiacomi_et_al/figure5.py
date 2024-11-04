# -*- coding: utf-8 -*-
"""
Created on Tue Jun  1 11:24:11 2021

@author: Tullio Degiacomi
"""

from netCDF4 import Dataset
import numpy as np

import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

import matplotlib.pyplot as plt
from matplotlib.cm import get_cmap,colors
import os
from wrf import (to_np, getvar, smooth2d, get_cartopy, cartopy_xlim,
                 cartopy_ylim, latlon_coords, interplevel, get_basemap, ALL_TIMES)


path_data = os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_200m")
ncfile = Dataset(path_data + '/SIM_200m_domain3')

# Get the latitude and longitude
lat = getvar(ncfile, "lat")
lon = getvar(ncfile, "lon")
TIMES = getvar(ncfile,"times",timeidx=ALL_TIMES)
HGT = getvar(ncfile,"HGT")
RAINNC = getvar(ncfile,"RAINNC",timeidx = ALL_TIMES) 
z = getvar(ncfile,"z",units="m")

# set dimensions of domain
imax = 975
jmax = 975
RAIN_H = np.arange((len(TIMES)-1)*imax*jmax).reshape((len(TIMES)-1, imax,jmax))

for i in range(len(TIMES)-2):
    RAIN_H[i]=RAINNC[i+1]-RAINNC[i]

cells1 = 400
X1 = np.arange(1.5,1201.5,3) 
Y1=X1


# definition of the domain for domain 2 and 3 in order to have KM on axis
i_start_2 = 163
i_start_3 = 41
j_start_3 = 51
start2 = (i_start_2-1)*3
cells2 = 225
cells3 = 390
cells3_x = 975
cells3_y = 975
X2 = np.arange(start2,start2+cells2) 
Y2 = X2


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

y_dominio3 = RAIN_H[5,448,179:804]
xf_dominio3 = fftfreq(N,dx)[:N//2]

y1 = y_dominio3-y_dominio3.mean()

plt.plot(X3[0:N]+69/2,y1)

yf_dominio3 = fft(y1)


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


path_data = os.path.abspath("/archive/tullio/SIMULAZIONI/SIM_200m")

ncfile = Dataset(path_data + '/SIM_200m_domain2')

# Get the latitude and longitude
lat = getvar(ncfile, "lat")
lon = getvar(ncfile, "lon")
TIMES = getvar(ncfile,"times",timeidx=ALL_TIMES)
HGT = getvar(ncfile,"HGT")
RAINNC = getvar(ncfile,"RAINNC",timeidx = ALL_TIMES) 

z = getvar(ncfile,"z",units="m")
# set dimensions of domain
imax = 225
jmax = 225
RAIN_H = np.arange((len(TIMES)-1)*imax*jmax).reshape((len(TIMES)-1, imax,jmax))

#%%
for i in range(len(TIMES)-2):
    RAIN_H[i]=RAINNC[i+1]-RAINNC[i]


#%% computing the fourier transform
from scipy.fft import fft, fftfreq
from scipy.signal import blackman

# number of sample points
N = 125
#sample spacing
dx = 1

y_dominio2 = RAIN_H[5,105,50:175]      
xf_dominio2 = fftfreq(N,dx)[:N//2]

y = y_dominio2-y_dominio2.mean()

plt.plot(X2[0:N]+50,y)

yf_dominio2 = fft(y)
# plotting Fourier
plt.rcParams["figure.figsize"] = (16,24)
plt.rcParams.update({'font.size': 14})
plt.tick_params(axis='both', which='major')
plt.subplot(2,1,1)
plt.plot(xf_dominio2,2.0/N*np.abs(yf_dominio2[0:N//2]),alpha=0.8,linewidth=2,label='Spectrum 1 km')
plt.plot(xf_dominio3[:62],0.4/N*np.abs(yf_dominio3[:62]),alpha=0.8,linewidth=2,label = 'Spectrum 200 m')
plt.xlim((0,0.5))
plt.xlabel('Frequency [km$^{-1}$]')
plt.ylabel('Amplitude')
#plt.title("Fourier spectrum rainbands hour 6",fontsize=18)
plt.legend(fontsize=12)
plt.xticks()
plt.yticks()
plt.ylim((0,7))
plt.text(0.007, 6.3, '(a)', fontsize=14)


x1 = np.arange(539, 664, 1)
x2 = np.arange(539, 664, 0.2)
plt.subplot(2,1,2)
plt.plot(x1,y_dominio2,linewidth=2,label='Hourly Prec. 1 km')
plt.plot(x2,y_dominio3,linewidth=2,label='Hourly Prec. 200 m',alpha=0.7)
plt.legend(fontsize=12,loc=1)
plt.ylabel('Hourly Prec. [mm h$^{-1}$]')
plt.xlabel('$\it{x}$ [km]')
plt.tick_params(axis='both', which='major')
plt.ylim((0,65))
plt.text(535, 57, '(b)',  fontsize=14)


plt.subplots_adjust(left=0.2,
                    bottom=0.12, 
                    right=0.95, 
                    top=0.95, 
                    wspace=0.1, 
                    hspace=0.4)

plt.savefig('figure5',dpi=600)
