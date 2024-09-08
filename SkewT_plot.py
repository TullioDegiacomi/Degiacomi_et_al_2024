# -*- coding: utf-8 -*-
"""
Created on Mon Aug 30 20:41:19 2021

@author: Admin
"""

# -*- coding: utf-8 -*-
"""
Created on Sun May  9 11:18:27 2021

@author: Admin
"""
#%%
# creating the skew-T logP diagram
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
import metpy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

import metpy.calc as mpcalc
from metpy.cbook import get_test_data
from metpy.plots import SkewT,Hodograph
from metpy.units import units
import matplotlib.gridspec as gridspec
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

# reading the dataset
path_output = os.path.abspath("C:/Users/Admin/Desktop/TESI/WRF_DATA/")
# stability plot with potential temperature and equivalent potential temp.
data = pd.read_csv('C:/Users/Admin/Desktop/TESI/WRF_DATA/Data/sondaggio WRF/analisi sondaggi/sounding_all_variables.txt',delim_whitespace=True)

# Change default to be better for skew-T

p = data['P'].values*units.hPa
T = data['Temp'].values*units.degC
Td = data['Dewp'].values*units.degC
wind_speed = data['Speed'].values*units('m/s')
wind_speed = wind_speed.to(units.knots)
wind_dir = data['Dir'].values*units.degrees
u, v = mpcalc.wind_components(wind_speed,wind_dir)

fig = plt.figure(figsize=(20, 14))
gs = gridspec.GridSpec(2, 2,hspace=0)
gs.update(wspace = 0, hspace = .25)
skew = SkewT(fig, rotation=45, subplot=gs[:, :1])


#skew=SkewT(fig)

# Plot the data using normal plotting functions, in this case using
# log scaling in Y, as dictated by the typical meteorological plot
skew.plot(p, T, 'r',linewidth=2.5)
skew.plot(p, Td, 'g',linewidth=2.5)
skew.ax.set_ylim(1000, 100)
skew.ax.set_xlim(-50, 30)

# Add the relevant special lines
skew.plot_dry_adiabats()
skew.plot_moist_adiabats()
skew.plot_mixing_lines()
skew.ax.set_ylim(1000, 100)

# Defining vertucal barb spacing
interval = np.arange(100,1000,43)*units('mbar')
# Get indexes of values closest to defined interval
ix = mpcalc.resample_nn_1d(p,interval)

skew.plot_barbs(p[ix],u[ix],v[ix],length=7.8,xloc=1)


# calculate LCL level as a black dot
lcl_pressure, lcl_temperature = mpcalc.lcl(p[0], T[0], Td[0])
skew.plot(lcl_pressure, lcl_temperature, marker="_", color='blue', markersize=30, markeredgewidth=3,label='LCL')

# Calculate full parcel profile and add to plot as black line
prof = mpcalc.parcel_profile(p, T[0], Td[0]).to('degC')
skew.plot(p, prof, 'k', linewidth=2.5)

#calculate LFC
lfc_p, lfc_t = mpcalc.lfc(p,T,Td,which = 'bottom')
skew.plot(lfc_p, lfc_t,marker="_", color='red', markersize=30, markeredgewidth=3,label = 'LFC')

#calculate EL
el_p, el_t = mpcalc.lfc(p,T,Td)
skew.plot(el_p, el_t, marker="_", color='orange', markersize=30, markeredgewidth=3,label = 'EL')

# Shade areas of CAPE and CIN
# choosing only he first km of variables to avoid shade_cin mistakes
p1 = p[:1005]
T1= T[:1005]
prof1 = prof[:1005]
skew.shade_cin(p1, T1, prof1)
skew.shade_cape(p1, T1, prof1)

plt.xlabel('Temperature (°C)',fontsize = 30)
plt.ylabel('Pressure (hPa)',fontsize = 30)
plt.xticks(fontsize = 28)
plt.yticks(fontsize=28)
plt.legend(fontsize=30,loc='best',bbox_to_anchor=(0.4, 0.25))

cape_cin = metpy.calc.cape_cin(p, T, Td, prof, which_lfc='bottom', which_el='top')
CAPE = cape_cin[0].magnitude
CIN = cape_cin[1].magnitude
hght = data['HeightMSL']

skew.ax.set_xlim(-25, 30)

wind_speed = data['Speed'].values*units('m/s')
wind_dir = data['Dir'].values*units.degrees
u, v = mpcalc.wind_components(wind_speed,wind_dir)

hght = hght[0:1010]
u_new = u[0:1010]
v_new = v[0:1010]

from matplotlib import ticker
#ax_hod = inset_axes(skew.ax, width=4, height=4, loc=3,bbox_to_anchor=(550,550,1,1))
ax = fig.add_subplot(gs[0, 1])
h = Hodograph(ax,component_range=50.)
h.add_grid(increment=20)
l=h.plot_colormapped(u_new, v_new, hght)
cb = plt.colorbar(l)
cb.set_label('m a.s.l.', labelpad=-80, y=1.08, rotation=0,fontsize=30)
plt.xlabel('u (m/s)',Fontsize = 30,labelpad=10)
plt.ylabel('v (m/s)',Fontsize = 30,labelpad=-10)
plt.xticks(Fontsize=28)
plt.yticks(Fontsize=28)
h.wind_vectors(u_new[::100],v_new[::100],scale = 1.,width = 0.3)
majors = [-40,-20,0, 20, 40]
ax.xaxis.set_major_locator(ticker.FixedLocator(majors))
minors = np.linspace(0, 1, 11)[1:-1]
ax.xaxis.set_minor_locator(ticker.FixedLocator(minors))




ax = fig.add_subplot(gs[1, 1])
ax.axis('off')
col_labels=['Values']
row_labels=['CAPE (J/kg)','CIN (J/kg)','EL (hPa)','LCL (hPa)','LFC (hPa)']
table_vals=[[704],[43],[233],[950],[876]]
# the rectangle is where I want to place the table
table = plt.table(cellText=table_vals,
                  colWidths = [0.1]*5,
                  cellLoc='center',
                  rowLabels=row_labels,
                  rowLoc='center',
                  colLabels=col_labels,
                  bbox = [0.47, 0.3, 0.35, 0.6],
                  loc='center')

table.set_fontsize(26)
table.scale(3.5, 3.5) 

plt.text(-1.13, 2.35, 'a)', style='italic', fontsize=44)
plt.show()


#%% making plot of vertical profiles
hght = data['HeightMSL']
thetae = data['thetae']
RH = data['RH']
speed = data['Speed']

fig=plt.figure()
ax=fig.add_subplot(111, label="1")
ax2=fig.add_subplot(111, label="2", frame_on=False)
ax3=fig.add_subplot(111, label="3", frame_on=False)
ax3 = ax.twiny()

ax.plot(thetae, hght, color="C0")
ax.set_xlabel("$\Theta_e$ (K)", color="C0")
ax.set_ylabel("Height (m)", color="C0")
ax.tick_params(axis='x', colors="C0")
ax.tick_params(axis='y', colors="C0")
ax.set_xlim((305,335))
ax.set_ylim((0,8000))

ax2.plot(RH, hght, color="C1")
ax2.xaxis.tick_top()
ax2.yaxis.tick_right()
ax2.set_ylim((0,8000))
ax2.set_xlabel('RH (%)', color="C1") 
ax2.set_ylabel('Height (m)', color="C1")       
ax2.xaxis.set_label_position('top') 
ax2.yaxis.set_label_position('right') 
ax2.tick_params(axis='x', colors="C1")
ax2.tick_params(axis='y', colors="C1")

ax3.plot(speed, hght, color="C1")
ax3.xaxis.tick_bottom()
ax3.yaxis.tick_right()
ax3.set_ylim((0,8000))
ax3.set_xlabel('Wind Speed (m/s)', color="C1",labelpad=40) 
ax3.set_ylabel('Height (m)', color="C1")       
ax3.xaxis.set_label_position('bottom') 
ax3.yaxis.set_label_position('right') 
ax3.tick_params(axis='x', colors="C1")
ax3.tick_params(axis='y', colors="C1")


#%%
import numpy as np
import matplotlib.pyplot as plt

fig = plt.figure(figsize=(14,12))
plt.rcParams.update({'font.size': 26})
ax1 = fig.add_subplot(111)
ax2 = fig.add_subplot(111,frame_on=False)
ax3=fig.add_subplot(111, label="3", frame_on=False)


# Add some extra space for the second axis at the bottom
fig.subplots_adjust(bottom=0.2)

ax1.plot(thetae,hght,linewidth=2.5,color='blue',alpha=0.6)
ax1.set_xlabel("$\Theta_e$ (K)",color='blue')
ax1.set_ylabel("Height (m)", color="black")
ax1.tick_params(axis='x',colors='blue')
ax1.tick_params(axis='y')
ax1.set_xlim((310,340))
ax1.set_ylim((0,8000))


ax2.plot(RH,hght,color='red',linestyle='dashed',dashes=(6, 4),linewidth=2.5)
# Move twinned axis ticks and label from top to bottom
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")

# Offset the twin axis below the host
ax2.spines["bottom"].set_position(("axes", -0.15))

# Turn on the frame for the twin axis, but then hide all 
# but the bottom spine
ax2.set_frame_on(True)
ax2.patch.set_visible(False)
ax2.set_xlabel("RH (%)",color='red')
ax2.tick_params(axis='x', colors="red")

ax2.set_ylim((0,8000))

ax3.plot(speed, hght, color="green",linestyle='dotted',linewidth=4.)
ax3.xaxis.tick_top()
ax3.set_yticks([],[])
ax3.set_ylim((0,8000))
ax3.set_xlabel('Wind Speed (m/s)', color="green",labelpad=20) 
ax3.set_ylabel('', color="green")       
ax3.xaxis.set_label_position('top') 
ax3.yaxis.set_label_position('right') 
ax3.tick_params(axis='x', colors="green")
ax3.tick_params(axis='y', colors="green")



plt.show()