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
path_output = os.path.abspath("/archive/tullio/scripts/")
# stability plot with potential temperature and equivalent potential temp.
data = pd.read_csv('sounding_all_variables.txt',delim_whitespace=True)

# Change default to be better for skew-T

hght = data.HeightMSL
p= data.P
T = data.Temp
Td = data.Td
wind_speed = data.Speed


p = p.values*units.hPa
T = T.values*units.degC
Td = Td.values*units.degC
wind_speed_ms = wind_speed.values*units('m/s')
wind_speed = wind_speed_ms.to(units.knots)
wind_dir = data['Dir'].values*units.degrees
u, v = mpcalc.wind_components(wind_speed_ms,wind_dir)


fig = plt.figure(figsize=(20, 14))
gs = gridspec.GridSpec(2, 2,hspace=0)
gs.update(wspace = 0, hspace = .25)
skew = SkewT(fig, rotation=45, subplot=gs[:, :1])


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
p1 = p[:2730]
p2 = p[:1000]
T1= T[:2730]
T2 = T[:1000]
prof1 = prof[:2730]
prof2 = prof[:1000]
skew.shade_cin(p2, T2, prof2)
skew.shade_cape(p1, T1, prof1)

plt.xlabel('Temperature [$^{\circ}$C]',fontsize = 30)
plt.ylabel('Pressure [hPa]',fontsize = 30)
plt.xticks(fontsize = 28)
plt.yticks(fontsize=28)
#plt.legend(fontsize=30,loc='best',bbox_to_anchor=(0.4, 0.25))

cape_cin = metpy.calc.cape_cin(p, T, Td, prof, which_lfc='bottom', which_el='top')
CAPE = cape_cin[0].magnitude
CIN = cape_cin[1].magnitude

skew.ax.set_xlim(-25, 30)


hght = hght[0:2730]
u_new = u[0:2730]
v_new = v[0:2730]

from matplotlib import ticker
#ax_hod = inset_axes(skew.ax, width=4, height=4, loc=3,bbox_to_anchor=(550,550,1,1))
ax = fig.add_subplot(gs[0, 1])
h = Hodograph(ax,component_range=50.)
h.add_grid(increment=20)
l=h.plot_colormapped(u_new, v_new, hght)
cb = plt.colorbar(l)
cb.ax.tick_params(labelsize=28)
cb.set_label('m MSL', labelpad=-80, y=1.08, rotation=0,fontsize=30)
plt.xlabel('$\it{u}$ [m s$^{-1}$]',fontsize = 30,labelpad=10)
plt.ylabel('$\it{v}$ [m s$^{-1}$]',fontsize = 30,labelpad=-10)
plt.xticks(fontsize=28)
plt.yticks(fontsize=28)
h.wind_vectors(u_new[::100],v_new[::100],scale = 1.,width = 0.3)
majors = [-40,-20,0, 20, 40]
ax.xaxis.set_major_locator(ticker.FixedLocator(majors))
minors = np.linspace(0, 1, 11)[1:-1]
ax.xaxis.set_minor_locator(ticker.FixedLocator(minors))


ax = fig.add_subplot(gs[1, 1])
ax.axis('off')
col_labels=['Values']
row_labels=['CAPE [J kg$^{-1}$]','CIN [J kg$^{-1}$]','EL [hPa]','LCL [hPa]','LFC [hPa]']
table_vals=[[876],[43],[233],[950],[876]]
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

plt.text(-1.13, 2.35, '(a)', fontsize=44)

plt.savefig('figure2a.png')



