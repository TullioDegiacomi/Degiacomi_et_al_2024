# -*- coding: utf-8 -*-
"""
Created on Tue Jan 18 10:17:06 2022

@author: Admin
"""
import pandas as pd
import numpy as np
from matplotlib import rc
rc('text.latex', preamble=r'\usepackage{color}')
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

data = pd.read_csv('C:/Users/Admin/Desktop/TESI/WRF_DATA/Immagini/Sounding_idealizzati/Idealized_Wind_Profile/CTRL_SIMULATION/Parameters_CTRL_sounding.txt',delimiter = '\t')
data2 = pd.read_csv('C:/Users/Admin/Desktop/TESI/WRF_DATA/Data/sondaggio WRF/analisi sondaggi/sounding_all_variables.txt',delim_whitespace=True)

#%%
theta = data.Theta
hght = data.Height
theta_e = data.theta_e
RH = data.RH
speed = data.v

hght2 = data2['HeightMSL']
theta2 = data2['theta']
thetae2 = data2['thetae']
RH2 = data2['RH']
speed2 = data2['Speed']



import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator,FixedLocator

fig = plt.figure(figsize=(14,13))
plt.rcParams.update({'font.size': 30})
ax1 = fig.add_subplot(111)
ax2 = fig.add_subplot(111,frame_on=False)
ax3=fig.add_subplot(111, label="3", frame_on=False)


# Add some extra space for the second axis at the bottom
fig.subplots_adjust(bottom=0.2)
#ax1.text(291, 7500, 'Transparent = Original', fontsize = 26)
#ax1.text(291, 7000, 'Thick = CTRL',fontsize = 26)

#ax1.plot(theta_e,hght,linewidth=4,color='blue',alpha=0.6)
#ax1.plot(thetae2,hght2,linewidth=2.5,color='blue',alpha=0.3)
ax1.plot(speed,hght,color='gray',linewidth=4)
ax1.plot(speed2,hght2,color='gray',linewidth=2.5,alpha=0.5)
ax1.set_xlabel("Wind Speed (m/s)",color='gray')
ax1.set_ylabel("Height (m)", color="black")
ax1.tick_params(axis='x',colors='gray')
ax1.tick_params(axis='y')
ax1.set_xlim((0,100))
ax1.set_ylim((0,8000))


#ax2.plot(RH,hght,color='red',linestyle='dashed',dashes=(10, 4),linewidth=2)
#ax2.plot(RH2,hght2,color='red',linestyle='dashed',dashes=(10, 4),linewidth=2)
ax2.plot(RH,hght,color='red',linewidth=4)
ax2.plot(RH2,hght2,color='red',linewidth=2.5,alpha=0.5)

# Move twinned axis ticks and label from top to bottom
ax2.xaxis.set_ticks_position("bottom")
ax2.xaxis.set_label_position("bottom")

# Offset the twin axis below the host
ax2.spines["bottom"].set_position(("axes", -0.15))

# Turn on the frame for the twin axis, but then hide all 
# but the bottom spine
ax2.set_frame_on(True)
ax2.patch.set_visible(False)
x_ticks = [0 + 20 * i for i in range((100 - 0) // 20 + 1)]
ax2.xaxis.set_major_locator(FixedLocator(x_ticks))
#ax2.xaxis.set_major_locator(MaxNLocator(integer=True, prune='both'))
#ax2.set_xticks([tick for tick in RH if tick >= 0])
ax2.text(55, -1950, 'RH (%)', color='red', ha='center')
ax2.set_xlabel("",color='red')
ax2.tick_params(axis='x', colors="red")
ax2.set_xlim((-100,110))

ax2.set_ylim((0,8000))

ax3.plot(theta, hght, color="green",linewidth=4)
ax3.plot(theta2, hght2, color="green",linewidth=2.5,alpha = 0.5)
ax3.plot(theta_e,hght,linewidth=4,color='blue',alpha=0.6)
ax3.plot(thetae2,hght2,linewidth=2.5,color='blue',alpha=0.3)
ax3.xaxis.tick_top()
#ax3.set_yticks([],[])
ax3.set_ylim((0,8000))
ax3.set_xlim((270,345))
#ax3.set_xlabel('$\Theta$, $\Theta_e$ (K)', color="green",labelpad=20) 

label_part1 = r'$\theta$ (K), '
label_part2 = r'$\theta_e$ (K)'
ax3.set_xlabel('')  # Clear the default x-label
ax3.text(0.44, 1.08, label_part1, color='green', transform=ax3.transAxes, ha='center')
ax3.text(0.54, 1.08, label_part2, color='blue', transform=ax3.transAxes, ha='center')
ax3.text(266, 8600, 'b)',style = 'italic', color='black', ha='center',fontsize=36)



ax3.set_ylabel('', color="black")       
ax3.xaxis.set_label_position('top') 
ax3.yaxis.set_label_position('right') 
ax3.tick_params(axis='x', colors="green")
ax3.tick_params(axis='y', colors="black")

ax3.text(305.2, 2500, '$N^{2}_{3}$ = 0.00008 s$^{-2}$', style='italic', fontsize=24,
          bbox={'facecolor': 'grey', 'alpha': 0.5})

ax3.text(302., 1400, '$N^{2}_{2}$ = 0.00019 s$^{-2}$', style='italic', fontsize=24,
          bbox={'facecolor': 'grey', 'alpha': 0.5})

ax3.text(298, 500, '$N^{2}_{1}$ = 0.00008 s$^{-2}$', style='italic', fontsize=24,
          bbox={'facecolor': 'grey', 'alpha': 0.5})


plt.show()

#%%
plt.plot(hght,speed)
plt.plot(hght2,speed2)