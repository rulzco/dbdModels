#!/usr/bin/env python
"""
Created on Sun May  9 16:00:07 2021

@author: ruloz
"""
# import readmesh function from fluidfoam package
import os
from fluidfoam import readmesh
from scipy.interpolate import griddata
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from matplotlib import ticker

mainRoute = os.getcwd() # get current directory route
case = os.path.basename(mainRoute) # case name
sol = mainRoute
x, y, z = readmesh(sol)

dir_list = os.listdir(mainRoute) # get list of directories inside case path

times = []
for i in range(len(dir_list)):
    if(dir_list[i].isdigit()): 
         times.append(dir_list[i]) # get time directories
    else:
        try:
            float(dir_list[i])
            times.append(dir_list[i])
        except:
            print("Not a float or int")
   # if(dir_list[i].isdigit()): 
   #      times.append(dir_list[i]) # get time directories
timename = max(times) # get latest time      

# import readvector and readscalar functions from fluidfoam package
from fluidfoam import readvector, readscalar, readfield

timename = '149'
vel = readvector(sol, timename, 'U')
F = readscalar(sol, '100', 'Force')
#bForce = readvector(sol, timename, 'bForce')

# Number of division for linear interpolation
ngridx = 300
ngridy = 200

# Interpolation grid dimensions
xinterpmin = -0.015
xinterpmax = 0.1
yinterpmin = 0
yinterpmax = 0.05

xinterpmin2 = -0.001
xinterpmax2 = 0.005
yinterpmin2 = 0
yinterpmax2 = 0.003

# Interpolation grid
xi = np.linspace(xinterpmin, xinterpmax, ngridx)
yi = np.linspace(yinterpmin, yinterpmax, ngridy)

xi2 = np.linspace(xinterpmin2, xinterpmax2, ngridx)
yi2 = np.linspace(yinterpmin2, yinterpmax2, ngridy)

# Structured grid creation
xinterp, yinterp = np.meshgrid(xi, yi)
xinterp2, yinterp2 = np.meshgrid(xi2, yi2)

# Interpolation of scalar fields and vector field components
F = griddata((x, y), F, (xinterp2, yinterp2), method='linear')
velx_i = griddata((x, y), vel[0, :], (xinterp, yinterp), method='linear')
vely_i = griddata((x, y), vel[1, :], (xinterp, yinterp), method='linear')
# bForcex = griddata((x, y), bForce[0, :], (xinterp2, yinterp2), method='linear')
# bForcey = griddata((x, y), bForce[1, :], (xinterp2, yinterp2), method='linear')

# Calculation of the streamline width as a function of the velociy magnitude
vel_i = np.sqrt(velx_i**2 + vely_i**2)
#lw = pow(vel_i, 1.5)/vel_i.max()

# bForce = np.sqrt(bForcex**2 + bForcey**2)

#NaN to zeros
vel_i = np.nan_to_num(vel_i)
F = np.nan_to_num(F)

te = 1.0e-4 # thickness of electrodes m  
le = 5.0e-4 # electrodes length m
# ye = 1.274e-4 # depth of emebded electrode m
# xe = 0.0 # distance between electrodes m

# Define plot parameters
cm = ['coolwarm','jet','rainbow','Blues','Greens','Greys','viridis']
plt.rcParams.update({'font.size':10})

# fig, ax = plt.subplots(1, figsize=(8, 5), dpi=90, facecolor='w', edgecolor='w') 
# plt.contourf(xi,yi,vel_i,100,cmap='jet')
# plt.xlim(-0.012, 0.03) #x-axis range
# plt.ylim(0.0, 0.01) #y-axis range
# plt.xlabel('x [m]')
# plt.ylabel('y [m]')
# # ax.set_aspect('equal')
# cbar = plt.colorbar(orientation='horizontal', shrink=1.0, fraction=0.050)
# cbar.set_label('U [m/s]')
# #plt.savefig(sol+'/'+case+'vel.png')


# fig, ax = plt.subplots(1, figsize=(8, 6), dpi=90, facecolor='w', edgecolor='w') 
# cs = ax.contourf(xi2,yi2,F,30,cmap=cm[0]) 
# # plt.plot([xinterpmin2,xinterpmax2],[0,0],'-',color='k') 
# #ax.contour(X, Y, F, 10, colors = "grey")
# rect1 = matplotlib.patches.Rectangle((-le, 0), 
                                      # le, te, 
                                      # color ='black')
# # rect2 = matplotlib.patches.Rectangle((xe, -(ye+te)), 
                                     # # (xe+le), te, 
                                      # # color ='black')  
# ax.add_patch( rect1 ) 
# # ax.add_patch( rect2 )
# #plt.locator_params(axis='y', nbins=6)
# plt.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
# plt.xlim(-0.0006, 0.0035) #límites del eje x
# plt.ylim(0, 0.003) #límites del eje y 
# # ax.set_aspect('equal')
# cbar = fig.colorbar(cs)
# cbar.set_label(label='F [N/m$^3$]', )
# plt.xlabel('x [m]');plt.ylabel('y [m]')
# #plt.savefig(sol+'/'+case+'force.png')


# # Single figure for velcoity and force

fig, (ax1, ax2) = plt.subplots(2, 1,figsize=(8, 6), dpi=90)
im = ax1.contourf(xi,yi,vel_i,100,cmap='jet')
ax1.title.set_text('Velocity field')
ax1.set_xlabel('x [m]')
ax1.set_ylabel('y [m]')
ax1.set_xlim([-0.012, 0.03])
ax1.set_ylim([0.0, 0.01])
cbar = fig.colorbar(im, ax=ax1)
# tick_locator = ticker.MaxNLocator(nbins=5)
# cbar.locator = tick_locator
# cbar.update_ticks()
cbar.set_label('U [m/s]')

im2 = ax2.contourf(xi2,yi2,F,100,cmap='jet')
# ax1.streamplot(xi,yi,velx_i,vely_i,color='k',
               # linewidth=0.6,density=[0.8, 1]) # stream lines 
# plt.plot([xinterpmin2,xinterpmax2],[0,0],'-',color='k')
ax2.title.set_text('Force')
ax2.set_xlabel('x [m]')
ax2.set_ylabel('y [m]') 
rect1 = matplotlib.patches.Rectangle((-le, 0), le, te, color ='black')
# rect2 = matplotlib.patches.Rectangle((xe, -(ye+te)), (xe+le), te, color ='black')  
ax2.add_patch( rect1 ) 
# ax2.add_patch( rect2 )
ax2.set_xlim([-0.0006, 0.0035])
ax2.set_ylim([0.0, 0.003])
ax2.ticklabel_format(axis="both", style="sci", scilimits=(0,0))
cbar2 = fig.colorbar(im2, ax=ax2)
# tick_locator = ticker.MaxNLocator(nbins=7)
# cbar2.locator = tick_locator
# cbar2.update_ticks()
cbar2.set_label('F [N/m$^3$]')
fig.tight_layout()
plt.savefig(sol+'/'+case+'.png')
plt.show() 

print('U max:{0:.4f} m/s'.format(np.amax(vel_i/5))) 
#print('U max:{0:.4f} m/s, Fbmax:{1:.4f} N/m^3'.format(np.amax(vel_i), np.amax(bForce))) 