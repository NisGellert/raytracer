import os
os.chdir(r"C:\Users\nige\DTU Space\PhD\Projects\Raytracing")
from IPython import get_ipython
get_ipython().magic('reset -sf') # Clear all
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import random
import copy
import math
import bisect
import linecache
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from functions_raytrace import *


'''
# ===============================================================
# Author: Nis C. Gellert
# DTU Space
# Last updated: September 2023
# The following scripts computes the on- and off-axis performance of HEX-P
# ===============================================================
# ----- Read before use: ----
# Generate the Optics Geometry file using "compute_telescope_geometry.py"


# ===============================================================

 
'''





# To Do: 
# Plot distribution of killedflag and reflection flag
# 

# 0,0.5,1,1.5,2,2.5, 3,3.5,4,4.5,5,5.5,6,6.5,6.75,7
min_r = 0.0595# [m], Mast radius in the optics, which is the minimum radius. Must be put in optics file! TO DO!

''' ... Define Source ...'''
Its = 300e3  # Number of photons aka. iterations
SrcRadius = 0.3 # [m], radius of the circle/source
Src_x = 0.00 # [m] center of the circle/source (x, y)
Src_y = 0 # [m] center of the circle/source (x, y)
Src_z = -0.15 # [m] Photons start 15 cm behind optics 
Src_oaa = -0.0 #.5 -4.0 # # Source Off-axis angle, x-plane [arcmin] 
Src_azimuth = 0.0 # Source Off-axis angle, y-plane [arcmin]
PlotResult = True # True/False, plot results
scatter =  1.352e-5  # Scatter [unitless]     # For HPD of 5 arcsec (5'') = 0.0833 arcmin (0.083''), scatter = 1.352e-5

''' ... Define Aperture ...'''
Aperture = True
dia_detector = 0.08 # [m] detector diameter
dia_aperture_1 = 0.126 # [m] Inner aperture_diameter
dia_outer_aperture_2 = dia_aperture_1+0.06  # [m] Outer aperture_diameter
dis_aperture = 2.0 # [m] aperture detector distance
thickness_aperture  = 0.00188 # [m] Thickness of aperture disk (From Nustar)
Max_reflection = 60 # 60 Maximum number of reflection before killed

''' ... Load Optics Geometry File...'''
fname = "C:/Users/nige/DTU Space/PhD/Projects/Raytracing/Telescope geometries/Geometry_HEXP_20230902.txt"
shells = np.loadtxt(fname)[:,0] # Shells
r_ut = np.loadtxt(fname)[:,1] # Radius to top of the upper shell, [m].
r_ub = np.loadtxt(fname)[:,2] # Radius to bottom of the upper shell, [m].
r_lt = np.loadtxt(fname)[:,3] # Radius to top of the lower shell, [m].
r_lb = np.loadtxt(fname)[:,4] # Radius to bottom of the lower shell, [m].
alpha = np.loadtxt(fname)[:,5]/1000 # Grazing incident angle, [rad]
z_ut = np.loadtxt(fname)[:,6] # z-distance to top of the upper shell, [m].
z_ub = np.loadtxt(fname)[:,7] # z-distance to bottom of the upper shell, [m].
z_lt = np.loadtxt(fname)[:,8] # z-distance to top of the lower shell, [m].
z_lb = np.loadtxt(fname)[:,9] # z-distance to bottom of the lower shell, [m].
x = linecache.getline(fname, 3).strip(); FL = float(x.split(":")[1]) # Focal length, [m]
x = linecache.getline(fname, 4).strip(); ML = float(x.split(":")[1]) # Mirror length, [m]
x = linecache.getline(fname, 5).strip(); N_shells = float(x.split(":")[1]) # Numbers of shells
x = linecache.getline(fname, 6).strip(); gap = float(x.split(":")[1]) # Optics gap, [m]
x = linecache.getline(fname, 8).strip(); ST = float(x.split(":")[1]) # Substrate thickness, [m]
#x = linecache.getline(fname, 9).strip(); r1_uo = float(x.split(":")[1]) # Radius to the center of the upper shell, [m].

''' --- Raytracer start --- '''
KilledFlag = np.empty(int(Its)); KilledFlag[:] = np.nan
ReflectionFlag = np.empty(int(Its)) ; ReflectionFlag[:] = np.nan
N_reflections = np.zeros(int(Its))

r = np.zeros(int(Its)) # radius of ray
P0 = np.empty([int(Its),7]); P0[:,:] = np.nan  # (x,y,z,r,alpha,alpha,alpha) Position P0, at source
P1 = np.empty([int(Its),7]); P1[:,:] = np.nan  # (x,y,z,r,alpha,alpha,alpha) Position P1, at first mirror intersection
P2 = np.empty([int(Its),7]); P2[:,:] = np.nan  # (x,y,z,r,alpha,alpha,alpha) Position P2, at second mirror intersection
P3 = np.empty([int(Its),7]); P3[:,:] = np.nan  # (x,y,z,r,alpha,alpha,alpha) Position P3, at focal point
ACP = np.empty([int(Its),7]); P3[:,:] = np.nan  # (x,y,z,r,alpha,alpha,alpha) Aperture Check point 

shell = np.empty(int(Its)); shell[:] = np.nan # Photon is in which shell
# All this is for off axis!! Understand this!!!
yawrad =float((Src_oaa/60.)*(np.pi/180)) # [rad] angle of the photon
pitchrad = float((Src_azimuth/60.)*(np.pi/180)) # [rad] angle of the photon
noaa = 0.0
kphi = np.random.rand()*2*np.pi 
px = np.sin(noaa)*np.cos(kphi)
py = np.sin(noaa)*np.sin(kphi)
pz = np.cos(noaa)
k1 = px*np.cos(yawrad)+pz*np.sin(yawrad) #;	k1 = (px*np.cos(yawrad)+pz*np.sin(yawrad))*np.cos(pitchrad) - py*np.sin(pitchrad) # Kristines code
k2 = py*np.cos(pitchrad) - np.sin(pitchrad)*(pz*np.cos(yawrad)-px*np.sin(yawrad)) # k2 = (px*np.cos(yawrad)+pz*np.sin(yawrad))*np.sin(pitchrad) + py*np.cos(pitchrad) # Kristines code
k3 = py*np.sin(pitchrad) + np.cos(pitchrad)*(pz*np.cos(yawrad)-px*np.sin(yawrad)) #k3 = -px*np.sin(yawrad) + pz*np.cos(yawrad) # Kristines code

i = -1
while i < int(Its)-1: #for i in range(int(Its)): 
    i = i + 1
    ReflectionFlag[i] = np.nan
    KilledFlag[i] = np.nan
    N_reflections[i] = 0
    ''' --- Position P0 --- '''
    ran_angle = 2 * math.pi * random.random() # random angle 
    P0[i,3] =  SrcRadius * math.sqrt(random.random()) # Random radius, [m] 
    r[i] = P0[i,3]
    P0[i,0] =  P0[i,3] * math.cos(ran_angle) + Src_x # x [m] 
    P0[i,1] =  P0[i,3] * math.sin(ran_angle) + Src_y # y [m] 
    P0[i,2] = Src_z # z  [m]
    # Rotate the [x,y,z] according to yaw and pitch
    P0[i,0] = P0[i,0]*np.cos(yawrad)+P0[i,2]*np.sin(yawrad)  # x [m] 
    P0[i,1] = P0[i,1]*np.cos(pitchrad) - np.sin(pitchrad)*(P0[i,2]*np.cos(yawrad)-P0[i,0]*np.sin(yawrad))# y [m] 
    P0[i,2] = P0[i,1]*np.sin(pitchrad) + np.cos(pitchrad)*(P0[i,2]*np.cos(yawrad)-P0[i,0]*np.sin(yawrad)) # z [m] 
    P0[i,3] = np.sqrt(P0[i,0]**2+P0[i,1]**2)
    P0[i,4] = k1
    P0[i,5] = k2
    P0[i,6] = k3
    
    
    ''' --- Position CP1  ---''' 
    # The photon enters the upper optics at z = 0
    z = 0
    CP1 = movephoton(P0[i],z)    
    if CP1[3] > max(r_ut): # KilledFlag, photon outside of optics, radius > max(r_ut).
        #KilledFlag[i] = -2
        i = i-1 # Repeat
        continue # Skip iteration  
    if CP1[3] < min(r_lb): ## KilledFlag, photon outside of optics, radius < min(r_lb). NOTE/OBS: Only for on-axis
        KilledFlag[i] = -2
        #i = i-1 # Repeat
        continue # Skip iteration   
    shell[i] = bisect.bisect(r_ut,CP1[3])+1 # What shell the photon is in given the radius of the photon and the radius of all shells
    if shell[i] != 1: # There is no mirror infront of inner shell
        if CP1[3] < (r_ut[int(shell[i]-2)] + ST): # KilledFlag, the photon hits the substrate of shell[i]-1
            KilledFlag[i] = -2
            #i = i-1 # Repeat
            continue # Skip iteration

    if alpha[int(shell[i]-1)] == 0: # Photon hits gap between modules
        i = i-1 # Repeat
        continue # Skip iteration
    '''# Position P1 and CP2'''

    P1_z = []
    apex1 = []
    apex1_backside = []
    P1_z_backside = 0
    apex1 = r_ut[int(shell[i]-1)] / np.tan(alpha[int(shell[i]-1)]) # [m]
    P1_z = cone_reflection_point(CP1,apex1, alpha[int(shell[i]-1)])  # At what distance does the photon hit the cone

    if alpha[int(shell[i]-2)] != 0:
        apex1_backside = (r_ut[int(shell[i]-2)]+ST) / np.tan(alpha[int(shell[i]-2)]) # [m]
        P1_z_backside = cone_reflection_point(CP1,apex1_backside,alpha[int(shell[i]-2)])
        
    if  0. < P1_z_backside < z_ub[int(shell[i]-2)] and P1_z > P1_z_backside:
        ReflectionFlag[i] =  '5' # Backside reflection
        N_reflections[i] = N_reflections[i]+1
        P1[i] = movephoton(CP1,P1_z_backside)   
        (angle,P1[i,4],P1[i,5],P1[i,6])=reflection(P1[i],alpha[int(shell[i]-2)],scatter)     # NIS HERE! Alpha! -alpha or 360*(np.pi/180)-aplha
        CP2 = movephoton(P1[i],z_ub[int(shell[i]-2)])
    
    elif P1_z > z_ub[int(shell[i]-1)]: # Missed first reflection
        ReflectionFlag[i] = '2' # Missed first reflection
        P1[i] = movephoton(CP1,0)
        CP2 = movephoton(CP1,z_ub[int(shell[i]-1)]) # Move to CP2
 
    else: 
        ReflectionFlag[i] = '1' # first reflection completed
        P1[i] = movephoton(CP1,P1_z)
        (angle,P1[i,4],P1[i,5],P1[i,6])=reflection(P1[i],alpha[int(shell[i]-1)],scatter)
        CP2 = movephoton(P1[i],z_ub[int(shell[i]-1)])
        N_reflections[i] = N_reflections[i]+1
                            
        
        
    '''# Position CP2. Check if photon hits backside'''  
    #CP2: If photon is not between the shell[i] r_ub and r_ub-1, then check how many reflections is needed to reach.
    r_type =[]
    if ReflectionFlag[i] == 5.0:
        r_type = 1 # Backside reflection
        
    if ReflectionFlag[i] == 1.0: 
        r_type = 0 # Primary mirror reflection
        
    while  N_reflections[i] < Max_reflection and (r_ub[int(shell[i]-1)] < CP2[3] or CP2[3] < r_ub[int(shell[i]-2)]+ST) and shell[i] !=1 and shell[i] !=57:   
        apex1 = r_ut[int(shell[i]-1)] / np.tan(alpha[int(shell[i]-1)]) # [m]
        P1_z = cone_reflection_point(P1[i],apex1,alpha[int(shell[i]-1)])  # At what distance does the photon hit the cone
        
        apex1_backside = (r_ut[int(shell[i]-2)]+ST) / np.tan(alpha[int(shell[i]-2)]) # [m]
        P1_z_backside = cone_reflection_point(P1[i],apex1_backside,alpha[int(shell[i]-2)])
            
        if r_type == 1 :#and (P1[i,2] < P1_z < z_ub[int(shell[i]-1)]): # Hit backside, now hits primary
            P1[i] = movephoton(P1[i],P1_z)
            (angle,P1[i,4],P1[i,5],P1[i,6])=reflection(P1[i],alpha[int(shell[i]-1)],scatter)
            CP2 = movephoton(P1[i],z_ub[int(shell[i]-1)])
            N_reflections[i] = N_reflections[i]+1
            r_type = 0 
            ReflectionFlag[i] =  '5' # Backside reflection
                     
        elif r_type == 0 :#and (P1[i,2] < P1_z_backside < z_ub[int(shell[i]-2)]): # Hit Primary, now hits backside
            P1[i] = movephoton(P1[i],P1_z_backside)
            (angle,P1[i,4],P1[i,5],P1[i,6])=reflection(P1[i],alpha[int(shell[i]-2)],scatter)
            CP2 = movephoton(P1[i],z_ub[int(shell[i]-1)])
            N_reflections[i] = N_reflections[i]+1
            r_type = 1
            ReflectionFlag[i] =  '5' # Backside reflection
            
        else:
            N_reflections[i] = Max_reflection
            KilledFlag[i] = 1 # photon killed at structure inside MM 
            ReflectionFlag[i] =  '55' # Backside reflection
            #i = i-1 # Repeat
            continue # Skip iteration  
            
    '''# Position CP3'''    
    CP3 = movephoton(CP2,z_lt[int(shell[i]-1)])
    if shell[i]==1: # Inner shell, photons hits center of optics
       if CP3[3] < min_r: # Need to define at top
            #KilledFlag[i] = 2
            #i = i-1 # Repeat
            continue # Skip iteration
    if shell[i]!=1:
        if r_lt[int(shell[i]-1)] < CP3[3] or CP3[3] < r_lt[int(shell[i]-2)]+ST: # Photon killed in the gap
            #KilledFlag[i] = 2
            #i = i-1 # Repeat
            continue # Skip iteration
    
    ''' --- Position P2 and CP4  ---''' 
    P2_z = []
    apex2 = r_lt[int(shell[i]-1)] / np.tan(3*alpha[int(shell[i]-1)]) + z_lt[int(shell[i]-1)] # Is it 3* or 4*?
    P2_z = cone_reflection_point(CP3,apex2,3*alpha[int(shell[i]-1)]) # Is it 3* or 4*?
    apex2_backside = []
    P2_z_backside = 0
    
    if alpha[int(shell[i]-2)] != 0:
        apex2_backside = (r_lt[int(shell[i]-2)]+ST) / np.tan(3*alpha[int(shell[i]-2)]) + z_lt[int(shell[i]-2)]  # [m]
        P2_z_backside = cone_reflection_point(CP3,apex2_backside,3*alpha[int(shell[i]-2)])

    if  z_lt[int(shell[i]-2)] < P2_z_backside < z_lb[int(shell[i]-2)] and P2_z > P2_z_backside:
        ReflectionFlag[i] =  '5' # Backside reflection
        N_reflections[i] = N_reflections[i]+1
        P2[i] = movephoton(CP3,P2_z_backside)   
        (angle,P2[i,4],P2[i,5],P2[i,6])=reflection(P2[i],3*alpha[int(shell[i]-2)],scatter)     # NIS HERE! Alpha! -alpha or 360*(np.pi/180)-aplha
        CP4 = movephoton(P2[i],z_lb[int(shell[i]-2)])    
    
    elif P2_z > z_lb[int(shell[i]-1)]: # Missed second reflection
        ReflectionFlag[i] = str(int(ReflectionFlag[i]))+'4'
        P2[i] = movephoton(CP3,z_lt[int(shell[i]-1)])
        CP4 = movephoton(CP3,z_lb[int(shell[i]-1)]) # Move to CP4
    # elif ReflectionFlag[i] == '2': # Second reflection completed and therefore lower bounce reflection
    #       
    else: 
        ReflectionFlag[i] = str(int(ReflectionFlag[i]))+'3' # Second reflection completed
        N_reflections[i] = N_reflections[i]+1
        P2[i] = movephoton(CP3,P2_z)
        (angle,P2[i,4],P2[i,5],P2[i,6])=reflection(P2[i],3*alpha[int(shell[i]-1)],scatter)# Is it 3* or 4*?
        CP4 = movephoton(P2[i],z_lb[int(shell[i]-1)])
        
    ''' CP4: Check if photon hits backside '''
    
    r_type =[]
    if ReflectionFlag[i] == 5.0:
        r_type = 1 # Backside reflection
        
    if ReflectionFlag[i] == 13 or ReflectionFlag[i] == 23 or ReflectionFlag[i] == 53: 
        r_type = 0 # Primary mirror reflection
        
        # If photons hits mirror before exiting: 
    if shell[i]==1 and CP4[3] < min_r: # Inner shell, photons hits center of optics
        KilledFlag[i] = 2
        #i = i-1 # Repeat
        continue # Skip iteration
   
            
    while  N_reflections[i] < Max_reflection and (r_lb[int(shell[i]-1)] < CP4[3] or CP4[3] < r_lb[int(shell[i]-2)]+ST) and shell[i] !=1 and shell[i] !=57:   
        apex2 = r_lt[int(shell[i]-1)] / np.tan(3*alpha[int(shell[i]-1)]) + z_lt[int(shell[i]-1)] # [m]
        P2_z = cone_reflection_point(P2[i],apex2,3*alpha[int(shell[i]-1)])  # At what distance does the photon hit the cone
        apex2_backside = (r_lt[int(shell[i]-2)]+ST) / np.tan(3*alpha[int(shell[i]-2)]) + z_lt[int(shell[i]-2)] # [m]
        P2_z_backside = cone_reflection_point(P2[i],apex2_backside,3*alpha[int(shell[i]-2)])
            
        if r_type == 1 :#and (P1[i,2] < P1_z < z_ub[int(shell[i]-1)]): # Hit backside, now hits primary
            P2[i] = movephoton(P2[i],P2_z)
            (angle,P2[i,4],P2[i,5],P2[i,6])=reflection(P2[i],3*alpha[int(shell[i]-1)],scatter)
            CP4 = movephoton(P2[i],z_lb[int(shell[i]-1)])
            N_reflections[i] = N_reflections[i]+1
            r_type = 0 
            ReflectionFlag[i] =  '5' # Backside reflection
                     
        elif r_type == 0 :#and (P1[i,2] < P1_z_backside < z_ub[int(shell[i]-2)]): # Hit Primary, now hits backside
            P2[i] = movephoton(P2[i],P2_z_backside)
            (angle,P2[i,4],P2[i,5],P2[i,6])=reflection(P2[i],3*alpha[int(shell[i]-2)],scatter)
            CP4 = movephoton(P2[i],z_lb[int(shell[i]-1)])
            N_reflections[i] = N_reflections[i]+1
            r_type = 1
            ReflectionFlag[i] =  '5' # Backside reflection
            
        else:
            N_reflections[i] = Max_reflection
            KilledFlag[i] = 1 # photon killed at structure inside MM 
            ReflectionFlag[i] =  '55' # Backside reflection
            #i = i-1 # Repeat
            continue # Skip iteration  
            
    if shell[i]!=1:
        if r_lb[int(shell[i]-1)] < CP4[3] or CP4[3] < r_lb[int(shell[i]-2)]+ST:   
            KilledFlag[i] = 3 # photon killed at the backside of the Secondary
            ReflectionFlag[i] = 0
            #i = i-1 # Repeat
            continue # Skip iteration
    if ReflectionFlag[i] == 53 or ReflectionFlag[i] == 54: 
        ReflectionFlag[i] = '5'
        
    ''' After Optic modules '''
    z_FL = z_ub[0] + (gap/2) + FL
    z_aperture = z_FL - dis_aperture
    ''' --- Position APC  ---''' 
    if Aperture: # If photon hits aperture
        #z_aperture = (FL + z_lt[int(shell[i]-1)])-dis_aperture
        ACP[i] = movephoton(CP4,z_aperture)
        if ACP[i,3] >= dia_aperture_1/2: # If photon radius > aperture radius
            KilledFlag[i] = 5
            ReflectionFlag[i] ='55'
            #i = i-1 # Repeat
            continue # Skip iteration
    
    ''' --- Position P3  ---'''    
    P3[i] = movephoton(CP4,z_FL)  
    if (abs(P3[i,0]) > dia_detector/2) or (abs(P3[i,1]) > dia_detector/2):  # If photon is within detector r = 0.04 m
        KilledFlag[i] = 6

      
if Aperture and PlotResult: #  Aperture plane 
    idx13 = np.where(ReflectionFlag == 13) # double bounce photons.
    idx23 = np.where(ReflectionFlag == 23) # lower single bounce photons.
    idx14 = np.where(ReflectionFlag == 14) # Upper single bounce photons.
    idx5 = np.where(ReflectionFlag == 5) # back reflection
    
    
    fig, (ax1) = plt.subplots(nrows=1, ncols=1,figsize=(6,6))  
    fig.tight_layout(pad=4.0)
    ax1.plot(ACP[idx23,0], ACP[idx23,1], 'r.',markersize = 2)
    ax1.plot(ACP[idx14,0], ACP[idx14,1], 'b.',markersize = 2 ) # label="Upper single bounce"
    ax1.plot(ACP[idx13,0], ACP[idx13,1], 'k.',markersize = 2 ) # label="Double bounce"
    ax1.plot(ACP[idx5,0], ACP[idx5,1], 'g.',markersize = 2 ) # label="Back reflection"
    ax1.add_patch(plt.Circle((0,0), dia_aperture_1/2, edgecolor = 'gray', fill=False, lw=3))

    red_patch = mpatches.Patch(color='red', label="Single lower bounce")
    blue_patch = mpatches.Patch(color='blue', label="Single upper bounce")
    black_patch = mpatches.Patch(color='black', label="Double bounce")
    green_patch = mpatches.Patch(color='green', label="Backside reflection")
    gray_patch = mpatches.Patch(color='gray', label="Aperture outline")
    plt.legend(handles=[red_patch, blue_patch, black_patch, green_patch,gray_patch], loc = "upper right")

    ax1.set_ylabel("y [m]",fontsize=14)
    ax1.set_xlabel("x [m]",fontsize=14) 
    ax1.set_xlim(-(dia_aperture_1-0.026),(dia_aperture_1-0.026))
    ax1.set_ylim(-(dia_aperture_1-0.026),(dia_aperture_1-0.026))
    
    

if PlotResult: #  Detector plane photons 
    idx13 = np.where(ReflectionFlag == 13) # double bounce photons.
    idx23 = np.where(ReflectionFlag == 23) # lower single bounce 
    idx14 = np.where(ReflectionFlag == 14) # Upper single bounce 
    idx5 = np.where(ReflectionFlag == 5) # Upper single bounce 
    fig, (ax1) = plt.subplots(nrows=1, ncols=1,figsize=(6,6))  
    fig.tight_layout(pad=4.0)
    ax1.plot(P3[idx23,0], P3[idx23,1], 'r.',markersize = 2) # label="Lower single bounce"
    ax1.plot(P3[idx14,0], P3[idx14,1], 'b.',markersize = 2 ) # label="Upper single bounce"
    ax1.plot(P3[idx13,0], P3[idx13,1], 'k.',markersize = 2 ) # label="Double bounce"
    ax1.plot(P3[idx5,0], P3[idx5,1], 'g.',markersize = 2 ) # label="Backside reflection"
    #ax1,axes()
    ax1.add_patch(plt.Rectangle((-0.04, -0.04), 0.08, 0.08, edgecolor = 'gray', fill=False, lw=3))
    #ax1.legend(handlelength=6)
    red_patch = mpatches.Patch(color='red', label="Single lower bounce")
    blue_patch = mpatches.Patch(color='blue', label="Single upper bounce")
    black_patch = mpatches.Patch(color='black', label="Double bounce")
    magenta_patch = mpatches.Patch(color='green', label="Backside reflection")
    green_patch = mpatches.Patch(color='grey', label="Detector outline")
    plt.legend(handles=[red_patch, blue_patch, black_patch,magenta_patch, green_patch], loc = "upper right")
    ax1.set_ylabel("y [m]",fontsize=14)
    ax1.set_xlabel("x [m]",fontsize=14) 
    ax1.set_xlim(-0.08,0.08)
    ax1.set_ylim(-0.08,0.08)
    # Black = double bounce = reflectionFlag 13
    # Blue = upper single bounce = ReflectionFlag 14
    # Red = lower single bounce =  23



if PlotResult: #  Detector plane photons within detector
    
    idx23_hit = np.logical_and((ReflectionFlag == 23), (KilledFlag!=6))
    idx14_hit = np.logical_and((ReflectionFlag == 14), (KilledFlag!=6))
    idx13_hit = np.logical_and((ReflectionFlag == 13), (KilledFlag!=6)) 
    idx5_hit = np.logical_and((ReflectionFlag == 5), (KilledFlag!=6)) 
    print("Off-axis angle = " + str(Src_oaa))
    print("Photon within the detector:")
    print("DB = " + str(len(P3[idx13_hit,0])))
    print("UB = " + str(len(P3[idx14_hit,0])))
    print("LB = " + str(len(P3[idx23_hit,0])))
    print("BR = " + str(len(P3[idx5_hit,0])))
    fig, (ax1) = plt.subplots(nrows=1, ncols=1,figsize=(6,6))  
    fig.tight_layout(pad=4.0)
    ax1.plot(P3[idx23_hit,0], P3[idx23_hit,1], 'r.',markersize = 2) 
    ax1.plot(P3[idx14_hit,0], P3[idx14_hit,1], 'b.',markersize = 2) 
    ax1.plot(P3[idx5_hit,0], P3[idx5_hit,1], 'g.',markersize = 2 ) # label="Backside reflection"
    ax1.plot(P3[idx13_hit,0], P3[idx13_hit,1], 'k.',markersize = 2 ) # label="Double bounce"
    ax1.add_patch(plt.Rectangle((-0.04, -0.04), 0.08, 0.08, edgecolor = 'gray', fill=False, lw=3))
    red_patch = mpatches.Patch(color='red', label="Single lower bounce")
    blue_patch = mpatches.Patch(color='blue', label="Single upper bounce")
    magenta_patch = mpatches.Patch(color='green', label="Backside reflection")
    black_patch = mpatches.Patch(color='black', label="Double bounce")
    green_patch = mpatches.Patch(color='gray', label="Detector outline")
    plt.legend(handles=[red_patch, blue_patch, black_patch,magenta_patch, green_patch], loc = "upper right")
    ax1.set_ylabel("y [m]",fontsize=14)
    ax1.set_xlabel("x [m]",fontsize=14) 
    ax1.set_xlim(-0.06,0.06)
    ax1.set_ylim(-0.06,0.06)

#%%
if PlotResult: # Plot Shell vs. counts, legend (DB, USB, LSB)
    #print("yay")
    DB_shells = shell[idx13_hit]
    LB_shells = shell[idx23_hit]
    UB_shells = shell[idx14_hit]
    BR_shells = shell[idx5_hit]
    
    fig, (ax10) = plt.subplots(nrows=1, ncols=1,figsize=(6,6))  
    fig.tight_layout(pad=4.0)
    ax10.grid()
    ax10.hist(DB_shells,bins=100,histtype='stepfilled', facecolor='none', edgecolor='k',label='Double bounce',linewidth=1)
    ax10.hist(LB_shells,bins=44,histtype='stepfilled', facecolor='none', edgecolor='r',label='Single lower bounce',linewidth=1)
    ax10.hist(UB_shells,bins=45,histtype='stepfilled', facecolor='none', edgecolor='b',label='Single upper bounce',linewidth=1)
    ax10.hist(BR_shells,bins=30,histtype='stepfilled', facecolor='none', edgecolor='g',label='Backside reflection',linewidth=1)

    
    ax10.set_xlabel("Shell",fontsize=14)
    ax10.set_ylabel("Photon counts on detector",fontsize=14) 
    
    ax10.set_xlim(0,101)
    #ax10.set_ylim(0.0,0.15)
    ax10.legend(handlelength=6,loc=('upper left'),fontsize=12)
    

  
if PlotResult: # Plot PSF of detector area
    data13 = {'Type': ["Double bounce"]*len(P3[idx13_hit,0]),
            'X': P3[idx13_hit,0],
            'Y': P3[idx13_hit,1]}
    data14 = {'Type': ["Upper single bounce"]*len(P3[idx14_hit,0]),
            'X': P3[idx14_hit,0],
            'Y': P3[idx14_hit,1]}
    data23 = {'Type': ["Lower single bounce"]*len(P3[idx23_hit,0]),
            'X': P3[idx23_hit,0],
            'Y': P3[idx23_hit,1]}
    df1 = pd.DataFrame(data13)
    df2 = pd.DataFrame(data14)
    df3 = pd.DataFrame(data23)
    
    frames1 = [df1]#, df2, df3]
    frames2 = [df2]
    frames3 = [df3]

    result1 = pd.concat(frames1)
    result2 = pd.concat(frames2)
    result3 = pd.concat(frames3)

    sns.jointplot(data=result1, x="X", y="Y", hue="Type")
    sns.jointplot(data=result2, x="X", y="Y", hue="Type")
    sns.jointplot(data=result3, x="X", y="Y", hue="Type")



#%%
#res = stats.cumfreq(P3[idx13,3][0], numbins=25)
if PlotResult: #  Get HPD of double bounce photons
    n_bins = 40
    fig, ((ax1,ax2),(ax3,ax4)) = plt.subplots(nrows=2, ncols=2,figsize=(4*2,4*2))  
    fig.tight_layout(pad=4.0)
    
    
    
    ax1.plot(P3[idx13_hit,0]*1000,P3[idx13_hit,1]*1000,'k.',label='DB photons')
    ax1.plot(np.median(P3[idx13_hit,0])*1000,np.median(P3[idx13_hit,1])*1000,'r.',label='median')
    ax1.plot(np.mean(P3[idx13_hit,0])*1000,np.mean(P3[idx13_hit,1])*1000,'b.',label='mean')
    ax1.set_xlabel("x [mm]")
    ax1.set_ylabel("y [mm]")
    #ax1.set_xlim(-3.2, 3.2)
    #ax1.set_ylim(-3.2, 3.2)
    ax1.grid()
    ax1.legend(handlelength=6)
    
    DB = [(P3[idx13_hit,0]-np.median(P3[idx13_hit,0])),(P3[idx13_hit,1]-np.median(P3[idx13_hit,1]))]
    ax2.plot(DB[0]*1000,DB[1]*1000,'k.',label='Centered aroud median')
    ax2.grid()
    ax2.set_xlabel("x [mm]")
    ax2.set_ylabel("y [mm]")
    ax2.set_xlim(-3.2, 3.2)
    ax2.set_ylim(-3.2, 3.2)
    ax2.legend(handlelength=6)
    
    DB_r=np.sqrt(DB[0]**2+DB[1]**2)
    ax3.hist(DB_r*1000, bins=n_bins)
    ax3.set_xlabel("Radus [mm]")
    ax3.set_ylabel("Double bounce photons")
    ax3.grid()
    
    
    values, base = np.histogram(DB_r, bins=n_bins*500)
    cumulative = np.cumsum(values)
    #ax4.hist(DB_r, bins=n_bins, density=True, histtype='step', cumulative=True, label='Normalized Cumulative')
    cumulative_norm = cumulative/max(cumulative)
    base_norm = base[:-1]
    ax4.plot(base_norm*1000, cumulative_norm, c='blue')
    v_0p5 = find_closest(cumulative_norm, 0.5) # or should it be 0.25?? becuase the radius is actually d?
    HPD = (math.atan(base_norm[v_0p5]/FL))*(180/np.pi)*3600 # HPD [arcsec]  /60 to get in arcmin
    #HPD = (math.atan(base_norm[v_0p5]/FL)*2)*3437.746771 # HPD [arcmin]
    #HPD = ((base_norm[v_0p5]*2*1000)/20000)*(180/np.pi)*3600/60 # KKM 
    ax4.plot(base_norm[v_0p5]*1000,cumulative_norm[v_0p5],'ro',label = ("HPD = "+ str(round(HPD,3)) + ' arcsec'))
    
    ax4.set_xlabel("Radus [mm]")
    #ax4.set_ylabel("Fraction of on-axis photons") 
    ax4.set_ylabel("Cumulative double bounce photons") 
    ax4.legend(handlelength=6)
    ax4.grid()
    
    '''
    
    ax1.hist(P3[idx13,3][0], bins=n_bins)
    ax1.set_xlabel("Radus [m]")
    ax1.set_ylabel("Double Bounce Photons") 
    values, base = np.histogram(P3[idx13,3][0], bins=n_bins*100)
    cumulative = np.cumsum(values)
    ax2.hist(P3[idx13,3][0], bins=n_bins, density=True, histtype='step', cumulative=True, label='Normalized Cumulative')
    cumulative_norm = cumulative/max(cumulative)
    base_norm = base[:-1]
    ax2.plot(base_norm, cumulative_norm, c='blue')
    v_0p5 = find_closest(cumulative_norm, 0.5) # or should it be 0.25?? becuase the radius is actually d?
    HPD = (math.atan(base_norm[v_0p5]/FL))*3437.746771 # HPD [arcmin]
    #HPD = (math.atan(base_norm[v_0p5]/FL)*2)*3437.746771 # HPD [arcmin]
    #HPD = ((base_norm[v_0p5]*2*1000)/20000)*(180/np.pi)*3600/60 # KKM 
    ax2.plot(base_norm[v_0p5],cumulative_norm[v_0p5],'ro',label = ("HPD = "+ str(round(HPD,2)) + ' arcmin'))
    
    ax2.set_xlabel("Radus [m]")
    ax2.set_ylabel("Fraction of on-axis photons") 
    ax2.legend(handlelength=6)
    '''
#%%
if PlotResult: # Plot 2D Raytrace, of Radius vs FL
    
    fig, (ax10) = plt.subplots(nrows=1, ncols=1,figsize=(6,6))  
    fig.tight_layout(pad=4.0)
    
    if len(P0[idx13,0][0])!=0:
        for i in range(len(P0[idx13,0][0])): # Plot Double Bounce
            X0 = P0[idx13,0][0][i]
            Y0 = P0[idx13,1][0][i]
            Z0 = P0[idx13,2][0][i]
            R0 = P0[idx13,3][0][i]
                
            
            X1 = P1[idx13,0][0][i]
            Y1 = P1[idx13,1][0][i]
            Z1 = P1[idx13,2][0][i]
            R1 = P1[idx13,3][0][i]
            
            X2 = P2[idx13,0][0][i]
            Y2 = P2[idx13,1][0][i]
            Z2 = P2[idx13,2][0][i]
            R2 = P2[idx13,3][0][i]
            
            X3 = P3[idx13,0][0][i]
            Y3 = P3[idx13,1][0][i]
            Z3 = P3[idx13,2][0][i]
            R3 = P3[idx13,3][0][i]
        
            r_DB = [R0,R1,R2,R3]
            z_DB = [-Z0,-Z1,-Z2,-Z3]
        
            #.plot(,DB, 'k-',label="Double Bounce",markersize = 2 )
            ax10.plot(r_DB,z_DB, 'k-',markersize = 2 )
        ax10.plot(r_DB,z_DB, 'k-',label="Double Bounce",markersize = 2 )
     
    for i in range(len(P0[idx23,0][0])): # Plot Single Lower Bounce
        X0 = P0[idx23,0][0][i]
        Y0 = P0[idx23,1][0][i]
        Z0 = P0[idx23,2][0][i]
        R0 = P0[idx23,3][0][i]
            
        
        X1 = P1[idx23,0][0][i]
        Y1 = P1[idx23,1][0][i]
        Z1 = P1[idx23,2][0][i]
        R1 = P1[idx23,3][0][i]
        
        X2 = P2[idx23,0][0][i]
        Y2 = P2[idx23,1][0][i]
        Z2 = P2[idx23,2][0][i]
        R2 = P2[idx23,3][0][i]
        
        X3 = P3[idx23,0][0][i]
        Y3 = P3[idx23,1][0][i]
        Z3 = P3[idx23,2][0][i]
        R3 = P3[idx23,3][0][i]
    
        r_LB = [R0,R2,-R3]
        z_LB = [-Z0,-Z2,-Z3]        
        
        ax10.plot(r_LB,z_LB, 'r--',markersize = 2 )
    ax10.plot(r_LB,z_LB, 'r--',label="Single Lower Bounce",markersize = 2 )
    
    if len(P0[idx14,0][0])!=0:
        for i in range(len(P0[idx14,0][0])): # Plot Single Upper Bounce
            X0 = P0[idx14,0][0][i]
            Y0 = P0[idx14,1][0][i]
            Z0 = P0[idx14,2][0][i]
            R0 = P0[idx14,3][0][i]
                
            
            X1 = P1[idx14,0][0][i]
            Y1 = P1[idx14,1][0][i]
            Z1 = P1[idx14,2][0][i]
            R1 = P1[idx14,3][0][i]
            
            X2 = P2[idx14,0][0][i]
            Y2 = P2[idx14,1][0][i]
            Z2 = P2[idx14,2][0][i]
            R2 = P2[idx14,3][0][i]
            
            X3 = P3[idx14,0][0][i]
            Y3 = P3[idx14,1][0][i]
            Z3 = P3[idx14,2][0][i]
            R3 = P3[idx14,3][0][i]
        
            r_UB = [R0,R1,R3]
            z_UB = [-Z0,-Z1,-Z3]        
            
            ax10.plot(r_UB,z_UB, 'b--',markersize = 2 )
        ax10.plot(r_UB,z_UB, 'b--',label="Single Upper Bounce",markersize = 2 )
    
    if len(P0[idx5,0][0])!=0:
        for i in range(len(P0[idx5,0][0])): # Plot Back reflectio
            X0 = P0[idx5,0][0][i]
            Y0 = P0[idx5,1][0][i]
            Z0 = P0[idx5,2][0][i]
            R0 = P0[idx5,3][0][i]
                
            
            X1 = P1[idx5,0][0][i]
            Y1 = P1[idx5,1][0][i]
            Z1 = P1[idx5,2][0][i]
            R1 = P1[idx5,3][0][i]
            
            X2 = P2[idx5,0][0][i]
            Y2 = P2[idx5,1][0][i]
            Z2 = P2[idx5,2][0][i]
            R2 = P2[idx5,3][0][i]
            
            X3 = P3[idx5,0][0][i]
            Y3 = P3[idx5,1][0][i]
            Z3 = P3[idx5,2][0][i]
            R3 = P3[idx5,3][0][i]
        
            r_BR = [R0,R1,R2,R3]
            z_BR = [-Z0,-Z1,-Z2,-Z3]
        
            #.plot(,DB, 'k-',label="Back reflectio",markersize = 2 )
            ax10.plot(r_BR,z_BR, 'g--',markersize = 2 )
        ax10.plot(r_BR,z_BR, 'g--',label="Back reflection",markersize = 2 )

    if Aperture:
        
        x = [dia_aperture_1/2,dia_aperture_1/2, dia_outer_aperture_2,dia_outer_aperture_2]
        y = [-z_aperture-thickness_aperture,-z_aperture,-z_aperture,-z_aperture-thickness_aperture]
        ax10.add_patch(patches.Polygon(xy=list(zip(x,y)), fill=True,color='m',label="Aperture disk"))
        
        x =[-dia_outer_aperture_2,-dia_outer_aperture_2,-dia_aperture_1/2,-dia_aperture_1/2]
        ax10.add_patch(patches.Polygon(xy=list(zip(x,y)), fill=True,color='m'))

    for i in range(len(r_ut)):
        rut = r_ut[i]
        rub = r_ub[i]
        zut = z_ut[i]
        zub = z_ub[i]
        
        rlt = r_lt[i]
        rlb = r_lb[i]
        zlt = z_lt[i]
        zlb = z_lb[i]
        
        #t_sub = 0.0005 # [m] substrate thickness
        
        xp = [rub,rut,rut+ST,rub+ST] # primary mirros
        yp = [-zub,-zut,-zut,-zub] # primary mirros
        
        xs = [rlb,rlt,rlt+ST,rlb+ST] # secondary mirrors
        ys = [-zlb,-zlt,-zlt,-zlb] # secondary mirrors
        
        ax10.add_patch(patches.Polygon(xy=list(zip(xp,yp)), fill=True,color='grey'))
        ax10.add_patch(patches.Polygon(xy=list(zip(xs,ys)), fill=True,color='grey'))
    ax10.add_patch(patches.Polygon(xy=list(zip(xp,yp)), fill=True,color='grey',label="Mirror"))   
    
    
    
    x_det = [-dia_detector/2,-dia_detector/2, dia_detector/2, dia_detector/2]
    y_det = [-z_FL-thickness_aperture,-z_FL,-z_FL, -z_FL-thickness_aperture]
    ax10.add_patch(patches.Polygon(xy=list(zip(x_det,y_det)), fill=True,color='gray',label="Detector"))   

    ax10.set_xlim(-0.080,0.26)
    ax10.set_ylim(-21.0,0.15)
    # ax10.set_xlim(0.082,0.1)
    # ax10.set_ylim(-0.9,0.1)
    
    ax10.set_xlabel("Radius [m]",fontsize=14)
    ax10.set_ylabel("Focal plane [m]",fontsize=14)
    ax10.legend(handlelength=2,loc=('upper left'),fontsize=8)
    
    
    
    

#%%

# Plot histrogram of reflectionflag
if PlotResult:
    #print("bla")
    binwidth=1
    plt.hist(ReflectionFlag, bins=range(int(min(ReflectionFlag)), int(max(ReflectionFlag)) + binwidth, binwidth))

#%% Plot 3D
if PlotResult: # Plot P0, P1, P2, P3
    fig, (ax1) = plt.subplots(nrows=1, ncols=1,figsize=(6,6))  
    fig.tight_layout(pad=4.0)
    ax1.plot(P0[:,0], P0[:,1], 'k.',label="P0",markersize = 2 )
    ax1.plot(P1[:,0], P1[:,1], 'b.',label="P1",markersize = 2 )
    ax1.plot(P2[:,0], P2[:,1], 'r.',label="P2",markersize = 2 )
    ax1.plot(P3[:,0], P3[:,1], 'm.',label="P3",markersize = 2 )
    ax1.legend(handlelength=6)
    ax1.set_ylabel("y [m]")
    ax1.set_xlabel("x [m]") 
    ax1.set_xlim(-0.3,0.3)
    ax1.set_ylim(-0.3,0.3)

#%%
if PlotResult: # 3D plot
    
    #fig, (ax2) = plt.subplots(nrows=1, ncols=1,figsize=(6,6))  
    #fig.tight_layout(pad=4.0)
    #ax2.tight_layout(pad=4.0) 
    ax = plt.axes(projection='3d')
    ax.plot3D(P0[:,0], P0[:,1], P0[:,2], 'k.')
    ax.plot3D(P1[:,0], P1[:,1], P1[:,2], 'b.')
    ax.plot3D(P2[:,0], P2[:,1], P2[:,2], 'r.')
    ax.plot3D(ACP[:,0], ACP[:,1], ACP[:,2], 'g,')
    ax.plot3D(P3[:,0], P3[:,1], P3[:,2], 'm,')
    
#%%
if PlotResult: # Detector histogram
    fig, (ax1) = plt.subplots(nrows=1, ncols=1,figsize=(6,6))  
    fig.tight_layout(pad=4.0)   
    ax1.hist(P3[:,1],bins=200)
    ax1.hist(P3[:,0],bins=200)
