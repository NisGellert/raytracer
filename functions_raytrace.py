import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import copy
import math
import bisect
import linecache
import matplotlib.patches as mpatches
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from scipy import stats
from scipy.interpolate import interp1d

'''
# ===============================================================
# Author: Nis C. Gellert
# DTU Space
# Last updated: February 2023
# The following script has the functions used in main_rautrace.py
# ===============================================================
# ----- Read before use: ----
# Generate the Optics Geometry file using "compute_telescope_geometry.py"

# ===============================================================
'''
 


def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

def movephoton(P, z):
    # P is the positionof the photon
    # z is the new absolute z position you want the photon to move to
    newP = copy.deepcopy(P) # New position
    newP[0] = P[0] + (z - P[2])*P[4]/P[6] # x
    newP[1] = P[1] + (z - P[2])*P[5]/P[6] # y
    newP[2] = z # z
    newP[3] = np.sqrt(newP[0]**2+newP[1]**2) # r 
    return newP
  

def intersect(P,shell,alpha,r_ut,ML):
    '''
    This is a rough numerical estimation of the intersection point
    '''
    # Calculate 100 intersections point. Intersection point is where the radius of the photon minus the radius of the plate at point-z is 0
    z=np.linspace(P[2],P[2]+ML,500)
    #z=np.linspace(0,ML,500)
    RX =[]; RY =[]; RZ =[]; Radii =[]; 
    for i in range(len(z)):
        #result.append(movephoton(P,z[i]))
        RX.append(movephoton(P,z[i])[0])
        RY.append(movephoton(P,z[i])[1])
        RZ.append(movephoton(P,z[i])[2])
        Radii.append(movephoton(P,z[i])[3])
    r_shell = r_ut - np.tan(alpha) * z 
    difference_array = np.absolute(r_shell-Radii)
    index = difference_array.argmin()
    if index == 0: # Missed first reflection
        ReflectionFlag = '2'
    elif index == len(z)-1: # Missed first reflection
        ReflectionFlag = '2'
    else: # Compledted first reflection
        ReflectionFlag = '1'
    return z[index], ReflectionFlag
    
def reflection(P,alpha,scatter):
    # Rotate P 
    newP = copy.deepcopy(P) # New position
    #photon - photon origin  and direction\n",
    #alpha - slope of mirror\n",
    xr=newP[0]
    yr=newP[1]
    zr=newP[2]
    dxk=newP[4]
    dyk=newP[5]
    dzk=newP[6]
    v_body = np.zeros(3)
    v_world = np.zeros(3)
    r0 = np.sqrt((xr**2) + (yr)**2)
    n3 = -np.sin(alpha)
    n_perp = np.sqrt(1. - (n3**2))
    n1 = -n_perp * xr/r0
    n2 = -n_perp * yr/r0
    dtheta = 0.0
    dz = alpha
    v_body[0]=np.cos(dz)*np.cos(dtheta)
    v_body[1]=np.cos(dz)*np.sin(dtheta)
    v_body[2]=-np.sin(dz)
    #calculate tranformation between World and Body\n",
    azimuth = np.arctan2(n2,n1)
    #transform into new normal vector\n",
    v_world[0] = np.cos(azimuth)*v_body[0] - np.sin(azimuth)*v_body[1]
    v_world[1] = np.sin(azimuth)*v_body[0] + np.cos(azimuth)*v_body[1] 
    v_world[2] = v_body[2]
    n1 = v_world[0]
    n2 = v_world[1]
    n3 = v_world[2]
    
    if scatter != 0:
        # carried over from ctrace, have not modified it yet
        #;;PERTURB NORMAL VECTOR */
        #;;first pick on a sphere, then squash point onto the plane perp. to normal */
        r_sphere = np.zeros(3)
        r_pancake = np.zeros(3)
        theta = np.arcsin(2*np.random.rand()-1)
        phi = 2*np.pi*np.random.rand()
        r_sphere[0] = np.sin(theta)*np.cos(phi)
        r_sphere[1] = np.sin(theta)*np.sin(phi)
        r_sphere[2] = np.cos(theta)
        dot = r_sphere[0]*n1 + r_sphere[1]*n2 + r_sphere[2]*n3
        r_pancake[0] = r_sphere[0] - dot*n1
        r_pancake[1] = r_sphere[1] - dot*n2
        r_pancake[2] = r_sphere[2] - dot*n3

        l = np.sqrt((r_pancake[0])**2 + (r_pancake[1])**2 + (r_pancake[2])**2) # idl
        r_pancake[0] = r_pancake[0]/l
        r_pancake[1] = r_pancake[1]/l
        r_pancake[2] = r_pancake[2]/l
        
        le = scatter*np.random.normal(loc=0.0, scale=1.0, size=None)
        r_pancake[0] *= le
        r_pancake[1] *= le
        r_pancake[2] *= le

        leng = np.sqrt(1 + (le**2))
        n1 = (n1 + r_pancake[0])/leng
        n2 = (n2 + r_pancake[1])/leng
        n3 = (n3 + r_pancake[2])/leng
    
    k_dot_n = n1*dxk + n2*dyk +n3*dzk
    dxk = dxk - 2*k_dot_n*n1
    dyk = dyk - 2*k_dot_n*n2
    dzk = dzk - 2*k_dot_n*n3
    #print(\"dx,dy,dz\",dxk,dyk,dzk)\n",
    return np.absolute(np.pi/2. - np.arccos(k_dot_n)),dxk,dyk,dzk
    
def cone_reflection_point(r,apex,alpha):
    #r = [xc,yc,zc]\n",
    #k = [dxc,dyc,dzc]\n",
    xc=r[0]
    yc=r[1]
    zc=r[2]
    dxc=r[4]
    dyc=r[5]
    dzc=r[6]
    A = ((dxc**2) + (dyc**2))/(dzc**2) -(np.tan(alpha)**2)
    B = 2. * ( xc*dxc + yc*dyc ) / dzc  - 2. * ( (dxc**2) + (dyc**2) ) * zc/ (dzc**2)  + 2. * apex * (np.tan(alpha)**2)
    C = (xc**2) + (yc**2) - (apex * np.tan(alpha))**2 - 2. * (xc*dxc + yc*dyc) * zc/dzc + ( (dxc**2) + (dyc**2) ) * (zc)**2 /(dzc)**2 
    #Z = ( -B + np.sqrt(np.absolute((B**2)-4.*A*C )) ) / ( 2.*A )# There can be two solutions to Z. Take minimum
    #Z2 = ( -B - np.sqrt(np.absolute((B**2)-4.*A*C )) ) / ( 2.*A )# Solution for anden grads løsning. There can be two solutions to Z. Take minimum

    Z1 = ( -B + np.sqrt(np.absolute((B**2)-4.*A*C )) ) / ( 2.*A )# There can be two solutions to Z. Take minimum
    Z2 = ( -B - np.sqrt(np.absolute((B**2)-4.*A*C )) ) / ( 2.*A )# Solution for anden grads løsning. There can be two solutions to Z. Take minimum
    
    Z = np.array([Z1,Z2])  
    if (Z1 < 0) and (Z2 < 0):
        result = 100
    else:  
        result = min(Z[Z >= 0])
    
    return result       
    #return Z