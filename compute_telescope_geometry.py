# ===============================================================
# Author: Nis C. Gellert
# DTU Space
# Last updated: September 2022
# The following scripts computes the optics geometry of HEX-P used for raytracing
# Equations based on "Instrument Design Documentr", Jason Koglin, 2010, V5.2
# ===============================================================


import os
os.chdir(r"C:\Users\nige\DTU Space\PhD\Projects\Raytracing")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import math
from datetime import date

# HEX-P telescope design
save_as_txt = [True,"Geometry_HEXP.txt"]
FL = 20 # Focal length, [m]
ML = 0.40 # Mirror length, [m]
N = 98 # Number of shells
rmax = 0.250# Rmax, [m]
ST = 0.0005 # Substrate thickness, [m] 
gap = 0.002 # Optics gap, [m]


f_gap = 1. # Variable gap parameter, [unitless]

r1_ub = 0.06 # Radius to the center of the upper shell, [m]
r1_uo = 0.06 # Radius to the center of the upper shell, [m]

r_uo = np.zeros(N) # Radius to center of the upper shell, Initialiser.
r_ut = np.zeros(N) # Radius to top of the upper shell, Initialiser.
r_ub = np.zeros(N) # Radius to bottom of the upper shell, Initialiser.
r_lt = np.zeros(N) # Radius to top of the lower shell, Initialiser.
r_lb = np.zeros(N) # Radius to bottom of the lower shell, Initialiser.


z_ut = np.zeros(N) # z to top of the upper shell, Initialiser.
z_ub = np.zeros(N) # z to bottom of the upper shell, Initialiser.
z_lt = np.zeros(N) # z to top of the lower shell, Initialiser.
z_lb = np.zeros(N) # z to bottom of the lower shell, Initialiser.



alpha = np.zeros(N) # Grazing incident angle, [rad]. Initialiser


r_uo[0] = r1_uo
alpha[0] = r_uo[0]/(4*FL)
r_ut[0] = r_uo[0] + alpha[0]*ML/2
r_ub[0] = r_ut[0] - alpha[0]*ML
r_lt[0] = r_ut[0] - alpha[0]*(ML+4*gap)
r_lb[0] = r_ut[0] - 4*alpha[0]*(ML+gap)








for i in range(N-1):  
    #theta[i+1] = (r_ut[i]+ST)/(4*FL+ML/2-ML*f_gap)*(180/np.pi)  # Shell incident angle [degree]
    r_ut[i+1] = (r_ut[i]+ST)*(1+((4*FL)/(ML*f_gap)+1/(2*f_gap)-1)**(-1)) # Shell radius [m]
    alpha[i+1] = (r_ut[i]+ST)/(4*FL+ML/2-ML*f_gap)  # Shell incident angle [rad]
    
    r_ub[i+1] = r_ut[i+1] - alpha[i+1] * ML  
    r_lt[i+1] = r_ut[i+1] - alpha[i+1] * (ML + 4*gap)  
    r_lb[i+1] = r_ut[i+1] - 4*alpha[i+1] * (ML + gap)  
    
    z_ut[i+1] = z_ut[i]
    z_ub[i+1] = ML*np.cos(alpha[i+1])
    z_lt[i+1] = z_lt[i]
    z_lb[i+1] = ML*(3*np.cos(alpha[i+1]))



#Txt file with 
#Shell r1 r2 r3 r4 

if save_as_txt[0]:
    fname = save_as_txt[1]
    file1 = open("Telescope geometries/"+fname,"w+") 
    today = date.today()
    file1.write("# Saved on %s\n\n" % today.strftime("%d/%m/%Y"))
    file1.write("# Focal length, [m]: %s \n" % FL)
    file1.write("# Mirror length, [m]: %s \n" % ML)
    file1.write("# Number of shells: %s \n" % N)
    file1.write("# Optics gap, [m]: %s \n" % gap)
    file1.write("# Variable gap parameter: %s \n" % f_gap)
    file1.write("# Substrate thickness, [m]: %s \n" % ST)
    file1.write("# Radius to the center of the upper shell, [m]: %s \n" % r1_uo )
    file1.write("# Shell, r_ut [m], r_ub [m], r_lt [m], r_lb [m], alpha [mrad], z_ut [m], z_ub [m], z_lt [m], z_lb [m] \n" )

    for i in range(int(N)):
         file1.write("%s %.9E %.9E %.9E %.9E %.9E %.9E %.9E %.9E %.9E\n" % (i+1, r_ut[i], r_ub[i], r_lt[i], r_lb[i], alpha[i]/1000, z_ut[i], z_ub[i], z_lt[i], z_lb[i]))
    file1.close()     

'''
# Plot radius of the shells
for i in range(len(r_ut)):
    plt.plot(r_ut[i],1.45,'o')
    plt.plot(r_ub[i],1.25,'o')
    plt.plot(r_lt[i],1.2,'o')
    plt.plot(r_lb[i],1,'o')
'''   
    
#%%
import os
os.chdir(r"C:\Users\nige\DTU Space\PhD\Projects\Raytracing")
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import random
import math
from datetime import date
#calculate_mirror_parameters_will, 'ND15_F20', 20000.,60.,250.,400.,0.5,1.,5.
# HEX-P telescope design

save_as_txt = [True,"Geometry_HEXP_20230902.txt"]
f=20000. # Focal length, [mm]
rmin=60. # Rmin, [mm]
rmax=250.# Rmax, [mm]
lshell=400. # Mirror length, [mm]
thickness=0.5 # Substrate thickness [mm]
shell_spacing=1. # [mm]
gap=5. # [mm] hey wait. or 2.5??


#;; Formula's
#;gamma=1.-shell_spacing*l/4./f
#;rshell[i] = rmax * gamma^i - (1.-gamma^i) / (1.-gamma) * t
#;alpha_i=R_i/4./f
#;; From Jason
#;ruo = 4*alpha*f		; radius center upper
#;rut = ruo+alpha*l/2.		; radius upper top
#;alpha = rut/(4*f+l/2.)
#;rub = rut-alpha*l		; radius upper bottom
#;rlt = rut-alpha*(l+4*g)	; radius lower top
#;rlb = rut-4*alpha*(l+g)	; radius lower bottom
#;; Build up the telescope from the inside out

i = 0
alpha = [0.0]
R1 = [rmin]
fgap = shell_spacing
t = thickness
L = lshell
F = f
g = gap
fov = 5./60.*(np.pi/180)
while R1[i] <= rmax:
    if i == 54:
        tempA = (t+R1[i]+L*fov)/(4.*F-L/2.)
        tempR = tempA*(4.*F+L/2.)
        R1.append(tempR)
        R1.append(tempR+12)
        
        alpha.append(tempA)
        alpha.append(tempA)
        
        #alpha.append(0.0)
        #alpha.append(0.0)
        i+=2	
	#alpha = [alpha,(t+R1[i]+L*fov)/(4.*F-L/2.)]
    alpha.append((t+R1[i]+L*fov)/(4.*F-L/2.)) # [rad]
    #print("what")
    #alpha[i]=(t+R1[i]+L*fov)/(4.*F-L/2.)
	#R1 = [R1,alpha[i+1]*(4.*F+L/2.)]
    R1.append(alpha[i+1]*(4.*F+L/2.)) # [mm]
    i+=1
	

R1 = np.array(R1)
alpha = np.array(alpha)

R2 = R1 - alpha*L
R3 = R1 - alpha*(L+4.*g)
R4 = R1 - 4.*alpha*(L+g)	
R0 = R1 - alpha*L/2.


r_ut = R1
r_ub = R2
r_lt = R3
r_lb = R4
#r1_uo = R0[0]


z_ut = np.zeros(len(R1)) # z to top of the upper shell, Initialiser.
z_ub = np.zeros(len(R1)) # z to bottom of the upper shell, Initialiser.
z_lt = np.zeros(len(R1)) # z to top of the lower shell, Initialiser.
z_lb = np.zeros(len(R1)) # z to bottom of the lower shell, Initialiser.

z_ut[:] = 0
z_ub[:] = lshell*np.cos(alpha)
z_lt[:] = (lshell+gap/2)*np.cos(alpha) + (gap/2)*np.cos(alpha*3)
#z_lt[:] = z_ub + gap + gap
z_lb[:] = (lshell+gap/2)*np.cos(alpha) + (gap/2)*np.cos(alpha*3) + lshell*(np.cos(alpha*3))

alpha[55] = 0
alpha[56] = 0
N = np.count_nonzero(alpha)
alpha = alpha[1:]
r_ut = r_ut[1:]
r_ub = r_ub[1:]
r_lt = r_lt[1:]
r_lb = r_lb[1:]
z_ut = z_ut[1:]
z_ub = z_ub[1:]
z_lt = z_lt[1:]
z_lb = z_lb[1:]
#nshell = n_elements(R0)
#print, [transpose(indgen(nshell)+1),transpose(R0[0:nshell-1]),transpose(alpha[0:nshell-1])]


	#;for i=0, nshell-1 do print, alpha[i]*1000, rut[i], rub[i], rlt[i], rlb[i]
	#openw, 1, 'HEXP_'+label+'_ctrace.dat'
	#printf, 1,';nshells TotalLength Gap Tsub FocalLength ngroups nEnergy nAngle'
	#printf, 1, nshell-1, l, g, t, f, ngrp, n_e, na, format='(i5,f6.1,f6.1,f6.3,f8.1, i4, i5, i5)'
	#for i=0, nshell - 1 do printf, 1, alpha[i]*1000, R1[i], R2[i], R3[i], R4[i]
	#close, 1

	#;; Get mirror group angles

#i=findgen(ngrp+1)	#; mirror group number, 1..n_group
#R_i=exp(alog(Rmin) + i/ngrp * (alog(Rmax)-alog(Rmin)))
#alpha_grp = R_i/(4.*F)

	#if not keyword_set(anglefile) then anglefile = 'HEXP_'+label+'_angle.dat'
	#openw,1,anglefile
	#printf,1,strtrim(string(ngrp),1)
	#printf,1,strtrim(string(transpose(alpha_grp)*1000),1)
	#close,1


if save_as_txt[0]:
    fname = save_as_txt[1]
    file1 = open("Telescope geometries/"+fname,"w+") 
    today = date.today()
    file1.write("# Saved on %s\n\n" % today.strftime("%d/%m/%Y"))
    file1.write("# Focal length, [m]: %s \n" % str(F/1000))
    file1.write("# Mirror length, [m]: %s \n" % str(lshell/1000))
    file1.write("# Number of shells: %s \n" % N)
    file1.write("# Optics gap, [m]: %s \n" % str(gap/1000))
    file1.write("# Shell spacing [m]: %s \n" % str(shell_spacing/1000))
    file1.write("# Substrate thickness, [m]: %s \n" % str(thickness/1000))
    #file1.write("# Radius to the center of the , [m]: %s \n" % r1_uo )
    file1.write("# Shell, r_ut [m], r_ub [m], r_lt [m], r_lb [m], alpha [mrad], z_ut [m], z_ub [m], z_lt [m], z_lb [m] \n" )

    for i in range(len(alpha)):
         file1.write("%s %.9E %.9E %.9E %.9E %.9E %.9E %.9E %.9E %.9E\n" % (i+1, r_ut[i]/1000, r_ub[i]/1000, r_lt[i]/1000, r_lb[i]/1000, alpha[i]*1000, z_ut[i]/1000, z_ub[i]/1000, z_lt[i]/1000, z_lb[i]/1000))
    file1.close()     
    
    
#%%

r_ut = []
r_ub = []
r_lt = []
r_lb = []
alpha =[]

# HEX-P telescope design
save_as_txt = [True,"Geometry_HEXP.txt"]
FL = 20 # Focal length, [m]
ML = 0.40 # Mirror length, [m]

r_ub.append(0.06) # # Radius to the center of the upper shell, [m]
r_ub_max = 0.250
rmax = 0.250# Rmax, [m]
ST = 0.0005 # Substrate thickness, [m] 
gap = 0.002 # Optics gap, [m]
alpha =[]
#f_gap = 1. # Variable gap parameter, [unitless]

#r_uo[0] = r1_uo




r_ut.append(r_ub[0]  + alpha[0]*ML)

#%%
alpha[0] = r_uo[0]/(4*FL)
r_ut[0] = r_uo[0] + alpha[0]*ML/2
r_ub[0] = r_ut[0] - alpha[0]*ML
r_lt[0] = r_ut[0] - alpha[0]*(ML+4*gap)
r_lb[0] = r_ut[0] - 4*alpha[0]*(ML+gap)

while r_ub[i] <= r_ub_max:
    #if i == 54:
    #    tempA = (t+R1[i]+L*fov)/(4.*F-L/2.)
    #    tempR = tempA*(4.*F+L/2.)
    #    R1.append(tempR)
    #    R1.append(tempR+12)
    #    alpha.append(0.0)
    #    alpha.append(0.0)
    #    i+=2	
    
	
    alpha.append((t+R1[i]+L*fov)/(4.*F-L/2.)) # KKM
        r_ut[i+1] = (r_ut[i]+ST)*(1+((4*FL)/(ML*f_gap)+1/(2*f_gap)-1)**(-1)) # Shell radius [m]
    alpha[i+1] = (r_ut[i]+ST)/(4*FL+ML/2-ML*f_gap) 
    
    #print("what")
    #alpha[i]=(t+R1[i]+L*fov)/(4.*F-L/2.)
	#R1 = [R1,alpha[i+1]*(4.*F+L/2.)]
    R1.append