#Usage:
# python (this_file.py) (Path to directory including OUTCARs)
#
import math
import numpy as np
import os, sys, glob
from random import shuffle

#####SETTINGS#####
PATH     = sys.argv[1]
c_length = 401
Energy_ini =  0.0
Energy_fin = 18.0
Estep    = 0.01
##################

def rrOUT(path,tag):
 out = open(path,'r')
 tmp =  out.readlines()
 Re, Im = np.array([]), np.array([])
 c=0
 for i in range(len(tmp)):
  if "REAL" in tmp[i].split():
   if c==tag:
    for n in range(c_length):
     Re = np.append(Re, map(float,tmp[i+3+n].split()))
  elif "IMAGINARY" in tmp[i].split(): 
   c+=1
   if c==tag:
    for n in range(c_length):
     Im = np.append(Im, map(float,tmp[i+3+n].split()))
 Re = np.transpose(Re.reshape((c_length,7)))
 Im = np.transpose(Im.reshape((c_length,7)))
 Re_avg = np.array([Re[0],(Re[1]+Re[2]+Re[3])/3])
 Im_avg = np.array([Im[0],(Im[1]+Im[2]+Im[3])/3])
 out.close()
 return [Re_avg, Im_avg]

#-------------------------------------------------#
outcars = glob.glob(PATH+'/*')
outcars.sort()
Energy  = np.linspace(Energy_ini, Energy_fin, Energy_fin/Estep)
eps1    = np.array([])
eps2    = np.array([])
for n in range(len(outcars)):
 tmpOUT = rrOUT(outcars[n],1)
 eps1   = np.append(eps1, np.interp(Energy, tmpOUT[0][0],tmpOUT[0][1]))
 eps2   = np.append(eps2, np.interp(Energy, tmpOUT[1][0],tmpOUT[1][1]))
eps1 = eps1.reshape((len(outcars),len(Energy)))
eps2 = eps2.reshape((len(outcars),len(Energy)))
eps1_avg = np.mean(np.transpose(eps1),axis=1)
eps2_avg = np.mean(np.transpose(eps2),axis=1)

output1 = open('ABScoff.dat','w')
output2 = open('sqrt_aE.dat','w')
c    =   3.0*10**(8)
hbar = 1.054*10**(-34)
qe   = 1.602*10**(-19)
const= (2**0.5)/(hbar/qe*c)/100
abscoff = const*Energy*((eps1_avg**2+eps2_avg**2)**0.5 - eps1_avg)**0.5
for i in range(len(Energy)):
 output1.write('%6.2f\t%6.2f\n' %(Energy[i],abscoff[i]))
 output2.write('%6.2f\t%6.2f\n' %(Energy[i],(abscoff[i]*Energy[i])**0.5))
output1.close()
output2.close()
