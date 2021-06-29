#Usage:
# python (this_file.py) 
#
import math
import numpy as np
import os, sys, glob

#####SETTINGS#####
PATH     = './outcars/*OUTCAR*'
c_length = 401
Energy_ini =  0.0
Energy_fin = 18.0
##################

def rrOUT(path):
 out = open(path,'r')
 tmp =  out.readlines()
 Re, Im = np.array([]), np.array([])
 for i in range(len(tmp)):
  if "REAL" in tmp[i].split():
   for n in range(c_length):
    Re = np.append(Re, map(float,tmp[i+3+n].split()))
  elif "IMAGINARY" in tmp[i].split(): 
   for n in range(c_length):
    Im = np.append(Im, map(float,tmp[i+3+n].split()))
 Re = np.transpose(Re.reshape((c_length,7)))
 Im = np.transpose(Im.reshape((c_length,7)))
 Re_avg = np.array([Re[0],(Re[1]+Re[2]+Re[3])/3])
 Im_avg = np.array([Im[0],(Im[1]+Im[2]+Im[3])/3])
 out.close()
 return [Re_avg, Im_avg]

def calcabscoeff(E, eps1, eps2):
 c    =   3.0*10**(8)
 hbar = 1.054*10**(-34)
 qe   = 1.602*10**(-19)
 const= (2**0.5)/(hbar/qe*c)/100
 return const*E*((eps1**2+eps2**2)**0.5 - eps1)**0.5

#-------------------------------------------------#
outcars = glob.glob(PATH)
outcars.sort()
Energy  = np.linspace(Energy_ini, Energy_fin, Energy_fin/0.01)
res     = np.array([])
for n in range(len(outcars)):
 eps1   = np.array([])
 eps2   = np.array([])
 tmpOUT = rrOUT(outcars[n])
 eps1   = np.append(eps1, np.interp(Energy, tmpOUT[0][0],tmpOUT[0][1]))
 eps2   = np.append(eps2, np.interp(Energy, tmpOUT[1][0],tmpOUT[1][1]))
 res    = np.append(res , calcabscoeff(Energy, eps1, eps2))

res     = res.reshape((len(outcars),len(Energy)))
res_avg = np.mean(np.transpose(res),axis=1)

output1 = open('ABScoff_test.dat','w')
output2 = open('sqrt_aE_test.dat','w')
for i in range(len(Energy)):
 output1.write('%6.2f\t%6.2f\n' %(Energy[i],res_avg[i]))
 output2.write('%6.2f\t%6.2f\n' %(Energy[i],(res_avg[i]*Energy[i])**0.5))
output1.close()
output2.close()
