import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
import matplotlib as mpl
mpl.rc('font',family = 'Arial')
mpl.rcParams['axes.linewidth']=2
mpl.rcParams['xtick.major.width']=2
mpl.rcParams['xtick.major.size']=5
mpl.rcParams['ytick.major.width']=2
mpl.rcParams['ytick.major.size']=5
fsize = 20
afont = {'fontname':'Arial', 'fontsize':fsize}

###
def readout(a):
 for l in range(len(a)):
  line = a[l].split()
  if 'free' in line:
   E  = float(line[-2])
  if 'volume' in line:
   cellV = float(line[-1])
 return E, cellV

def readlog(a):
 for l in range(len(a)):
  if 'TotEng' in a[l].split():
   lastE = float(a[l+1].split()[-2])
 return lastE

###

isoEDFT=np.array([])
isoEDFTf=np.array([])
isoENNP=np.array([])
isoENNPf=np.array([])
cellV = np.array([])
cellVf = np.array([])
for i in range(5,15):
 num = i+1
 tmp  = open('/data/2move/Calc/GST/NN_GeTe/c-GeTe/hexagonal/sample_%d/OUTCAR'%num,'r')
 tmpf  = open('/data/2move/Calc/GST/NN_GeTe/c-GeTe/FCC/sample_%d/OUTCAR'%num,'r')
 a  = tmp.readlines()
 e, v = readout(a)
 isoEDFT = np.append(isoEDFT, e)
 cellV = np.append(cellV, v)
 af  = tmpf.readlines()
 e, v = readout(af)
 isoEDFT = np.append(isoEDFT, e)
 cellVf = np.append(cellVf, v)

 tmp2 = open('./calc_lmp/log_hex_%d'%num,'r')
 tmp2f = open('./calc_lmp/log_fcc_%d'%num,'r')
 a2 = tmp2.readlines()
 a2f = tmp2f.readlines()
 isoENNP = np.append(isoENNP, readlog(a2))
 isoENNPf = np.append(isoENNPf, readlog(a2f))

###
###Bulk modulus
from scipy.optimize import curve_fit
from scipy.optimize import leastsq
def BMeos(params, V):
 k = 2.0/3.0
 E0, V0, B0, Bp0 = params
 return E0+9.0*V0*B0/16.0*(Bp0*((V0/V)**k-1)**3+(6-4*(V0/V)**k)*((V0/V)**k-1)**2)

def obj(params, E, V):
 return E - BMeos(params,V)

p0 = [-24, 170, 40, 2]
#plsq  = leastsq(obj, p0, args=(isoEDFT, cellV[0]))
plsq  = leastsq(obj, p0, args=(isoENNP, cellV))
V_fit = np.linspace(cellV[0], cellV[-1], 100)
E = BMeos(plsq[0], V_fit)
fig, ax = plt.subplots(num=None, figsize=(8,6), dpi=80, facecolor='w',edgecolor='k')
plt.plot(V_fit, E,'r')
pisoENNP,  = ax.plot(cellV,isoENNP,'ko')
prop_leg={'loc':1,'frameon':False,'fontsize':fsize}
ax.tick_params(direction='in',labelsize=fsize)
plt.ylabel("Total energy (eV)",**afont)
plt.xlabel("Volume (angstrom^3)",**afont)
plt.tight_layout()

print('Bulk modulus (energy): %.1f (%.4f)'%(plsq[0][2]*160.21766208, plsq[0][0]))

p0 = [-32, 180, 40, 2]
plsqf  = leastsq(obj, p0, args=(isoENNPf, cellVf))
Vf_fit = np.linspace(cellVf[0], cellVf[-1], 100)
E = BMeos(plsqf[0], Vf_fit)
fig, ax = plt.subplots(num=None, figsize=(8,6), dpi=80, facecolor='w',edgecolor='k')
plt.plot(Vf_fit, E,'r')
pisoENNPf,  = ax.plot(cellVf,isoENNPf,'ko')
prop_leg={'loc':1,'frameon':False,'fontsize':fsize}
ax.tick_params(direction='in',labelsize=fsize)
plt.ylabel("Total energy (eV)",**afont)
plt.xlabel("Volume (angstrom^3)",**afont)
plt.tight_layout()

print('Bulk modulus (energy): %.1f (%.4f)'%(plsqf[0][2]*160.21766208, plsqf[0][0]))
print('Ediff(HEX-FCC) (eV/atom) : %.3f' %( plsq[0][0]/6 - plsqf[0][0]/8))

plt.show()
