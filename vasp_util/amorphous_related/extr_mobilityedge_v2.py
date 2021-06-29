import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
####################################################################
#python (*.py) splitdos.dat VB_range CB_range
#xlim, ylim <-- setting in code
####################################################################

#####SETTING#####
vb_ran= float(sys.argv[2]) # 0.45
cb_ran= float(sys.argv[3]) # 0.45
vb_expinf = -0.3
cb_expinf = 0.3
inf_ran = 0.05

xlim  = [-2,2]
ylim  = [0,160]
inipts = 21   # The initial points of x for fitting
#################

class read:
 def __init__(self):
  self.data = np.array([])
  self.ini  = 0
  self.end  = 0
 def get(self, filename):
  with open(filename,'r') as ff:
   tmp = ff.readlines()
   for i in range(len(tmp)):
    if i == 0: pass
    else:
     self.data = np.append(self.data, map(float,tmp[i].split())) 
   self.data = np.transpose(self.data.reshape(len(tmp)-1,3))
  return self.data

 def getpara(self,ini,end,fit):
  #fit == 1: VB fitting
  #fit == 2: CB fitting
  self.ini = int((ini-self.data[0][0])/(self.data[0][1]-self.data[0][0])) 
  self.end = int((end-self.data[0][0])/(self.data[0][1]-self.data[0][0]))
  self.para = np.polyfit(self.data[0][self.ini:self.end+1], self.data[fit][self.ini:self.end+1], 1)
  return self.para

 def poly1d(self,x):
  return self.para[1]+x*self.para[0]

 def R2(self,fit):
  #fit == 1: VB fitting
  #fit == 2: CB fitting
  self.xpts = self.data[0][self.ini:self.end+1]
  self.ypts = self.data[fit][self.ini:self.end+1]
  self.ybar = np.mean(self.ypts)
  self.SStot = np.sum((self.ypts-self.ybar)**2)
  self.SSres = np.sum((self.poly1d(self.xpts)-self.ypts)**2)
  return 1-(self.SSres/self.SStot)
####################################################################
def deriv(x,y):
 drv=np.zeros((2,len(x)))
 h = (x[1]-x[0])
 tot_gau = np.zeros(len(x))
 sigma = 0.1
 if len(x) == len(y):
  for i in range(len(x)):
   tot_gau = tot_gau + y[i]*mlab.normpdf(x,x[i],sigma)
  y = tot_gau/len(x)
  drv[0] = y[:]
  for i in range(len(x)):
   if i==0:
    drv[1][i] = (y[2]-2*y[1]+y[0])/h**2
   elif i==len(x)-1:
    drv[1][i] = drv[1][i-1]
   else:
    drv[1][i] = (y[i+1]-2*y[i]+y[i-1])/h**2
  return drv
 else: return 0

def findroot(xp,yp,a,b):
 def f(pt,xp,yp):
  return np.interp(pt,xp,yp)
 m = (a+b)/2.0
 c = 0
 while(np.absolute(b-a) > 0.001 and c < 1000):
  if   f(a,xp,yp)*f(m,xp,yp) <= 0: b=m
  else: a=m
  m=(a+b)/2.0
  c+=1
 return m
####################################################################
avg = read()
avg.get(sys.argv[1])

step = avg.data[0][1]-avg.data[0][0]
zero_idx = int(np.absolute(avg.data[0][0])/step)
arrini = zero_idx+int(xlim[0]/step)
arrend = zero_idx+int(xlim[1]/step)
idxlist = range(arrini,arrend)

###CB-fit###
cx2 = avg.data[2][idxlist]
yp=deriv(avg.data[0][idxlist],cx2**0.5)
cbR = int(findroot(avg.data[0][idxlist],yp[1],cb_expinf-inf_ran,cb_expinf+inf_ran)/step)
xini = np.linspace(avg.data[0][zero_idx + cbR]-inf_ran,avg.data[0][zero_idx + cbR]+inf_ran,inipts)
gau_cb = yp

r2array = np.zeros(inipts)
for iii in range(len(xini)):
 avg.getpara(xini[iii],xini[iii]+cb_ran,2)
 r2array[iii] = avg.R2(2)
cbpara=[avg.para[0],avg.para[1]]
cini = np.abs(cbpara[1]/cbpara[0])
cbx=np.linspace(cini,cini+2,1001)
cfitini = xini[r2array.argmax()]
print('CB-fit results')
print('Initial x pos :\t%.2f'%cfitini)
print('fitting range :\t%.2f'%cb_ran)
print('R^2 value     :\t%.4f'%r2array[r2array.argmax()])
print('mobility edge :\t%.2f'%cini)


###VB-fit###
vx2 = avg.data[1][idxlist]
yp=deriv(avg.data[0][idxlist],vx2**0.5)
vbR = int(findroot(avg.data[0][idxlist],yp[1],vb_expinf-inf_ran,vb_expinf+inf_ran)/step)
xini = np.linspace(avg.data[0][zero_idx+vbR]-inf_ran,avg.data[0][zero_idx+vbR]+inf_ran,inipts)
gau_vb = yp

for iii in range(len(xini)):
 avg.getpara(xini[iii]-vb_ran,xini[iii],1)
 r2array[iii] = avg.R2(1)
vbpara=[avg.para[0],avg.para[1]]
vini = -np.abs(vbpara[1]/vbpara[0])
vbx=np.linspace(vini,vini-2,1001)
vfitini = xini[r2array.argmax()]
print('VB-fit results')
print('Initial x pos :\t%.2f'%vfitini)
print('fitting range :\t%.2f'%vb_ran)
print('R^2 value     :\t%.4f'%r2array[r2array.argmax()])
print('mobility edge :\t%.2f'%vini)
print('mobility gap  :\t%.2f'%(cini-vini))

###plot###
fig, ax = plt.subplots()
prop_fit={'linestyle':'--','lw':1}
axes=plt.gca()
axes.set_xlim(xlim)
axes.set_ylim(ylim)
vb,  = ax.plot(avg.data[0],(avg.data[1]**2)**0.25,color='red')
vb_fit = ax.plot(vbx, (vbpara[1]+(vbx)*vbpara[0])**0.5,':r',**prop_fit)
cb,  = ax.plot(avg.data[0],avg.data[2]**0.5,color='blue')
cb_fit = ax.plot(cbx, (np.absolute(cbpara[1]+(cbx)*cbpara[0]))**0.5,':b',**prop_fit)

ax.axvspan(vfitini,vfitini-vb_ran, alpha=0.2, color='red')
ax.axvspan(cfitini,cfitini+cb_ran, alpha=0.2, color='blue')
prop_leg={'loc':2,'fontsize':'small','frameon':False}
plt.ylabel("DOS (states/eV)")
plt.xlabel("Fermi energy (eV)")

fig, (ax1,ax2) = plt.subplots(2,1)
v_gau, = ax1.plot(avg.data[0][idxlist],gau_vb[0]*4,'r')
c_gau, = ax1.plot(avg.data[0][idxlist],gau_cb[0]*4,'b')

ax2.plot(avg.data[0][idxlist],gau_vb[1],'r')
ax2.plot(avg.data[0][idxlist],gau_cb[1],'b')
ax2.axvspan(vb_expinf-inf_ran,vb_expinf+inf_ran, alpha=0.2, color='red')
ax2.axvspan(cb_expinf-inf_ran,cb_expinf+inf_ran, alpha=0.2, color='blue')

plt.show()
