# Use: python ~.py sqrt_aE.dat 3.5 6.0
import matplotlib.pyplot as plt
import numpy as np
import sys
import matplotlib as mpl
mpl.rc('font',family = 'Arial')


fff = open(sys.argv[1],'r')
xini = float(sys.argv[2])
xend = float(sys.argv[3])

tmp = fff.readlines()

data=np.array([])

for i in range(len(tmp)):
 data = np.append(data, map(float, tmp[i].split()))
data =np.transpose(data.reshape(len(tmp),2))

for i in range(len(tmp)):
 if xini <= data[0][i]: 
  idx_ini = i
  break
idx_end = i + int((xend-xini)/(data[0][1]-data[0][0]))

para = np.polyfit(data[0][idx_ini:idx_end+1], data[1][idx_ini:idx_end+1], 1)

print para
print -(para[1]/para[0])
#fig, ax = plt.subplots()
#ax.plot(data[0],data[1])


afont15 = {'fontname':'Arial', 'fontsize' : 15}

axes=plt.gca()
axes.set_xlim([0,xend])
axes.set_ylim([0,data[1][idx_end]])
plt.xticks(**afont15)
plt.yticks(**afont15)
plt.xlabel(r'$h \nu $ (eV)', **afont15)
plt.ylabel(r'$(\alpha h \nu)^{1/2}$ $(cm^{-1}eV)^{1/2}$', **afont15)
plt.plot(data[0],data[1])
plt.plot(data[0], para[1]+data[0]*para[0])
plt.tight_layout()
plt.show()



