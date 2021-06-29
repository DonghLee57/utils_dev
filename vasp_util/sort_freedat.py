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


dat = open(sys.argv[1],'r')

tmp = dat.readlines()
sort = []
for i in range(len(tmp)):
 rmlen =len('outcars_M2/out_snap')
 step = int(tmp[i].split()[0][rmlen:-1])
 sort.append([step,float(tmp[i].split()[-2])])
sort.sort(key=lambda x: x[0])
tr_sort = map(list, zip(*sort))
plt.plot(tr_sort[0],tr_sort[1],'rd')
plt.ylabel("Total energy (eV)",**afont)
plt.xlabel("Ionic step",**afont)
plt.tight_layout()
plt.show()
