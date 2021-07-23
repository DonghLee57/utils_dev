#python3 code
import sys
import numpy as np
#np.set_printoptions(threshold=sys.maxsize)
import multiprocessing as mp
import random
#####SET params#####
fr = 200
Natom = 5616
cutR = 3.42
nproc = mp.cpu_count()
#nproc = 4
fname = sys.argv[1]
####################

#####Functions #####
def get_fr_dat(fname=fname, frame=fr, N=Natom):
 stl = frame*(N+9)
 pos = []
 lat = np.zeros((3,3))
 with open(fname,'r') as o:
  for i in range(stl): o.readline()
  for i in range(5):   o.readline()
  for i in range(3):
   lat[i][i] = float(o.readline().split()[-1])
  o.readline()
  for n in range(N):
   pos.append(list(map(float,o.readline().split())))
 return np.array(pos), lat

def get_2Ddist(a, b, lat):
 minD = 999
 for i in range(-1,2):
  for j in range(-1,2):
   tmp = b + i*lat[0] + j*lat[1]
   calD = np.linalg.norm(a-tmp)
   if calD < minD:
    minD = calD
 return minD

def get_3Ddist(a, b, lat):
 minD = 999
 for i in range(-1,2):
  for j in range(-1,2):
   for k in range(-1,2):
    tmp = b + i*lat[0] + j*lat[1] + k*lat[2]
    calD = np.linalg.norm(a-tmp)
    if calD < minD:
     minD = calD
 return minD

def MatDist(x):
 nlist = []
 for m in range(len(POS)):
  if m > x:  
   nlist.append(get_2Ddist(POS[x],POS[m],lat))
#   nlist.append(get_3Ddist(POS[x],POS[m],lat))
  else: 
   nlist.append(0)
 return nlist

####################
pos, lat = get_fr_dat()

### Parallellization
# Mask for Ta atoms
mask = []
for n in range(Natom):
 if np.abs(pos[n][1]-1.0) < 0.01:
  mask.append(True)
 else: mask.append(False)

# Make Bond Matrix between Ta atoms 
POS = pos[mask,2:]
pool = mp.Pool(nproc)
neighbors = pool.map(MatDist,range(len(POS)))
pool.close()
pool.join()
neighbors = np.array(neighbors) + np.array(neighbors).T

# Find centers of divid stars
cTa = []
for i in range(len(neighbors)):
 if np.less(neighbors[i][:],cutR).tolist().count(True) == 7 : #Full coordination = 7
  cTa.append(i)

# Remove close centers
rm_list = []
for i in cTa:
 for j in cTa:
  if i != j:
   if neighbors[i][j] <= 3*0.97*cutR:
    rm_list.append(i) 
    break
for i in rm_list:
 cTa.remove(i)

# Clustering as David stars
#CMAP - Hot [0, 0.3-1] scale
CMAP = list(range(len(POS)))
random.shuffle(CMAP)
CMAP = (np.array(CMAP) - min(CMAP))/(max(CMAP)-min(CMAP))*0.7+0.3
CMAP[0] = 0

# Method1 - Clustering the shortest 12 atoms
gr_list = np.zeros(len(POS))
c_nei = neighbors.copy()
for i in cTa:
 c_nei[i].sort()
 for ind in range(13):
  gr_list[np.where(neighbors[i] == c_nei[i][ind])] = i

###test output###
t = np.arange(len(POS))
random.shuffle(t)
o = open('output_%d.lammps' %fr,'w')
o.write("ITEM: TIMESTEP\n")
o.write("%d\n" %fr)
o.write("ITEM: NUMBER OF ATOMS\n")
o.write("%d\n" %Natom)
o.write("ITEM: BOX BOUNDS pp pp pp\n")
for i in range(3):
 o.write("0.0 %f\n" %lat[i][i])
o.write("ITEM: ATOMS id type x y z c\n")
"""
#Color center atoms
c=-1
for n in range(Natom):
 if mask[n]:
  c+=1
  if c in cTa:
   o.write("%d  %d  %f %f %f 2\n" %(int(pos[n][0]),int(pos[n][1]),pos[n][2],pos[n][3],pos[n][4]))
  elif c in rm_list:
   o.write("%d  %d  %f %f %f 1\n" %(int(pos[n][0]),int(pos[n][1]),pos[n][2],pos[n][3],pos[n][4]))
  else:
   o.write("%d  %d  %f %f %f 0\n" %(int(pos[n][0]),int(pos[n][1]),pos[n][2],pos[n][3],pos[n][4]))
 else:
  o.write("%d  %d  %f %f %f 0\n" %(int(pos[n][0]),int(pos[n][1]),pos[n][2],pos[n][3],pos[n][4]))
o.close()
"""
#Color clusters
c=-1
for n in range(Natom):
 if mask[n]:
  c+=1
  o.write("%d  %d  %14.8f %14.8f %14.8f %.2f\n" %(int(pos[n][0]),int(pos[n][1]),pos[n][2],pos[n][3],pos[n][4],CMAP[int(gr_list[c])]))
 else:
  o.write("%d  %d  %f %f %f 0\n" %(int(pos[n][0]),int(pos[n][1]),pos[n][2],pos[n][3],pos[n][4]))
o.close()
