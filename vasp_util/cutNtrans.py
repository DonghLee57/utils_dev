#
# python ~.py POSCAR1 POSCAR2
#
import sys
import numpy as np
###INPUT###
pos = sys.argv[1]
system = ['Ge','Te']
cutR = 7.0 #for test
###########
def getCenter(pos):
 return int(pos.split('_')[-1])-1

def readPOSCAR(file_loc):
 lat = np.array([])
 ntypes = []
 ename = []
 Datoms = np.array([])
 Catoms = np.array([])
 f = open(file_loc,'r')
 tmp = f.readlines()
 lat = np.append(lat,map(float,tmp[2].split()))
 lat = np.append(lat,map(float,tmp[3].split()))
 lat = np.append(lat,map(float,tmp[4].split()))
 lat = np.reshape(lat,(3,3))
 ename = tmp[5].split()
 ntypes = map(int,tmp[6].split())
 if tmp[7][0] == 'S':
  if 'D' == tmp[8][0]:
   for i in range(sum(ntypes)):
    Datoms = np.append(Datoms,map(float,tmp[9+i].split()[:3]))
   Datoms = np.reshape(Datoms,(sum(ntypes),3))
   Catoms = np.matmul(Datoms,lat)
  elif 'C' == tmp[8][0]:
   for i in range(sum(ntypes)):
    Catoms = np.append(Catoms,map(float,tmp[9+i].split()[:3]))
   Catoms = np.reshape(Catoms,(sum(ntypes),3))
   ilat = np.linalg.inv(lat)
   Datoms = np.matmul(Catoms,ilat)
 if tmp[7][0] != 'S':
  if 'D' == tmp[7][0]:
   for i in range(sum(ntypes)):
    Datoms = np.append(Datoms,map(float,tmp[8+i].split()[:3]))
   Datoms = np.reshape(Datoms,(sum(ntypes),3))
   Catoms = np.matmul(Datoms,lat)
  elif 'C' == tmp[7][0]:
   for i in range(sum(ntypes)):
    Catoms = np.append(Catoms,map(float,tmp[8+i].split()[:3]))
   Catoms = np.reshape(Catoms,(sum(ntypes),3))
   ilat = np.linalg.inv(lat)
   Datoms = np.matmul(Catoms,ilat)
 return lat, ntypes, ename, Datoms, Catoms

def getDIS(lat,pos1,pos2):#Position arguments should be Cartensian.
 real_dist = 9999
 for a in range(-1,2):
  for b in range(-1,2):
   for c in range(-1,2):
    trans_pos = pos2 + a*lat[0] + b*lat[1] + c*lat[2]
    dist = sum((pos1-trans_pos)**2)**0.5
    if dist < real_dist:
     real_dist = dist
     close_pos = trans_pos[:]
 return real_dist, close_pos

def neighbor_in_R(lat,positions,center,R,ntypes,ename,transition=np.array([0,0,0])):
 #positions: cartensian coordinates
 #center: Atom number (0,1,2,3,...)
 def find_element(AN,ntypes,ename):
  for n in range(len(ntypes)):
   if AN < sum(ntypes[:n])+ntypes[n]:   
    return ename[n]
 neighbors = np.array([])
 neighbors = np.append(neighbors, (find_element(center,ntypes,ename)))
 neighbors = np.append(neighbors, positions[center]+transition)
 neighbors = np.append(neighbors, '0.0')
 for i in range(len(positions)):
  if i != center:
   r, closeP = getDIS(lat,positions[center],positions[i])
  if r < R:
    neighbors = np.append(neighbors, find_element(i,ntypes,ename))
    neighbors = np.append(neighbors, closeP+transition)
    neighbors = np.append(neighbors, r)
 neighbors = np.reshape(neighbors, (len(neighbors)/5,5))
# test = neighbors[:]
# neighbors = neighbors[neighbors[:,-1].argsort(kind='mergesort')]
 return neighbors #[...,['Element','x','y','z','r'],...]

###
cent = getCenter(pos)

dic_nei = {}
nei_size_list = []
cmat_eig = []

tmp = readPOSCAR(pos) #lat, ntypes, ename, Datoms, Catoms
centP = np.matmul(np.array([0.5,0.5,0.5]),tmp[0])
transV = centP - tmp[-1][cent]

nei = neighbor_in_R(tmp[0],tmp[-1],cent,cutR,tmp[1],tmp[2],transV)
#print nei
ntypes = [0]*len(tmp[2])
for i in range(len(tmp[2])):
 for j in range(len(nei)):
  if tmp[2][i] == nei[j][0]:
   ntypes[i] += 1
###
out = open(sys.argv[1]+'_center','w')
out.write('translated system\n')
out.write(' 1.0\n')
for i in range(3):
 out.write('  %12.16f  %12.16f %12.16f\n'%(tmp[0][i][0],tmp[0][i][1],tmp[0][i][2]))
for i in range(len(tmp[2])):
 if i == len(tmp[2])-1:
  out.write('  %s\n'%(tmp[2][-1]))
 else:
  out.write('  %s'%(tmp[2][i]))
for i in range(len(tmp[2])):
 if i == len(tmp[2])-1:
  out.write('  %d\n'%(ntypes[-1]))
 else:
  out.write('  %d'%(ntypes[i]))
out.write('Selective dynamics\n')
out.write('Cartesian\n')
for i in range(len(nei)):
 if nei[i][0] == nei[0][0]: 
  out.write(' %12.16s  %12.16s %12.16s T T T\n' %(nei[i][1],nei[i][2],nei[i][3]))
for i in range(len(nei)):
 if nei[i][0] != nei[0][0]: 
  out.write(' %12.16s  %12.16s %12.16s T T T\n' %(nei[i][1],nei[i][2],nei[i][3]))
