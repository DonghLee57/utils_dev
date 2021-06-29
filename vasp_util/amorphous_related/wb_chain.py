import sys
import numpy as np

def readPOSCAR(file_loc):
 lat = np.array([])
 ntypes = []
 Datoms = np.array([])
 Catoms = np.array([])
 f = open(file_loc,'r')
 tmp = f.readlines()
 lat = np.append(lat,map(float,tmp[2].split()))
 lat = np.append(lat,map(float,tmp[3].split()))
 lat = np.append(lat,map(float,tmp[4].split()))
 lat = np.reshape(lat,(3,3))
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
 else:
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
 return lat, ntypes, Datoms, Catoms

def getDIS(lat,pos1,pos2):#Position arguments should be Cartensian.
 real_dist = 9999
 for a in range(-1,2):
  for b in range(-1,2):
   for c in range(-1,2):
    trans_pos = pos2 + a*lat[0] + b*lat[1] + c*lat[2]
    dist = sum((pos1-trans_pos)**2)**0.5
    if dist < real_dist:
     real_dist = dist 
 return real_dist

lat, ntypes, Datoms, Catoms = readPOSCAR(sys.argv[1])
natoms = sum(ntypes)

###Generate lists of wrong bond pairs and atoms
# wbp : list of wrong bond pairs
# wba : list of wrong bond atoms
wbp = []
wba = []
AtomL = range(48) #Select atoms
Rc = 3.2
for i in AtomL:
 for j in AtomL:
  if j > i:
   if getDIS(lat,Catoms[i],Catoms[j]) < Rc:
    wbp.append([i,j])
    wba.append(i)
    wba.append(j)
wba = list(set(wba))

###Check wrong bond chain
wbpernode={}
for i in range(len(wba)):
 idx_list=[]
 cluster =[]
 for idx in range(len(wbp)):
  if wba[i] in wbp[idx]:
   idx_list.append(idx)
 for j in reversed(idx_list):
  cluster += wbp[j]
 if cluster != []:
  final_wbp = list(set(cluster))
  wbpernode[wba[i]] = final_wbp[:]
  wbpernode[wba[i]].remove(wba[i])
#  for k in reversed(wbpernode[wba[i]]):
#   if k < wba[i]:
#    wbpernode[wba[i]].remove(k)
#  print wba[i], final_wbp, wbpernode[wba[i]]
#print wba

#Get list of atoms in a wrong bond chain
def c_size(node, trace=[]):
#global var: wbpernode
 res = []
 tmp = []
 trace.append(node)
 for i in wbpernode[node]:
  if i not in trace:
   tmp.append(i)
 if len(tmp) == 0: pass
 else:
  for i in tmp:
   res.append(c_size(i,trace))
  return trace 

lll = []
for i in wba:
 tmp = c_size(i,[])
 tmp.sort()
 lll.append(tmp)
lll = list(set([tuple(set(lll)) for lll in lll]))

length = []
for i in lll:
 length.append(len(i))
length.sort(reverse=True)

out = sys.argv[2] + '   '
for i in length:
 out += str(i) + '   '
print out
