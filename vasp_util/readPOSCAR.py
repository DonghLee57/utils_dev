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
 if tmp[7][0] == 'S': tmp.pop(7)
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

def getANG(lat,pos1,pos2,pos3):#Position arguments should be Cartensian.
 #pos1 : Center atom
 #pos2,3 : Neighboring atoms
 a = getDIS(lat,pos2,pos3)
 b = getDIS(lat,pos1,pos2)
 c = getDIS(lat,pos1,pos3)
 #Unit of output : degree
 return np.arccos((b**2+c**2-a**2)/(2*b*c))*180/np.pi


"""
lat, ntypes, Datoms, Catoms = readPOSCAR(sys.argv[1])
natoms = sum(ntypes)
for i in range(natoms):
 for j in range(natoms):
  if j > i:
   dist = getDIS(lat,Catoms[i],Catoms[j])
"""
