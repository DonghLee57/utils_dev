# 
# Calculate the local bond order parameter q_l(i)
# Ref.
#   https://lammps.sandia.gov/doc/compute_orientorder_atom.html
#   J. Chem. Phys. 129, 114707 (2008). eq (3)
#   J. Chem. Phys. 138, 044501 (2013). eq (1) & Table I.
#
# Python ~.py POSCAR.vasp
# output: LAMMPS dump style
#
import sys
import numpy as np
from scipy.special import sph_harm

###INPUT###
pos  = sys.argv[3]
cutR = float(sys.argv[1])
l = int(sys.argv[2])                 #Free parameter in spherical harmonics
###########
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

def neighbor_in_R(lat,positions,center,R,ntypes,ename):
 #positions: cartensian coordinates
 #center: Atom number (0,1,2,3,...)
 #neighbors[0]: itself position
 def find_element(AN,ntypes,ename):
  for n in range(len(ntypes)):
   if AN < sum(ntypes[:n])+ntypes[n]:
    return ename[n]
 neighbors = np.array([])
 neighbors = np.append(neighbors, (find_element(center,ntypes,ename)))
 neighbors = np.append(neighbors, center)
 neighbors = np.append(neighbors, positions[center])
 neighbors = np.append(neighbors, '0.0')
 for i in range(len(positions)):
  if i != center:
   r, closeP = getDIS(lat,positions[center],positions[i])
   if r < R:
    neighbors = np.append(neighbors, find_element(i,ntypes,ename))
    neighbors = np.append(neighbors, i)
    neighbors = np.append(neighbors, closeP)
    neighbors = np.append(neighbors, r)
 neighbors = np.reshape(neighbors, (len(neighbors)/6,6))
 return neighbors #[...,['Element','index','x','y','z','r'],...]

def getAngle(v1,v2): #radian
 return  np.arccos(np.dot(v1,v2)/np.linalg.norm(v1)/np.linalg.norm(v2))

def q(m,nei_atoms,l=l):
 # Should be -l <= m <= l
 # For the given vector, the polar angles (theta, phi) is defined in polar coordinates.
 x = np.array([1.0,0])
 z = np.array([0,0,1.0])
 q = 0
  for i in range(1,len(nei_atoms)):
  vec = nei_atoms[i][2:5].astype(float) - nei_atoms[0][2:5].astype(float)
  phi = getAngle(vec,z)
  #When phi = 0 or pi, theta is undefined...any arbitrary value is assigned.
  if phi < 10**(-10) or phi - np.pi < 10**(-10): theta = 0.0
  else: theta = getAngle(vec[:2], x)
  q += sph_harm(m,l,theta,phi)
 q /= float(len(nei_atoms))
 return q

def avgq(q, nei_atoms):
 qbar = 0
 for i in range(len(nei_atoms)):
  qbar += q[int(nei_atoms[i][1])]
 qbar /= float(len(nei_atoms))
 return qbar

def avgQ(center, qbar, nei_atoms,l=l):
 #avgQ4 of crystalline atoms is greater than 0.9 with the cutoff distance of 3.6
 #k runs over the neighboring atoms including the i-th atom itself
 Q = 0

 for i in range(len(nei_atoms)):
  if int(nei_atoms[i][1]) != center:
   S = 0
   term_i = np.array([])
   term_j = np.array([])
   for m in range(-l,l+1):
    S += np.dot(qbar[m][center],np.conj(qbar[m][i]))
    term_i = np.append(term_i, np.absolute(qbar[m][center]))
    term_j = np.append(term_j, np.absolute(qbar[m][i]))
   a = (np.linalg.norm(term_i)*np.linalg.norm(term_j))
   Q += S/a
 Q /= float(len(nei_atoms))-1
 return Q
 
 
ts = 0
tmp = readPOSCAR(pos) #lat, ntypes, ename, Datoms, Catoms
NA = sum(tmp[1])
dic_nei = {}
dic_q = {}
dic_qbar = {}
for m in range(-l,l+1):
 dic_q[m] = {}
 dic_qbar[m] = {}

dic_Qbar = {}
for i in range(len(tmp[3])):
 dic_nei[i] = neighbor_in_R(tmp[0],tmp[4],i,cutR,tmp[1],tmp[2])
 for m in range(-l,l+1):
  dic_q[m][i]= q(m,dic_nei[i],l)

for i in range(len(tmp[3])):
 for m in range(-l,l+1):
  dic_qbar[m][i] = avgq(dic_q[m],dic_nei[i])

for i in range(len(tmp[3])):
 dic_Qbar[i] = avgQ(i, dic_qbar,dic_nei[i])
path =sys.argv[3]
out_name = path.split('/')[0]+'/Q'+str(l)+'/'+path.split('/')[1]+'.lammpstrj'
print out_name
TRJ = open(out_name,'w')
TRJ = open('Cry_ref.lammpstrj','w')
TRJ.write('ITEM: TIMESTEP\n')
TRJ.write('%d\n' %ts)
TRJ.write('ITEM: NUMBER OF ATOMS\n')
TRJ.write('%d\n' %NA)
TRJ.write('ITEM: BOX BOUNDS xy xz yz pp pp pp\n')
#change lat
for a in range(3):
 TRJ.write('0.0000000000000000e+00 %f 0.0000000000000000e+00\n' %tmp[0][0][0])
TRJ.write('ITEM: ATOMS id type x y z Qreal Qimage\n')
for n in range(NA):
 if n < NA/2: e =1
 else: e = 2
 trATOM = n
 TRJ.write('%d   %d   '%(trATOM+1,e))
 TRJ.write('%f   %f   %f' %(tmp[4][trATOM][0], tmp[4][trATOM][1], tmp[4][trATOM][2]))
 TRJ.write('   %.4f  %.4fi' %(dic_Qbar[n].real,dic_Qbar[n].imag))
 if dic_Qbar[n] > 0.9: #0.9 is the value in REF.
  TRJ.write('   1.0  0.0  0.0')
 else:
  TRJ.write('   1.0  1.0  1.0')
 if n != NA-1:
  TRJ.write('\n')
