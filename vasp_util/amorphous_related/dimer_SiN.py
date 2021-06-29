#Usage:
# python (this_file.py) POSCAR

import sys
import math


#####SETTINGS#####
selectE = 0	 # In POSCAR, 0 -> 1st atom speices
Rc_Si   = float(sys.argv[2])     # Unit : Angstrom
Rc_N    = float(sys.argv[3])	 # Unit : Angstrom
##################

def getDIST(a,b,lat):
 #a-fixed coordinate
 #b-translational coordinate
 real_dis=999
 for i in range(-2,3):
  for j in range(-2,3):
   for k in range(-2,3):
    tX=i*lat[0][0]+j*lat[1][0]+k*lat[2][0]+b[0]
    tY=i*lat[0][1]+j*lat[1][1]+k*lat[2][1]+b[1]
    tZ=i*lat[0][2]+j*lat[1][2]+k*lat[2][2]+b[2]
    transP=[tX, tY, tZ]
    dis=((transP[0]-a[0])**2+(transP[1]-a[1])**2+(transP[2]-a[2])**2)**0.5
    if ((dis < real_dis) and (dis !=0)):
     real_dis=dis
 return real_dis

#####Open POSCAR
p=open(sys.argv[1],'r')
tot_pos=[]
for i in enumerate(open(sys.argv[1],'r')):
 numL = i[0]
for i in range(numL+1):
 tot_pos.append(p.readline().split())
p.close()
if tot_pos[7][0][0] =='S':
 tot_pos.pop(7)

#####Cartesian coordinate converting part
#lattice
lat=[tot_pos[2],tot_pos[3],tot_pos[4]]
for i in range(3):
 lat[i] = map(float, lat[i])
#The number of atom
num_atom_type=[]
for i in range(len(tot_pos[6])):
 num_atom_type.append(tot_pos[6][i])
numTYPE=len(num_atom_type)
num_atom_type=map(int, num_atom_type)
numATOM=sum(num_atom_type)
#remove Selective dynamics T/F
for i in range(numATOM):
 if 'T' in tot_pos[8+i] or 'F' in tot_pos[8+i]:
  tot_pos[8+i].pop(-1)
  tot_pos[8+i].pop(-1)
  tot_pos[8+i].pop(-1)
 tot_pos[8+i] = map(float, tot_pos[8+i])
#convert all coordinates to Cartensian
Atom_CP=[]
if tot_pos[7][0][0] =='C':
 for i in range(numATOM):
  Atom_CP.append([i+1, tot_pos[8+i][0],tot_pos[8+i][1],tot_pos[8+i][2]])
elif tot_pos[7][0][0] =='D':
 for i in range(numATOM):
  Atom_CP.append([i+1,tot_pos[8+i][0]*lat[0][0]+tot_pos[8+i][1]*lat[1][0]+tot_pos[8+i][2]*lat[2][0],tot_pos[8+i][0]*lat[0][1]+tot_pos[8+i][1]*lat[1][1]+tot_pos[8+i][2]*lat[2][1],tot_pos[8+i][0]*lat[0][2]+tot_pos[8+i][1]*lat[1][2]+tot_pos[8+i][2]*lat[2][2]])

#Calculation
pairs = []
atoms = []
for i in range(num_atom_type[selectE]):
 for j in range(num_atom_type[selectE]):
  if j > i:
   if getDIST(Atom_CP[i][1:],Atom_CP[j][1:],lat) < Rc_Si:
    pairs.append((i,j))
    atoms.append(i)
    atoms.append(j)
atoms = list(set(atoms))

print pairs
print atoms

searchspace = []
neighbors = []
for i in range(len(atoms)):
 tmp_nei = []
 for j in range(num_atom_type[0],numATOM):
  if getDIST(Atom_CP[atoms[i]][1:],Atom_CP[j][1:],lat) < Rc_N:
   tmp_nei.append(j)
 neighbors.append(tmp_nei)

dimers=[]
for i in range(len(pairs)):
# print pairs[i][0], atoms.index(pairs[i][0]), neighbors[atoms.index(pairs[i][0])]
 c=0
 for a in neighbors[atoms.index(pairs[i][0])]:
  if a in neighbors[atoms.index(pairs[i][1])]: c += 1
 if c < 1 :
  dimers.append(pairs[i])

print dimers
# atoms.index(pairs[i][1])
# print atoms[i], neighbors[i]
