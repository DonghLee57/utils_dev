import numpy as np
import sys, glob

directory = glob.glob('./*')
directory.sort()
directory.remove('./LOP.py')
directory.remove('./input_R_A')
directory.pop()

def classify(n,tmp):
 #n   - bond number
 #tmp - bonded atoms
 Bool =  False
 for i in range(n):
  if 'C' in tmp[i].split()[0]:
   Bool = True
 if Bool: return True  # 1 - bonded to Dopants
 else:    return False # 2 - not

def dat_len(n,tmp,save):
 #n   - bond number
 #tmp - bonded atoms
 for i in range(n): save = np.append(save,float(tmp[i].split()[-1]))
 return save

def dat_ang(n,tmp,save):
 #n   - bond number
 #tmp - bonded atoms
 for i in range(n*(n-1)/2):
  angle = float(tmp[i].split()[-1])
  if angle <= 140:
   save = np.append(save,angle)
 return save

### 1 - bonded to Dopants 
### 2 - not
Ge1, Ge2 = np.array([]), np.array([])     # Bond number
GeTe1, GeTe2 = np.array([]), np.array([]) # Ge-Te bond length
aGe1, aGe2 = np.array([]), np.array([])   # angles around Ge

#test =  [directory[0]]
#for s in test:
for s in directory:
 tmp = open('%s/bond'%s,'r')
 dat = tmp.readlines()
 for l in range(len(dat)):
  line = dat[l].split()
  if 'bond' in line and 'number' in line and 'is' in line:
   if 'Ge' in dat[l-1].split():
    natoms = int(line[-1])
    if classify(natoms, dat[l+2:l+2+natoms]):
     Ge1   = np.append(Ge1, natoms)
     GeTe1 = dat_len(natoms,dat[l+2:l+2+natoms],GeTe1)
     aGe1  = dat_ang(natoms,dat[l+4+natoms:l+4+natoms*(natoms+1)/2],aGe1)
    else:
     Ge2   = np.append(Ge2, natoms)
     GeTe2 = dat_len(natoms,dat[l+2:l+2+natoms],GeTe2)
     aGe2  = dat_ang(natoms,dat[l+4+natoms:l+4+natoms*(natoms+1)/2],aGe2)
 tmp.close()

ttt = open('GeTe2.dat','w')
for i in range(len(GeTe2)):
 ttt.write('%6.3f\n' %GeTe2[i])
ttt.close()
print("Data for Ge atoms bonded to Carbon")
print("averaged CN     :  %6.2f"%np.mean(Ge1))
print("averaged Ge-Te  :  %6.2f"%np.mean(GeTe1))
print("averaged angles :  %6.2f"%np.mean(aGe1))
print("Data for Ge atoms not bonded to Carbon")
print("averaged CN     :  %6.2f"%np.mean(Ge2))
print("averaged Ge-Te  :  %6.2f"%np.mean(GeTe2))
print("averaged angles :  %6.2f"%np.mean(aGe2))



