import sys
import glob

outcars = glob.glob('./out*')
outcars.sort()
EE=[]
for i in range(len(outcars)):
 OUT = open(outcars[i],'r')
 out = OUT.readlines()
 for n in range(len(out)):
  if 'NIONS'   in out[n].split(): Natom = int(out[n].split()[-1])
  elif 'TOTEN' in out[n].split(): Energy= float(out[n].split()[-2])

 EE.append([Energy, Energy/Natom,0])

for i in range(len(EE)-1): EE[i][2] = EE[i][1] - EE[-1][1]

res = open('energy.dat','w')
for i in range(len(EE)):
 res.write('%6s %14.6f %8.4f %8.4f\n' %(outcars[i].split('_')[-1], EE[i][0], EE[i][1], EE[i][2]))
