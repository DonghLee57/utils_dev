import sys
import glob

outcars = glob.glob('./out*')
outcars.sort()
P = []
for i in range(len(outcars)):
 OUT = open(outcars[i],'r')
 out = OUT.readlines()
 for n in range(len(out)):
  if 'NIONS'      in out[n].split(): Natom = int(out[n].split()[-1])
  elif 'pressure' in out[n].split(): press = float(out[n].split()[-7])

 P.append([press, 0])

for i in range(len(P)-1): P[i][1] = P[i][0] - P[-1][0]

res = open('press.dat','w')
for i in range(len(P)):
 res.write('%6s %10.4f %10.4f\n' %(outcars[i].split('_')[-1], P[i][0], P[i][1]))
