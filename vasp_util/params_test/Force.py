import sys
import math

def diff(f1,f2):
 mf1 = (f1[0]**2+f1[1]**2+f1[2]**2)**0.5
 mf2 = (f2[0]**2+f2[1]**2+f2[2]**2)**0.5
 dot = f1[0]*f2[0]+f1[1]*f2[1]+f1[2]*f2[2]
 angle = math.acos(round(dot/mf1/mf2,4))/math.pi*180
 print('%8.4f  %8.4f  %8.4f  %8.2f'%(round(mf1,4),round(mf2,4),round(mf2-mf1,4),round(angle,2)))

OUT1= open(sys.argv[1],'r')
OUT2= open(sys.argv[2],'r')
out1= OUT1.readlines()
out2= OUT2.readlines()

for i in range(len(out1)):
 if 'TOTAL-FORCE' in out1[i].split(): l1 = i
 elif 'NIONS'     in out1[i].split(): Natom = int(out1[i].split()[-1])

for i in range(len(out2)):
 if 'TOTAL-FORCE' in out2[i].split(): l2 = i

f1, f2 = [], []
for i in range(Natom):
 f1.append(map(float,out1[l1+2+i].split()[3:]))
 f2.append(map(float,out2[l2+2+i].split()[3:]))

for i in range(Natom): diff(f1[i],f2[i])
