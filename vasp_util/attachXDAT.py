#Usauge:
# python (this_file.py) XDATCAR1 XDATCAR2
# XDATCAR2 - Its head part, above the 'Direct configuration', have to be removed.
#Output:
# xdat_attached
output='xdat_attached'

import sys, os

xdat1=sys.argv[1]
xdat2=sys.argv[2]

###attach files
os.system('cat %s %s > %s' %(xdat1, xdat2, 'tp'))

###Open output_file
out=open('tp','r')
tot=[]
for i in enumerate(open('tp','r')):
 line=i[0]
for i in range(line+1):
 tot.append(out.readline().split())
out.close()
os.system('rm tp')

###Necessary Values
numATOM=sum(map(int,tot[6]))
head=7
Iter=(line+1-head)/(numATOM+1)

###Write output again
op=open(output,'w')
c=0
for i in range(line+1):
 for l in range(len(tot[i])):
  if l == len(tot[i])-1:
   if i%((numATOM+1)*c+7) == 0:
    c=c+1
    op.write('%s\n' %(str(c)))
   else:
    op.write('%s\n' %(tot[i][-1]))
  else:
   op.write('%s ' %(tot[i][l]))
