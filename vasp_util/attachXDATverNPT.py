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
tot = out.readlines()
print tot[0]
print tot[1]
print tot[2]
print tot[3]
print tot[4]
#os.system('rm tp')
###Necessary Values
numATOM=sum(map(int,tot[6].split()))
Iter=len(tot)/(numATOM+8)

###Write output again
op=open(output,'w')
#for i in range(1):
for i in range(Iter):
 for h  in range(7):
  op.write('%s'%tot[h+(numATOM+8)*(i)])
 op.write('Direct configuration = %d\n'%(i+1)) 
 for n in range(numATOM):
  op.write('%s' %tot[n+8+(numATOM+8)*(i)])
#  print h+(numATOM+8)*(i-1),8+(numATOM+8)*(i)
op.close()
