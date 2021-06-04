#Python2
import sys
import numpy as np

XDAT = open(sys.argv[1] , 'r').readlines()
header = ''
for i in range(7):
 header += XDAT[i]
posis = XDAT[7:]

tot_atnum = np.sum(np.array(XDAT[6].replace('\n','').strip().split()).astype(np.int))

i = int(sys.argv[3])
fil = open(sys.argv[2], 'w')
fil.write(header)
#fil.write('Selective dynamics\n')
fil.write(posis[(i-1)*(tot_atnum+1)])
for j in range((i-1)*(tot_atnum+1) + 1, (i)*(tot_atnum+1)):
 fil.write(posis[j])
