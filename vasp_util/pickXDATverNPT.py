import sys, os
import numpy as np

XDAT = open(sys.argv[1] , 'r').readlines()
NA   = np.sum(np.array(XDAT[6].replace('\n','').strip().split()).astype(np.int))
frame = int(sys.argv[3])

header = ''
for i in range(7):
 header += XDAT[i+(NA+8)*(frame-1)]
line = (NA+8)*(frame-1)
posis = XDAT[line+8 :line+NA+8]

fil = open(sys.argv[2], 'w')
fil.write(header)
fil.write('Selective dynamics\n')
fil.write('%s'%(XDAT[i+1+(NA+8)*(frame-1)]))
for i in range(NA):
 fil.write(posis[i])
fil.close()

