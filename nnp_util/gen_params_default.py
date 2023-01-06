# python (.py) N > params_X
#
# N: the number of elements
#
import sys
params = ''
ntypes = int(sys.argv[1])

Rc = 6.0
#G2
eta = [0.003214,0.035711,0.071421,0.124987,0.214264,0.357106,0.714213,1.428426]
for i in range(ntypes):
    for e in eta:
        params += f'{2}{i+1:>2}{0:>2}{Rc:>6.1f}{e:>10.6f}{0.0:>6.1f}{0.0:>6.1f}\n'

#G4
eta = [0.000357,0.028569,0.089277]
zeta = [1,2,4]
lam = [1,-1]
for i in range(ntypes):
    for j in range(ntypes):
        if j >= i:
            for l in lam:
                for z in zeta:
                    for e in eta:
                        params += f'{4}{i+1:>2}{j+1:>2}{Rc:>6.1f}{e:>10.6f}{z:>6.1f}{l:>6.1f}\n'

print(params[:-1])                                                                              
