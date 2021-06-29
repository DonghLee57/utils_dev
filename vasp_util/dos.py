#######################################################################
#                                                                     #
#         Extracting TDOS and PDOS from DOSCAR of VASP result         #
#                                                                     #
#      e.g. python dos.py DOSCAR out_file_prefix > pdos_sum.dat       #
#                                                                     #
#######################################################################

import sys

try:
    DOSCAR = sys.argv[1]; filePrefix = sys.argv[2]
except:
    print "Usage:", sys.argv[0], "DOSCAR out_file_prefix"; sys.exit(1)
    
#######################################################################
# FILL THE FOLLOWINGS MANUALLY !                                      #
tdosCalc = 1               # 1 with tdos calc. / 0 without tdos calc.
spinCalc = 0               # 1 for spin polar. / 0 for non spin polar.
targetAtomList = [0,1,2,3,4]              #list of target atoms (NOTE: index - 1) / [] for all
targetOrbital = []     # list of target orbitals (VASP convention + u/d) / [] for all
#######################################################################

dosData = open(str(DOSCAR),'r')

# Reading headers: 1st line gives total number of atoms (numAtom)
#                  last line gives title of the simulation (simTitle)

dosReadLine = dosData.readline().split()
numAtom = map(lambda x:int(x), dosReadLine)[0]

dosData.readline()
dosData.readline()
dosData.readline()

simTitle = str(dosData.readline())

#######################################################################
# Necessary lists and dictionaries to be used                         #
## List for energy scale value
energyScaleValue = []
## List for PDOS information archive
atomOrbitalArchive = [[] for i in xrange(numAtom)]
## Orbital dictionary
nonpolarOrbitalDic = {'s':0, \
              'py':1, 'pz':2, 'px':3, \
              'dxy':4, 'dyz':5, 'dz2':6, 'dxz':7, 'dx2':8} #, \
             # 'f1':9, 'f2':10, 'f3':11, 'f4':12, 'f5':13, 'f6':14, 'f7':15}
orbitalDic = {'su':0, 'sd':1, \
              'pyu':2, 'pyd':3, 'pzu':4, 'pzd':5, 'pxu':6, 'pxd':7,\
              'dxyu':8, 'dyzu':10, 'dz2u':12, 'dxzu':14, 'dx2u':16, \
              'dxyd':9, 'dyzd':11, 'dz2d':13, 'dxzd':15, 'dx2d':17}
              #, \
              #'f1u':18, 'f2u':20, 'f3u':22, 'f4u':24, 'f5u':26, 'f6u':28, 'f7u':30, \
              #'f1d':19, 'f2d':21, 'f3d':23, 'f4d':25, 'f5d':27, 'f6d':29, 'f7d':31}
## Target Atom List
if (not len(targetAtomList)):
    targetAtomList = [i for i in xrange(numAtom)]
targetAtomList.sort()
## Target Orbital List
if (not len(targetOrbital)):
    if (spinCalc):
        targetOrbitalList = orbitalDic.values()
    else:
        targetOrbitalList = nonpolarOrbitalDic.values()
else:
    if (orbitalDic.has_key(targetOrbital[0])):
        targetOrbitalList = [orbitalDic.get(i) for i in targetOrbital]
    else:
        targetOrbitalList = [nonpolarOrbitalDic.get(i) for i in targetOrbital]
targetOrbitalList.sort()
#######################################################################

# TDOS file generation when tdosCalc = 1
#

if (tdosCalc):
    tdosFile = open(str(filePrefix)+'_tdos.dat','w')
    numEnergyLevel = int(dosData.readline().split()[2])
    for member in range(numEnergyLevel):
        tdosLine = dosData.readline()
        energyScaleValue.append(float(tdosLine.split()[0]))
        tdosFile.write(tdosLine)
    tdosFile.close()
else:
    numEnergyLevel = 3000
    energyScaleValue = [0 for i in range(numEnergyLevel)]
# PDOS manipulation
#

## Filling up the PDOS-related lists: atomOrbitalArchive[atom][energy]

for atomMember in range(numAtom):
    dosData.readline()
    for member in range(numEnergyLevel):
        pdosTemp = map(lambda x:float(x),dosData.readline().split())
        atomOrbitalArchive[atomMember].append(pdosTemp[1:])
        atomOrbitalArchive[atomMember][member].append(reduce(lambda x,y:x+y,atomOrbitalArchive[atomMember][member]))
dosData.close()

## Spin summation for the targets: spinSum[energy][spin]

spinSum = [[0]*len(targetOrbitalList) for i in xrange(numEnergyLevel)]

for atomTarget in targetAtomList:
    for member in range(numEnergyLevel):
        targetLine = [atomOrbitalArchive[atomTarget][member][i] for i in targetOrbitalList]
        spinSum[member] = map(lambda x,y:x+y, spinSum[member], targetLine)

for i in range(numEnergyLevel):
    print energyScaleValue[i],'\t',
    for j in range(len(spinSum[i])):
        print spinSum[i][j],'\t',
    print





