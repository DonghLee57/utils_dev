#########################################################
#               gen_cry_varient_cell.py                 #
#                                                       #
#       CODE for generating training set using          #
#       crystal structure.                              #
#       In this code, lattice parameter of              #
#       input crystal structure is slightly changed.    #
#                                                       #
#       USAGE: python code_name.py [input structure]    #
#########################################################

import sys
import Cellinfo
import numpy as np

# Function for generating rotation matrix
def genRotMatrix(axis, angle):
	rot_matrix = np.zeros([3,3])
	rot_matrix[0][0] = np.cos(angle) + axis[0]**2 *(1-np.cos(angle))
	rot_matrix[0][1] = axis[0]*axis[1]*(1-np.cos(angle)) - axis[2]*np.sin(angle)
	rot_matrix[0][2] = axis[0]*axis[2]*(1-np.cos(angle)) + axis[1]*np.sin(angle)
	rot_matrix[1][0] = axis[1]*axis[0]*(1-np.cos(angle)) + axis[2]*np.sin(angle)
	rot_matrix[1][1] = np.cos(angle) + axis[1]**2 *(1-np.cos(angle))
	rot_matrix[1][2] = axis[1]*axis[2]*(1-np.cos(angle)) - axis[0]*np.sin(angle)
	rot_matrix[2][0] = axis[2]*axis[0]*(1-np.cos(angle)) + axis[1]*np.sin(angle)
	rot_matrix[2][1] = axis[2]*axis[1]*(1-np.cos(angle)) + axis[0]*np.sin(angle)
	rot_matrix[2][2] = np.cos(angle) + axis[2]**2 *(1-np.cos(angle))

	return rot_matrix

# training POSCAR index
st_idx = 1

# get input POSCAR to Poscar class in Cellinfo.py
POS = Cellinfo.Poscar(sys.argv[1])
# change POSCAR type to Direct(if not)
# Because in this code, we change lattice parameter,
# direct coordinate is needed.
POS.changeCandD('Direct')
# dummy Poscar class variable 
tPOS = Cellinfo.Poscar(sys.argv[1])

# Lattice scaling ratio list
lat_scale_list = np.linspace(0.95, 1.05, 20)
# Lattice expansion and compression
for i in lat_scale_list:
    tPOS.lparam[0] *= i
    tPOS.lparam[1] *= i
    tPOS.lparam[2] *= i

    tPOS.writePOS('POSCAR_T' + str(st_idx))

    tPOS.lparam[0] = POS.lparam[0]
    tPOS.lparam[1] = POS.lparam[1]
    tPOS.lparam[2] = POS.lparam[2]

    st_idx += 1
print st_idx-1
# orthorhombic strain (xx, yy, zz / volume conserving)
ost_scale_list = np.linspace(-0.05, 0.05, 21)
print len(ost_scale_list), st_idx
for i in ost_scale_list:
    for j in range(1):
        tPOS.lparam[j] *= 1+i
        tPOS.lparam[(j+1)%3] *= np.sqrt(1/(1+i))
        tPOS.lparam[(j+2)%3] *= np.sqrt(1/(1+i))

        tPOS.writePOS('POSCAR_T' + str(st_idx))

        tPOS.lparam[0] = POS.lparam[0]
        tPOS.lparam[1] = POS.lparam[1]
        tPOS.lparam[2] = POS.lparam[2]

        st_idx += 1

# shear strain (xy, yz, zx / volume conserving)
sst_scale_list = np.linspace(0., 0.05, 21)[1:]
print len(sst_scale_list), st_idx
for i in sst_scale_list:
    tPOS.lparam[0] += i*tPOS.lparam[1] # xy
    tPOS.writePOS('POSCAR_T' + str(st_idx))
    tPOS.lparam[0] = POS.lparam[0]
    st_idx += 1
"""

    tPOS.lparam[0] += i*tPOS.lparam[2] # xz
    tPOS.writePOS('POSCAR_T' + str(st_idx))
    tPOS.lparam[0] = POS.lparam[0]
    st_idx += 1

    tPOS.lparam[1] += i*tPOS.lparam[2] # yz
    tPOS.writePOS('POSCAR_T' + str(st_idx))
    tPOS.lparam[1] = POS.lparam[1]
    st_idx += 1

# Reciprocal lattice for angular distortion
reci_matrix = np.zeros([3,3])
reci_matrix[0] = np.cross(POS.lparam[1],POS.lparam[2]) / np.dot(POS.lparam[0], np.cross(POS.lparam[1],POS.lparam[2]))
reci_matrix[1] = np.cross(POS.lparam[2],POS.lparam[0]) / np.dot(POS.lparam[1], np.cross(POS.lparam[2],POS.lparam[0]))
reci_matrix[2] = np.cross(POS.lparam[0],POS.lparam[1]) / np.dot(POS.lparam[2], np.cross(POS.lparam[0],POS.lparam[1]))
# Normalize
reci_matrix[0] /= np.linalg.norm(reci_matrix[0])
reci_matrix[1] /= np.linalg.norm(reci_matrix[1])
reci_matrix[2] /= np.linalg.norm(reci_matrix[2])

# Angular scaling ratio list
ang_scale_list = np.linspace(-2*np.pi/45, 2*np.pi/45, 10)

for i in ang_scale_list:
	for j in ang_scale_list:
		for k in ang_scale_list:
			tPOS.lparam[0] = np.dot(genRotMatrix(reci_matrix[2], i), tPOS.lparam[0].reshape([3,1])).transpose()
			tPOS.lparam[0] = np.dot(genRotMatrix(reci_matrix[1], j), tPOS.lparam[0].reshape([3,1])).transpose()
			tPOS.lparam[1] = np.dot(genRotMatrix(reci_matrix[0], k), tPOS.lparam[1].reshape([3,1])).transpose()

			tPOS.writePOS('POSCAR_T' + str(st_idx))
			# Reset lattice parameter
			tPOS.lparam[0] = POS.lparam[0]
			tPOS.lparam[1] = POS.lparam[1]
			tPOS.lparam[2] = POS.lparam[2]

			st_idx += 1
"""

