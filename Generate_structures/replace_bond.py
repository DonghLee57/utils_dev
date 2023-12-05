import sys
import numpy as np
from ase.io import read, write
from ase import Atoms

#
IDX4REP = np.array([0,1,2,3])
CIDX_REF= -1
IDX4RM  = 0
CIDX_MOL= -1

#
def main():
    REF_STR = read(sys.argv[1],format='vasp')

    # Replace the assinged bonds
    NEW = Atoms('X',cell=REF_STR.cell, pbc=[1,1,1])
    NEW += REF_STR[CIDX_REF]
    for bidx, bond in enumerate(IDX4REP):
        SUB_MOL = read(sys.argv[2],format='vasp')
        SUB_MOL.translate(-SUB_MOL.get_center_of_mass())
        AX_m = SUB_MOL.positions[CIDX_MOL] - SUB_MOL.positions[IDX4RM]
        del SUB_MOL[IDX4RM]
        AX_b = REF_STR.positions[CIDX_REF] - REF_STR.positions[bond]
        AX_r = np.cross(AX_m, AX_b)
        theta = np.arccos(np.round(np.dot(AX_m,AX_b)/np.linalg.norm(AX_m)/np.linalg.norm(AX_b),4))
        if theta < 1E-4:
            SUB_MOL.positions = (-1)*(SUB_MOL.positions)
        else:
            SUB_MOL.positions = (-1)*rot_coords(SUB_MOL.positions, -theta, AX_r)
        SUB_MOL.translate(REF_STR.positions[bond])
        NEW += SUB_MOL.copy()

    # Write the output       
    del NEW[0]
    NEW = NEW[NEW.numbers.argsort()]
    write('POSCAR_replace.vasp',images=NEW, format='vasp')

#
def rot_coords(coords, theta, axis):
    # theta: radian
    # axis: arbitarary axis to rotate
    ux, uy, uz = axis/np.linalg.norm(axis)
    R = np.array([[np.cos(theta)+ux**2*(1-np.cos(theta)),       ux*uy*(1-np.cos(theta))-uz*np.sin(theta),    ux*uz*(1-np.cos(theta))+uy*np.sin(theta)],\
                  [uy*ux*(1-np.cos(theta))+uz*np.sin(theta),    np.cos(theta)+uy**2*(1-np.cos(theta)),       uy*uz*(1-np.cos(theta))-ux*np.sin(theta)],\
                  [uz*ux*(1-np.cos(theta))-uy*np.sin(theta),    uz*uy*(1-np.cos(theta))+ux*np.sin(theta),    np.cos(theta)+uz**2*(1-np.cos(theta))]])
    return np.matmul(coords,R)

def vec_ang(v1, v2, exp='rad'):
    # exp: 'rad' or 'deg'
    angle = np.arccos(np.round(np.fabs(np.dot(v1,v2))/np.linalg.norm(v1)/np.linalg.norm(v2),4))
    if exp == 'rad':
        return angle
    elif exp == 'deg':
        return angle*180/np.pi

#
if __name__ == "__main__":
    main()
