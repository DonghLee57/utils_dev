myLAMMPS = '/path/lammps/python'
import sys
sys.path.append(myLAMMPS)

import lammps
import ase.io
from ase import Atoms

def main():
    input_str = sys.argv[1]
    str4lmp = sys.argv[2]
    nnp = sys.argv[3]
    
    info = readPOSCAR(input_str)
    poscar = ase.io.read(input_str,index=':',format='vasp')

    ase.io.write(str4lmp, images=poscar, format='lammps-data')
    write_lmp_script(str4lmp,nnp,info[1])

    lmp = lammps.lammps()
    lmp.file('py_script.in')
    E0 = get_lmp_log('energy')
    lmp.close()

def write_lmp_script(str_name, pot_name, elements):
    e_dict={'H':1, 'N':7, 'O':8, 'Cl':17,'Ti':22,'Mo':42}
    mass_dict={'H':1, 'N':14.001, 'O':15.999, 'Cl':35.453,'Ti':47.880}
    o=open(pot_name,'r').readline()
    pot_elements = o.split()[1:]
    o=open('py_script.in','w')
    header='units           metal\n'
    header+='boundary        p p p\n'
    header+='atom_style      atomic\n'
    header+='box tilt        large\n'
    header+='read_data       %s\n' %str_name
    header+='\npair_style      nn\n'

    elements.sort(key=lambda x: e_dict[x])
    for i in range(len(elements)):
        pot_elements.remove(elements[i])
    coeff_ele = elements + pot_elements

    header+='pair_coeff      * * %s ' %pot_name
    header+=' '.join(coeff_ele)
    header+='\n\n'
    for idx, e in enumerate(coeff_ele):
        header+=f'mass\t{idx+1} {mass_dict[e]}\n'

    header+='\nrun 0'
    o.write(header)
    return 1


def get_lmp_log(prop):
    o=enumerate(open('log.lammps'))
    for idx, line in o:
        tmp = line.split()
        if len(tmp) > 0:
            if tmp[0] =='Step':
                step, T, PE = list(map(float,next(o)[-1].split()[:3]))
    if prop.lower() == 'energy': return PE

##########################
if __name__ == "__main__":
    main()
