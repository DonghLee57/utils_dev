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
    insert_mol = sys.argv[4]
    
    # For calculating energy difference, energies of the reactant, product molecules should calculated.
    mol_list = glob.glob('./molecules/*')
    mol_E = {}
    for idx, m in enumerate(mol_list):
        mol = ase.io.read(m,index=':',format='vasp')
        mol_E[m.split('_')[-1]] = get_lmp_E(mol,nnp)

    base = ase.io.read(input_str,index=':',format='vasp')
    E_old = get_lmp_E(base, nnp)

    # Scan the active sites for precursors or reactants
    # Insert the molecule as the product form
    base_new  = insert_molecule(base, insert_mol)
    E_new = get_lmp_E(base_new,nnp)

    # Reaction: surface + insert_mol -> (surface+insert_mol*) + product
    
    
def write_lmp_script(str_name, pot_name, elements):
    periodic_table = [
    'H', 'He',
    'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
    'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
    'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn','Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
    'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd','In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
    'Cs', 'Ba',    'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
               'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi','Po', 'At', 'Rn',
    'Fr', 'Ra',    'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk','Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
               'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Nh', 'Fl', 'Mc','Lv', 'Ts', 'Og']
    e_dict ={}
    for Z,symbol in enumerate(periodic_table):
        e_dict[symbol] = Z+1
    # In the present ASE, num(Cl) < num(N)...
    e_dict['Cl'] = 6
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

def get_lmp_E(poscar, nnp, str4lmp='str.dat'):
    pos = ase.io.read(poscar,index=':',format='vasp')
    types = list(set(pos[0].get_chemical_symbols()))
    ase.io.write(str4lmp, images=pos, format='lammps-data')
    write_lmp_script(str4lmp,nnp,types)
    lmp = lammps.lammps()
    lmp.file('py_script.in')
    E = get_lmp_log('energy')
    lmp.close()
    return E

def insert_molecule(base, mol_name):
    lat= base[0].cell
    ilat= np.linalg.inv(lat)
    types = base[0].get_chemical_symbols()

    # Search the active sites on surface
    ### crt_h...arbitrary now...
    crt_h = 0.4
    candidates = []
    for idx, pos in enumerate(base[0].positions):
        if np.matmul(pos,ilat)[2] > crt_h:
            if types[idx] == 'X':
                candidates.append(idx)

    if len(candidates) > 0 :
        S = np.random.randint(0,len(candidates))
        ### cutoff setting
        cut = [1.]*len(base[0].positions)
        nei =  neighborlist.NeighborList(cut, self_interaction=False,bothways=True)
        nei.update(base[0])
        nei_idx = nei.get_neighbors(candidates[S])[0]
        vec = np.array([0,0,0])
        for idx, i in enumerate(nei_idx):
            vec = vec + get_closePOS(lat, base[0].positions[candidates[S]], base[0].positions[i])-base[0].positions[candidates[S]]
        vec = vec/np.linalg.norm(vec)#*(-1)*2.2
        #print(vec)

        # insert molecule
        mol = ase.io.read('./molecules/POSCAR_'+mol_name,index=':',format='vasp')
        mol[0].positions = mol[0].positions - np.mean(mol[0].positions,axis=0)
        ### how to set universal molecule axis??
        mol_axis = mol[0].positions[1] - mol[0].positions[0]
        theta = vec_ang(vec, mol_axis)
        rot_axis = np.cross(vec, mol_axis)
        ### theta, coordinates sign...
        mol[0].positions = (-1)*rot_coords(mol[0].positions, -theta, rot_axis)

        types = mol[0].get_chemical_symbols()
        for idx, e in enumerate(types):
            if idx != 1:
             base[0].append(e)
             base[0].positions[-1] = mol[0].positions[idx] + base[0].positions[candidates[S]] + vec*(-1)*2.2
        types = list(set(base[0].get_chemical_symbols()))

    else:
        print('There is no atom near the surface to satisfy the conidtion. Stop the program.')
        sys.exit()

    return base

##########################
if __name__ == "__main__":
    main()
