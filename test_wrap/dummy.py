# Write LAMMPS structure file with dummy
pot_type = ['Ti','Cl','N','H']
nelements = len(pot_type)

mol = ase.io.read(gpath + 'TiCl4_mol.poscar',index=':',format='vasp')
types = list(set(mol[0].get_chemical_symbols()))
ndummy = 4 - len(types)
for idx, p in enumerate(pot_type):
  if not p in types:
    mol[0].append(p)

ase.io.write(gpath + '/TiCl4_mol.lammps',images=mol,format='lammps-data')
tmp = open(gpath+'/TiCl4_mol.lammps','r').readlines()
for i in range(len(tmp)-2):
  line = tmp[i].split()
  if 'atoms' in line: print(f"{int(line[0])-2} {line[1]}")
  else: print(tmp[i])
