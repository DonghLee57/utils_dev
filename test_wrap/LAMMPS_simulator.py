import numpy as np
import sys, subprocess
from ase.io import read, write
from ase.data import atomic_numbers, atomic_names, atomic_masses, chemical_symbols

# Modify the path provided below to suit your environment.
LAMMPS = 'lmp_binary'
MPIRUN = 'mpirun'
NNCORE = sys.argv[1]

def main():
  POT = sys.argv[2]
  pot_type = 'NNP'
  # If not NNP,
  symbols = [['A'],['B']]
  
  if pot_type == 'NNP':
    REF_FILE = './POSCAR_generated.vasp'
    symbols = get_sorted_symbols(REF_FILE)
    tmp = read(REF_FILE, format='vasp')
    write('test.lammps', images=tmp, format='lammps-data')
    
  my = SIMULATOR(POT, pot_type)
  #def INIT_SCRIPT(self, FILE, symbols, unit='metal', style='atomic', pbc=[1,1,1]):
  #def REPLICATE(self, nx, ny, nz):
  #def MINIMIZE(self, cell=False):
  #def ANNEAL(self, T, TIME=500, ensemble='nvt'):
  #def QUENCH(self, T_init, T_end=300.0, Qrate=10.0, ensemble='nvt'):
  my.INIT_SCRIPT('test.lammps', symbols,'metal','atomic',[1,1,1])
  my.MINIMIZE()
  my.REPLICATE(2,2,2)
  my.ANNEAL(1000, 10)
  my.QUENCH(1000)
  my.MINIMIZE(cell=True)

  #def RUN_SCRIPT(self, SCRIPT=self.SCRIPT):
  my.RUN_SCRIPT(my.SCRIPT)


###
class SIMULATOR:
  def __init__(self, POT, pot_type):
    self.TIMESTEP = 0.001
    self.DUMPSTEP = 1000
    self.POT = POT
    self.pot_type = pot_type
    if self.pot_type == 'NNP':
      self.symbols = subprocess.check_output(['head',f'{self.POT}','-n','1'],universal_newlines=True).split()[1:]
      self.NTYPES = len(self.symbols)
    elif self.pot_type in ['tersoff','sw'] :
      self.symbols = None
    self.SCRIPT = 'script.in'
    self.EXEC = 0
    self.LOG = 'simulator.log'
    self.create_vel = True
    self.write_dump = True
    self.run_MD = True

  def INIT_SCRIPT(self, FILE, symbols, unit='metal', style='atomic', pbc=[1,1,1]):
    if self.EXEC == 0:
      with open(self.LOG,'w') as o: o.write(f'{self.EXEC:4d} : Initialize "{SCRIPT}"\n')
    else:
      with open(self.LOG,'w') as o: o.write(f'{self.EXEC:4d} : You tried to initialize the script again. You must initialize it once at first.\n')
      sys.exit()
    header  = f'processors\t* * * grid numa\n'
    header  = f'units\t\t{unit}\n'
    header += f'atom_style\t{style}\n'
    PBC = []
    for _, b in enumerate(pbc):
        if b: PBC.append('p')
        else: PBC.append('f')
    header += f'boundary\t{" ".join(PBC)}\n'
    header += f'box\t\ttilt large\n'
    header += f'read_data\t{FILE}\n\n'
    header += set_pot(self.POT, self.pot_type, symbols, self.symbols)
    if self.pot_type == 'NNP': NNP_prepare(FILE, self.NTYPES)
    header += f'thermo\t\t{self.DUMPSTEP} \n'
    header +=  'thermo_modify\tlost ignore flush yes\n\n'
    with open(self.SCRIPT,'w') as o: o.write(header)

  def RUN_SCRIPT(self, SCRIPT=None):
    self.EXEC += 1
    if SCRIPT==None: SCRIPT = self.SCRIPT
    res = subprocess.check_output([MPIRUN, '-n', NNCORE, LAMMPS, '-in', SCRIPT],universal_newlines=True)
    with open('0_log.x','w') as o: o.write(res)
    if self.EXEC == 1:
      with open(self.LOG,'w') as o: o.write(f'{self.EXEC:4d} : LAMMPS running with "{SCRIPT}"\n')
    else:
      with open(self.LOG,'a') as o: o.write(f'\n{self.EXEC:4d} : LAMMPS running with "{SCRIPT}"\n')

  def REPLICATE(self, nx, ny, nz):
    self.EXEC += 1
    with open(self.SCRIPT,'a') as o: o.write(f'replicate\t\t{nx} {ny} {nz}\n')
    if self.EXEC == 1:
      with open(self.LOG,'w') as o: o.write(f'{self.EXEC:4d} : The structure is replicated along the x, y, and z axes by a factor of {nx}, {ny}, and {nz}, respectively\n')
    else:
      with open(self.LOG,'a') as o: o.write(f'\n{self.EXEC:4d} : The structure is replicated along the x, y, and z axes by a factor of {nx}, {ny}, and {nz}, respectively\n')

  def MINIMIZE(self, cell=False):
    self.EXEC += 1
    with open(self.SCRIPT,'a') as o:
      if cell: o.write('fix\t\t1 all box/relax iso 0.0 vmax 0.001 fixedpoint 0 0 0\n')
      o.write('minimize 1e-6 1e-10 10000 100000\n')
      if cell: o.write('unfix\t\t1\n')
      o.write(f'write_data\t{self.EXEC}_MINIMIZE.lammps\n\n')
    if self.EXEC == 1:
      with open(self.LOG,'w') as o: o.write(f'{self.EXEC:4d} : Minimize the structure.\n')
    else:
      with open(self.LOG,'a') as o: o.write(f'\n{self.EXEC:4d} : Minimize the structure.\n')

  def ANNEAL(self, T, TIME=500, ensemble='nvt'):
    self.EXEC += 1
    It = int(TIME/self.TIMESTEP)
    with open(self.SCRIPT,'a') as o:
      if self.create_vel:
        o.write(f'velocity\tall create {T} {np.random.randint(low=1,high=1000)} dist gaussian\n')
        self.create_vel = False
      if ensemble =='nvt':   o.write(f'fix\t\tTHERMOSTAT all nvt temp {T} {T} 1\n')
      elif ensemble =='npt': o.write(f'fix\t\tTHERMOSTAT all npt temp {T} {T} 1 iso 0.0 0.0 1000.0 fixedpoint 0 0 0\n')
      if self.write_dump:
        o.write(f'dump\t\tDUMP all custom {self.DUMPSTEP} dump.lammpstrj id type x y z\n')
        o.write('dump_modify\tDUMP sort id\n\n')
        self.write_dump = False
      if self.run_MD:
        o.write(f'timestep\t\t{self.TIMESTEP:.4f}\n')
        self.run_MD = False
      o.write(f'run\t\t{It}\n')
      o.write(f'unfix\t\tTHERMOSTAT\n')
      o.write(f'write_data\t{self.EXEC}_ANNEAL.lammps\n\n')
    if self.EXEC == 1:
      with open(self.LOG,'w') as o: o.write(f'{self.EXEC:4d} : Annealing at {T} K during {TIME} ps.\n')
    else:
      with open(self.LOG,'a') as o: o.write(f'{self.EXEC:4d} : Annealing at {T} K during {TIME} ps.\n')

  def QUENCH(self, T_init, T_end=300.0, Qrate=10.0, ensemble='nvt'):
    self.EXEC += 1
    It = int(np.fabs(T_init-T_end)/(Qrate*self.TIMESTEP))
    with open(self.SCRIPT,'a') as o:
      if self.create_vel:
        o.write(f'velocity\tall create {T} {np.random.randint(low=1,high=1000)} dist gaussian\n')
        self.create_vel = False
      if self.write_dump:
        o.write(f'dump\t\tDUMP all custom {self.DUMPSTEP} dump.lammpstrj id type x y z\n')
        o.write('dump_modify\t\tDUMP sort id\n\n')
        self.write_dump = False
      if ensemble =='nvt':   o.write(f'fix\t\tTHERMOSTAT all nvt temp {T_init} {T_end} 1\n')
      elif ensemble =='npt': o.write(f'fix\t\tTHERMOSTAT all npt temp {T_init} {T_end} 1 iso 0.0 0.0 1000.0 fixedpoint 0 0 0\n')
      if self.run_MD:
        o.write(f'timestep\t\t{self.TIMESTEP:.4f}\n')
        self.run_MD = False
      o.write(f'run\t\t{It}\n')
      o.write(f'unfix\t\tTHERMOSTAT\n')
      o.write(f'write_data\t{self.EXEC}_QUENCH.lammps\n\n')
    if self.EXEC == 1:
      with open(self.LOG,'w') as o: o.write(f'{self.EXEC:4d} : Quenching with {Qrate} K/ps from {T_init} to {T_end}.\n')
    else:
      with open(self.LOG,'a') as o: o.write(f'{self.EXEC:4d} : Quenching with {Qrate} K/ps from {T_init} to {T_end}.\n')

  def DEL_OVERLAP(self, distance, g1=None, g2=None):
    self.EXEC +=1
    if g1 == None: g1 = 'all'
    if g2 == None: g2 = 'all'
    with open(self.SCRIPT,'a') as o:
      o.write('delete_atoms\t\toverlap {distance:.3f} {g1} {g2}\n')
    if self.EXEC == 1:
        with open(self.LOG,'w') as o: o.write(f'{self.EXEC:4d} : Delete overlapping atoms within {distance:.3f} Angstrom.')
    else:
        with open(self.LOG,'w') as o: o.write(f'\n{self.EXEC:4d} : Delete overlapping atoms within {distance:.3f} Angstrom.')

  def MK_BLOCK_GROUP(self, name, x=None, y=None, z=None):
    self.EXEC +=1
    r_cmd = f'region\t\tr{name} block'
    if x != None: r_cmd += f' {x[0]} {x[1]} '
    else:         r_cmd += f' EDGE EDGE '
    if y != None: r_cmd += f' {y[0]} {y[1]} '
    else:         r_cmd += f' EDGE EDGE '
    if z != None: r_cmd += f' {z[0]} {z[1]} '
    else:         r_cmd += f' EDGE EDGE '
    with open(self.SCRIPT,'a') as o:
      o.write(r_cmd) 
      o.write(f'group\t\t{name} region r{name}\n')
    if self.EXEC == 1:
        with open(self.LOG,'w') as o: o.write(f'{self.EXEC:4d} : Group atoms in the region as {name}.')
    else:
        with open(self.LOG,'w') as o: o.write(f'\n{self.EXEC:4d} : Group atoms in the region as {name}.')

  def NEB(self, FILE, symbols, NEBSTEP):
    self.EXEC += 1
    """
    # Check running
    with open(self.SCRIPT,'w') as o:
      o.write(f'read_data\t{FILE}\n\n')
      o.write(set_pot(self.POT, self.pot_type, symbols, self.symbols))
      o.write(f'variable\t\tu uloop {NEBSTEP}\n')
      o.write(f'reset_timestep 0\n')
      o.write(f'timestep {self.TIMESTEP:.4f}\n')
      o.write(f'dump\t\tDUMP all custom {self.DUMPSTEP} dump_$u.lammpstrj id type x y z\n')
      o.write('dump_modify\tDUMP sort id\n\n')
      o.write('### quickmin fire cg sd hftn\n')
      o.write('min_style\t\tquickmin\n')
      o.write('neb\t\t0.0 0.01 100 100 10 each coords_$i.lammps\n')
    """
    
def C2K(T):
  return T+273.15
  
def get_sorted_symbols(FILE):
  sample = read(FILE,format='vasp')
  write('test.lammps',images=sample,format='lammps-data')
  compare = read('test.lammps',style='atomic',format='lammps-data')
  vasp_symbols = list(set(sample.get_chemical_symbols()))
  lammps_symbols = list(set(compare.get_chemical_symbols()))
  sorted_symbols = []
  for i, item in enumerate(lammps_symbols):
    sorted_symbols.append([vasp_symbols[i],atomic_numbers[item]])
  sorted_symbols.sort()
  return sorted_symbols

# Matching the number of atom types with the number of elements in NNP
def NNP_prepare(FILE, NTYPES):
  with open(FILE,'r') as o: tmp=o.readlines()
  change    = tmp[3].split()
  change[0] = f'{NTYPES}'
  tmp[3] = ' '.join(change) + '\n'
  with open(FILE,'w') as o:
    for i, item in enumerate(tmp):
      o.write(item)

# Write the part related to potential and atom types in the LAMMPS script
def set_pot(POT, pot_type, symbols, all_symbols):
  lines, elements, masses, count, bins = '', '', '', 0, []
  for i, item in enumerate(symbols):
    bins.append(item[0])
    elements += f' {item[0]}'
    masses += f'mass\t\t{i+1} {atomic_masses[atomic_numbers[item[0]]]}\n'
  if pot_type == 'NNP':
    pot_type = 'nn'
    #pot_type = 'nn/intel' #if SIMD compile
    for i, item in enumerate(all_symbols):
      if not item in bins:
        count    += 1
        elements += f' {item}'
        masses   += f'mass\t\t{count+len(symbols)} {atomic_masses[atomic_numbers[item]]}\n'
  #elif pot_type == 'tersoff':
  #  pot_type = 'tersoff'
  lines += f'pair_style  {pot_type}\npair_coeff  * * "{POT}" {elements}\n'
  lines += masses
  return lines
      
###
if __name__ == "__main__":
  main()
