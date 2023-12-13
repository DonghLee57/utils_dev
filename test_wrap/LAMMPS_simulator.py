import numpy as np
import os, sys, subprocess
from ase.io import read, write
from ase.data import atomic_numbers, atomic_names, atomic_masses, chemical_symbols

# Modify the path provided below to suit your environment.
LAMMPS = 'lmp_binary'
MPIRUN = 'mpirun'
NNCORE = sys.argv[1]

def main():
  POT = sys.argv[2]
  pot_type = 'NNP'
  
  if pot_type == 'NNP':
    REF_FILE = './POSCAR'
    symbols = get_sorted_symbols(REF_FILE)
    tmp = read(REF_FILE, format='vasp')
    write('test.lammps', images=tmp, format='lammps-data')
    
  my = SIMULATOR(POT, pot_type)
  # def by_script(self, SCRIPT):
  #my.by_script('test.in')

  # def ANNEAL(self, FILE, symbols, T, TIME=500, ensemble='nvt', min_init=False, min_end=False, out=None):
  #my.ANNEAL()

  # def MQ(self, FILE, symbols, T_melt, T_init, T_end=300.0, Qrate=10.0, ensemble='nvt', min_init=False, min_end=False, out=None):
  #my.MQ()

  # def NEB(self, FILE, symbols, NEBSTEP):
  #my.NEB()

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

  def by_script(self, SCRIPT):
    res = subprocess.check_output([MPIRUN,'-n',NNCORE,LAMMPS,'-in',SCRIPT],universal_newlines=True)
    with open(f'log.x','w') as o: o.write(res)
    return 0

  def ANNEAL(self, FILE, symbols, T, TIME=500, ensemble='nvt', min_init=False, min_end=False, out=None):
    if self.pot_type == 'NNP': NNP_prepare(FILE, self.NTYPES)
    It = int(TIME/self.TIMESTEP)
    with open(self.SCRIPT,'w') as o:
      o.write('processors\t* * * grid numa\n')
      o.write('units\t\tmetal\n')
      o.write('atom_style\tatomic\n')
      o.write('boundary\tp p p\n')
      o.write('box\t\ttilt large\n')
      o.write(f'read_data\t{FILE}\n\n')
      o.write(set_pot(self.POT, self.pot_type, symbols, self.symbols))
      o.write(f'thermo\t\t{self.DUMPSTEP} \n')
      o.write('thermo_modify\tlost ignore flush yes\n\n')
      if min_init:
        o.write('fix 1 all box/relax iso 0.0 vmax 0.001 fixedpoint 0 0 0\n')
        o.write('minimize 1e-6 1e-10 10000 100000\n')
        if out == None: o.write(f'write_data\tmin_init.lammps\n')
        else:           o.write(f'write_data\tmin_init_{out}.lammps\n')
      o.write(f'velocity\tall create {T} {np.random.randint(low=1,high=100)} dist gaussian\n')
      if ensemble =='nvt':   o.write(f'fix\t\tTHERMOSTAT all nvt temp {T} {T} 1\n')
      elif ensemble =='npt': o.write(f'fix\t\tTHERMOSTAT all npt temp {T} {T} 1 iso 0.0 0.0 1000.0 fixedpoint 0 0 0\n')
      o.write(f'timestep {self.TIMESTEP:.4f}\n')
      o.write(f'dump\t\tDUMP all custom {self.DUMPSTEP} dump_anneal_{T:d}.lammpstrj id type x y z\n')
      o.write('dump_modify\tDUMP sort id\n\n')
      o.write(f'run\t\t{It}\n')
      o.write(f'write_data\t{T:d}.lammps\n')
      if min_end:
        o.write('unfix THERMOSTAT\n')
        o.write('fix 1 all box/relax iso 0.0 vmax 0.001 fixedpoint 0 0 0\n')
        o.write('minimize 1e-6 1e-10 10000 100000\n')
        if out == None: o.write(f'write_data\tmin_end.lammps\n')
        else:           o.write(f'write_data\tmin_end_{out}.lammps\n')
      res = subprocess.check_output([MPIRUN,'-n',NNCORE,LAMMPS,'-in',self.SCRIPT],universal_newlines=True)
      with open(f'log_anneal_{T:d}.x','w') as o: o.write(res)

  def MQ(self, FILE, symbols, T_melt, T_init, T_end=300.0, Qrate=10.0, ensemble='nvt', min_init=False, min_end=False, out=None):
    if self.pot_type == 'NNP': NNP_prepare(FILE, self.NTYPES)
    It = int((T_init-T_end)/(Qrate*self.TIMESTEP))
    with open(self.SCRIPT,'w') as o:
      o.write('processors\t* * * grid numa\n')
      o.write('units\t\tmetal\n')
      o.write('atom_style\tatomic\n')
      o.write('boundary\tp p p\n')
      o.write('box\t\ttilt large\n')
      o.write(f'read_data\t{FILE}\n\n')
      o.write(set_pot(self.POT, self.pot_type, symbols, self.symbols))
      o.write(f'thermo\t\t{self.DUMPSTEP} \n')
      o.write('thermo_modify\tlost ignore flush yes\n\n')
      if min_init:
        o.write('fix 1 all box/relax iso 0.0 vmax 0.001 fixedpoint 0 0 0\n')
        o.write('minimize 1e-6 1e-10 10000 100000\n')
        if out == None: o.write(f'write_data\tmin_init.lammps\n')
        else:           o.write(f'write_data\tmin_init_{out}.lammps\n')
      o.write(f'velocity\tall create {T_melt} {np.random.randint(low=1,high=100)} dist gaussian\n')
      if ensemble =='nvt':   o.write(f'fix\t\tTHERMOSTAT all nvt temp {T_melt} {T_melt} 1\n')
      elif ensemble =='npt': o.write(f'fix\t\tTHERMOSTAT all npt temp {T_melt} {T_melt} 1 iso 0.0 0.0 1000.0 fixedpoint 0 0 0\n')
      o.write(f'dump\t\tDUMP all custom {self.DUMPSTEP} dump_mq.lammpstrj id type x y z\n')
      o.write('dump_modify\tDUMP sort id\n\n')
      o.write(f'timestep {self.TIMESTEP:.4f}\n')
      o.write(f'run\t\t50000\n')
      o.write('unfix THERMOSTAT\n')
      if ensemble =='nvt':   o.write(f'fix\t\tTHERMOSTAT all nvt temp {T_init} {T_end} 1\n')
      elif ensemble =='npt': o.write(f'fix\t\tTHERMOSTAT all npt temp {T_init} {T_end} 1 iso 0.0 0.0 1000.0 fixedpoint 0 0 0\n')
      o.write(f'run\t\t{It}\n')
      if min_end:
        o.write('unfix THERMOSTAT\n')
        o.write('fix 1 all box/relax iso 0.0 vmax 0.001 fixedpoint 0 0 0\n')
        o.write('minimize 1e-6 1e-10 10000 100000\n')
        if out == None: o.write(f'write_data\tmin_end.lammps\n')
        else:           o.write(f'write_data\tmin_end_{out}.lammps\n')
      res = subprocess.check_output([MPIRUN,'-n',NNCORE,LAMMPS,'-in',self.SCRIPT],universal_newlines=True)
      with open(f'log_mq.x','w') as o: o.write(res)
  
  def NEB(self, FILE, symbols, NEBSTEP):
    # Check running
    if self.pot_type == 'NNP': NNP_prepare(FILE, self.NTYPES)
    with open(self.SCRIPT,'w') as o:
      o.write('processors\t* * * grid numa\n')
      o.write('units\t\tmetal\n')
      o.write('atom_style\tatomic\n')
      o.write('boundary\tp p p\n')
      o.write('box\t\ttilt large\n')
      o.write(f'read_data\t{FILE}\n\n')
      o.write(self.set_pot('nn', self.POT, symbols, self.symbols))
      o.write(f'thermo\t\t{self.DUMPSTEP} \n')
      o.write('thermo_modify\tlost ignore flush yes\n\n')
      o.write(f'variable\t\tu uloop {NEBSTEP}\n')
      o.write(f'reset_timestep 0\n')
      o.write(f'timestep {self.TIMESTEP:.4f}\n')
      o.write(f'dump\t\tDUMP all custom {self.DUMPSTEP} dump_$u.lammpstrj id type x y z\n')
      o.write('dump_modify\tDUMP sort id\n\n')
      o.write('### quickmin fire cg sd hftn\n')
      o.write('min_style\t\tquickmin\n')
      o.write('neb\t\t0.0 0.01 100 100 10 each coords_$i.lammps\n')
      res = subprocess.check_output([MPIRUN,'-n',NNCORE,LAMMPS,'-in',self.SCRIPT],universal_newlines=True)
      with open(f'log_NEB.x','w') as o: o.write(res)
        
###
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
  if pot_type == 'NNP': pot_type = 'nn'
  lines, elements, masses, count, bins = '', '', '', 0, []
  for i, item in enumerate(symbols):
    bins.append(item[0])
    elements += f' {item[0]}'
    masses += f'mass\t\t{i+1} {atomic_masses[atomic_numbers[item[0]]]}\n'
  for i, item in enumerate(all_symbols):
    if not item in bins:
      count    += 1
      elements += f' {item}'
      masses   += f'mass\t\t{count+len(symbols)} {atomic_masses[atomic_numbers[item]]}\n'
  lines += f'pair_style  {pot_type}\npair_coeff  * * "{pot}" {elements}\n'
  lines += masses
  return lines
      
###
if __name__ == "__main__":
    main()
