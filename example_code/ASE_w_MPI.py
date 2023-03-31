import ase
from mpi4py import MPI
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size= comm.Get_size()

test_in  = 'POSCAR_in'
test_out = 'POSCAR_out'

# !mpirun -np $NN python (*.py)
#
# Case 1 - Automatic broadcast in ASE
#
structure = ase.io.read(test_in,format='vasp')
ase.io.write(test_out,images=structure,format='vasp')


# Case 2 - Manual parallelization
# parallel = True or False
if rank==0:
  structure = ase.io.read(test_in,format='vasp',parallel=False)
  ase.io.write(test_out,images=structure,format='vasp',parallel=False)
else:
  pass
