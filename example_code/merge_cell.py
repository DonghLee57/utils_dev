# Python3.x
# Reference - J. Phys.: Condens. Matter 29 (2017) 185901
# Algorithm
# 1. Construction of the surface supercells
# 2. Rotational alignment
# 3. Supercell matching
# 4. Acceptance of the interface supercell

import sys
import numpy as np
import ase.io
import ase.build

NMAX, MMAX = 6,6
[str_bot,str_top] = ase.io.read(sys.argv[1:3], format='vasp')
[bot_hkl, top_hkl] = [sys.argv[3], sys.argv[4]]

# 1. Construction of the surface supercells
