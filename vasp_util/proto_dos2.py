# POSCAR & DOSCAR & OUTCAR
# PROCAR for IPR (should be checked)

import sys, subprocess
import numpy as np
from ase.io import read
import matplotlib.pyplot as plt
import re

# ================================================================
# User settings (can be moved to command line arguments or config)
# ================================================================
energy_rng = [-3, 3]
plot_pdos = 1
plot_ldos = 0
plot_ipr = 0
norbit = 9 # s, py, pz, px, dxy, dyz, dz2, dxz, x2-y2
fs = 12

def main():
    symbols, ntypes = get_poscar_header('POSCAR')
    my_dos = DOS(doscar='DOSCAR', norbit=norbit)
    my_dos.ISPIN = int(subprocess.check_output('grep ISPIN OUTCAR'.split(),universal_newlines=True).split()[2])
    energy = np.linspace(my_dos.EMIN, my_dos.EMAX, my_dos.NEDOS)

    if plot_ipr:
        my_ipr = IPR(procar='PROCAR')

    # Get reference potential from OUTCAR
    reference_potential = get_reference_potential('OUTCAR', 'POSCAR')
    print("Reference potential:", reference_potential)

    # Plotting
    fig, ax = plt.subplots()

    # Local DOS
    if plot_ldos:
        poscar = read('POSCAR')
        z_coords = poscar.positions.T[2]
        site_dos = np.sum(my_dos.DOS[:, :, :], axis=2).T
        z_grid, e_grid = np.meshgrid(z_coords, energy)
        ldos_norm = site_dos / np.max(site_dos)
        cmesh = ax.pcolormesh(z_grid, e_grid - my_dos.fermi,
                              ldos_norm, cmap='rainbow', shading='auto')
        fig.colorbar(cmesh, label="LDOS")
        ax.set_xlabel("Direction (Å)", fontsize=fs)
        ax.set_ylabel("Energy (eV)", fontsize=fs)
        ax.set_title("Local Density of States", fontsize=fs)

    # Partial DOS
    elif plot_pdos:
        ATOMS = {}
        for sym, count in zip(symbols, ntypes):
            idx = sum(ntypes[:ntypes.index(count)])
            if sym not in ATOMS.keys():
                ATOMS[sym] = list(range(idx, idx+count))
            else:
                ATOMS[sym] += list(range(idx, idx+count))
        if my_dos.ISPIN == 2:
            up = np.arange(1,norbit*my_dos.ISPIN+1,2)
            down = np.arange(2,norbit*my_dos.ISPIN+1,2)
            for sid, sym in enumerate(symbols):
                pdos = np.sum(np.sum(my_dos.DOS[ATOMS[sym]].T[up],axis=0),axis=1)
                ax.plot(energy-my_dos.fermi, pdos)
                pdos = np.sum(np.sum(my_dos.DOS[ATOMS[sym]].T[down],axis=0),axis=1)
                ax.plot(energy-my_dos.fermi, -pdos,c=f'C{sid}', label=sym)
        else:
            up = np.arange(1,norbit*my_dos.ISPIN+1)
            for sid, sym in enumerate(symbols):
                pdos = np.sum(np.sum(my_dos.DOS[ATOMS[sym]].T[up],axis=0),axis=1)
                ax.plot(energy-my_dos.fermi, pdos, c=f'C{sid}', label='Total DOS')
            ax.set_ylim(bottom=0)

    # Total DOS
    else:
        ax.plot(energy-my_dos.fermi, my_dos.TDOS[:,1],'k')
        if my_dos.ISPIN == 2:
            ax.plot(energy-my_dos.fermi, -my_dos.TDOS[:,2],'k')
        else:
            ax.set_ylim(bottom=0)

    # Plot IPR
    if plot_ipr:
        ax_r = ax.twinx()
        stem_container = ax_r.stem(my_ipr.energies - my_dos.fermi, my_ipr.ipr_values, markerfmt=' ')
        stem_container.stemlines.set_alpha(0.2)
        ax_r.set_ylabel('IPR', fontsize=fs)
        ax_r.set_ylim(bottom=0)

    ax.axvline(0,c='gray', ls='--')
    ax.set_xlim(energy_rng)
    ax.set_xlabel(r'E-E$_f$ (eV)',fontsize=fs)
    ax.set_ylabel(r'DOS (a.u.)',fontsize=fs)
    ax.set_yticks([])
    ax.legend(fontsize=fs)
    plt.tight_layout()
    plt.show()

# ================================================================
# Core classes & functions
# ================================================================
class DOS:
    def __init__(self, doscar='DOSCAR', norbit=9):
        self.NIONS = 0
        self.EMAX, self.EMIN, self.NEDOS, self.fermi = 0, 0, 0, 0
        self.ISPIN = 1
        self.NORBIT= norbit
        self.DOS = []
        self.TDOS= []
        self.read(doscar)

    def read(self, DOSCAR):
        with open(DOSCAR,'r') as o: tmp = enumerate(o.readlines())
        self.NIONS = int(next(tmp)[1].split()[0])
        for i in range(4): next(tmp)
        self.EMAX, self.EMIN, self.NEDOS, self.fermi, _ = map(float,next(tmp)[1].split())
        self.NEDOS = int(self.NEDOS)
        for e in range(self.NEDOS):
            _, dat = next(tmp)
            if len(self.TDOS) == 0:
                if len(dat.split()) <=3: self.ISPIN = 1
                else: self.ISPIN = 2
                self.TDOS = np.zeros((self.NEDOS, self.ISPIN*2+1))
                self.DOS = np.zeros((self.NIONS, self.NEDOS, self.ISPIN*self.NORBIT+1))
            else:
                self.TDOS[e] += np.array(list(map(float,dat.split())))
        n = 0
        for idx, line in tmp:
            for e in range(self.NEDOS):
                _, dat = next(tmp)
                self.DOS[n][e] += np.array(list(map(float,dat.split())))
            n+=1
        return 0

class IPR:
    def __init__(self, procar='PROCAR'):
        self.kpoints = 0
        self.bands = 0
        self.ions = 0
        self.energies = []
        self.ipr_values = []
        self.read(procar)

    def read(self, procar):
        with open(procar, 'r') as f:
            lines = f.readlines()
        # Read header information
        header = lines[1].split()
        self.kpoints = int(header[3])
        self.bands = int(header[7])
        self.ions = int(header[-1])
        data = []
        for i in range(len(lines)):
            if 'band ' in lines[i]:
                sum_sq = 0
                line = lines[i].split()
                energy = float(line[4])
                # Get normalization
                norm_line = lines[i+3+self.ions].split()
                norm = float(norm_line[-1])
                # Calculate IPR
                for j in range(self.ions):
                    line = lines[i+j+3].split()
                    proj_value = float(line[-1])
                    sum_sq += proj_value**2
                # Normalize IPR
                if norm > 0:
                    ipr = sum_sq/norm**2
                else:
                    ipr = 0
                data.append([energy, ipr])
        # Sort by energy and separate into arrays
        data.sort(key=lambda x: x[0])
        self.energies = np.array([x[0] for x in data])
        self.ipr_values = np.array([x[1] for x in data])

def get_poscar_header(poscar):
    try:
        with open(poscar, 'r') as o:
            tmp = o.readlines()
        symbols = tmp[5].split()
        ntypes = list(map(int, tmp[6].split()))
        return symbols, ntypes
    except IndexError:
        print('Please ensure you provide a POSCAR file in the correct format.')
        print('Command >> python (this.py) (YOUR_POSCAR)\n')
        print('##### Example of a POSCAR file format #####')
        print('Comment')
        print('Scaling factor')
        print(' Lattice[0][0] Lattice[0][1] Lattice[0][2]')
        print(' Lattice[1][0] Lattice[1][1] Lattice[1][2]')
        print(' Lattice[2][0] Lattice[2][1] Lattice[2][2]')
        print('Speicies names <-- This line is required.')
        print('Ions per species <-- This line is required.')
        print('Selective dynamics <-- optional')
        print(' Ion positions...')


def get_reference_potential(ref_element, outcar='OUTCAR', poscar='POSCAR'):
    """
    Extracts the reference potential from the 'average (electrostatic) potential at core' section of OUTCAR.
    Uses grep to extract the relevant lines based on the total number of atoms from POSCAR,
    considering data is displayed in 5 columns per line.
    Parses the extracted section to obtain the potential values.
    Returns the first potential value (can be changed to average or other reference as needed).

    Args:
        outcar (str): Path to the OUTCAR file.
        poscar (str): Path to the POSCAR file.

    Returns:
        float: The reference potential value.

    Raises:
        ValueError: If no potential values are found.
    """
    # Get total number of atoms from POSCAR
    _, ntypes = get_poscar_header(poscar)
    total_atoms = sum(ntypes)
    # Calculate number of data lines needed: ceil(total_atoms / 5)
    lines_needed = (total_atoms + 4) // 5
    # Add a few extra lines for header and blank lines, but ensure at least lines_needed
    lines_to_extract = lines_needed + 2

    # Use grep to extract the section
    command = f"grep -A {lines_to_extract} 'average (electrostatic) potential at core' {outcar}"
    outcar_snippet = subprocess.check_output(command, shell=True, universal_newlines=True)

    # Parse the extracted section for potential values
    lines = outcar_snippet.split('\n')
    potentials = []
    start_parsing = False

    for line in lines:
        if 'average (electrostatic) potential at core' in line:
            start_parsing = True
            continue
        if start_parsing:
            if not line.strip() or 'test charge' in line:
                continue
            pairs = re.findall(r'(\d+)\s+(-?\d+\.\d+)', line)
            for _, potential in pairs:
                potentials.append(float(potential))
            if len(pairs) == 0:
                break
    potentials = np.array(potentials)
    selected = np.where( np.array(read(poscar).get_chemical_symbols()) == ref_element )[0]
    if len(potentials) > 0:
        return np.mean(potentials[selected])
    else:
        raise ValueError("No potential values found in OUTCAR section.")

#--------------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
