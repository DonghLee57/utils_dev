import sys, subprocess
import numpy as np
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

# ================================================================
# Core functions (unchanged)
# ================================================================
class DOS:
    # ... (keep the same as before)
    pass

class IPR:
    # ... (keep the same as before)
    pass

def get_poscar_header(poscar):
    # ... (keep the same as before)
    pass

def get_reference_potential(outcar='OUTCAR', poscar='POSCAR'):
    # ... (keep the same as before)
    pass

# ================================================================
# New: initialization and plotting functions
# ================================================================
def initialize_dos_pdos(path):
    """Initialize DOS and partial DOS data for a given path."""
    poscar = f"{path}/POSCAR"
    doscar = f"{path}/DOSCAR"
    outcar = f"{path}/OUTCAR"

    symbols, ntypes = get_poscar_header(poscar)
    my_dos = DOS(doscar=doscar, norbit=norbit)
    try:
        my_dos.ISPIN = int(subprocess.check_output(['grep', 'ISPIN', outcar], universal_newlines=True).split()[2])
    except:
        my_dos.ISPIN = 1
    energy = np.linspace(my_dos.EMIN, my_dos.EMAX, my_dos.NEDOS)

    # Get reference potential
    ref_pot = get_reference_potential(outcar, poscar)
    print(f"Reference potential for {path}: {ref_pot}")

    if plot_ipr:
        procar = f"{path}/PROCAR"
        my_ipr = IPR(procar=procar)
    else:
        my_ipr = None

    return symbols, ntypes, my_dos, energy, my_ipr

def plot_dos(ax, symbols, ntypes, my_dos, energy, label=None):
    """Plot DOS/PDOS/LDOS for a given set of data."""
    ATOMS = {}
    for sym, count in zip(symbols, ntypes):
        idx = sum(ntypes[:ntypes.index(count)])
        if sym not in ATOMS.keys():
            ATOMS[sym] = list(range(idx, idx+count))
        else:
            ATOMS[sym] += list(range(idx, idx+count))

    # Local DOS
    if plot_ldos:
        # NOTE: read('POSCAR') is not defined; you may need to implement or remove this part
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
        if my_dos.ISPIN == 2:
            up = np.arange(1, norbit*my_dos.ISPIN+1, 2)
            down = np.arange(2, norbit*my_dos.ISPIN+1, 2)
            for sid, sym in enumerate(symbols):
                pdos = np.sum(np.sum(my_dos.DOS[ATOMS[sym]].T[up], axis=0), axis=1)
                ax.plot(energy-my_dos.fermi, pdos, label=f"{label}_{sym}_up" if label else f"{sym}_up")
                pdos = np.sum(np.sum(my_dos.DOS[ATOMS[sym]].T[down], axis=0), axis=1)
                ax.plot(energy-my_dos.fermi, -pdos, label=f"{label}_{sym}_down" if label else f"{sym}_down")
        else:
            up = np.arange(1, norbit*my_dos.ISPIN+1)
            for sid, sym in enumerate(symbols):
                pdos = np.sum(np.sum(my_dos.DOS[ATOMS[sym]].T[up], axis=0), axis=1)
                ax.plot(energy-my_dos.fermi, pdos, label=f"{label}_{sym}" if label else sym)
            ax.set_ylim(bottom=0)

    # Total DOS
    else:
        ax.plot(energy-my_dos.fermi, my_dos.TDOS[:,1], label=label)
        if my_dos.ISPIN == 2:
            ax.plot(energy-my_dos.fermi, -my_dos.TDOS[:,2], label=f"{label}_down")
        else:
            ax.set_ylim(bottom=0)

# ================================================================
# Main function (revised for multiple paths)
# ================================================================
def main(paths=None):
    if paths is None:
        paths = ['.']  # default: current directory

    fig, ax = plt.subplots()

    for i, path in enumerate(paths):
        try:
            symbols, ntypes, my_dos, energy, my_ipr = initialize_dos_pdos(path)
            plot_dos(ax, symbols, ntypes, my_dos, energy, label=f"path{i+1}")
        except Exception as e:
            print(f"Error processing {path}: {e}")
            continue

    ax.axvline(0, c='gray', ls='--')
    ax.set_xlim(energy_rng)
    ax.set_xlabel(r'E-E$_f$ (eV)', fontsize=fs)
    ax.set_ylabel(r'DOS (a.u.)', fontsize=fs)
    ax.set_yticks([])
    ax.legend(fontsize=fs)
    plt.tight_layout()
    plt.show()

# ================================================================
# Example usage
# ================================================================
if __name__ == "__main__":
    # Example: python dos_ipr.py path1 path2 path3
    if len(sys.argv) > 1:
        paths = sys.argv[1:]
    else:
        paths = ['.']
    main(paths)
