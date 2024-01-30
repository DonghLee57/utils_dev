import sys, os, glob
import numpy as np
import itertools
from ase.io import read, write
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d


def main():
    my = StructureAnalysis('POSCAR', 'vasp')
    #my = StructureAnalysis('XDATCAR', 'vasp-xdatcar',index=slice(0,1000,10), compress=True)
    #my = StructureAnalysis('my.lammps', 'lammps-data')
    #my = StructureAnalysis('dump.lammpstrj', 'lammps-dump-text', index=slice(0,1000,10), compress=True)
    
    rdf  = my.calculate_rdf(6)
    np.savetxt('rdf.out', rdf, fmt='%.4f')

    prdf = my.calculate_prdf('Si','Si',6)
    np.savetxt('prdf.out', prdf, fmt='%.4f')

    
    adf  = my.calculate_adf()
    return 0

###
class StructureAnalysis:
    def __init__(self, filename, fileformat, index='-1', compress=None):
        self.structure = read(filename, index=index, format=fileformat)
        if compress: write('INPUT_STR.vasp',images=self.structure, format='vasp-xdatcar')

    def calculate_rdf(self, rmax, dr=0.02):
        bins = np.arange(dr/2, rmax+dr/2, dr)
        rdf = np.zeros(len(bins)-1)
        if ( str(type(self.structure[0])) == "<class 'ase.atom.Atom'>" ):
            atoms = self.structure.copy()
            if rmax > atoms.get_cell().diagonal().min() / 2:
                print('WARNING: The input maximum radius is over the half the smallest cell dimension.')
            nions = atoms.get_global_number_of_atoms()
            dist = np.zeros((nions, nions))
            for i in range(nions):
                dist[i] = atoms.get_distances(i, range(nions), mic=True)
            res, bin_edges = np.histogram(dist, bins=bins)
            rdf += res / ( (nions**2 / atoms.get_volume()) * 4 * np.pi * dr * bin_edges[:-1]**2 )
        elif ( str(type(self.structure[0])) == "<class 'ase.atoms.Atoms'>" ):
            nimg = len(self.structure)
            for n, atoms in enumerate(self.structure):
                if ( n == 0 ) and (rmax > atoms.get_cell().diagonal().min() / 2) :
                    print('WARNING: The input maximum radius is over the half the smallest cell dimension.')
                nions = atoms.get_global_number_of_atoms()
                dist = np.zeros((nions, nions))
                for i in range(nions):
                    dist[i] = atoms.get_distances(i, range(nions), mic=True)
                res, bin_edges = np.histogram(dist, bins=bins)
                rdf += res / ( (nions**2 / atoms.get_volume()) * 4 * np.pi * dr * bin_edges[:-1]**2 )
            rdf /= nimg
        return [bin_edges[:-1], rdf]

    def calculate_prdf(self, elemA, elemB, rmax,  dr=0.02):
        bins = np.arange(dr/2, rmax+dr/2, dr)
        rdf = np.zeros(len(bins)-1)
        if ( str(type(self.structure[0])) == "<class 'ase.atom.Atom'>" ):
            atoms = self.structure.copy()
            if rmax > atoms.get_cell().diagonal().min() / 2:
                print('WARNING: The input maximum radius is over the half the smallest cell dimension.')
            idA = np.where( np.array(atoms.get_chemical_symbols()) == elemA )[0]
            nelemA = len(idA)
            idB = np.where( np.array(atoms.get_chemical_symbols()) == elemB )[0]
            nelemB = len(idB)
            dist = np.zeros((nelemA, nelemB))
            for i, idx in enumerate(idA):
                dist[i] = atoms.get_distances(idx, idB, mic=True)
            res, bin_edges = np.histogram(dist, bins=bins)
            rdf += res / (nelemA * nelemB / atoms.get_volume() * 4 * np.pi * dr * bin_edges[:-1]**2)
        elif ( str(type(self.structure[0])) == "<class 'ase.atoms.Atoms'>" ):
            nimg = len(self.structure)
            for n, atoms in enumerate(self.structure):
                if ( n == 0 ) and (rmax > atoms.get_cell().diagonal().min() / 2) :
                    print('WARNING: The input maximum radius is over the half the smallest cell dimension.')
                idA = np.where( np.array(atoms.get_chemical_symbols()) == elemA )[0]
                nelemA = len(idA)
                idB = np.where( np.array(atoms.get_chemical_symbols()) == elemB )[0]
                nelemB = len(idB)
                dist = np.zeros((nelemA, nelemB))
                for i, idx in enumerate(idA):
                    dist[i] = atoms.get_distances(idx, idB, mic=True)
                res, bin_edges = np.histogram(dist, bins=bins)
                rdf += res / (nelemA * nelemB / atoms.get_volume() * 4 * np.pi * dr * bin_edges[:-1]**2)
            rdf /= nimg
        return [bin_edges[:-1], rdf]
    
    def calculate_adf(self, bins=100, r_max=None):
        nimg = len(self.structure)
        rdf = np.zeros(bins)
        if r_max is None:
            r_max = atoms.get_cell().diagonal().min() / 2  # Half the smallest cell dimension
        bins = np.arange(dr/2, r_max, dr)
        angles = []
        for i in range(len(atoms)):
            neighbors_i = atoms.get_neighbors(i, r_max)[0]
            for j in neighbors_i:
                if j.index <= i:  # Avoid double counting
                    continue
                neighbors_j = atoms.get_neighbors(j.index, r_max)[0]
                for k in neighbors_j:
                    if k.index <= j.index:
                        continue
                    # Calculate angle i-j-k
                    vec_ij = atoms.positions[j.index] - atoms.positions[i]
                    vec_jk = atoms.positions[k.index] - atoms.positions[j.index]
                    angle = np.arccos(np.dot(vec_ij, vec_jk) / (np.linalg.norm(vec_ij) * np.linalg.norm(vec_jk)))
                    angle = np.degrees(angle)  # Convert to degrees
                    angles.append(angle)
        # Create a histogram of angles
        adf, bin_edges = np.histogram(angles, bins=bins, range=(0, 180), density=True)
        return bin_edges[1:], adf     

###
if __name__ == "__main__":
    main()
