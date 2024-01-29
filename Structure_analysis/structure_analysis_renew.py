import sys, os, glob
import numpy as np
import itertools
from ase.io import read
from mpi4py import MPI
comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d


def main():
    my = StructureAnalysis('POSCAR', 'vasp')
    #my = StructureAnalysis('XDATCAR', 'vasp-xdatcar',index=slice(0,1000,10))
    #my = StructureAnalysis('my.lammps', 'lammps-data')
    #my = StructureAnalysis('dump.lammpstrj', 'lammps-dump-text', index=slice(0,1000,10))
    
    rdf  = my.calculate_rdf()
    prdf = my.calculate_prdf()
    adf  = my.calculate_adf()
    return 0

###
class StructureAnalysis:
    def __init__(self, filename, fileformat, index='-1'):
        self.structure = read(filename, index=index, format=fileformat)

    def calculate_rdf(self, dr=0.02, r_max=None):
        nimg = len(self.structure)
        rdf = np.zeros(bins)
        for n, atoms in enumerate(self.structure):
            if r_max is None:
                r_max = atoms.get_cell().diagonal().min() / 2  # Half the smallest cell dimension
            bins = np.arange(dr/2, r_max, dr)
            distances = []
            for i in range(len(atoms)):
                for j in range(i + 1, len(atoms)):
                    # Compute distance between atoms i and j
                    distance = atoms.get_distance(i, j, mic=True)
                    if distance < r_max:
                        distances.append(distance)
            res, bin_edges = np.histogram(distances, bins=bins)
            rdf += res / ( (len(atoms)**2 / atoms.get_volume()) * 4 * np.pi * dr * bin_edges[:-1]**2 )
        return bin_edges[1:], rdf

    def calculate_prdf(self, element1, element2,  dr=0.02, r_max=None):
        nimg = len(self.structure)
        rdf = np.zeros(bins)
        if r_max is None:
            r_max = atoms.get_cell().diagonal().min() / 2  # Half the smallest cell dimension
        bins = np.arange(dr/2, r_max, dr)
        distances = []
        for n, atoms in enumerate(self.structure):
        # indentation
        for i in range(len(atoms)):
            if atoms[i].symbol != element1:
                continue
            for j in range(len(atoms)):
                if i != j and atoms[j].symbol == element2:
                    # Compute distance between atoms i and j
                    distance = atoms.get_distance(i, j, mic=True)
                    if distance < r_max:
                        distances.append(distance)
        res, bin_edges = np.histogram(distances, bins=bins)
        n_element1 = len([atom for atom in atoms if atom.symbol == element1])
        n_element2 = len([atom for atom in atoms if atom.symbol == element2])
        rdf += res / (n_element1 * n_element2 / atoms.get_volume() * 4 * np.pi * dr * bin_edges[:-1]**2)
        # indentation
        return bin_edges[1:], rdf
    

    def calculate_adf(self, bins=100, r_max=None):
        if r_max is None:
            r_max = atoms.get_cell().diagonal().min() / 2  # Half the smallest cell dimension
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
