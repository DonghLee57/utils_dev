#import cProfile
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

    prdf = my.calculate_prdf(['Si','Si'],6)
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

    def calculate_prdf(self, targets, rmax,  dr=0.02):
        # targets: ["Element 1", "Element 2"]
        # ex) ["Si", "O"]
        [elemA, elemB] = targets
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
    
    def calculate_adf(self, targets, cutoff, angle_bins=np.arange(0,180.1,2), expr='degree'):
        # To facilitate a clear bond angle analysis, the target triplet with cutoff distances should be provided.
        # targets: ["Center", "Neighbor_1", "Neighbor_2"]
        # ex) ["Si", "Si", "O"]
        # cutoff: [r_1, r_2]
        # ex) [2.6, 2.0]
        if expr != 'degree': angle_bins = angle_bins*np.pi/180
        if ( str(type(self.structure[0])) == "<class 'ase.atom.Atom'>" ):
            atoms = self.structure.copy()
            symbols = np.array(atoms.get_chemical_symbols())
            cidx = np.where( targets[0] == symbols )[0]
            nidx = np.where( targets[1] == symbols )[0]
            midx = np.where( targets[2] == symbols )[0]
            theta = []
            psize = len(cidx)//size
            if targets[1] == targets[2]:
                if rank == size-1:
                    for c in range(rank*psize,len(cidx)):
                        for n in range(len(nidx)):
                            vec1 = atoms.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                            dist1 = np.linalg.norm(vec1)
                            if dist1 < cutoff[0]:
                                for m in range(len(midx)):
                                    if m > n:
                                        vec2 = atoms.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                        dist2 = np.linalg.norm(vec2)
                                        if dist2 < cutoff[1]:
                                            theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
                else:
                    for c in range(rank*psize,(rank+1)*psize):
                        for n in range(len(nidx)):
                            vec1 = atoms.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                            dist1 = np.linalg.norm(vec1)
                            if dist1 < cutoff[0]:
                                for m in range(len(midx)):
                                    if m > n:
                                        vec2 = atoms.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                        dist2 = np.linalg.norm(vec2)
                                        if dist2 < cutoff[1]:
                                            theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
            else:
                if rank == size-1:
                    for c in range(rank*psize,len(cidx)):
                        for n in range(len(nidx)):
                            vec1 = atoms.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                            dist1 = np.linalg.norm(vec1)
                            if dist1 < cutoff[0]:
                                for m in range(len(midx)):
                                    if m != n:
                                        vec2 = atoms.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                        dist2 = np.linalg.norm(vec2)
                                        if dist2 < cutoff[1]:
                                            theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
                else:
                    for c in range(rank*psize,(rank+1)*psize):
                        for n in range(len(nidx)):
                            vec1 = atoms.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                            dist1 = np.linalg.norm(vec1)
                            if dist1 < cutoff[0]:
                                for m in range(len(midx)):
                                    if m != n:
                                        vec2 = atoms.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                        dist2 = np.linalg.norm(vec2)
                                        if dist2 < cutoff[1]:
                                            theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
        elif ( str(type(self.structure[0])) == "<class 'ase.atoms.Atoms'>" ):
            nimg = len(self.structure)
            for n, atoms in enumerate(self.structure):
                symbols = np.array(atoms.get_chemical_symbols())
                cidx = np.where( targets[0] == symbols )[0]
                nidx = np.where( targets[1] == symbols )[0]
                midx = np.where( targets[2] == symbols )[0]
                theta = []
                psize = len(cidx)//size
                if targets[1] == targets[2]:
                    if rank == size-1:
                        for c in range(rank*psize,len(cidx)):
                            for n in range(len(nidx)):
                                vec1 = atoms.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                                dist1 = np.linalg.norm(vec1)
                                if dist1 < cutoff[0]:
                                    for m in range(len(midx)):
                                        if m > n:
                                            vec2 = atoms.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                            dist2 = np.linalg.norm(vec2)
                                            if dist2 < cutoff[1]:
                                                theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
                    else:
                        for c in range(rank*psize,(rank+1)*psize):
                            for n in range(len(nidx)):
                                vec1 = atoms.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                                dist1 = np.linalg.norm(vec1)
                                if dist1 < cutoff[0]:
                                    for m in range(len(midx)):
                                        if m > n:
                                            vec2 = atoms.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                            dist2 = np.linalg.norm(vec2)
                                            if dist2 < cutoff[1]:
                                                theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
                else:
                    if rank == size-1:
                        for c in range(rank*psize,len(cidx)):
                            for n in range(len(nidx)):
                                vec1 = atoms.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                                dist1 = np.linalg.norm(vec1)
                                if dist1 < cutoff[0]:
                                    for m in range(len(midx)):
                                        if m != n:
                                            vec2 = atoms.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                            dist2 = np.linalg.norm(vec2)
                                            if dist2 < cutoff[1]:
                                                theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
                    else:
                        for c in range(rank*psize,(rank+1)*psize):
                            for n in range(len(nidx)):
                                vec1 = atoms.get_distances(cidx[c], nidx[n], mic=True, vector=True)
                                dist1 = np.linalg.norm(vec1)
                                if dist1 < cutoff[0]:
                                    for m in range(len(midx)):
                                        if m != n:
                                            vec2 = atoms.get_distances(cidx[c], midx[m], mic=True, vector=True)
                                            dist2 = np.linalg.norm(vec2)
                                            if dist2 < cutoff[1]:
                                                theta.append(np.arccos(np.round(np.dot(vec1, vec2.T)[0]/dist1/dist2,6)))
        theta = comm.reduce(theta,op=MPI.SUM,root=0)
        theta = comm.bcast(theta, root=0)
        if expr == 'degree': theta = np.array(theta)*180/np.pi
        res, bin_edges = np.histogram(theta, bins=angle_bins, density=True)
        return [bin_edges[:-1], res]
        
    def unwrapping(self, selected_atoms=None, selected_elem=None):
        if ( str(type(self.structure[0])) == "<class 'ase.atoms.Atoms'>" ):
            atoms = self.structures.copy()
            nimg = len(self.structure)
            unwrap_traj = np.zeros((nimg, len(targets), 3))
            if selected_elem !=  None:
                targets = np.where( selected_elem == np.array(self.structure[0].get_chemical_symbols()) )[0]
            if selected_atoms != None:
                targets = selected_atoms
            img0 = atoms.get_positions()[targets]
            unwrap_traj[0] = img0.copy()
            lattice = img0.cell.copy()
            old = img0.get_scaled_positions()
            psize = len(targets)//size
            if rank == size-1:
                parts = np.arange(rank*psize,len(targets))
                for idx in range(1, nimg):
                    new = atoms[idx].get_scaled_positions()
                    diff = new[targets[parts]] - old[targets[parts]]
                    unwrap_traj[idx][parts] = unwrap_traj[idx-1][parts] + (np.round(-diff)+diff)@lattice
                    old = atoms[idx].get_scaled_positions()
            else:
                parts = np.arange(rank*psize,(rank+1)*psize)
                for idx in range(1, nimg):
                    new = atoms[idx].get_scaled_positions()
                    diff = new[targets[parts]] - old[targets[parts]]
                    unwrap_traj[idx][parts] = unwrap_traj[idx-1][parts] + (np.round(-diff)+diff)@lattice
                    old = atoms[idx].get_scaled_positions()
            unwrap_traj = comm.reduce(unwrap_traj, op=MPI.SUM, root=0)
            unwrap_traj = comm.bcast(unwrap_traj, root=0)
            return unwrap_traj
        else:
            print("WARNING: Cannot unwrap atomic coordinates in a single structure.")
            return None

###
if __name__ == "__main__":
    main()
