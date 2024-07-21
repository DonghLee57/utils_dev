from mpi4py import MPI
import sys
import numpy as np
from ase.io import read, write
from ase.neighborlist import neighbor_list
import torch
from torch_geometric.data import Data
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
threshold=10

class UnionFind:
    def __init__(self, size):
        self.root = np.arange(size, dtype=int)
        self.rank = np.ones(size, dtype=int)

    def find(self, x):
        if self.root[x] != x:
            self.root[x] = self.find(self.root[x])  # 경로 압축
        return self.root[x]

    def union(self, x, y):
        rootX = self.find(x)
        rootY = self.find(y)

        if rootX != rootY:
            if self.rank[rootX] > self.rank[rootY]:
                self.root[rootY] = rootX
            elif self.rank[rootX] < self.rank[rootY]:
                self.root[rootX] = rootY
            else:
                self.root[rootY] = rootX
                self.rank[rootX] += 1

def process_graph(edges, num_nodes):
    uf = UnionFind(num_nodes)
    for x, y in edges:
        uf.union(x, y)
    return uf

def create_graph_object(filename):
    atoms = read(filename)
    i, j, distances = neighbor_list('ijd', atoms, cutoff=3.0)
    node_features = torch.tensor(atoms.numbers, dtype=torch.float)
    edge_index = torch.tensor([i, j], dtype=torch.long)
    graph = Data(x=node_features, edge_index=edge_index)
    return graph, atoms
  
def main():
    """
    # example
    from ase.build import fcc111
    slab = fcc111('Al', size=(4,4,3), vacuum=20.0)
    zmax = max(slab.positions.T[2])

    slab.append('H')
    slab.positions[-1] = np.array([0,0, zmax + 5])
    slab.append('H')
    slab.positions[-1] = np.array([1,0, zmax + 5])
    write('test.vasp',images=slab,format='vasp')
    """

    # 파일명 지정 및 그래프 객체 생성
    filename = sys.argv[1]

    """
    # ver.1
    graph, atoms = create_graph_object(filename)
    num_nodes = graph.num_nodes  
    edges = graph.edge_index.t().tolist()
    num_components = process_graph(edges, num_nodes)
    print("Number of connected components:", num_components)
    local_edges = np.array_split(edges, size)[rank]
    local_uf = process_graph(local_edges, num_nodes)
    #
    """
    # ver.2
    if rank == 0:
        graph, atoms = create_graph_object(filename)
        num_nodes = graph.num_nodes  
        edges = graph.edge_index.t().tolist()
        local_edges = np.array_split(edges, size)[rank]
    else:
        local_edges = None
        atoms = None
        num_nodes = None
    num_nodes = comm.bcast(num_nodes, root=0)
    local_edges = comm.scatter(local_edges, root=0)
    local_uf = process_graph(local_edges, num_nodes)
    #
  
    global_roots = np.empty(num_nodes, dtype=int)
    comm.Allreduce(local_uf.root, global_roots, op=MPI.MAX)
  
    if rank == 0:
        final_uf = UnionFind(num_nodes)
        final_uf.root = global_roots
        unique_components = len(set(final_uf.find(x) for x in range(num_nodes)))
        print("Number of connected components:", unique_components)
        nodes_to_remove = []
        for component in unique_components:
            component_size = len(component)
            print("Component size:", component_size)
            if component_size <= threshold:
                print("Small component atoms:", atoms[np.array(list(component))].symbols)
                nodes_to_remove.extend(component)

        if nodes_to_remove:
            atoms = atoms[[i for i in range(len(atoms)) if i not in nodes_to_remove]]
    write('after.vasp',images=atoms, format='vasp')

if __name__ == "__main__":
    main()
