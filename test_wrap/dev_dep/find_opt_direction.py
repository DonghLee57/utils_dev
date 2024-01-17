import numpy as np
from ase import Atoms
from ase.io import read, write
import matplotlib.pyplot as plt

def fibonacci_sphere(samples:int=1, radius:float=1):
    """
    Generates points on a sphere using Fibonacci grid.

    :param samples: Number of points to generate.
    :param radius: Radius of the sphere.
    :return: Array of points on the sphere.
    """
    points = np.zeros((samples, 3))
    phi = np.pi * (3. - np.sqrt(5.))  # Golden angle

    for i in range(samples):
        y =  1 - (i / float(samples - 1)) * 2  # y goes from 1 to -1
        R = np.sqrt( 1 - y**2)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * R
        z = np.sin(theta) * R

        points[i] = np.array([x, y, z])*radius
    return np.array(points)


def find_optimal_passivation_direction(IDX, Obj, radius=1.0, samples=1, min_distance=0.5):
    """
    Find the optimal direction for passivation on an atom.
    :param atom: Atom to passivate.
    :param atoms: ASE Atoms object representing the structure.
    :param radius: Radius for the Fibonacci grid.
    :param samples: Number of points to sample on the sphere.
    :param min_distance: Minimum allowed distance from other atoms.
    :return: Optimal direction vector for passivation, None if no valid direction found.
    """
    fibonacci_points = fibonacci_sphere(samples, radius)
    optimal_direction = np.zeros(3)

    for ii, point in enumerate(fibonacci_points):
        test_position = Obj.positions[IDX] + point
        test_atom = Atoms(['X'], positions=[test_position])
        test_structure = Obj + test_atom

        distances = test_structure.get_distances(-1, range(test_structure.get_global_number_of_atoms()-1), mic=True)
        if len(np.where( distances < min_distance )[0]) == 0:
            optimal_direction += test_structure.get_distance(IDX, -1, mic=True, vector=True)
            
    return optimal_direction/np.linalg.norm(optimal_direction)

tmp = read('POSCAR_generated.vasp')
ID = np.argmax(tmp.positions.T[2])

X = find_optimal_passivation_direction(ID, tmp, 2, 30, 1.5)
tmp.append('X')
tmp.positions[-1] = tmp.positions[ID] + X

write('test.vasp',images=tmp,format='vasp')
