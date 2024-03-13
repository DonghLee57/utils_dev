# Under testing
import numpy as np
from ase import Atoms
from ase.build import bulk, make_supercell
from scipy.spatial import Voronoi, voronoi_plot_2d
from scipy.spatial.transform import Rotation as R
import matplotlib.pyplot as plt

def create_polycrystal(seed_points, crystal_structure, cell_size, pad_width):
    """
    Generates a polycrystalline structure by placing randomly rotated and translated
    instances of a base crystal structure at Voronoi seed points.
    
    Parameters:
    - seed_points: Coordinates of the seed points for the Voronoi diagram, representing grain centers.
    - crystal_structure: ASE Atoms object of the base crystal structure for each grain.
    - cell_size: The size of the simulation cubic cell (assuming a cubic cell for simplicity).
    - pad_width: The padding width for trimming the final structure.
    
    Returns:
    - polycrystal: An ASE Atoms object representing the polycrystalline structure.
    """
    # Assuming periodic boundary conditions, extend seed points
    extended_points = np.concatenate([seed_points + disp for disp in np.diag([cell_size]*3)])
    vor = Voronoi(extended_points)

    polycrystal = None

    for point in seed_points:
        # Create a copy of the base crystal structure
        grain = crystal_structure.copy()

        # Apply random rotation
        rotation = R.random().as_matrix()
        grain.rotate(rotation, center=(0, 0, 0))

        # Apply translation to the Voronoi seed point
        grain.translate(point)

        # Merge the grain into the polycrystal
        if polycrystal is None:
            polycrystal = grain
        else:
            polycrystal += grain

    # Trim the structure to the desired size with padding
    # This is a simplified placeholder. In practice, you would check atoms' positions
    # and remove those outside the desired cubic cell minus the pad width.
    min_corner = np.array([0, 0, 0]) + pad_width
    max_corner = np.array([cell_size, cell_size, cell_size]) - pad_width
    trimmed_atoms = [atom for atom in polycrystal if np.all(atom.position >= min_corner) and np.all(atom.position <= max_corner)]
    polycrystal = Atoms(trimmed_atoms)

    return polycrystal

# Example parameters and usage
seed_points = np.random.rand(5, 3) * 10  # Generate random seed points within a 10x10x10 cube
base_crystal_structure = bulk('Si', 'diamond', a=5.43)  # Silicon crystal as an example
cell_size = 10
pad_width = 1

polycrystal = create_polycrystal(seed_points, base_crystal_structure, cell_size, pad_width)
print(f"Generated a polycrystalline structure with {len(polycrystal)} atoms.")
