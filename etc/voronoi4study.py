import numpy as np
from scipy.spatial import Voronoi
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Line3DCollection

def plot_voronoi_edges_3d(vor, highlight_point_index):
    """
    Plots the edges of a 3D Voronoi diagram and highlights the edges belonging to a selected point.
    
    Parameters:
    - vor: Voronoi object.
    - highlight_point_index: Index of the Voronoi point to highlight the planes for.
    """
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    # Regular edges
    edge_points = []
    # Highlighted edges for the selected point
    highlighted_edge_points = []

    for ridge_idx, ridge in enumerate(vor.ridge_points):
        if -1 not in vor.ridge_vertices[ridge_idx]:  # Ignore ridges that extend to infinity
            p0, p1 = ridge
            if p0 == highlight_point_index or p1 == highlight_point_index:
                # Ridge belongs to the selected point
                highlighted_edge_points.append(vor.vertices[vor.ridge_vertices[ridge_idx]])
            else:
                # Regular ridge
                edge_points.append(vor.vertices[vor.ridge_vertices[ridge_idx]])

    # Regular edges
    if edge_points:
        edge_collection = Line3DCollection(edge_points, colors='k', linewidths=0.2)
        ax.add_collection3d(edge_collection)

    # Highlighted edges
    if highlighted_edge_points:
        highlighted_edge_collection = Line3DCollection(highlighted_edge_points, colors='r', linewidths=1.0)
        ax.add_collection3d(highlighted_edge_collection)

    # Plotting the Voronoi points
    ax.scatter(vor.points[:,0], vor.points[:,1], vor.points[:,2], s=5, c='b')
    # Highlight the selected point
    ax.scatter(vor.points[highlight_point_index,0], vor.points[highlight_point_index,1], vor.points[highlight_point_index,2], s=50, c='r')

    ###
    test_points = np.array([[5, 5, 5], [2,2,2], [8, 8, 8], vor.points[0]+np.array([0.2,0.2,0.2])])  # Define multiple test points
    ax.scatter(test_points.T[0], test_points.T[1], test_points.T[2], s=50, c='g')

    # Setting the axes limits to frame all points and edges
    min_corner = vor.vertices.min(axis=0) - 1.0
    max_corner = vor.vertices.max(axis=0) + 1.0
    ax.set_xlim(min_corner[0], max_corner[0])
    ax.set_ylim(min_corner[1], max_corner[1])
    ax.set_zlim(min_corner[2], max_corner[2])

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.title('3D Voronoi Diagram with Highlighted Planes')
    plt.show()

def create_periodic_voronoi(points, cell_size):
    """
    Generate a Voronoi diagram with periodic boundary conditions.

    Parameters:
    - points: An array of points within the original cell.
    - cell_size: The size of the periodic cell along each axis (assuming a cubic cell for simplicity).

    Returns:
    - vor: A Voronoi object containing the vertices and regions of the diagram respecting PBCs.
    """
    # Duplicate points in all directions to simulate periodicity
    extended_points = []
    for dx in [-cell_size, 0, cell_size]:
        for dy in [-cell_size, 0, cell_size]:
            for dz in [-cell_size, 0, cell_size]:
                if dx == dy == dz == 0:
                    # Skip the original cell in the center
                    continue
                displacement = np.array([dx, dy, dz])
                extended_points.extend(points + displacement)

    # Include original points
    extended_points = np.array(extended_points)
    extended_points = np.vstack((points, extended_points))

    # Compute the Voronoi diagram for the extended points
    return extended_points

def are_points_in_voronoi_cell(vor, cell_index, test_points):
    """
    Determines if given points are within the Voronoi cell of a specified generating point, efficiently for multiple points.
    
    Parameters:
    - vor: A Voronoi object from scipy.spatial.
    - cell_index: The index of the generating point for the Voronoi cell of interest.
    - test_points: An array of points to test, each as a numpy array [x, y, z].
    
    Returns:
    - is_within: A boolean array where True indicates that the corresponding test_point is within the Voronoi cell of the generating point at cell_index.
    """
    # Compute squared distances from each test point to each generating point to avoid sqrt for efficiency
    distances_squared = np.sum((vor.points[:, np.newaxis, :] - test_points) ** 2, axis=2)
    
    # Find the indices of the closest generating points to each test point
    closest_indices = np.argmin(distances_squared, axis=0)
    
    # Check if the closest generating points are the one of interest
    is_within = closest_indices == cell_index
    
    return is_within

# Continuing from the provided example
np.random.seed(42)
points = np.random.rand(10, 3) * 10  # Adjusted to generate more points for better visualization

epts = create_periodic_voronoi(points, 10)
vor = Voronoi(epts)

# Choose a point index to highlight. This should be less than len(points)
highlight_point_index = 0  # For example, highlighting the first point and its planes


# Example usage with multiple test points
test_points = np.array([[5, 5, 5], [2,2,2], [8, 8, 8], vor.points[0]+np.array([0.2,0.2,0.2])])  # Define multiple test points
cell_index = 0  # Index of the Voronoi cell/generating point to check against

is_within_array = are_points_in_voronoi_cell(vor, cell_index, test_points)
for i, is_within in enumerate(is_within_array):
    print(f"Is test point {test_points[i]} within the Voronoi cell of point index {cell_index}? {is_within}")

plot_voronoi_edges_3d(vor, highlight_point_index)
