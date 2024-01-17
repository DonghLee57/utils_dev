import numpy as np
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
        y =  radius * (1 - (i / float(samples - 1)) * 2)  # y goes from radius to -radius
        R = np.sqrt(radius**2 - y**2)  # radius at y

        theta = phi * i  # golden angle increment

        x = np.cos(theta) * R
        z = np.sin(theta) * R

        points[i]= [x, y, z]
    return np.array(points)

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

fs = fibonacci_sphere(50, 10)

ax.set_box_aspect((np.ptp(fs.T[0]), np.ptp(fs.T[1]), np.ptp(fs.T[2])))
ax.scatter(fs.T[0], fs.T[1], fs.T[2])

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.tight_layout()
plt.show()
