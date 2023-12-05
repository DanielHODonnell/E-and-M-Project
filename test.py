import numpy as np
import matplotlib.pyplot as plt


q = 1.0
sphere_radius = 2.0
sphere_center = (0.0, 0.0)


slice_radius = 4.0
resolution = 20


x, y = np.meshgrid(np.linspace(-slice_radius, slice_radius, resolution), np.linspace(-slice_radius, slice_radius, resolution))


r = np.sqrt((x - sphere_center[0]) ** 2 + (y - sphere_center[1]) ** 2)
E = q / (4 * np.pi * r**2)


Ex = E * (x - sphere_center[0]) / r
Ey = E * (y - sphere_center[1]) / r


plt.quiver(x, y, Ex, Ey, scale=50, pivot='middle', angles='xy', color='b')


plt.scatter(sphere_center[0], sphere_center[1], color='red', label='Charge')
circle = plt.Circle(sphere_center, sphere_radius, color='gray', fill=False, label='Sphere')
plt.gca().add_patch(circle)


plt.xlim(-slice_radius, slice_radius)
plt.ylim(-slice_radius, slice_radius)
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Electric Field Vectors in a 2D Slice')
plt.legend()
plt.grid(True)


plt.savefig('electric_field_vectors_2d.png')
