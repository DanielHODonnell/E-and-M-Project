import numpy as np

# Part A

"""
Define constants
"""
e = 1.602176634e-19
r = 0.09572e-9
theta = 104.5 * np.pi / 180

"""
Define the unit vectors pointing from the center of the molecule to the atoms
"""
r1 = np.array([0, 0, 0])
r2 = np.array([r * np.cos(theta/2), r * np.sin(theta/2), 0])
r3 = np.array([-r * np.cos(theta/2), -r * np.sin(theta/2), 0])

"""
Calculate the dipole moment
"""
p = e * (2 * r2 - r1 - r3)

"""
Calculate the quadrupole moment
"""
Q = np.zeros((3, 3))
Q[0, 0] = e/2 * (r1.dot(r1) - 1/3)
Q[1, 1] = e/2 * (r2.dot(r2) - 1/3)
Q[2, 2] = e/2 * (r3.dot(r3) - 1/3)
Q[0, 1] = e/2 * r1.dot(r2)
Q[1, 0] = Q[0, 1]
Q[1, 2] = e/2 * r2.dot(r3)
Q[2, 1] = Q[1, 2]

"""
Define the unit vector pointing from the center of the system to the point where the potential is being evaluated
"""
r_hat = np.array([1, 0, 0])

"""
Calculate the first three non-zero terms of the multipole expansion of the electrostatic potential
"""
result = e * r_hat.dot(p) / (4 * np.pi * r**3 * np.linalg.norm(r_hat)) \
      - e**2 * (3 * r_hat.dot(r1)**2 - 1) * (r_hat.dot(r2)**2 - 1) / (4 * np.pi * r**5 * np.linalg.norm(r_hat)**5) \
      - 6 * e**2 * r_hat.dot(r1) * r_hat.dot(r2) / (4 * np.pi * r**5 * np.linalg.norm(r_hat)**5)

print(f"Part A result: {result}")

# Part B & C

"""
Calculate the electric field at the geometric center of the water molecule using only the dipole term
&
Calculate the actual electric field at the geometric center of the water molecule
"""
E_dipole = p / (4 * np.pi * r**3)
E_actual = np.zeros(3)
E_actual += e * (r2 - r1) / (4 * np.pi * r**3)
E_actual += e * (r3 - r2) / (4 * np.pi * r**3)

"""
Calculating the error of the electric field
"""
error_E = np.linalg.norm(E_actual - E_dipole)

"""
Calculate the potential at a distance of 10 nm from the center of the water molecule in a direction perpendicular 
to the dipole moment using only the dipole term
"""
potential_dipole = p.dot(np.array([1, 0, 0])) / (4 * np.pi * 0.01**3)

"""
Calculate the actual potential at a distance of 10 nm from the center of the water molecule in a 
direction perpendicular to the dipole moment
"""
potential_actual = 0
for i in range(3):
    potential_actual += e * (r1[i] - r2[i]) / (4 * np.pi * 0.01**3)
    potential_actual += e * (r3[i] - r2[i]) / (4 * np.pi * 0.01**3)

"""
Calculate the error in the potential
"""
error_potential = potential_actual - potential_dipole

print(f"Error in electric field magnitude (Part B): {error_E}")
print(f"Error in potential (Part C): {error_potential}")
