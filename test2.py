import numpy as np
import matplotlib.pyplot as plt


p = 1.0
R = 1.0


theta = np.linspace(0, np.pi, 100)


sigma = (p / (2 * np.pi * R)) * np.cos(theta)


plt.figure(figsize=(8, 6))
plt.plot(theta, sigma, label='Surface Charge Density')
plt.xlabel('Polar Angle (Theta)')
plt.ylabel('Surface Charge Density (sigma)')
plt.title('Surface Charge Density on a Conducting Sphere due to a Dipole')
plt.grid()
plt.legend()

plt.savefig('surface_charge_density.png')
