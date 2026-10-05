import numpy             as np
import matplotlib.pyplot as plt

ld_dns = np.genfromtxt('dns_j2/lodi.out', delimiter=',', dtype=np.float64)

ld_fe2 = np.genfromtxt('fe2_j2/lodi.out', delimiter=',', dtype=np.float64)

fig, ax = plt.subplots()

p1 = ax.plot(ld_dns[:,1], ld_dns[:,2], 'or' , label='dns')
p2 = ax.plot(ld_fe2[:,1], ld_fe2[:,2], '--b', label='fe2')

ax.legend()
ax.grid()

ax.set_xlabel('Disp. [mm]')
ax.set_ylabel('Force [kN]')

plt.show()

