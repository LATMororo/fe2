import numpy             as np
import matplotlib.pyplot as plt

ld  = np.genfromtxt('lodi.out', delimiter=',', dtype=np.float64)
ld2 = np.genfromtxt('lodi_fe2.out', delimiter=',', dtype=np.float64)

fig, ax = plt.subplots()

p1 = ax.plot(ld[:,1], ld[:,2], 'or' , label='prnn')
p2 = ax.plot(ld[:,1], ld[:,2], 'b' , label='fe2')

ax.legend()
ax.grid()

ax.set_xlabel('Disp. [mm]')
ax.set_ylabel('Force [kN]')

plt.show()

