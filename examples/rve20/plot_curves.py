
import numpy             as np
import matplotlib.pyplot as plt

# NOTE: experimental results given in N. Chawla, C. Andres, J.W. Jones. Cyclic
#       stress-strain behavior of particle reinforced metal material composites.
#       Script Materialia, vol. 38, no. 10, pp. 1595-1600.

exper = np.genfromtxt('experiment_vf20.dat', delimiter=',', dtype=np.float64)
ld1   = np.genfromtxt('lodi.out', delimiter=',', dtype=np.float64)

fig, ax = plt.subplots()

p1 = ax.plot(exper[:,0], exper[:,1], 'or', label='exper')
p2 = ax.plot(ld1[:,1], ld1[:,2], '--b', label='current')

ax.legend()
ax.grid()

ax.set_xlabel('Engn. strain [-]')
ax.set_ylabel('Stress [MPa]')

plt.show()
