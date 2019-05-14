import numpy as np
import matplotlib.pyplot as plt

penalty = np.loadtxt('p450/penalty.dat')
penalty_ref = np.loadtxt('../bead_modeling_flat/test.txt')

diff = penalty[:,3][:10000] - penalty_ref[:,3][:10000]

plt.plot( penalty[:,3][:10000] )
plt.show()
print(diff.max(), diff.min() )
