import numpy as np
import matplotlib.pyplot as plt

penalty = np.loadtxt('penalty.dat')

plt.semilogy( penalty[:,0], penalty[:,6] )

#plt.plot( penalty[:,0], penalty[:,4] )
plt.show()
