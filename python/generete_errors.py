import numpy as np
import matplotlib.pyplot as plt

toy = np.loadtxt('toy_intensity.dat')
errs = np.loadtxt('3a4nd2.rad')
zero = np.loadtxt('../BENCHMARKS/only_prot1/intensities/0.dat')

ref = errs[:,2]/errs[:,1]

for i in range(len(toy)):
    #toy[i,1] = np.random.normal( toy[i,1], toy[i,1] * 0.02 )
    toy[i,2] = toy[i,1] * ref[i]


#plt.loglog( toy[:,0], toy[:,1], linewidth = 0 )
plt.errorbar( toy[1:,0], toy[1:,1], toy[1:,2], marker = 'o' )
plt.loglog( zero[1:,0], zero[1:,1] )

chi2 = np.sum( ( toy[:,1] - zero[:,1])**2 / toy[:,2]**2 )
print(chi2)

np.savetxt('toy_int.dat', toy)

plt.show()
