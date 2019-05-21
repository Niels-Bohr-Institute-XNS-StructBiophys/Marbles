import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.collections import LineCollection

penalty = np.loadtxt('p450/penalty.dat')
npasses = 80

# plt.figure(0)
# plt.semilogy( penalty[:,0]/1000, penalty[:,2] )
# plt.xlabel('Steps/1000')
# plt.ylabel(r'$\chi^2$')
# plt.show()
#
# plt.figure(1)
# plt.plot( penalty[:,0]/1000, penalty[:,3] )
# plt.xlabel('Steps/1000')
# plt.ylabel('Type Penalty')
# plt.show()
#
# plt.figure(2)
# plt.semilogy( penalty[:,0]/1000, penalty[:,4] )
# plt.xlabel('Steps/1000')
# plt.ylabel('Histogram Penalty')
# plt.show()
#
# plt.figure(3)
# plt.semilogy( penalty[:,0]/1000, penalty[:,5] )
# plt.xlabel('Steps/1000')
# plt.ylabel('Connect Penalty')
# plt.show()

plt.figure(0)
lw = 1.3
plt.semilogy( penalty[:,0]/1000, penalty[:,2], label = r'$\chi^2$', linewidth = lw, c = 'tab:blue' )
plt.semilogy( penalty[:,0]/1000, penalty[:,3], label = 'Type', linewidth = lw, c = 'tab:red' )
plt.semilogy( penalty[:,0]/1000, penalty[:,4], label = 'Histogram', linewidth = lw, c = 'tab:green' )
plt.semilogy( penalty[:,0]/1000, penalty[:,5], label = 'Connect', linewidth = lw, c = 'tab:grey' )
plt.semilogy( penalty[:,0]/1000, penalty[:,6], label = 'Total', linewidth = lw, c = 'k' )
plt.legend()
plt.xlabel('Steps/1000')
plt.ylabel('Penalty')
plt.tight_layout()
plt.savefig('penalty.pdf', format = 'pdf')

plt.figure(1)
cm_subsection = np.linspace(0, 1, npasses)

colors = [ cm.coolwarm(x) for x in cm_subsection ]

for i, color in enumerate(colors):
    intensity = np.loadtxt( f'p450/intensities/{i}.dat' )
    plt.loglog( intensity[1:,0], intensity[1:,1], color=color )

#plt.colorbar()
plt.show()

plt.figure(2)
intensity = np.loadtxt( f'p450/intensities/{npasses-1}.dat' )

rad = np.loadtxt('3a4nd2.rad')
plt.loglog( intensity[1:,0], intensity[1:,1], color='k', label = 'Fit' )
plt.loglog( rad[1:,0], rad[1:,1], marker = 'o', c = 'tab:red', markersize = 2, linewidth=0, label = 'Experimental data' )
plt.fill_between( rad[1:,0], rad[1:,1] - rad[1:,2], rad[1:,1] + rad[1:,2], alpha = 0.4, linewidth=0, color = 'tab:red', label = 'Error' )

plt.xlabel(r'q [$\AA^{-1}$]')
plt.ylabel(r'Intensity [cm$^{-1}$]')
plt.legend()
plt.tight_layout()
plt.savefig('fit.pdf', format = 'pdf')
#plt.show()
