import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
matplotlib.rc('font', size=18)

if len(sys.argv) > 1:
	xlabel, ylabel, title = sys.argv[1:]
else:
	import os
	os.chdir('../..')
	xlabel, ylabel, title = 'Time (ns)', 'Energy (MeV)', 'Original neutron spectrum'

X = np.genfromtxt('output/{}_x.csv'.format(title), delimiter=',')
Y = np.genfromtxt('output/{}_y.csv'.format(title), delimiter=',')
Z = np.genfromtxt('output/{}_z.csv'.format(title), delimiter=',')

if 'euteron' in title:
	Y *= 8/9

plt.pcolormesh(X, Y, Z, cmap='plasma', norm=matplotlib.colors.SymLogNorm(vmin=0, vmax=Z.max(), linthresh=Z.max()/100, linscale=1))
plt.xlabel(xlabel, fontsize=18)
plt.ylabel(ylabel, fontsize=18)
plt.gca().xaxis.set_tick_params(labelsize=18)
plt.gca().yaxis.set_tick_params(labelsize=18)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=18) 
plt.title(" ", fontsize=18)

plt.tight_layout()
plt.show()
