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
	xlabel, ylabel, title = 'Streak direction (cm)', 'Slit direction (cm)', 'Camera 0 image'

X = np.genfromtxt('output/{}_x.csv'.format(title), delimiter=',')
Y = np.genfromtxt('output/{}_y.csv'.format(title), delimiter=',')
Z = np.genfromtxt('output/{}_z.csv'.format(title), delimiter=',')

print(Z)
first_nonzero = (Z.sum(axis=1) > 0).argmax()
last_nonzero = Z.shape[0] - (Z.sum(axis=1) > 0)[::-1].argmax()
maximum = Z.max(where=Z != 0, initial=-np.inf)
minimum = Z.min(where=Z != 0, initial=np.inf)

if maximum / minimum > 5e3:
	norm = matplotlib.colors.SymLogNorm(vmin=0, vmax=maximum, linthresh=max(1, maximum/3e3), linscale=1/np.log(10))
elif maximum / minimum > 5e1:
	norm = matplotlib.colors.LogNorm(vmin=minimum, vmax=maximum)
else:
	norm = matplotlib.colors.Normalize(vmax=maximum)
plt.pcolormesh(X, Y, Z, cmap='plasma', norm=norm, rasterized=True)
plt.xlabel(xlabel, fontsize=18)
plt.ylabel(ylabel, fontsize=18)
plt.ylim(Y[first_nonzero], Y[last_nonzero])
plt.gca().xaxis.set_tick_params(labelsize=18)
plt.gca().yaxis.set_tick_params(labelsize=18)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=18) 
plt.title(title, fontsize=18)

plt.tight_layout()
plt.savefig(f"output/{title}.png", dpi=300)
plt.savefig(f"output/{title}.eps")
plt.show()
