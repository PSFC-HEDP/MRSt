import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys

from cmap import CUSTOM_CMAP

matplotlib.rc('font', size=18)

if len(sys.argv) > 1:
	xlabel, ylabel, title = sys.argv[1:]
else:
	import os
	os.chdir('../..')
	xlabel, ylabel, title = "Time (ns)", "Energy (MeV)", "Synthetic signal distribution"
	# xlabel, ylabel, title = 'x (cm)', 'y (cm)', 'Camera 1 image'

X = np.genfromtxt('output/{}_x.csv'.format(title), delimiter=',')
Y = np.genfromtxt('output/{}_y.csv'.format(title), delimiter=',')
Z = np.genfromtxt('output/{}_z.csv'.format(title), delimiter=',')

first_nonzero = (Z.sum(axis=1) > 0).argmax()
last_nonzero = Z.shape[0] - (Z.sum(axis=1) > 0)[::-1].argmax()
maximum = Z.max(where=Z != 0, initial=-np.inf)
minimum = Z.min(where=Z != 0, initial=np.inf)

if maximum / minimum > 5e1:
	norm = matplotlib.colors.SymLogNorm(vmin=0, vmax=maximum, linthresh=max(2*minimum, maximum/1e2), linscale=1/np.log(10))
else:
	norm = matplotlib.colors.Normalize(vmax=maximum)
plt.pcolormesh(X, Y, Z, cmap=CUSTOM_CMAP["coffee"], norm=norm, rasterized=True)
plt.xlabel(xlabel, fontsize=18)
plt.ylabel(ylabel, fontsize=18)
plt.ylim(Y[first_nonzero], Y[last_nonzero])
plt.gca().xaxis.set_tick_params(labelsize=18)
plt.gca().yaxis.set_tick_params(labelsize=18)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=18)
cbar.set_label("Counts per bin")
plt.title(title, fontsize=18)

if "image" in title:
	plt.arrow(np.min(X)*.72 + np.max(X)*.28, np.min(Y)*.85 + np.max(Y)*.15, 0, Y.ptp()*.2, width=.03, color="w", head_width=.15)
	plt.text(np.min(X)*.95 + np.max(X)*.05, np.min(Y)*.95 + np.max(Y)*.05, "Sweep direction", verticalalignment="bottom", color="w")

with h5py.File(f"output/{title}.h5", "w") as f:
	for name, value in [("x", X), ("y", Y), ("z", Z)]:
		dataset = f.create_dataset(name, value.shape)
		dataset[...] = value

plt.tight_layout()
plt.savefig(f"output/{title}.png", dpi=300, transparent=True)
plt.savefig(f"output/{title}.eps")
plt.show()
