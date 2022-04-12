import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import sys
matplotlib.rc('font', size=16)

def true_max(*args):
	return max([np.max(arg, where=~np.isnan(arg), initial=-np.inf) for arg in args])

def true_min(*args):
	return min([np.min(arg, where=~np.isnan(arg), initial=np.inf) for arg in args])

if len(sys.argv) <= 1:
	import os
	os.chdir('../..')
	# xlabel, ylabel, titleA, titleB = "Time (ns)", "Energy (MeV)", "Synthetic signal distribution", "Detector background"
	xlabel, ylabel, titleA, titleB = "Time (ns)", "Energy (MeV)", "Corrected signal distribution", "Synthetic deuteron spectrum"
else:
	xlabel, ylabel, titleA, titleB = sys.argv[1:]

X = np.genfromtxt('output/{}_x.csv'.format(titleA), delimiter=',')
Y = np.genfromtxt('output/{}_y.csv'.format(titleA), delimiter=',')
Z = [
	np.genfromtxt('output/{}_z.csv'.format(titleA), delimiter=','),
	np.genfromtxt('output/{}_z.csv'.format(titleB), delimiter=','),
]

Y = (Y[1:] + Y[:-1])/2

if "(ns)" in xlabel:
	x0 = X[np.argmax(np.sum(Z[0], axis=0))]
	X = (X - x0)*1000
	xlabel = xlabel.replace("ns", "ps")

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)

axis = plt.axes([0.20, 0.05, 0.60, 0.03])
slider = Slider = Slider(axis, xlabel, X[0], X[-1], valinit=X[len(X)//2])
def update(*args):
	x = slider.val
	j = np.minimum(np.digitize(x, X) - 1, Z[0].shape[1] - 1)
	ax.clear()
	ax.plot(Y, Z[0][:, j], '-', label=titleA)
	ax.plot(Y, Z[1][:, j], '--', label=titleB)
	# ax.legend()
	if Z[1][:, j].max()/Z[1][:, j].min(initial=np.inf, where=Z[1][:, j] > 0) > 1e3:
		ax.set_yscale('symlog', linthresh=max(1, true_max(Z[0][:, j], Z[1][:, j])/100), linscale=1/np.log(10))
	ax.set_xlabel(ylabel)
	ax.set_title("Slice comparison of {} & {}".format(titleA, titleB))
slider.on_changed(update)

update()
plt.show()
