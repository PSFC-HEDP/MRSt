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
	xlabel, ylabel, titles = "Time (ns)", "Energy (MeV)", ["Corrected signal distribution", "Fit signal distribution", "True spectrum"]
else:
	xlabel, ylabel, *titles = sys.argv[1:]

X = np.genfromtxt('output/{}_x.csv'.format(titles[0]), delimiter=',')
Y = np.genfromtxt('output/{}_y.csv'.format(titles[0]), delimiter=',')
Zs = np.array([
	np.genfromtxt('output/{}_z.csv'.format(title), delimiter=',')
	for title in titles
])

Y = (Y[1:] + Y[:-1])/2

if "(ns)" in xlabel:
	x0 = X[np.argmax(np.sum(Zs[-1, :, :], axis=0))]
	X = (X - x0)*1000
	xlabel = xlabel.replace("ns", "ps")

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)

axis = plt.axes([0.20, 0.05, 0.60, 0.03])
slider = Slider = Slider(axis, xlabel, X[0], X[-1], valinit=X[len(X)//2])
def update(*args):
	x = slider.val
	j = np.minimum(np.digitize(x, X) - 1, Zs.shape[2] - 1)
	ax.clear()
	if Zs[0, :, j].max()/Zs[0, :, j].min(initial=np.inf, where=Zs[1, :, j] > 0) > 1e3:
		ax.set_yscale('symlog', linthresh=max(1, Zs[0, :, j].max()/100), linscale=1/np.log(10))
	limits = None
	for i in range(Zs.shape[0]):
		fmt = ["C0-", "C1--", "C2--"]
		ax.plot(Y, Zs[i, :, j], fmt[i], label=titles[i])
		if limits is None:
			limits = ax.axis()
	ax.axis(limits)
	ax.set_ylim(0, None)
	ax.set_xlabel(ylabel)
	ax.set_title(f"Slice comparison of {' & '.join(titles)}")
slider.on_changed(update)

update()
plt.show()
