import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import sys

if len(sys.argv) <= 1:
	import os
	os.chdir('../..')
	print(os.getcwd())
	xlabel, ylabel, titleA, titleB = 'Energy', 'Count', 'Original neutron spectrum', 'Fitted neutron spectrum'
else:
	xlabel, ylabel, titleA, titleB = sys.argv[1:]

X = np.genfromtxt('working/{}_x.csv'.format(titleA), delimiter=',')
Y = np.genfromtxt('working/{}_y.csv'.format(titleA), delimiter=',')
Z = (
	np.genfromtxt('working/{}_z.csv'.format(titleA), delimiter=','),
	np.genfromtxt('working/{}_z.csv'.format(titleB), delimiter=','),
)

if 'euteron' in titleA:
	Y *= 8/9

X = (X[1:] + X[:-1])/2
Y = (Y[1:] + Y[:-1])/2

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)

axis = plt.axes([0.20, 0.05, 0.60, 0.03])
slider = Slider = Slider(axis, 'Time', X[0], X[-1], valinit=X[len(X)//2])
def update(*args):
	t = slider.val
	j = int(round(np.interp(t, X, np.arange(len(X)))))
	ax.clear()
	ax.plot(Y, Z[0][:, j], '-', label=titleA)
	ax.plot(Y, Z[1][:, j], '--', label=titleB)
	ax.legend()
	if np.max([Z[0][:, j], Z[1][:, j]]) > 0:
		ax.set_yscale('symlog', linthresh=np.max([Z[0][:, j], Z[1][:, j]])/100, linscale=1)
	ax.set_xlabel(ylabel)
	ax.set_title("Slice comparison of {} & {}".format(titleA, titleB))
slider.on_changed(update)

update()
plt.show()
