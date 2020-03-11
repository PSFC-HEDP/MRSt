import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import sys

xlabel, ylabel, titleA, titleB = sys.argv[1:]

X = np.genfromtxt('working/{}_x.csv'.format(titleA), delimiter=',')
Y = np.genfromtxt('working/{}_y.csv'.format(titleA), delimiter=',')
Z = (
	np.genfromtxt('working/{}_z.csv'.format(titleA), delimiter=','),
	np.genfromtxt('working/{}_z.csv'.format(titleB), delimiter=','),
)

fig, ax = plt.subplots()
plt.subplots_adjust(bottom=0.2)

axis = plt.axes([0.20, 0.05, 0.60, 0.03])
slider = Slider = Slider(axis, 'Time', 0, len(X)-2, valinit=len(X)//2)
def update(*args):
	t = slider.val
	ax.clear()
	ax.plot((Y[1:] + Y[:-1])/2, Z[0][:, int(t)], '-', label=titleA)
	ax.plot((Y[1:] + Y[:-1])/2, Z[1][:, int(t)], '--', label=titleB)
	ax.legend()
	ax.set_yscale('symlog', linthresh=np.max(Z)/100, linscale=1)
	ax.set_xlabel(ylabel)
	ax.set_title("Slice comparison of {} & {}".format(titleA, titleB))
slider.on_changed(update)

update()
plt.show()
