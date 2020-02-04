import numpy as np
import matplotlib.pyplot as plt
import sys

xlabel, ylabels, title, n = sys.argv[1:]
ylabels = ylabels.split('\n')
n = int(n)

X = np.genfromtxt('working/{}_x.csv'.format(title), delimiter=',')
Ys = [np.genfromtxt('working/{}_y_{}.csv'.format(title, i), delimiter=',') for i in range(n)]

fig, host_ax = plt.subplots()
fig.subplots_adjust(right=0.37*(n-1))
axes = [host_ax]
plots = []
for i in range(n):
	if i > 0:
		axes.append(axes[0].twinx())
	if i > 1:
		axes[i].spines['right'].set_position(('axes', 1.2))
		axes[i].set_frame_on(True)
		axes[i].patch.set_visible(False)
		for sp in axes[i].spines.values(): sp.set_visible(False)
		axes[i].spines['right'].set_visible(True)
	plots.append(axes[i].plot(X, Ys[i], label=ylabels[i], color='C'+str(i))[0])
	axes[i].set_ylabel(ylabels[i])

axes[0].set_xlabel(xlabel)

axes[0].legend(plots, [p.get_label() for p in plots])

plt.show()
