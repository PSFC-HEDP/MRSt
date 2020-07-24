import numpy as np
import matplotlib.pyplot as plt
import sys

# import os
# os.chdir('../..')
# print(os.getcwd())
# xlabel, ylabels, title, n = 'sate', 'Ti\nρR\nYn\na2', 'data', 4
xlabel, ylabels, title, n = sys.argv[1:]

ylabels = ylabels.split('\n')
n = int(n)

X = np.genfromtxt('working/{}_x.csv'.format(title), delimiter=',')
Ys = [np.genfromtxt('working/{}_y_{}.csv'.format(title, i), delimiter=',') for i in range(n)]
Δs = [np.genfromtxt('working/{}_err_{}.csv'.format(title, i), delimiter=',') for i in range(n)]

fig, host_ax = plt.subplots()
fig.subplots_adjust(right=1 - (0.12*(n-1)))
axes = [host_ax]
plots = []
for i in range(n):
	if i > 0:
		axes.append(axes[0].twinx())
		axes[i].spines['right'].set_position(('axes', 1 + (i-1)*.18))
	if i > 1:
		axes[i].set_frame_on(True)
		axes[i].patch.set_visible(False)
		for sp in axes[i].spines.values(): sp.set_visible(False)
		axes[i].spines['right'].set_visible(True)

	rainge = {'Y':(0,None), 'T':(0,16), 'ρ':(0,2), 'V':(-100,100), 'a':(-1, 1)}[ylabels[i][0]]
	Ys[i][np.isnan(Δs[i])] = np.nan
	plots.append(axes[i].plot(X, Ys[i], label=ylabels[i], color='C'+str(i))[0])
	axes[i].fill_between(X, Ys[i] - Δs[i], Ys[i] + Δs[i], color='C'+str(i), alpha=0.3)
	axes[i].set_ylabel(ylabels[i])
	axes[i].set_ylim(*rainge)

	if ylabels[i].startswith('Y'):
		Ymax = Ys[i].max(initial=0, where=np.isfinite(Ys[i]))
		lims = np.min(np.where(Ys[i]/Ymax >= 1e-3, X, np.inf)), np.max(np.where(Ys[i]/Ymax >= 1e-3, X, -np.inf))
		if not all(np.isfinite(lims)):
			lims = X[0], X[-1]

axes[0].set_xlabel(xlabel)
axes[0].set_xlim(*lims)

axes[0].legend(plots, [p.get_label() for p in plots])

plt.tight_layout()
plt.show()
