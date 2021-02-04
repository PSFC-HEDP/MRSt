import numpy as np
import matplotlib.pyplot as plt
import sys

if len(sys.argv) <= 1:
	import os
	os.chdir('../..')
	print(os.getcwd())
	xlabel, ylabels, title, n = 'Time (ns)', 'Ti (keV)\nρR (g/cm^2)\nYn (10^15/ns)\nVi (km/s)', 'data', 4
else:
	xlabel, ylabels, title, n = sys.argv[1:]

ylabels = ylabels.split('\n')
n = int(n)

XA = np.loadtxt('working/{}_x.csv'.format(title), delimiter=',')
YAs = [np.loadtxt('working/{}_y_{}.csv'.format(title, i), delimiter=',') for i in range(n)]
ΔAs = [np.loadtxt('working/{}_err_{}.csv'.format(title, i), delimiter=',') for i in range(n)]

if True:
	data = np.loadtxt('data/Yn-rR-Ti_150327_16p26 - Yn-rR-Ti_150327_16p26.csv', delimiter=',', skiprows=1)
	XB = data[:,0]
	YBs = [data[:,2], data[:,4] + data[:,5], data[:,1], np.zeros(XB.shape)]
	YBs[2] *= np.sum(YAs[2]*(XA[1] - XA[0]))/np.sum(YBs[2]*(XB[1] - XB[0]))

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

	rainge = {'Y':(0,None), 'T':(0,16), 'ρ':(0,1.5), 'V':(-100,100), 'a':(-1, 1)}[ylabels[i][0]]
	YAs[i][np.isnan(ΔAs[i])] = np.nan
	plots.append(axes[i].plot(XA, YAs[i], label=ylabels[i], color=f'C{i}')[0])
	if True:     axes[i].plot(XB, YBs[i], '--', color=f'C{i}')[0]
	axes[i].fill_between(XA, YAs[i] - ΔAs[i], YAs[i] + ΔAs[i], color='C'+str(i), alpha=0.3)
	axes[i].set_ylabel(ylabels[i])
	axes[i].set_ylim(*rainge)

	if ylabels[i].startswith('Y'):
		Ymax = YAs[i].max(initial=0, where=np.isfinite(YAs[i]))
		lims = np.min(np.where(YAs[i]/Ymax >= 1e-3, XA, np.inf)), np.max(np.where(YAs[i]/Ymax >= 1e-3, XA, -np.inf))
		if not all(np.isfinite(lims)):
			lims = XA[0], XA[-1]

axes[0].set_xlabel(xlabel)
axes[0].set_xlim(*lims)

axes[0].legend(plots, [p.get_label() for p in plots])

plt.tight_layout()
plt.show()
