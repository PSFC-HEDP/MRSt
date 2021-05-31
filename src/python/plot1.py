import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
matplotlib.rc('font', size=18)

if len(sys.argv) <= 1:
	import os
	os.chdir('../..')
	print(os.getcwd())
	xlabel, ylabels, title, answer, n = 'Time (ns)', 'Ti (keV)\nρR (g/cm^2)\nYn (10^15/ns)', 'data', 'marginal', 3
else:
	xlabel, ylabels, title, answer, n = sys.argv[1:]

ylabels = ylabels.split('\n')
n = int(n)

XA = np.loadtxt(f'output/{title}_x.csv', delimiter=',')
YAs = [np.loadtxt(f'output/{title}_y_{i}.csv', delimiter=',') for i in range(n)]
ΔAs = [np.loadtxt(f'output/{title}_err_{i}.csv', delimiter=',') for i in range(n)]

if answer != '-':
	try:
		data = np.loadtxt(f'input/trajectories {answer}.csv', delimiter=',', skiprows=1) # get the true curves
		XB = data[:,0]
		YBs = [data[:,4], data[:,3], data[:,1], np.zeros(XB.shape)] # extract the relevant info from them
		YBs[2] *= (0.1e6/1e-6)/(1e15*14.1e6*1.6e-19/1e-9)
		while np.sum(YAs[2]*np.gradient(XA)) < np.sum(YBs[2]*(XB[1] - XB[0]))/3:
			YBs[2] /= 10
		# YBs[2] *= np.sum(YAs[2]*(XA[1] - XA[0]))/np.sum(YBs[2]*(XB[1] - XB[0])) # normalize the yield curves to account for any magnitude discrepancy
	except IOError:
		pass

fig, host_ax = plt.subplots(figsize=(8,5))
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

	rainge = {'Y':(0,None), 'T':(0,10), 'ρ':(0,1.5), 'V':(-100,100), 'a':(-1, 1)}[ylabels[i][0]]
	YAs[i][np.isnan(ΔAs[i])] = np.nan
	plots.append(axes[i].plot(XA, YAs[i], '-o', label=ylabels[i], color=f'C{i}')[0])
	if True:     axes[i].plot(XB, YBs[i], '--', color=f'C{i}')[0]
	axes[i].fill_between(XA, YAs[i] - ΔAs[i], YAs[i] + ΔAs[i], color='C'+str(i), alpha=0.3)
	axes[i].set_ylabel(ylabels[i])
	axes[i].set_ylim(*rainge)

	if ylabels[i].startswith('Y'):
		Ymax = YAs[i].max(initial=0, where=np.isfinite(YAs[i]))
		lims = np.min(XA[YAs[i]/Ymax >= 1e-3]), np.max(XA[YAs[i]/Ymax >= 1e-3])
		if not all(np.isfinite(lims)):
			lims = XA[0], XA[-1]

axes[0].set_xlabel(xlabel)
axes[0].set_xlim(*lims)

axes[0].legend(plots, [p.get_label() for p in plots])

plt.tight_layout()

# fig, axis = plt.subplots()
# for i in range(n):
# 	if ylabels[i][0] == 'T':
# 		i_temp = i
# 	elif ylabels[i][0] == 'ρ':
# 		i_dens = i
# 	elif ylabels[i][0] == 'Y':
# 		valid = YAs[i]/Ymax >= 1e-3
# axis.errorbar(x=YAs[i_dens][valid], xerr=ΔAs[i_dens][valid], y=YAs[i_temp][valid], yerr=ΔAs[i_temp][valid], fmt='-.')
# axis.plot(YBs[i_dens], YBs[i_temp], '--k')
# axis.set_xlabel("ρR (g/cm^2)")
# axis.set_ylabel("Ti (keV)")

plt.show()
