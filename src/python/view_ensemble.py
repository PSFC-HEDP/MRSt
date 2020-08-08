import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'serif', 'font.size': 10})

X_LABEL = "Yield factor"

Y_LABELS = [
	("Total yield", 4.802e17),
	("Burn-average ρR (g/cm^2)", .692), ("Burn-average Ti (keV)", 10.205),
	("Bang time (ns)", 16.364), ("Max \u03C1R (ns)", 16.302),
	("Burn width (ps)", 66.705),
	("Burn skew", -1.142), ("Burn kurtosis", 6.496),
	("dρR/dt at BT (mg/cm^2/(100ps))", -214), ("dTi/dt at BT (keV/(100ps))", 3.8),
	("Burn-average vi (km/s)", -.722), ("dvi/dt at BT (km/s/(100ps))", -114.2)
]

COLUMNS = 2
SIZE = (8.0, 10.5)
MARGIN = dict(bottom=.06, top=.94, left=.12, right=.99, wspace=.35, hspace=.05)
# COLUMNS = 3
# SIZE = (16, 9)
# MARGIN = dict(bottom=.06, top=.94, left=.06, right=.99, wspace=.26, hspace=.05)


def text_wrap(s):
	if len(s) > 14:
		i = len(s)//2
		for j in range(i):
			if s[i+j] == ' ':
				return s[:i+j] + '\n' + s[i+j+1:]
			elif s[i-j] == ' ':
				return s[:i-j] + '\n' + s[i-j+1:]
	return s


simulations = pd.read_csv('../../working/ensemble.csv', na_values=["Infinity"])
for suf in ["", " error"]:
	simulations["Total yield"+suf] = simulations["Total yield (10^15)"+suf]*1e15
	simulations["Burn width (ps)"+suf] = simulations["Burn width (ns)"+suf]/1e-3
	simulations["dρR/dt at BT (mg/cm^2/(100ps))"+suf] = simulations["dρR/dt at BT (g/cm^2/ns)"+suf]/1e-2
	simulations["dTi/dt at BT (keV/(100ps))"+suf] = simulations["dTi/dt at BT (keV/ns)"+suf]/1e1
	simulations["dvi/dt at BT (km/s/(100ps))"+suf] = simulations["dvi/dt at BT (km/s/ns)"+suf]/1e1

fig, axs = plt.subplots((len(Y_LABELS) + COLUMNS-1)//COLUMNS, COLUMNS, figsize=SIZE)
fig.subplots_adjust(**MARGIN)

for i, (axis, true) in enumerate(Y_LABELS):
	if 'keV' in axis:      yFactor = simulations["Temperature factor"]
	elif 'g/cm^2' in axis: yFactor = simulations["Down-scatter factor"]
	elif 'yield' in axis:  yFactor = simulations["Yield factor"]
	else:                  yFactor = np.ones(len(simulations.index))
	order = np.argsort(simulations[X_LABEL])
	order = order[np.isfinite(simulations["Total yield"].values[order])].values

	ax = axs[i//COLUMNS,i%COLUMNS]
	valid = np.logical_not(np.isnan(simulations[axis+" error"]))
	ax.scatter(simulations[X_LABEL][valid], simulations[axis][valid], s=10, zorder=1, label="Based on fit to synthetic data")
	ax.set_ylim(auto=False)
	ax.errorbar(simulations[X_LABEL][valid], simulations[axis][valid], yerr=simulations[axis+" error"][valid], elinewidth=1, linestyle='none')
	if simulations[axis].min() > 0 and simulations[axis].max()/simulations[axis].min() >= 100:
		ax.set_yscale('log')
		ax.set_ylim(simulations[axis].min()/1.5, simulations[axis].max()*1.5)
	ax.plot(simulations[X_LABEL][order], yFactor[order]*true, 'C1--', zorder=0, label="Based on original data")
	if 'ield' in X_LABEL:
		ax.set_xscale('log')
	# ax.legend()
	if i//COLUMNS == axs.shape[0]-1:
		ax.set_xlabel(X_LABEL)
	elif i//COLUMNS == 0:
		ax.set_xlabel(X_LABEL)
		ax.xaxis.set_label_position('top')
		ax.xaxis.tick_top()
	else:
		ax.xaxis.set_visible(False)
	ax.set_ylabel(text_wrap(axis))

plt.savefig('../../working/big-plot.eps')
plt.savefig('../../working/big-plot.png')
plt.show()
