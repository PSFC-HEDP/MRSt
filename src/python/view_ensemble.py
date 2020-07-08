import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

X_LABEL = "Yield factor"

Y_LABELS = [
	("Bang time (ns)", 16.364), ("Max ρR (ns)", 16.310), ("Max dρR/dt (ns)", 16.243),
	("Ti at BT (keV)", 10.552), ("ρR at BT (g/cm^2)", 1.069), ("vi at BT (μm/ns)", 10.595),
	("dTi/dt at BT (keV/ns)", 107.44), ("dρR/dt at BT (g/cm^2/ns)", -8.05),
	("dvi/dt at BT (μm/ns^2)", -1130), ("Max ρR (g/cm^2)", 1.426),
	("Total yield", 5.37e17), ("Burn mean (ns)", 16.357),
	("Burn width (ns)", .0666), ("Burn skew", -1.11), ("Burn kurtosis", 6.32)
]

COLUMNS = 3

simulations = pd.read_csv('../../working/ensemble.csv')
simulations["Total yield"] = simulations["Total yield (10^15)"]*1e15

fig, axs = plt.subplots((len(Y_LABELS) + COLUMNS-1)//COLUMNS, COLUMNS)
fig.subplots_adjust(hspace=0)

for i, (axis, true) in enumerate(Y_LABELS):
	if 'keV' in axis:      yFactor = simulations["Temperature factor"]
	elif 'g/cm^2' in axis: yFactor = simulations["Down-scatter factor"]
	elif 'yield' in axis:  yFactor = simulations["Yield factor"]
	else:                  yFactor = np.ones(len(simulations.index))
	order = np.argsort(simulations[X_LABEL])
	order = order[np.isfinite(simulations["Total yield"].values[order])].values

	title = axis
	title = title[0].lower() + title[1:]
	if '(' in title:
		title = title[:title.index('(')-1]

	ax = axs[i//COLUMNS,i%COLUMNS]
	ax.scatter(simulations[X_LABEL], simulations[axis], s=10, zorder=1, label="Based on fit to synthetic data")
	if simulations[axis].min() > 0 and simulations[axis].max()/simulations[axis].min() >= 1000:
		ax.set_yscale('log')
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
	ax.set_ylabel(axis)

plt.show()
