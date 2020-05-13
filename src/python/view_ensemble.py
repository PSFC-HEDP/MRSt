import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

X_LABEL = "Yield factor"

simulations = pd.read_csv('../../working/ensemble.csv')
simulations["Total yield"] = simulations["Total yield (10^15)"]*1e15

for axis, true in [
		("Bang time (ns)", 16.362), ("Max ρR (ns)", 16.309), ("Max dρR/dt (ns)", 16.257),
		("Ti at BT (keV)", 10.293), ("ρR at BT (g/cm^2)", 0.7), ("vi at BT (μm/ns)", -54.795),
		("dTi/dt at BT (keV/ns)", 111.827), ("dρR/dt at BT (g/cm^2/ns)", -3.478),
		("dvi/dt at BT (μm/ns^2)", -660.832), ("Max ρR (g/cm^2)", 0.877),
		("Total yield", 4.25e19), ("Burn mean (ns)", 16.357),
		("Burn width (ns)", .0657), ("Burn skew", -1.15), ("Burn kurtosis", 6.45)
		# ("Max ρR (ns)", 16.309), ("Max ρR (g/cm^2)", 0.877),# ("Max dρR/dt (ns)", 16.257), ("ρR at BT (g/cm^2)", 0.7),
		]:

	if 'keV' in axis:      yFactor = simulations["Temperature factor"]
	elif 'g/cm^2' in axis: yFactor = simulations["Density factor"]
	elif 'yield' in axis:  yFactor = simulations["Yield factor"]
	else:                  yFactor = np.ones(len(simulations.index))
	order = np.argsort(simulations[X_LABEL])
	order = order[np.isfinite(simulations["Total yield"].values[order])].values

	title = axis
	title = title[0].lower() + title[1:]
	if '(' in title:
		title = title[:title.index('(')-1]

	plt.figure()
	plt.plot(simulations[X_LABEL], simulations[axis], 'o', label="Based on fit to synthetic data")
	if simulations[axis].min() > 0 and simulations[axis].max()/simulations[axis].min() >= 100:
		plt.yscale('log')
	plt.plot(simulations[X_LABEL][order], yFactor[order]*true, '--', label="Based on original data")
	if 'ield' in X_LABEL:
		plt.xscale('log')
	plt.legend()
	plt.xlabel(X_LABEL)
	plt.ylabel(axis)
	plt.title("Variation in {} measurement with varying spectra".format(title))
	plt.tight_layout()

plt.show()
