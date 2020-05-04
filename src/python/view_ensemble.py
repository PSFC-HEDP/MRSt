import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

X_LABEL = "Yield factor"

simulations = pd.read_csv('../../working/yield.csv')
simulations["Total yield"] = simulations["Total yield (10^15)"]*1e15
simulations["Temperature factor"] = simulations["Temperature factor"]**2

for axis, true in [
		("Bang time (ns)", 16.362), ("Max ρR (ns)", 16.309), ("Max dρR/dt (ns)", 16.257),
		("Ti at BT (keV)", 10.293), ("ρR at BT (g/cm^2)", 0.7), ("vi at BT (μm/ns)", -54.795),
		("dTi/dt at BT (keV/ns)", 111.827), ("dρR/dt at BT (g/cm^2/ns)", -3.478),
		("dvi/dt at BT (μm/ns^2)", -660.832), ("Max ρR (g/cm^2)", 0.877),
		("Total yield", 4.25e19), ("Burn mean (ns)", 16.357),
		("Burn width (ns)", .0657), ("Burn skew", -1.15), ("Burn kurtosis", 6.45)
		]:
		
	title = axis
	title = title[0].lower() + title[1:]
	if '(' in title:
		title = title[:title.index('(')-1]

	plt.figure()
	plt.plot(simulations[X_LABEL], simulations[axis], 'o', label="Based on fit to synthetic data")
	if simulations[axis].min() > 0 and simulations[axis].max()/simulations[axis].min() >= 100:
		plt.yscale('log')
	else:
		plt.plot([simulations[X_LABEL].min(), simulations[X_LABEL].max()], [true, true], '--', label="Based on original data")
	if 'ield' in X_LABEL:
		plt.xscale('log')
	plt.legend()
	plt.xlabel(X_LABEL)
	plt.ylabel(axis)
	plt.title("Variation in {} measurement with varying yield".format(title))
	plt.tight_layout()

plt.show()
