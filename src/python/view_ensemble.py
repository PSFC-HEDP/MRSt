import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

simulations = pd.read_csv('../../working/ensemble.csv')
# print(simulations)

for axis in [
		"Bang time (ns)", "Max ρR (ns)", "Max dρR/dt (ns)", "Ti at BT (keV)", "ρR at BT (g/cm^2)",
		"vi at BT (μm/ns)", "dTi/dt at BT (keV/ns)", "dρR/dt at BT (g/cm^2/ns)",
		"dvi/dt at BT (μm/ns^2)", "Max ρR (g/cm^2)", "Total yield (10^15)", "Burn mean (ns)",
		"Burn width (ns)", "Burn skew", "Burn kurtosis"]:
	plt.figure()
	plt.xscale('log')
	if simulations[axis][0] > 0 and simulations[axis].max()/simulations[axis].min() >= 100:
		plt.yscale('log')
	plt.scatter(simulations['Yield factor'], simulations[axis])
	plt.xlabel("Yield factor")
	plt.ylabel(axis)
	plt.tight_layout()

plt.show()
