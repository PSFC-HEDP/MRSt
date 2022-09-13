# view_spectrums.py – plot a bunch of the same inferred trajectories

import numpy as np
from matplotlib import pyplot as plt


FILENAME = '../../output/spectrums_medium_driftt_10c_15ps_20_2022-09-12'
figsize = (8, 4)

times = np.loadtxt(f"{FILENAME}_time.csv", delimiter=",")
yields = np.loadtxt(f"{FILENAME}_yield.csv", delimiter=",")
temperatures = np.loadtxt(f"{FILENAME}_temperature.csv", delimiter=",")
densities = np.loadtxt(f"{FILENAME}_density.csv", delimiter=",")

t0 = times[np.unravel_index(np.argmax(yields), yields.shape)]
times = (times - t0)*1e3

num_to_plot = 10

for i, (quantity, name) in enumerate([(yields*1e15, "Burn (ns^-1)"), (temperatures, "Ion temperature (keV)"), (densities*1e3, "ρR (mg/cm^2")]):
	plt.figure(figsize=figsize)
	for j in range(num_to_plot):
		plt.plot(times[j, :], quantity[j, :], f"C{i}", linewidth=.5)
	plt.xlim(np.min(times), np.max(times))
	plt.ylim(0, None)
	plt.xlabel("Time (ps)")
	plt.ylabel(name)

plt.show()
