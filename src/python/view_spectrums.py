# view_spectrums.py – plot a bunch of the same inferred trajectories

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


FILENAME = '../../output/spectrums_medium_driftt_5c_15ps_20_base_2022-09-28'
figsize = (8, 3.5)

times = np.loadtxt(f"{FILENAME}_time.csv", delimiter=",")
yields = np.loadtxt(f"{FILENAME}_yield.csv", delimiter=",")
temperatures = np.loadtxt(f"{FILENAME}_temperature.csv", delimiter=",")
densities = np.loadtxt(f"{FILENAME}_density.csv", delimiter=",")

times = times*1e3

num_to_plot = 15

anser = pd.read_csv("../../input/scan/base trajectories.csv", delimiter=',') # get the true curves

if "time(us)" in anser.columns:
	anser["time(ps)"] = anser["time(us)"]*1e6
anser["Total rhor (mg/cm^2)"] = anser["Total rhor (gm/cm^2)"]*1e3
anser["burnrate (ns^-1)"] = anser["burnrate (us^-1)"]*1e-3
# anser["burnrate(us^-1)"] *= np.sum(yields*np.gradient(times))/np.sum(anser["burnrate(us^-1)"]*np.gradient(anser["time(ps)"]))

for i, (quantity, name, other_name) in enumerate([(yields*1e15, "Burn (ns^-1)", "burnrate (ns^-1)"), (temperatures, "Ion temperature (keV)", "tion(keV)"), (densities*1e3, "ρR (mg/cm^2)", "Total rhor (mg/cm^2)")]):
	plt.figure(figsize=figsize)
	plt.plot(anser["time(ps)"], anser[other_name], f"k--", linewidth=1)
	for j in range(num_to_plot):
		plt.plot(times[j, :], quantity[j, :], f"C{i}", linewidth=.5)
	plt.xlim(np.min(times), np.max(times))
	plt.ylim(0, None)
	plt.xlabel("Time (ps)")
	plt.ylabel(name)
	plt.tight_layout()
	plt.savefig(f"{FILENAME[:-11]}_{['burns', 'temps', 'rhors'][i]}.png", dpi=150)

plt.show()
