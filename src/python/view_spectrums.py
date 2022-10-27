# view_spectrums.py – plot a bunch of the same inferred trajectories

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

plt.rcParams.update({'font.family': 'sans', 'font.size': 12})


FILENAMES = [
	'../../output/spectrums_medium_p_5c_200_2022-10-25',
]
figsize = (8, 3.5)
num_to_plot = 15

anser = pd.read_csv("../../input/scan base trajectories.csv", delimiter=',') # get the true curves

if "time(us)" in anser.columns:
	anser["time(ps)"] = anser["time(us)"]*1e6
anser["Total rhor (mg/cm^2)"] = anser["Total rhor (gm/cm^2)"]*1e3
anser["burnrate (ns^-1)"] = anser["burnrate (us^-1)"]*1e-3
# anser["burnrate(us^-1)"] *= np.sum(yields*np.gradient(times))/np.sum(anser["burnrate(us^-1)"]*np.gradient(anser["time(ps)"]))

figures = {}

for k, filename in enumerate(FILENAMES):
	times = np.loadtxt(f"{filename}_time.csv", delimiter=",")
	yields = np.loadtxt(f"{filename}_yield.csv", delimiter=",")
	temperatures = np.loadtxt(f"{filename}_temperature.csv", delimiter=",")
	densities = np.loadtxt(f"{filename}_density.csv", delimiter=",")

	times = times*1e3

	for i, (quantity, name, other_name) in enumerate([(yields*1e15, "Burn (ns^-1)", "burnrate (ns^-1)"), (temperatures, "Ion temperature (keV)", "tion(keV)"), (densities*1e3, "ρR (mg/cm^2)", "Total rhor (mg/cm^2)")]):
		if name in figures:
			fig, ax = figures[name]
		else:
			fig, ax = plt.subplots(figsize=figsize)
			figures[name] = fig, ax
		if len(FILENAMES) == 1:
			colors = f"k--", f"C{i}"
		else:
			colors = f"C{k}--", f"C{k}"
		for j in range(num_to_plot):
			ax.plot(times[j, :], quantity[j, :], colors[1], linewidth=.5)
		if len(FILENAMES) == 1:
			ax.plot(anser["time(ps)"], anser[other_name], colors[0], linewidth=1)
		ax.set_xlim(-200, 150)
		ax.set_ylim(0, None)
		ax.set_xlabel("Time (ps)")
		ax.set_ylabel(name)
		ax.set_title(filename[23:-15])
		fig.tight_layout()
		fig.savefig(f"{FILENAMES[0][:-11]}_{['burns', 'temps', 'rhors'][i]}.png", dpi=150)

plt.show()
