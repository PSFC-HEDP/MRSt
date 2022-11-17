# view_spectrums.py – plot a bunch of the same inferred trajectories

import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

plt.rcParams.update({'font.family': 'sans', 'font.size': 12})

FILENAMES = [
	# '../../output/spectrums_medium_5c_100_2022-11-09',
	'../../output/spectrums_medium_5c_100_scan-base_2022-11-09',
	'../../output/spectrums_medium_5c_100_scan-p2_2022-11-09',
	'../../output/spectrums_medium_5c_scan-p2-p1_100_2022-11-09',
	'../../output/spectrums_medium_5c_scan-p2-p1-p4_100_2022-11-16',
	'../../output/spectrums_medium_5c_scan-p2-p1-p4-p2_100_2022-11-10',
	'../../output/spectrums_medium_5c_scan-p2-p1-p4-burnoff_100_2022-11-16',
]
FANCY_COLORS = ["#E34B66", "#D3A94C", "#037F58", "#2AA4E9", "#88379D", "#594047"]
figsize = (8, 4)
num_to_plot = 15

anser = pd.read_csv("../../input/scan base trajectories.csv", delimiter=',') # get the true curves

if "time(us)" in anser.columns:
	anser["time(ps)"] = anser["time(us)"]*1e6
anser["Total rhor (g/cm^2)"] = anser["Total rhor (gm/cm^2)"]
anser["burnrate (ns^-1)"] = anser["burnrate (us^-1)"]*1e-3
# anser["burnrate(us^-1)"] *= np.sum(yields*np.gradient(times))/np.sum(anser["burnrate(us^-1)"]*np.gradient(anser["time(ps)"]))

figures = {}

for k, filename in enumerate(FILENAMES):
	times = np.loadtxt(f"{filename}_time.csv", delimiter=",")
	yields = np.loadtxt(f"{filename}_yield.csv", delimiter=",")
	temperatures = np.loadtxt(f"{filename}_temperature.csv", delimiter=",")
	densities = np.loadtxt(f"{filename}_density.csv", delimiter=",")

	times = times*1e3

	for i, (quantity, name, other_name) in enumerate([(yields*1e15, "Burn (ns^-1)", "burnrate (ns^-1)"), (temperatures, "Ion temperature (keV)", "tion(keV)"), (densities, "ρR (mg/cm^2)", "Total rhor (g/cm^2)")]):
		if name in figures:
			fig, ax = figures[name]
		else:
			fig, ax = plt.subplots(figsize=figsize)
			figures[name] = fig, ax
		if len(FILENAMES) == 1:
			colors = f"--k", f"C{i}"
		else:
			colors = f"--{FANCY_COLORS[k]}", f"{FANCY_COLORS[k]}"
		for j in range(num_to_plot):
			ax.plot(times[j, :], quantity[j, :], colors[1], linewidth=.5, alpha=.5)
		if len(FILENAMES) == 1:
			ax.plot(anser["time(ps)"], anser[other_name], colors[0], linewidth=1)
		ax.set_xlim(-100, 100)
		ax.set_ylim(None, None)
		ax.set_xlabel("Time (ps)")
		ax.set_ylabel(name)
		ax.set_title(filename[23:-15])
		fig.tight_layout()
		fig.savefig(f"{FILENAMES[0][:-11]}_{['burns', 'temps', 'rhors'][i]}.png", dpi=300, transparent=True)

plt.show()
