import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

data = pd.read_csv("../../output/comparison_high_2slits_400um_0c_1x_200_2022-06-14.csv")

for title in data.columns:
	if "(10^15)" in title:
		data[title.replace(" (10^15)", "")] = data[title]*1e15
	elif "(ns)" in title:
		data[title.replace("ns", "ps")] = data[title]*1e3

simulations = ["P1, P2", "Base case", "P1, extra P2, P4",
               "P1, P2, P4", "P2", "P1, P2, P4, burn off"]

for x, y in [("Total yield", "Burn width (ps)"), ("Burn skewness", "Burn width (ps)"), ("Ti at BT (keV)", "dTi/dt at BT (keV/ns)"), ("Ti at stagnation (keV)", "dTi/dt at stagnation (g/cm^2/ns)")]:
	plt.figure()
	x_means, y_means = [], []
	for i, sim in enumerate(simulations):
		here = data["Yield"] == i
		plt.scatter(data[here][x], data[here][y], c=f"C{i}", s=5)
		x_means.append(np.mean(data[here][x]))
		y_means.append(np.mean(data[here][y]))
		plt.text(x_means[-1], y_means[-1], sim)
	if np.min(x_means) > 0 and np.max(x_means) > np.min(x_means)*30:
		plt.xscale("log")
	if np.min(y_means) > 0 and np.max(y_means) > np.min(y_means)*30:
		plt.yscale("log")
	plt.xlabel(x)
	plt.ylabel(y)

plt.show()
