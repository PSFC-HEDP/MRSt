import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

filename = "../../output/comparison_high_maximu_400um_0c_1x_200_2022-07-07.csv"
data = pd.read_csv(filename)

for title in data.columns:
	if "(10^15)" in title:
		data[title.replace(" (10^15)", "")] = data[title]*1e15
	elif "(ns)" in title:
		data[title.replace("ns", "ps")] = data[title]*1e3
	elif "/ns" in title:
		data[title.replace("/ns", "/(100ps)")] = data[title]*1e-1

simulations = [
	"Base case", "P2", "P1, P2", "P1, P2, P4", "P1, extra P2, P4",
	# "P1, P2, P4, burn off",
]
order = [1, 4, 0, 3, 2, 5]

# ground_truth = {
# 	"Base case":            [4.194e+17, 64.536, -0.385, 3.548, 14.559, 23.903, 0.887, -1.529],
# 	"P2":                   [1.749e+17, 81.188, -0.465, 3.512, 11.160, 2.098, 0.805, -0.900],
# 	"P1, P2":               [1.642e+17, 99.101, -0.501, 3.499, 10.593, 3.082, 0.864, -1.213],
# 	"P1, P2, P4":           [8.922e+16, 93.162, -0.882, 5.546, 9.485, 0.774, 0.945, -0.776],
# 	"P1, extra P2, P4":     [3.04e+16, 84.297, -0.761, 4.316, 7.628, 0.656, 1.208, -0.691],
# 	"P1, P2, P4, burn off": [2.246e+15, 98.651, -0.376, 4.091, 4.511, -2.467, 1.379, 1.007],
# }
ground_truth = {
	"Base case":            [421.3687017e15, 82.91398085,-0.2900906493,3.295152373,14.4519683,6.350179054,0.8569425781,-.9141120023],
	"P2":                   [174.6411255e15,102.1353122,-0.18426077,3.380619635,10.85817745,1.177295346,0.856886018,-.7814462052],
	"P1, P2":               [164.0954218e15,107.4330766,-0.1588381968,3.168865761,10.47252618,.6080325039,0.856005883,-.7120849961],
	"P1, P2, P4":           [88.71657896e15,102.0059732,-0.5901396969,4.288099759,9.413961466,1.208433783,0.9520658982,-.7003857988],
	"P1, extra P2, P4":     [30.41117331e15, 99.92889588,-0.7871359762,4.387469235,7.56604614,-.4128174512,1.20743159,-.5336200009],
	"P1, P2, P4, burn off": [2.245892065e15,124.3532818,-0.4010667207,3.985329307,4.581542467,-1.662217588,1.329866331,.70173167],
}
short_header = ["Total yield", "Burn width (ps)", "Burn skewness", "Burn kurtosis", "Ti at BT (keV)", "dTi/dt at BT (keV/(100ps))", "ρR at BT (g/cm^2)", "dρR/dt at BT (g/cm^2/(100ps))"]

for x, y in [("Total yield", "Burn width (ps)"),
             ("Burn skewness", "Burn width (ps)"),
             ("Burn skewness", "Burn kurtosis"),
             ("Ti at BT (keV)", "dTi/dt at BT (keV/(100ps))"),
             # ("Ti at stagnation (keV)", "dTi/dt at stagnation (keV/(100ps))"),
	]:
	plt.figure()
	x_means, y_means = [], []
	for i, sim in enumerate(simulations):
		here = data["Yield"] == order[i]
		plt.scatter(data[here][x], data[here][y], c=f"C{i}", s=4, label=sim)
		x_means.append(np.mean(data[here][x]))
		y_means.append(np.mean(data[here][y]))
		true_values = ground_truth[sim]
		plt.scatter(true_values[short_header.index(x)], true_values[short_header.index(y)],
			s=50, facecolor="white", edgecolor=f"C{i}", linewidth=2, zorder=10)
	plt.legend()
	if np.min(x_means) > 0 and np.max(x_means) > np.min(x_means)*30:
		plt.xscale("log")
	if np.min(y_means) > 0 and np.max(y_means) > np.min(y_means)*30:
		plt.yscale("log")
	plt.xlabel(x)
	plt.ylabel(y)
	plt.title(filename[24:-19])
	plt.tight_layout()
	plt.savefig("../../output/compare_{}_{}_{}".format(
		re.sub(r'[/ ]', '', re.sub(r'\(.*\)', '', x)),
		re.sub(r'[/ ]', '', re.sub(r'\(.*\)', '', y)),
		filename[24:-19],
	), dpi=300)

plt.show()