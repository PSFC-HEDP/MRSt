import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams.update({'font.family': 'sans', 'font.size': 12})


# filename = "../../output/comparison_perfect_perfect_15ps_.001_20_2022-12-09.csv"
filename = "../../output/comparison_medium_drifttub_0.001_1000_2023-03-09.csv"
CASES = ["014", "008", "023", "017", "021"]
data = pd.read_csv(filename)

for title in data.columns:
	if "(10^15)" in title:
		data[title.replace(" (10^15)", "")] = data[title]*1e15
	elif "(ns)" in title:
		data[title.replace("ns", "ps")] = data[title]*1e3
	elif "/ns" in title:
		data[title.replace("/ns", "/(100ps)")] = data[title]*1e-1

case_names = [
	"014", "008", "023", "017", "021",
]

ground_truth_min = ground_truth_max = {
	# based on the simulations
	"014": [8.29E+17, 68.10963735071582,-0.2900906493,3.295152373,14.4519683,6.350179054,0.8569425781,-.778,9.970616694,6.236026801,7.365489006,5.221992371,8.477012035,5.540807458,-82.70139249],
	"008": [5.62E+17, 81.57827724303338,-0.18426077,3.380619635,10.85817745,1.177295346,0.856886018,-.7814462052,8.960380402,3.400291759,6.978207048,4.202804277,7.874669961,4.256110662,-87.44125714],
	"023": [5.33E+17, 100.86974049543558,-0.1588381968,3.168865761,10.47252618,.6080325039,0.856005883,-.7120849961,9.015436804,3.10824185,7.061600627,3.981685591,7.484156813,4.064857076,-93.12940659],
	"017": [3.72E+17, 91.64947774226279,-0.5901396969,4.288099759,9.413961466,1.208433783,0.9520658982,-.7003857988,7.599155994,3.602491165,6.062004003,2.431564928,7.496507264,3.58146492,-58.44907688],
	"021": [2.01E+17, 82.80405354958481,-0.7871359762,4.387469235,7.56604614,-.4128174512,1.20743159,-.5336200009,6.913116212,1.990664254,5.71167482,1.591588135,7.053548379,1.845198033,-43.27693469],
}
# ground_truth_min = {
# 	# based on the model (min)
# 	"Base case":            [421.3687017e15, 68.10963735071582,-0.2900906493,3.295152373,14.4519683,6.350179054,0.8569425781,-.778,9.970616694,6.236026801,7.365489006,5.221992371,8.477012035,5.540807458,-82.70139249],
# 	"P2":                   [174.6411255e15, 81.57827724303338,-0.18426077,3.380619635,10.85817745,1.177295346,0.856886018,-.7814462052,8.960380402,3.400291759,6.978207048,4.202804277,7.874669961,4.256110662,-87.44125714],
# 	"P1, P2":               [164.0954218e15,100.86974049543558,-0.1588381968,3.168865761,10.47252618,.6080325039,0.856005883,-.7120849961,9.015436804,3.10824185,7.061600627,3.981685591,7.484156813,4.064857076,-93.12940659],
# 	"P1, P2, P4":           [88.71657896e15, 91.64947774226279,-0.5901396969,4.288099759,9.413961466,1.208433783,0.9520658982,-.7003857988,7.599155994,3.602491165,6.062004003,2.431564928,7.496507264,3.58146492,-58.44907688],
# 	"P1, extra P2, P4":     [30.41117331e15, 82.80405354958481,-0.7871359762,4.387469235,7.56604614,-.4128174512,1.20743159,-.5336200009,6.913116212,1.990664254,5.71167482,1.591588135,7.053548379,1.845198033,-43.27693469],
# 	"P1, P2, P4, burn off": [2.245892065e15,99.0216593136282,-0.4010667207,3.985329307,4.581542467,-1.662217588,1.329866331,.70173167,4.957506998,-0.519537713,5.174677867,0.240952608,2.825312798,-2.304634913,66.16427696],
# }
# ground_truth_max = {
# 	# based on the model (max)
# 	"Base case":            [421.3687017e15, 71.64003666274049,-0.2900906493,3.295152373,14.4519683,6.350179054,0.8569425781,-.778,9.970616694,6.236026801,7.365489006,5.221992371,8.477012035,5.540807458,-82.70139249],
# 	"P2":                   [174.6411255e15, 86.69302292763639,-0.18426077,3.380619635,10.85817745,1.177295346,0.856886018,-.7814462052,8.960380402,3.400291759,6.978207048,4.202804277,7.874669961,4.256110662,-87.44125714],
# 	"P1, P2":               [164.0954218e15,102.39040272449171,-0.1588381968,3.168865761,10.47252618,.6080325039,0.856005883,-.7120849961,9.015436804,3.10824185,7.061600627,3.981685591,7.484156813,4.064857076,-93.12940659],
# 	"P1, P2, P4":           [88.71657896e15, 93.21363880627495,-0.5901396969,4.288099759,9.413961466,1.208433783,0.9520658982,-.7003857988,7.599155994,3.602491165,6.062004003,2.431564928,7.496507264,3.58146492,-58.44907688],
# 	"P1, extra P2, P4":     [30.41117331e15, 86.64189694403553,-0.7871359762,4.387469235,7.56604614,-.4128174512,1.20743159,-.5336200009,6.913116212,1.990664254,5.71167482,1.591588135,7.053548379,1.845198033,-43.27693469],
# 	"P1, P2, P4, burn off": [2.245892065e15,103.96399027238204,-0.4010667207,3.985329307,4.581542467,-1.662217588,1.329866331,.70173167,4.957506998,-0.519537713,5.174677867,0.240952608,2.825312798,-2.304634913,66.16427696],
# }
# based on the model (seed 1)
short_header = ["total yield ()", "burn width (ps)", "burn skewness ()", "burn kurtosis ()",
                "Ti at BT (keV)", "dTi/dt at BT (keV/(100 ps))", "ρR at BT (g/cm^2)",
                "dρR/dt at BT (g/cm^2/(100 ps))", "Ti at BT-50ps (keV)", "dTi/dt at BT-50ps (keV/(100 ps))",
                "Ti at BT-100ps (keV)", "dTi/dt at BT-100ps (keV/(100 ps))", "Ti at compression (keV)",
                "dTi/dt at compression (keV/(100 ps))", "peak compression (ps)"]

# check that this is done correctly
for key in ground_truth_min:
	assert len(ground_truth_min[key]) == len(short_header)
	assert key in ground_truth_max
	assert len(ground_truth_max[key]) == len(short_header)

case_numbers = data["case"].unique()
case_yields = [-np.mean(data["total yield ()"][data["case"] == case]) for case in case_numbers]
sorted_cases = np.argsort(case_yields)
order = {}
for i, sim in enumerate(case_names):
	order[sim] = case_numbers[sorted_cases[i]]

requirements = {"total yield ()": -.05, "burn width (ps)": 7,
                "Burn skewness ()": .3, "burn kurtosis ()": 3,
                "Ti at BT (keV)": -.07, "dTi/dt at BT (keV/(100 ps))": 1.9,
                "Ti at stagnation (keV)": -.07, "dTi/dt at stagnation (keV/(100 ps))": 1.9,
                "ρR at BT (keV)": -.07, "dρR/dt at BT (g/cm^2/(100 ps))": .060
                }

COLORS = ["#E34B66", "#D3A94C", "#037F58", "#2AA4E9", "#88379D", "#594047"]

for axes in [
	# ("total yield ()", "burn width (ps)"),
	("burn skewness ()", "burn kurtosis ()"),
	("Ti at BT (keV)", "dTi/dt at BT (keV/(100 ps))"),
	# ("Ti at compression (keV)", "dTi/dt at compression (keV/(100 ps))"),
	# ("Ti at BT-50ps (keV)", "dTi/dt at BT-50ps (keV/(100 ps))"),
	# ("Ti at BT-100ps (keV)", "dTi/dt at BT-100ps (keV/(100 ps))"),
	("ρR at BT (g/cm^2)", "dρR/dt at BT (g/cm^2/(100 ps))"),
	# ("dρR/dt at BT (g/cm^2/(100 ps))", "dTi/dt at BT (keV/(100 ps))"),
	# ("total yield ()",),
	("burn width (ps)",),
	# ("burn skewness ()",),
	# ("burn kurtosis ()",),
	("peak compression (ps)",),
	# ("dTi/dt at BT (keV/(100 ps))",),
	# ("dρR/dt at BT (g/cm^2/(100 ps))",),
]:
	plt.figure()

	if len(axes) == 2:
		x, y = axes
		x_means, y_means = [], []
		for i, sim in enumerate(CASES):
			here = data["case"] == order[sim]
			plt.scatter(data[here][x], data[here][y], c=COLORS[i],
			            s=8, alpha=1, edgecolors="none", label=sim, zorder=10)
			x_means.append(np.mean(data[here][x]))
			y_means.append(np.mean(data[here][y]))
			if x in short_header and y in short_header:
				x_value = (ground_truth_min[sim][short_header.index(x)] + ground_truth_max[sim][short_header.index(x)])/2
				# x_requ = requirements[x] if requirements[x] > 0 else -requirements[x]*x_value
				y_value = (ground_truth_min[sim][short_header.index(y)] + ground_truth_max[sim][short_header.index(y)])/2
				# y_requ = requirements[y] if requirements[y] > 0 else -requirements[y]*y_value
				# if SIMPLE:
				# 	x_requ, y_requ = 0, 0
				x_variation = (ground_truth_max[sim][short_header.index(x)] - ground_truth_min[sim][short_header.index(x)])/2
				y_variation = (ground_truth_max[sim][short_header.index(y)] - ground_truth_min[sim][short_header.index(y)])/2
				# plt.errorbar(x=x_value, y=y_value, xerr=x_variation, yerr=y_variation,
				#              marker="o", markersize=6, markeredgewidth=2, markerfacecolor="white", markeredgecolor=COLORS[i],
				#              linewidth=1, zorder=20)
			else:
				print(f"warning: I'm missing the ground truth for {x} and/or {y}")
		plt.legend()
		if np.min(data[x]) > 0 and np.max(data[x]) > np.min(x_means)*50:
			plt.xscale("log")
		if np.min(data[y]) > 0 and np.max(data[y]) > np.min(y_means)*50:
			plt.yscale("log")
		plt.xlabel(x)
		plt.ylabel(y)
		output_filename = "../../output/compare_{}_{}_{}".format(
			re.sub(r'[/ ]', '', re.sub(r'\(.*\)', '', x)),
			re.sub(r'[/ ]', '', re.sub(r'\(.*\)', '', y)),
			filename[24:-19],
		)

	else:
		quantity, = axes
		for i, sim in enumerate(CASES):
			here = data["case"] == order[sim]
			values = data[here][quantity]
			bar_heights, bar_edges = np.histogram(values, bins=18, density=True)
			x = np.repeat(bar_edges, 2)
			y = np.concatenate([[0], np.repeat(bar_heights, 2), [0]])
			plt.fill_between(x, 0, y, color=COLORS[i], label=sim, alpha=0.3, zorder=2*i)
			plt.plot(x, y, COLORS[i], zorder=2*i+1)
			if quantity in short_header:
				x_value = (ground_truth_min[sim][short_header.index(quantity)] + ground_truth_max[sim][short_header.index(quantity)])/2
				# plt.axvline(x_value, color=COLORS[i], linestyle="dashed")
			else:
				print(f"warning: I'm missing the ground truth for {quantity}")
		plt.ylim(0, None)
		plt.legend()
		plt.xlabel(quantity)
		output_filename = "../../output/compare_{}_{}".format(
			re.sub(r'[/ ]', '', re.sub(r'\(.*\)', '', quantity)),
			filename[24:-19],
		)

	plt.grid()
	plt.title(filename[24:-19])
	plt.tight_layout()
	plt.savefig(output_filename.replace(".0", "p0"), dpi=300, transparent=True)

plt.show()
