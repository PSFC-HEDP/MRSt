import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
plt.rcParams.update({'font.family': 'serif', 'font.size': 10})

COLUMNS = 2
SIZE = (17.0, 10.5)
MARGIN = dict(bottom=.06, top=.94, left=.06, right=.99, wspace=.36, hspace=.05)
BIN_WIDTH = 0.3 # in bels
# COLUMNS = 3
# SIZE = (16, 9)
# MARGIN = dict(bottom=.06, top=.94, left=.06, right=.99, wspace=.26, hspace=.05)
FILENAME = '../../working/ensemble_1_2_2.0_2_1000_2020-09-06.csv'

X_LABEL = "True yield"

Y_LABELS = [
	(None, 0, 0, 0, 0, False), ("Total yield", 2e14, 4.7838e17, 9e17, 5e-2, True),
	("Burn-average ρR (g/cm^2)", 0.3, .68998, 1.1, 7e-2, True), ("Burn-average Ti (keV)", 7.8, 10.0186, 12.2, 7e-2, True),
	("Bang time (ns)", 16.339, 16.3632, 16.391, 1e-2, False),	("Burn width (ps)", 47, 66.711, 83, 7, False),
	("Burn skew", -2.1, -1.1572, 0.11, 3e-1, False), ("Burn kurtosis", -0.5, 6.5079, 15.5, 3, False),
	("dρR/dt at BT (mg/cm^2/(100ps))", -850, -196.9, 450, 60, False), ("dTi/dt at BT (keV/(100ps))", -0.2, 6.733, 12.2, 1.9, False),
	("Burn-average vi (km/s)", -15.2, 1.83, 15.2, 20, False), ("dvi/dt at BT (km/s/(100ps))", -170, -87.3, 20, 8, False)
]


def text_wrap(s):
	if len(s) > 14:
		i = len(s)//2
		for j in range(i):
			if s[i+j] == ' ':
				return s[:i+j] + '\n' + s[i+j+1:]
			elif s[i-j] == ' ':
				return s[:i-j] + '\n' + s[i-j+1:]
	return s


simulations = pd.read_csv(FILENAME, na_values=["Infinity"])
simulations = simulations[simulations["Yield factor"] != 0]
simulations["True yield"] = simulations["Yield factor"]*Y_LABELS[1][2]
for suf in ["", " error"]:
	simulations["Total yield"+suf] = simulations["Total yield (10^15)"+suf]*1e15
	simulations["Bang time (ps)"+suf] = simulations["Bang time (ns)"+suf]/1e-3
	simulations["Burn width (ps)"+suf] = simulations["Burn width (ns)"+suf]/1e-3
	simulations["dρR/dt at BT (mg/cm^2/(100ps))"+suf] = simulations["dρR/dt at BT (g/cm^2/ns)"+suf]/1e-2
	simulations["dTi/dt at BT (keV/(100ps))"+suf] = simulations["dTi/dt at BT (keV/ns)"+suf]/1e1
	simulations["dvi/dt at BT (km/s/(100ps))"+suf] = simulations["dvi/dt at BT (km/s/ns)"+suf]/1e1

fig, axs = plt.subplots((len(Y_LABELS) + COLUMNS-1)//COLUMNS, COLUMNS*2, figsize=SIZE)
fig.subplots_adjust(**MARGIN)
bins = np.geomspace(simulations[X_LABEL].min(), simulations[X_LABEL].max(), max(3, 1 + int(np.ptp(np.log10(simulations[X_LABEL]))/BIN_WIDTH)))

for i, (axis, y_min, y_true, y_max, presis, percent) in enumerate(Y_LABELS):
	ax = axs[i//COLUMNS, i%COLUMNS * 2 + 0]

	ax.set_xlim(simulations[X_LABEL].min()/1.2, simulations[X_LABEL].max()*1.2)
	if 'ield' in X_LABEL:
		ax.set_xscale('log')
	if i//COLUMNS == axs.shape[0]-1:
		ax.set_xlabel(X_LABEL)
	elif i//COLUMNS == 0:
		ax.set_xlabel(X_LABEL)
		ax.xaxis.set_label_position('top')
		ax.xaxis.tick_top()
	else:
		ax.xaxis.set_visible(False)

	if axis is None:
		ax.plot([], [], 'C1--', label="Original data")
		ax.scatter([], [], label="Fit to synthetic data")
		ax.legend()
		ax.yaxis.set_visible(False)
		continue

	if 'keV' in axis:      y_factor = simulations["Temperature factor"]
	elif 'g/cm^2' in axis: y_factor = simulations["Down-scatter factor"]
	elif 'yield' in axis:  y_factor = simulations["Yield factor"]
	else:                  y_factor = np.ones(len(simulations.index))
	order = np.argsort(simulations[X_LABEL])
	order = order[np.isfinite(simulations["Total yield"].values[order])].values

	ax = axs[i//COLUMNS, i%COLUMNS * 2 + 0]
	valid = np.logical_not(np.isnan(simulations[axis+" error"]))
	if not percent:
		ax.fill_between(simulations[X_LABEL][order],
			y_factor[order]*y_true - presis, y_factor[order]*y_true + presis, color='C1', alpha=0.3)
	else:
		ax.fill_between(simulations[X_LABEL][order],
			y_factor[order]*y_true*(1 - presis), y_factor[order]*y_true*(1 + presis), color='C1', alpha=0.3)
	ax.plot(simulations[X_LABEL][order], y_factor[order]*y_true, 'C1--', zorder=0, label="Based on original data")
	ax.scatter(simulations[X_LABEL][valid], simulations[axis][valid], s=2, zorder=1, label="Based on fit to synthetic data")
	# ax.errorbar(simulations[X_LABEL][valid], simulations[axis][valid], yerr=simulations[axis+" error"][valid], elinewidth=1, linestyle='none')
	if y_min > 0 and y_max/y_min >= 10:
		ax.set_yscale('log')
	ax.set_ylim(y_min, y_max)
	ax.set_ylabel(text_wrap(axis))

for i, (axis, y_min, y_true, y_max, presis, percent) in enumerate(Y_LABELS):
	ax = axs[i//COLUMNS, i%COLUMNS * 2 + 1]

	ax.set_xlim(simulations[X_LABEL].min()/1.2, simulations[X_LABEL].max()*1.2)
	if 'ield' in X_LABEL:
		ax.set_xscale('log')
	if i//COLUMNS == axs.shape[0]-1:
		ax.set_xlabel(X_LABEL)
	elif i//COLUMNS == 0:
		ax.set_xlabel(X_LABEL)
		ax.xaxis.set_label_position('top')
		ax.xaxis.tick_top()
	else:
		ax.xaxis.set_visible(False)

	if axis is None:
		ax.plot([], [], 'C1--', label="Required accuracy")
		ax.plot([], [], 'C0-', label="Standard deviation from actuality")
		ax.plot([], [], 'C2--', label="Reported error bar size")
		ax.legend()
		ax.yaxis.set_visible(False)
		continue

	if 'keV' in axis:      y_factor = simulations["Temperature factor"]
	elif 'g/cm^2' in axis: y_factor = simulations["Down-scatter factor"]
	elif 'yield' in axis:  y_factor = simulations["Yield factor"]
	else:                  y_factor = np.ones(len(simulations.index))

	if presis <= 1e-2 and "(n" in axis:
		axis = axis.replace("(n", "(p")
		presis *= 1e3
		y_true *= 1e3

	stds, errs = [], []
	for j in range(1, len(bins)):
		stds.append(
			np.sqrt(np.mean(np.square((simulations[axis] - y_factor*y_true)[(simulations[X_LABEL] >= bins[j-1]) & (simulations[X_LABEL] < bins[j])]))))
		errs.append(
			np.mean(simulations[axis+" error"][(simulations[X_LABEL] >= bins[j-1]) & (simulations[X_LABEL] < bins[j])]))
	if "/dt" in axis:
		print(f"σ = {stds[-1]}")
	bin_centers = np.sqrt(bins[1:]*bins[:-1])
	y_factor = np.interp(bin_centers, simulations[X_LABEL][order], y_factor[order])
	if not percent:
		ax.plot(bin_centers, stds, 'C0-', label="Standard deviation from actuality")
		ax.plot(bin_centers, y_factor*presis, 'C1--', label="Required accuracy")
		ax.plot(bin_centers, errs, 'C2--', label="Reported error bar size")
	else:
		ax.plot(bin_centers, stds/(y_factor*y_true), 'C0-', label="Standard deviation from actuality")
		ax.plot(bin_centers, np.full(bin_centers.shape, presis), 'C1--', label="Required accuracy")
		ax.plot(bin_centers, errs/(y_factor*y_true), 'C2--', label="Reported error bar size")
	ax.set_yscale('log')
	if 'ield' in X_LABEL:
		ax.set_xscale('log')
	ax.set_ylim(presis*3.2e-2, presis*5)
	if not percent:
		if "(" in axis:
			ax.set_ylabel(text_wrap(axis[:axis.index("(")].replace(" at BT", "") + "error " + axis[axis.index("("):]))
		else:
			ax.set_ylabel(text_wrap(axis+" error"))
	else:
		ax.set_ylabel(text_wrap(re.sub(r'(\([^)]+\))?$', '', axis) + " error"))

config = '-'+FILENAME[23:38] if len(FILENAME) > 23 else ''
fig.savefig('../../working/mrst{}.eps'.format(config))
fig.savefig('../../working/mrst{}.png'.format(config))

plt.show()
