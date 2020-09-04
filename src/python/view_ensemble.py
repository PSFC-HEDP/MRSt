import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
plt.rcParams.update({'font.family': 'serif', 'font.size': 10})

X_LABEL = "Yield factor"

Y_LABELS = [
	("Total yield", 2e14, 4.79839e17, 9e18, 5e-2, True), ("Total yield", 2e14, 4.7838e17, 9e18, 5e-2, True),
	("Burn-average ρR (g/cm^2)", 0.3, .68998, 1.1, 7e-2, True), ("Burn-average Ti (keV)", 8.8, 10.0186, 11.2, 7e-2, True),
	("Bang time (ns)", 16.349, 16.3632, 16.381, 1e-2, False),	("Burn width (ps)", 58, 66.711, 77, 7, False),
	("Burn skew", -2.1, -1.1572, 0.11, 3e-1, False), ("Burn kurtosis", 2.8, 6.5079, 10.2, 3, False),
	("dρR/dt at BT (mg/cm^2/(100ps))", -850, -196.9, 450, 60, False), ("dTi/dt at BT (keV/(100ps))", -0.2, 6.733, 12.2, 1.9, False),
	("Burn-average vi (km/s)", -20.2, 1.83, 20.2, 10, False), ("dvi/dt at BT (km/s/(100ps))", -170, -87.3, 20, 8, False)
]

COLUMNS = 2
SIZE = (8.0, 10.5)
MARGIN = dict(bottom=.06, top=.94, left=.12, right=.99, wspace=.35, hspace=.05)
BIN_WIDTH = 0.3 # in bels
# COLUMNS = 3
# SIZE = (16, 9)
# MARGIN = dict(bottom=.06, top=.94, left=.06, right=.99, wspace=.26, hspace=.05)
FILENAME = '../../working/ensemble_6_10_6.0_2_1000_2020-09-02.csv'


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
for suf in ["", " error"]:
	simulations["Total yield"+suf] = simulations["Total yield (10^15)"+suf]*1e15
	simulations["Burn width (ps)"+suf] = simulations["Burn width (ns)"+suf]/1e-3
	simulations["dρR/dt at BT (mg/cm^2/(100ps))"+suf] = simulations["dρR/dt at BT (g/cm^2/ns)"+suf]/1e-2
	simulations["dTi/dt at BT (keV/(100ps))"+suf] = simulations["dTi/dt at BT (keV/ns)"+suf]/1e1
	simulations["dvi/dt at BT (km/s/(100ps))"+suf] = simulations["dvi/dt at BT (km/s/ns)"+suf]/1e1

fig_p, axs_p = plt.subplots((len(Y_LABELS) + COLUMNS-1)//COLUMNS, COLUMNS, figsize=SIZE)
fig_p.subplots_adjust(**MARGIN)

bins = np.geomspace(simulations[X_LABEL].min(), simulations[X_LABEL].max(), max(3, 1 + int(np.ptp(np.log10(simulations[X_LABEL]))/BIN_WIDTH)))
fig_w, axs_w = plt.subplots((len(Y_LABELS) + COLUMNS-1)//COLUMNS, COLUMNS, figsize=SIZE)
fig_w.subplots_adjust(**MARGIN)

for i, (axis, y_min, y_true, y_max, presis, percent) in enumerate(Y_LABELS):
	if 'keV' in axis:      y_factor = simulations["Temperature factor"]
	elif 'g/cm^2' in axis: y_factor = simulations["Down-scatter factor"]
	elif 'yield' in axis:  y_factor = simulations["Yield factor"]
	else:                  y_factor = np.ones(len(simulations.index))
	order = np.argsort(simulations[X_LABEL])
	order = order[np.isfinite(simulations["Total yield"].values[order])].values

	ax = axs_p[i//COLUMNS,i%COLUMNS]
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
	if 'ield' in X_LABEL:
		ax.set_xscale('log')
	ax.set_ylim(y_min, y_max)
	if i == 0:
		ax.legend()
	if i//COLUMNS == axs_p.shape[0]-1:
		ax.set_xlabel(X_LABEL)
	elif i//COLUMNS == 0:
		ax.set_xlabel(X_LABEL)
		ax.xaxis.set_label_position('top')
		ax.xaxis.tick_top()
	else:
		ax.xaxis.set_visible(False)
	ax.set_ylabel(text_wrap(axis))

	ax = axs_w[i//COLUMNS,i%COLUMNS]
	stds, errs = [], []
	for j in range(1, len(bins)):
		stds.append(
			np.sqrt(np.mean(np.square((simulations[axis] - y_factor*y_true)[(simulations[X_LABEL] >= bins[j-1]) & (simulations[X_LABEL] < bins[j])]))))
		errs.append(
			np.mean(simulations[axis+" error"][(simulations[X_LABEL] >= bins[j-1]) & (simulations[X_LABEL] < bins[j])]))
	if "/dt" in axis:
		print(f"σ = {stds[-1]}")
	bin_centers = np.sqrt(bins[1:]*bins[:-1])
	y_factor = np.interp(bin_centers, simulations[X_LABEL], y_factor)
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
	if i == 0:
		ax.legend()
	if i//COLUMNS == axs_p.shape[0]-1:
		ax.set_xlabel(X_LABEL)
	elif i//COLUMNS == 0:
		ax.set_xlabel(X_LABEL)
		ax.xaxis.set_label_position('top')
		ax.xaxis.tick_top()
	else:
		ax.xaxis.set_visible(False)
	if not percent:
		ax.set_ylabel(text_wrap(axis))
	else:
		ax.set_ylabel(text_wrap(re.sub(r'(\([^)]+\))?$', '', axis)))

config = '-'+FILENAME[23:38] if len(FILENAME) > 23 else ''
fig_p.savefig('../../working/mrst-scatter{}.eps'.format(config))
fig_p.savefig('../../working/mrst-scatter{}.png'.format(config))
fig_w.savefig('../../working/mrst-errors{}.eps'.format(config))
fig_w.savefig('../../working/mrst-errors{}.png'.format(config))

plt.show()
