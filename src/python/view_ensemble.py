import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams.update({'font.family': 'serif', 'font.size': 10})

X_LABEL = "Yield factor"

Y_LABELS = [
	("Total yield", 2e14, 4.7838e17, 9e18, 1e13), ("Total yield", 2e14, 4.7838e17, 9e18, 1e13),
	("Burn-average ρR (g/cm^2)", 0.3, .6900, 1.1, 3e-3), ("Burn-average Ti (keV)", 8.8, 10.0198, 11.2, 3e-2),
	("Bang time (ns)", 16.349, 16.3645, 16.371, 3e-4),	("Burn width (ps)", 58, 66.732, 77, 1e-1),
	("Burn skew", -2.1, -1.1469, 0.11, 1e-2), ("Burn kurtosis", 3.8, 6.5166, 9.2, 1e-1),
	("dρR/dt at BT (mg/cm^2/(100ps))", -850, -282.4, 450, 3e1), ("dTi/dt at BT (keV/(100ps))", -0.2, 6.467, 12.2, 3e-1),
	("Burn-average vi (km/s)", -20.2, 1.734, 20.2, 3e-1), ("dvi/dt at BT (km/s/(100ps))", -170, -72.19, 20, 3e0)
]

COLUMNS = 2
SIZE = (8.0, 10.5)
MARGIN = dict(bottom=.06, top=.94, left=.12, right=.99, wspace=.35, hspace=.05)
N_BINS = 10
# COLUMNS = 3
# SIZE = (16, 9)
# MARGIN = dict(bottom=.06, top=.94, left=.06, right=.99, wspace=.26, hspace=.05)
FILENAME = '../../working/ensemble_t_1000_2020-08-24.csv'


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

bins = np.geomspace(simulations[X_LABEL].min(), simulations[X_LABEL].max(), N_BINS+1)
fig_w, axs_w = plt.subplots((len(Y_LABELS) + COLUMNS-1)//COLUMNS, COLUMNS, figsize=SIZE)
fig_w.subplots_adjust(**MARGIN)

for i, (axis, y_min, y_true, y_max, presis) in enumerate(Y_LABELS):
	if 'keV' in axis:      yFactor = simulations["Temperature factor"]
	elif 'g/cm^2' in axis: yFactor = simulations["Down-scatter factor"]
	elif 'yield' in axis:  yFactor = simulations["Yield factor"]
	else:                  yFactor = np.ones(len(simulations.index))
	order = np.argsort(simulations[X_LABEL])
	order = order[np.isfinite(simulations["Total yield"].values[order])].values

	ax = axs_p[i//COLUMNS,i%COLUMNS]
	valid = np.logical_not(np.isnan(simulations[axis+" error"]))
	ax.plot(simulations[X_LABEL][order], yFactor[order]*y_true, 'C1--', zorder=0, label="Based on original data")
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
	for j in range(N_BINS):
		stds.append(
			np.sqrt(np.mean(np.square((simulations[axis] - yFactor*y_true)[(simulations[X_LABEL] >= bins[j]) & (simulations[X_LABEL] < bins[j+1])]))))
		errs.append(
			np.mean(simulations[axis+" error"][(simulations[X_LABEL] >= bins[j]) & (simulations[X_LABEL] < bins[j+1])]))
	ax.plot(np.sqrt(bins[1:]*bins[:-1]), stds, 'C0-', label="Standard deviation from actuality")
	ax.plot(np.sqrt(bins[1:]*bins[:-1]), errs, 'C1--', label="Reported error bar size")
	ax.set_yscale('log')
	if 'ield' in X_LABEL:
		ax.set_xscale('log')
	ax.set_ylim(presis*0.8, presis*1.25e2)
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

config = '-'+FILENAME[23] if len(FILENAME) > 23 else ''
fig_p.savefig('../../working/mrst-scatter{}.eps'.format(config))
fig_p.savefig('../../working/mrst-scatter{}.png'.format(config))
fig_w.savefig('../../working/mrst-errors{}.eps'.format(config))
fig_w.savefig('../../working/mrst-errors{}.png'.format(config))

plt.show()
