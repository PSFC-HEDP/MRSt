import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import re
import sys
plt.rcParams.update({'font.family': 'sans', 'font.size': 11})
# plt.rcParams.update({'font.family': 'serif', 'font.size': 12})
import warnings
warnings.filterwarnings("ignore")

PLOT_ENVELOPE = False
PLOT_THEORETICAL_ERROR_BARS = False

# INCLUDE_ERRORS = False
# COLUMNS = 3
# SIZE = (12, 3.5/1)
# MARGIN = dict(bottom=.15, top=.97, left=.07, right=.99, wspace=.40, hspace=.05)

# INCLUDE_ERRORS = True
# COLUMNS = 2
# SIZE = (16, 7.5/4)
# MARGIN = dict(bottom=.07, top=.93, left=.06, right=.99, wspace=.35, hspace=.05)

# INCLUDE_ERRORS = True
# COLUMNS = 1
# SIZE = (8, 7/4)
# MARGIN = dict(bottom=.08, top=.92, left=.10, right=.97, wspace=.36, hspace=.05)

# INCLUDE_ERRORS = False
# COLUMNS = 2
# SIZE = (7.5, 9.0/5)
# MARGIN = dict(bottom=.06, top=.94, left=.12, right=.99, wspace=.30, hspace=.05)

# INCLUDE_ERRORS = False
# COLUMNS = 2
# SIZE = (10, 6/4)
# MARGIN = dict(bottom=.11, top=.89, left=.13, right=.99, wspace=.30, hspace=.05)

# INCLUDE_ERRORS = False
# COLUMNS = 1
# SIZE = (4.0, 5/3)
# MARGIN = dict(bottom=.09, top=.91, left=.18, right=.95, hspace=.05)

INCLUDE_ERRORS = False
INCLUDE_HISTOGRAMS = True
COLUMNS = 1
SIZE = (5.2, 5/3)
MARGIN = dict(bottom=.10, top=.90, left=.15, right=.98, hspace=.05, wspace=.04)

# INCLUDE_ERRORS = True
# COLUMNS = 1
# SIZE = (7, 5/3)
# MARGIN = dict(bottom=.10, top=.90, left=.11, right=.97, wspace=.37, hspace=.05)

assert not (INCLUDE_ERRORS and INCLUDE_HISTOGRAMS)


if len(sys.argv) <= 1:
	# FILENAME = '../../output/ensemble_high_2slits_400um_0c_15ps_200_2022-09-02.csv'
	FILENAME = '../../output/ensemble_medium_driftt_0c_15ps_200_2022-09-02.csv'
else:
	FILENAME = '../../output/'+sys.argv[1]
BIN_WIDTH = 0.3 # in bels
REFERENCE_YIELDS = [3e16, 3e17, 3e18]

X_LABEL = "Yield"

Y_LABELS = [
	("Burn width (ps)", 39, 56, 76, 7, False),
	# ("Burn skewness", -1.6, -.698, -0.1, 3e-1, False),
	# ("Burn kurtosis", -0.5, 4.7, 10.5, 3, False),
	("Ti at BT (keV)", 6.3, 7.56, 8.7, 5e-2, True),
	# ("Ti at stagnation (keV)", 3.8, 5.856, 7.2, 5e-2, True),
	("dTi/dt at BT (keV/100ps)", -2.2, 1.4, 4.4, 1.4, False),
	# ("ρR at BT (g/cm^2)", 0.77, 0.978, 1.13, 7e-2, True),
	# ("ρR at stagnation (g/cm^2)", 1.27, 1.416, 1.63, 7e-2, True),
	# ("dρR/dt at BT (g/cm^2/100ps)", -2.1, -1.1, 1.1, 0.95, False),
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

def hide_ticks(axis):
	"""
		hide the ticks on an axis without hiding the entire
		axis (this way the grid won't go away)
	"""
	for tick in axis.get_major_ticks() + axis.get_minor_ticks():
		tick.tick1line.set_visible(False)
		tick.tick2line.set_visible(False)
		tick.label1.set_visible(False)
		tick.label2.set_visible(False)

def rolling_average(y, n, bessel_correction=False):
	""" rolling average of y, using n neighbors in each direction """
	if n > (y.size-1)//2:
		raise IndexError("n is too big")
	c = -1 if bessel_correction else 0
	Σy = np.concatenate([[0], np.cumsum(y)])
	return np.concatenate([
		Σy[n+1:2*n+1]/(np.arange(n+1, 2*n+1) + c),
		(Σy[2*n+1:] - Σy[:-2*n-1])/(2*n+1 + c),
		(Σy[-1] - Σy[-2*n-1:-n-1])/(np.arange(2*n, n, -1) + c),
	])

def smooth_average(y, bessel_correction=False):
	""" double rolling average that automatically handles nans """
	valid = np.isfinite(y)
	output = y[:]
	output[valid] = rolling_average(
		rolling_average(
			y[valid],
			n=min((np.sum(valid)-1)//2, 36),
			bessel_correction=bessel_correction),
		n=min((np.sum(valid)-1)//2, 36))
	output[~valid] = (output[np.roll(~valid, 1)] + output[np.roll(~valid, -1)])/2
	return output

def symlog(x):
	x0 = np.median(x)
	return np.where(x > x0, x0*(np.log(x/x0) + 1), x)

def symexp(y):
	y0 = np.median(y)
	return np.where(y > y0, y0*np.exp(y/y0 - 1), x)


try:
	simulations = pd.read_csv(FILENAME, na_values=["Infinity"])
except FileNotFoundError:
	print("Ese archivo no existe.")
	quit()

if "Yield factor" in simulations:
	simulations["Yield"] = simulations["Yield factor"]*4.7838e17
simulations = simulations[simulations["Yield"] != 0]
for parameter in simulations:
	if '(ns)' in parameter:
		simulations[parameter.replace('(ns)', '(ps)')] = simulations[parameter]*1e-9/1e-12
	if '/ns' in parameter:
		simulations[parameter.replace('/ns', '/100ps')] = simulations[parameter]/1e-9*100e-12
	if ' (10^15)' in parameter:
		simulations[parameter.replace(' (10^15)', '')] = simulations[parameter]*1e15
for parameter in simulations:
	if '(g/' in parameter:
		simulations[parameter.replace('(g/', '(mg/')] = simulations[parameter]/1e-3

num_rows = (len(Y_LABELS) + COLUMNS - 1)//COLUMNS
if INCLUDE_ERRORS:
	num_columns = COLUMNS*2
	gridspec_keywords = dict()
elif INCLUDE_HISTOGRAMS:
	num_columns = COLUMNS*2
	gridspec_keywords = dict(width_ratios=[2, 1]*COLUMNS)
else:
	num_columns = COLUMNS
	gridspec_keywords = dict()
fig, axs = plt.subplots(num_rows, num_columns, figsize=(SIZE[0], SIZE[1]*num_rows), gridspec_kw=gridspec_keywords, squeeze=False)

fig.subplots_adjust(**MARGIN)

print("               reference yields:", end='')
for reference_yield in REFERENCE_YIELDS:
	print(f"  {reference_yield:< 17.2e}", end='')
print()

# first do the scatter plot
for i, (axis, y_min, y_original, y_max, presis, percent) in enumerate(Y_LABELS): # iterate through each desired plot
	if axis is None: continue
	if INCLUDE_ERRORS:
		ax = axs[i//COLUMNS, i%COLUMNS * 2 + 0] # identify the corresponding axis object
	else:
		ax = axs[i//COLUMNS, i%COLUMNS]

	x = simulations[X_LABEL].values # and the corresponding data
	y = simulations[axis].values
	ɛ = simulations[axis+" error"].values

	ax.set_xlim(2e16, 5e18) # set up the x axis
	if 'ield' in X_LABEL:
		ax.set_xscale('log')
	if i//COLUMNS == axs.shape[0]-1:
		ax.set_xlabel(X_LABEL)
	elif i//COLUMNS == 0:
		ax.set_xlabel(X_LABEL)
		ax.xaxis.set_label_position('top')
		ax.xaxis.tick_top()
	else:
		hide_ticks(ax.xaxis)
	ax.grid()

	if axis is None: # add a legend if there is space
		ax.plot([], [], 'C1--', label="Original data")
		ax.scatter([], [], label="Fit to synthetic data")
		ax.legend()
		ax.yaxis.set_visible(False)
		continue

	if 'yield' in axis:  y_true = simulations["Yield"].values # automatically determine the relevant scaling factor
	else:                y_true = y_original*np.ones(len(simulations.index))

	order = np.argsort(x) # get some useful indices of the data
	order = order[np.isfinite(simulations["Total yield"].values[order])]

	μ = smooth_average(y[order])
	σ = np.sqrt(smooth_average((y[order] - μ)**2, bessel_correction=False))

	if presis is not None:
		if not percent:
			shaded_region = (y_true[order] - presis, y_true[order] + presis)
		else:
			shaded_region = (y_true[order]*(1 - presis), y_true[order]*(1 + presis))
		ax.fill_between(x[order], shaded_region[0], shaded_region[1], color='#F7DFC8')
	else:
		shaded_region = (y_true[order], y_true[order])

	ax.plot(x[order], y_true[order], 'C1-', zorder=1, label="Based on original data")
	ax.scatter(x[order], y[order], s=1, zorder=2, label="Based on fit to synthetic data")
	if PLOT_ENVELOPE:
		ax.plot(x[order], μ + σ, 'C0-', linewidth=1, zorder=1, label="1σ variation")
		ax.plot(x[order], μ - σ, 'C0-', linewidth=1, zorder=1)
	if y_min > 0 and y_max/y_min > 10:
		ax.set_yscale('log')
	ax.set_ylim(y_min, y_max)
	ax.set_ylabel(text_wrap(axis.replace("^2", "²").replace("Ti", "Tᵢ")))

	# set up the histogram axes
	if INCLUDE_HISTOGRAMS:
		ax = axs[i//COLUMNS, i%COLUMNS * 2 + 1]
		ax.set_xlim(0, len(REFERENCE_YIELDS))
		x_tick_labels = [re.sub(r'([0-9.]+)e\+([0-9]+)', r"$\1\\times10^{\2}$", f"{Yn:.1g}") for Yn in REFERENCE_YIELDS]
		ax.set_ylim(y_min, y_max)
		ax.axhline(y_true[0], color='C1', zorder=1.5)
		ax.fill_between(np.linspace(0, len(REFERENCE_YIELDS) + 1, order.size),
		                shaded_region[0], shaded_region[1], color='#F7DFC8')
		if i//COLUMNS == axs.shape[0]-1:
			ax.set_xticks(np.arange(len(REFERENCE_YIELDS)), labels=x_tick_labels, rotation=-30, ha="left", rotation_mode='anchor')
			ax.tick_params(axis='x', which='major', pad=-0)
		elif i//COLUMNS == 0:
			ax.set_xticks(np.arange(len(REFERENCE_YIELDS)), labels=x_tick_labels, rotation=30, ha="left", rotation_mode='anchor')
			ax.tick_params(axis='x', which='major', pad=-0)
			ax.xaxis.set_label_position('top')
			ax.xaxis.tick_top()
		else:
			hide_ticks(ax.xaxis)
		hide_ticks(ax.yaxis)
		ax.grid()

	dy = y - y_true
	print(f"{axis:>25s} error:", end='')
	for j, reference_yield in enumerate(REFERENCE_YIELDS):
		at_yield = (~np.isnan(y)) & (np.absolute(np.log10(x/reference_yield)) <= 0.15)
		error_at_yield = np.sqrt(np.mean(np.square(dy[at_yield])))
		number_at_yield = np.sum(at_yield)
		print(f"  {error_at_yield:7.4f} ± {error_at_yield*np.sqrt(2/number_at_yield):7.4f}", end='')

		if INCLUDE_HISTOGRAMS:
			_, _, patches = ax.hist(
				y[at_yield], bottom=j, color="C0", zorder=2,
				bins=np.linspace(max(y_min, min(shaded_region[0].min(), y[at_yield].min())),
				                 min(y_max, max(shaded_region[1].max(), y[at_yield].max())),
				                 12),
				orientation="horizontal")
			peak = np.max(patches.datavalues)
			for rectangle in patches.patches:
				rectangle.set_width(rectangle.get_width()/peak*0.8)
			deviation = np.std(y[at_yield])
			mean = np.mean(y[at_yield])
			for sign in [-1, 1]:
				ax.arrow(j + .85, mean,
				         0, sign*deviation,
				         length_includes_head=True, zorder=2,
				         head_width=.2, head_length=0, color="#000")
			# ax.plot([j + 0.50, j + 0.50], [mean - deviation, mean + deviation], "k", zorder=5)
			# ax.plot([j + 0.42, j + 0.50, j + 0.58], [mean - deviation + arrow_hite, mean - deviation, mean - deviation + arrow_hite], "k", zorder=5)
			# ax.plot([j + 0.42, j + 0.50, j + 0.58], [mean + deviation - arrow_hite, mean + deviation, mean + deviation - arrow_hite], "k", zorder=5)

	print()

	if INCLUDE_ERRORS:
		ax = axs[i//COLUMNS, i%COLUMNS * 2 + 1]

		ax.set_xlim(10**round(np.log10(x.min())),
		            10**round(np.log10(x.max()))) # set up the x axis
		if 'ield' in X_LABEL:
			ax.set_xscale('log')
		if i//COLUMNS == axs.shape[0]-1:
			ax.set_xlabel(X_LABEL)
		elif i//COLUMNS == 0:
			ax.set_xlabel(X_LABEL)
			ax.xaxis.set_label_position('top')
			ax.xaxis.tick_top()
		else:
			hide_ticks(ax.xaxis)
		ax.grid()

		if axis is None:
			ax.plot([], [], 'C1--', label="Required accuracy")
			ax.plot([], [], 'C0-', label="Standard deviation from actuality")
			if PLOT_THEORETICAL_ERROR_BARS:
				ax.plot([], [], 'C3--', label="Reported error bar size")
			ax.legend()
			ax.yaxis.set_visible(False)
			continue

		stds = np.sqrt(smooth_average(np.square(dy[order]), bessel_correction=True))
		errs = smooth_average(ɛ[order])
		ax.plot(x[order], stds, 'C0-', label="Standard deviation from actuality")
		if presis is not None:
			if percent:
				presis *= y_true[order]
			ax.plot(x[order], np.full(order.shape, presis), 'C1--', label="Required accuracy")
		if PLOT_THEORETICAL_ERROR_BARS:
			ax.plot(x[order], errs, 'C3--', label="Reported error bar size")
		ax.set_yscale('log')

		# figure out the best y limits
		y_min = np.min(presis)*1e-1
		y_max = np.max(presis)*1e+1
		if np.median(stds) < y_min:
			y_max *= np.median(stds)/y_min
			y_min = np.median(stds)
		if np.min(presis) > y_max:
			y_min *= np.min(presis)/y_max
			y_max = np.min(presis)
		ax.set_ylim(y_min, y_max)

		# ax.grid(which='major', axis='y')
		if 'ield' in X_LABEL:
			ax.set_xscale('log')
		if "(" in axis:
			ax.set_ylabel(text_wrap(axis[:axis.index("(")].replace(" at BT", "").replace("^2", "²") + "error " + axis[axis.index("("):]))
		else:
			ax.set_ylabel(text_wrap(axis+" error"))

filename = os.path.splitext(FILENAME)[0]
fig.savefig(f'{filename}.eps', dpi=300)
fig.savefig(f'{filename}.png', dpi=300)

plt.show()
