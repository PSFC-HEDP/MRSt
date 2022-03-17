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

# INCLUDE_ERRORS = False
# COLUMNS = 3
# SIZE = (12, 3.5)
# MARGIN = dict(bottom=.15, top=.97, left=.07, right=.99, wspace=.40, hspace=.05)
# INCLUDE_ERRORS = True
# COLUMNS = 2
# SIZE = (16, 7.5)
# MARGIN = dict(bottom=.07, top=.93, left=.06, right=.99, wspace=.35, hspace=.05)
INCLUDE_ERRORS = True
COLUMNS = 1
SIZE = (8, 5)
MARGIN = dict(bottom=.10, top=.90, left=.13, right=.99, wspace=.35, hspace=.05)
# INCLUDE_ERRORS = False
# COLUMNS = 2
# SIZE = (7.5, 9.0)
# MARGIN = dict(bottom=.06, top=.94, left=.12, right=.99, wspace=.30, hspace=.05)
# INCLUDE_ERRORS = False
# COLUMNS = 2
# SIZE = (10, 6)
# MARGIN = dict(bottom=.11, top=.89, left=.13, right=.99, wspace=.30, hspace=.05)
# INCLUDE_ERRORS = False
# COLUMNS = 1
# SIZE = (4.5, 6)
# MARGIN = dict(bottom=.08, top=.92, left=.19, right=.99, wspace=.25, hspace=.05)
# INCLUDE_ERRORS = False
# COLUMNS = 1
# SIZE = (6, 4)
# MARGIN = dict(bottom=.15, top=.95, left=.15, right=.95)


if len(sys.argv) <= 1:
	FILENAME = '../../output/ensemble_high_30ps_2022-03-16.csv'
	# FILENAME = '../../output with varying response/ensemble_5_0_1000_og_2022-02-09.csv'
else:
	FILENAME = '../../output/'+sys.argv[1]
BIN_WIDTH = 0.3 # in bels
REFERENCE_YIELDS = [1e16, 1e17, 1e18, 1e19]

X_LABEL = "Yield"

Y_LABELS = [
	# ("Burn-average Ti (keV)", 5.25, 7.14320, 9.75, 5e-2, True),
	# ("Ti at BT (keV)", 5.25, 7.51554, 9.75, 5e-2, True),
	# ("dTi/dt at BT (keV/100ps)", -3.5, 1.09569, 3.5, 1.9, False),
	# ("Burn width (ps)", 56, 69.5532, 87, 7, False),
	# ("Burn-average Ti (keV)", 5.25, 7.52, 9.75, 5e-2, True),
	("Burn width (ps)", 56, 68, 87, 7, False),
	("Ti at BT (keV)", 5.25, 7.6, 9.75, 5e-2, True),
	("dTi/dt at BT (keV/100ps)", -5.2, 1.4, 5.2, 1.5, False),
	("dρR/dt at BT (mg/cm^2/100ps)", -5.2, 1.4, 5.2, 1.5, False),
	("ρR at BT (mg/cm^2)", -5.2, 1.4, 5.2, 1.5, False),
	# ("Burn-average Ti (keV)", 3.25, 4.4, 5.75, 5e-2, True),
	# ("Ti at BT (keV)", 2.25, 4.4, 6.75, 5e-2, True),
	# ("dTi/dt at BT (keV/100ps)", -3.5, -1.6, 3.5, 1.2, False),
	# ("Burn width (ps)", 100, 126, 150, 7, False),
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
			n=min((np.sum(valid)-1)//2, 72),
			bessel_correction=bessel_correction),
		n=min((np.sum(valid)-1)//2, 72))
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

if INCLUDE_ERRORS:
	fig, axs = plt.subplots((len(Y_LABELS) + COLUMNS-1)//COLUMNS, COLUMNS*2, figsize=SIZE)
else:
	fig, axs = plt.subplots((len(Y_LABELS) + COLUMNS-1)//COLUMNS, COLUMNS, figsize=SIZE)
try:
	if len(axs.shape) == 1: # force the axis matrix to be 2D since Matplotlib apparently automatically reduces any array where one of the dimensions is 1
		axs = axs[:, np.newaxis]
except AttributeError:
	axs = np.array([[axs]])
# if axs.shape[1] != COLUMNS: # and force it to be correct since Matplotlib apparently automatically transposes any array with one row‽‽
# 	axs = axs.T

fig.subplots_adjust(**MARGIN)

for i, (axis, y_min, y_original, y_max, presis, percent) in enumerate(Y_LABELS): # iterate through each desired plot
	if axis is None: continue
	if INCLUDE_ERRORS:
		ax = axs[i//COLUMNS, i%COLUMNS * 2 + 0] # identify the corresponding axis object
	else:
		ax = axs[i//COLUMNS, i%COLUMNS]

	x = simulations[X_LABEL].values # and the corresponding data
	y = simulations[axis].values
	ɛ = simulations[axis+" error"].values

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
		ax.xaxis.set_visible(False)

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

	if not percent: # plot the actual stuff
		ax.fill_between(x[order],
			y_true[order] - presis, y_true[order] + presis, color='#F7DFC8')
	else:
		ax.fill_between(x[order],
			y_true[order]*(1 - presis), y_true[order]*(1 + presis), color='#F7DFC8')
	ax.plot(x[order], y_true[order], 'C1-', zorder=1, label="Based on original data")
	ax.scatter(x[order], y[order], s=1, c=np.arange(x.size)[order], zorder=2, label="Based on fit to synthetic data")
	# ax.plot(x[order], μ + σ, 'C0-', linewidth=1, zorder=1, label="1σ variation")
	# ax.plot(x[order], μ - σ, 'C0-', linewidth=1, zorder=1)
	if y_min > 0 and y_max/y_min > 10:
		ax.set_yscale('log')
	ax.set_ylim(y_min, y_max)
	ax.set_ylabel(text_wrap(axis.replace("^2", "²")))

if INCLUDE_ERRORS:
	for i, (axis, y_min, y_original, y_max, presis, percent) in enumerate(Y_LABELS):
		if axis is None: continue
		ax = axs[i//COLUMNS, i%COLUMNS * 2 + 1]

		x = simulations[X_LABEL].values # and the corresponding data
		y = simulations[axis].values
		ɛ = simulations[axis+" error"].values

		ax.set_xlim(x.min(), x.max())
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
			ax.plot([], [], 'C3--', label="Reported error bar size")
			ax.legend()
			ax.yaxis.set_visible(False)
			continue

		if 'yield' in axis:  y_true = simulations["Yield"].values
		else:                y_true = y_original*np.ones(len(simulations.index))

		if percent:
			presis *= y_true[order]
		# if presis <= 1e-2 and "(n" in axis:
		# 	axis = axis.replace("(n", "(p")
		# 	presis *= 1e3
		# 	y_true *= 1e3

		dy = y - y_true
		# if percent:
		# 	dy /= y_true
		# print(f"{axis} is actually {np.mean(y[~np.isnan(y)])}")
		for reference_yield in REFERENCE_YIELDS:
			print(f"reference yield: {reference_yield:.2e}")
			at_yield = (~np.isnan(y)) & (np.absolute(np.log10(x/reference_yield)) <= 0.15)
			error_at_yield = np.sqrt(np.mean(np.square(dy[at_yield])))
			number_at_yield = np.sum(at_yield)
			print(f"{axis} error:    {error_at_yield:.4f} ± {error_at_yield*np.sqrt(2/number_at_yield):.4f}")

		stds = np.sqrt(smooth_average(np.square(dy[order]), bessel_correction=True))
		errs = smooth_average(ɛ[order])
		if not percent:
			ax.plot(x[order], stds, 'C0-', label="Standard deviation from actuality")
			ax.plot(x[order], np.full(order.shape, presis), 'C1--', label="Required accuracy")
			ax.plot(x[order], errs, 'C3--', label="Reported error bar size")
		else:
			ax.plot(x[order], stds, 'C0-', label="Standard deviation from actuality")
			ax.plot(x[order], np.full(order.shape, presis), 'C1--', label="Required accuracy")
			ax.plot(x[order], errs, 'C3--', label="Reported error bar size")
		# if np.max(errs)/np.max(errs) > 10:
		ax.set_yscale('log')
		ax.set_ylim(np.min(presis)*3.2e-2, np.max(presis)*5)
		# else:
		# 	ax.set_ylim(0, None)
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
