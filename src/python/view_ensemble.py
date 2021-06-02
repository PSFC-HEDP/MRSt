import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import re
import sys
plt.rcParams.update({'font.family': 'sans', 'font.size': 12})
# plt.rcParams.update({'font.family': 'serif', 'font.size': 12})
import warnings
warnings.filterwarnings("ignore")

# INCLUDE_ERRORS = False
# COLUMNS = 3
# SIZE = (12, 3.5)
# MARGIN = dict(bottom=.15, top=.97, left=.07, right=.99, wspace=.40, hspace=.05)
# INCLUDE_ERRORS = True
# COLUMNS = 2
# SIZE = (16, 7)
# MARGIN = dict(bottom=.07, top=.93, left=.06, right=.99, wspace=.42, hspace=.05)
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
# SIZE = (5, 7)
# MARGIN = dict(bottom=.08, top=.92, left=.19, right=.99, wspace=.25, hspace=.05)
# INCLUDE_ERRORS = False
# COLUMNS = 1
# SIZE = (6, 4)
# MARGIN = dict(bottom=.15, top=.95, left=.15, right=.95)


if len(sys.argv) <= 1:
	# FILENAME = '../../output/ensemble-solenoid.csv'
	FILENAME = '../../output/ensemble_4_9_5_2_1500_2021-06-01.csv'
else:
	FILENAME = '../../output/'+sys.argv[1]
BIN_WIDTH = 0.3 # in bels
REFERENCE_YIELD = 1e16

X_LABEL = "Yield"

Y_LABELS = [
	# (None, 0, 0, 0, 0, False),
	# ("Total yield", 2e14, 4.36508e17, 9e17, 5e-2, True),
	# ("Bang time (ns)", 16.248, 16.2544, 16.262, 0, False),#1e-2, False),
	# ("Burn width (ps)", 53, 67.7, 82, 7, False),
	# ("Burn skewness", -1.6, -.698, -0.1, 3e-1, False),
	# ("Burn kurtosis", -0.5, 4.7, 10.5, 3, False),
	("ρR at BT (g/cm^2)", .880, .977, 1.220, 0, False),
	# ("Burn-average ρR (g/cm^2)", 0.78, .9959, 1.22, 7e-2, True),
	# ("dρR/dt at BT (g/cm^2/100ps)", -0.9, -.78, -.35, .060, False),
	("Ti at BT (keV)", 5.7, 7.43, 8.3, 0, False),#7.63848, 11, 0, False),
	# ("Burn-average Ti (keV)", 5.4, 7.144, 8.6, 7e-2, True),
	# ("dTi/dt at BT (keV/100ps)", -2.3, 2.08, 6.3, 1.9, False),
	# ("Bang time (ns)", 16.238, 16.2544, 16.272, 1e-2, False), ("Burn width (ps)", 56, 73.06, 87, 7, False),
	# ("Burn skewness", -1.3, -.546, 0.1, 3e-1, False), ("Burn kurtosis", -0.5, 4.08, 10.5, 3, False),
	# ("Burn-average vi (km/s)", -15.2, 0, 15.2, 20, False), ("dvi/dt at BT (km/s/100ps)", -110, 0, 110, 8, False),
	# ("ρR at stagnation (g/cm^2)", 0.9, 1.39, 1.7, 0, False),
	# ("Peak Ti (keV)", 5, 8.073, 13, 0, False),	#("Energy confinement time (ps)", 0, 114.8, 200, 0, False),
	# ("dTi/dt at BT (keV/(100ps))", -2.5, 2.1, 6.3, 1.9, False),
	# ("d^2Ti/dt^2 at BT (keV/ns^2)", -2200, -400, 1200, 400, False),#("dTi/dt at stagnation (keV/(100ps))", -20, 30.20, 40, 1.9, False),
	("d^2V/dt^2/V at BT (1/ns^2)", -60, 616, 1100, 50, False),
	# ("Stagnation - BT (ps)", -160, -55.5 , 40, 0, False)
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
		rolling_average(y[valid], n=min((np.sum(valid)-1)//2, 72), bessel_correction=bessel_correction), 6)
	return output


try:
	simulations = pd.read_csv(FILENAME, na_values=["Infinity"])
except FileNotFoundError:
	print("Ese archivo no existe.")
	quit()

simulations = simulations[simulations["Yield factor"] != 0]
simulations["Yield"] = simulations["Yield factor"]*4.7838e17
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

for i, (axis, y_min, y_true, y_max, presis, percent) in enumerate(Y_LABELS): # iterate through each desired plot
	if axis is None: continue
	if INCLUDE_ERRORS:
		ax = axs[i//COLUMNS, i%COLUMNS * 2 + 0] # identify the corresponding axis object
	else:
		ax = axs[i//COLUMNS, i%COLUMNS]

	x = simulations[X_LABEL].values # and the corresponding data
	y = simulations[axis].values
	ɛ = simulations[axis+" error"].values

	ax.set_xlim(x.min(), x.max()) # set up the x axis
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

	if 'keV' in axis:      y_factor = simulations["Temperature factor"].values # automatically determine the relevant scaling factor
	elif 'g/cm^2' in axis: y_factor = simulations["Down-scatter factor"].values
	elif 'yield' in axis:  y_factor = simulations["Yield factor"].values
	else:                  y_factor = np.ones(len(simulations.index))

	order = np.argsort(x) # get some useful indices of the data
	order = order[np.isfinite(simulations["Total yield"].values[order])]

	μ = smooth_average(y[order])
	σ = np.sqrt(smooth_average((y[order] - μ)**2, bessel_correction=False))

	if not percent: # plot the actual stuff
		ax.fill_between(x[order],
			y_factor[order]*y_true - presis, y_factor[order]*y_true + presis, color='#F7DFC8')
	else:
		ax.fill_between(x[order],
			y_factor[order]*y_true*(1 - presis), y_factor[order]*y_true*(1 + presis), color='#F7DFC8')
	ax.plot(x[order], y_factor[order]*y_true, 'C1-', zorder=1, label="Based on original data")
	ax.scatter(x[order], y[order], s=1, zorder=2, label="Based on fit to synthetic data")
	ax.plot(x[order], μ + σ, 'C3-', zorder=1, label="1σ variation")
	ax.plot(x[order], μ - σ, 'C3-', zorder=1)
	if y_min > 0 and y_max/y_min > 10:
		ax.set_yscale('log')
	ax.set_ylim(y_min, y_max)
	ax.set_ylabel(text_wrap(axis.replace("^2", "²")))

if INCLUDE_ERRORS:
	for i, (axis, y_min, y_true, y_max, presis, percent) in enumerate(Y_LABELS):
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

		if 'keV' in axis:      y_factor = simulations["Temperature factor"].values
		elif 'g/cm^2' in axis: y_factor = simulations["Down-scatter factor"].values
		elif 'yield' in axis:  y_factor = simulations["Yield factor"].values
		else:                  y_factor = np.ones(len(simulations.index))

		# if presis <= 1e-2 and "(n" in axis:
		# 	axis = axis.replace("(n", "(p")
		# 	presis *= 1e3
		# 	y_true *= 1e3

		dy = y - y_factor*y_true
		if percent:
			dy /= y_factor*y_true
		error_at_yield = np.sqrt(np.mean(np.square(dy[np.absolute(np.log10(x/REFERENCE_YIELD)) <= 0.15])))
		number_at_yield = np.sum(np.absolute(np.log10(x/REFERENCE_YIELD)) <= 0.15)
		print(f"{axis} error:    {error_at_yield} ± {error_at_yield*np.sqrt(2/number_at_yield)}")

		stds = np.sqrt(smooth_average(np.square(dy[order]), bessel_correction=True))
		errs = smooth_average(ɛ[order])
		if not percent:
			ax.plot(x[order], stds, 'C0-', label="Standard deviation from actuality")
			ax.plot(x[order], y_factor[order]*presis, 'C1--', label="Required accuracy")
			ax.plot(x[order], errs, 'C3--', label="Reported error bar size")
		else:
			ax.plot(x[order], stds, 'C0-', label="Standard deviation from actuality")
			ax.plot(x[order], np.full(bin_centers.shape, presis), 'C1--', label="Required accuracy")
			ax.plot(x[order], errs/(y_factor[order]*y_true), 'C3--', label="Reported error bar size")
		# if np.max(errs)/np.max(errs) > 10:
		ax.set_yscale('log')
		ax.set_ylim(presis*3.2e-2, presis*5)
		# else:
		# 	ax.set_ylim(0, None)
		if 'ield' in X_LABEL:
			ax.set_xscale('log')
		if not percent:
			if "(" in axis:
				ax.set_ylabel(text_wrap(axis[:axis.index("(")].replace(" at BT", "").replace("^2", "²") + "error " + axis[axis.index("("):]))
			else:
				ax.set_ylabel(text_wrap(axis+" error"))
		else:
			ax.set_ylabel(text_wrap(re.sub(r'(\([^)]+\))?$', '', axis) + " error"))

config = '-'+FILENAME[23:36] if len(FILENAME) > 23 else ''
fig.savefig('../../output/mrst{}.eps'.format(config), dpi=300)
fig.savefig('../../output/mrst{}.png'.format(config), dpi=300)

plt.show()
