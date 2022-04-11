import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
matplotlib.rc('font', size=18)

if len(sys.argv) <= 1:
	import os
	os.chdir('../..')
	print(os.getcwd())
	# xlabel, ylabels, title, answer, n = 'Time (ns)', 'Yn (10^15/ns)\nTi (keV)\nρR (g/cm^2)', 'Trajectories', 'marginal', 3
	xlabel, ylabels, title, answer, n = 'Energy (MeV)', 'Deuterons\nDeuterons\nSignal', 'Integrated spectra', '-', 3
else:
	xlabel, ylabels, title, answer, n = sys.argv[1:]

ylabels = ylabels.split('\n')
n_curves = int(n)
assert n_curves == len(ylabels)
n_plots = len(set(ylabels))
assert n_plots <= n_curves

# load the A data from disc, where Java should have put it
XA = np.loadtxt(f'output/{title}_x.csv', delimiter=',')
YAs = [np.loadtxt(f'output/{title}_y_{i}.csv', delimiter=',') for i in range(n_curves)]
ΔAs = [np.loadtxt(f'output/{title}_err_{i}.csv', delimiter=',') for i in range(n_curves)]

# if an implosion name was given, load that as the B data
if answer != '-':
	try:
		data = np.loadtxt(f'input/trajectories {answer}.csv', delimiter=',', skiprows=1) # get the true curves
		XB = data[:,0]
		YBs = [data[:,1], data[:,4], data[:,3], np.zeros(XB.shape)] # extract the relevant info from them
		YBs[0] *= (0.1e6/1e-6)/(1e15*14.1e6*1.6e-19/1e-9)

		# while np.sum(YAs[0]*np.gradient(XA)) < np.sum(YBs[0]*(XB[1] - XB[0]))/3: # normalize the yield curves to account for any magnitude discrepancy
		# 	YBs[0] /= 10
		YBs[0] *= np.sum(YAs[0]*np.gradient(XA))/np.sum(YBs[0]*np.gradient(XB))

	except IOError:
		print(f"didn't find {answer}")
		XB, YBs = None, None

else:
	XB, YBs = None, None

# switch to a better coordinate system in x
if "(ns)" in xlabel:
	if XB is None:
		x0 = XA[np.argmax(YAs[0])]
	else:
		x0 = XB[np.argmax(YBs[0])]
	XA = (XA - x0)*1000
	if XB is not None:
		XB = (XB - x0)*1000
	xlabel = xlabel.replace("ns", "ps")

fig, host_ax = plt.subplots(figsize=(9,5))
fig.subplots_adjust(right=1 - (0.12*(n_plots-1)))
axes = [host_ax]
plots = []
j = 0
for i in range(n_curves):
	if i > 0 and ylabels[i] not in ylabels[:i]:
		j += 1

	if j > 0:
		axes.append(axes[0].twinx())
		axes[j].spines['right'].set_position(('axes', 1 + (j-1)*.18))
	if j > 1:
		axes[j].set_frame_on(True)
		axes[j].patch.set_visible(False)
		for sp in axes[j].spines.values(): sp.set_visible(False)
		axes[j].spines['right'].set_visible(True)

	rainge = {
		'Y':(0, None),
		'T':(0, 10),
		'ρ':(0, 2.0),
		'V':(-100, 100),
		'a':(-1, 1),
		'D':(0, None),
		'S':(0, None),
	}.get(ylabels[i][0], (None, None))
	YAs[i][np.isnan(ΔAs[i])] = np.nan
	plots.append(axes[j].plot(XA, YAs[i], '-o', label=ylabels[i], color=f'C{i}')[0])
	if XB is not None:
		axes[j].plot(XB, YBs[i], '--', color=f'C{i}')[0]
	axes[j].fill_between(XA, YAs[i] - ΔAs[i], YAs[i] + ΔAs[i], color='C'+str(i), alpha=0.3)
	axes[j].set_ylabel(ylabels[i])
	axes[j].set_ylim(*rainge)
	# if "MeV" in xlabel:
	# 	axes[j].set_yscale('symlog', linthresh=max(1, YAs[i].max()/1e3), linscale=1/np.log(10))

	if ylabels[i].startswith('Y'):
		Ymax = YAs[i].max(initial=0, where=np.isfinite(YAs[i]))
		lims = np.min(XA[YAs[i]/Ymax >= 1e-3]), np.max(XA[YAs[i]/Ymax >= 1e-3])
		if not all(np.isfinite(lims)):
			lims = XA[0], XA[-1]
		axes[0].set_xlim(*lims)
axes[0].set_xlabel(xlabel.replace("ns", "ps"))

# axes[0].legend(plots, [p.get_label() for p in plots])

plt.tight_layout()

# fig, axis = plt.subplots()
# for i in range(n_curves):
# 	if ylabels[i][0] == 'T':
# 		i_temp = i
# 	elif ylabels[i][0] == 'ρ':
# 		i_dens = i
# 	elif ylabels[i][0] == 'Y':
# 		valid = YAs[i]/Ymax >= 1e-3
# axis.errorbar(x=YAs[i_dens][valid], xerr=ΔAs[i_dens][valid], y=YAs[i_temp][valid], yerr=ΔAs[i_temp][valid], fmt='-.')
# axis.plot(YBs[i_dens], YBs[i_temp], '--k')
# axis.set_xlabel("ρR (g/cm^2)")
# axis.set_ylabel("Ti (keV)")

plt.savefig("inferred_trajectories.png", dpi=300)

plt.show()
