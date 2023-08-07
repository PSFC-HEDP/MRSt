import h5py
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sys
matplotlib.rc('font', size=18)

if len(sys.argv) <= 1:
	import os
	os.chdir('../..')
	print(os.getcwd())
	xlabel, ylabels, title, n = 'Time (ns)', 'Ti (keV)\nρR (g/cm^2)\nYn (10^15/ns)', 'data', 3
	xlabel, ylabels, title, answer, n = 'Time (ns)', 'Yn (10^15/ns)\nTi (keV)\nρR(g/cm^2)', 'Trajectories', 'input/rby 001 {}.csv', 3
	# xlabel, ylabels, title, answer, n = 'Energy (MeV)', 'Deuterons\nDeuterons\nSignal', 'Integrated spectra', '-', 3
else:
	xlabel, ylabels, title, answer, n = sys.argv[1:]

# collect the basic arguments
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
	filename = answer.format("trajectories")
	try:
		with open(filename, "r") as f:
			data_header = f.readline().split(",")
	except IOError:
		print(f"didn't find {filename}")
		XB, YBs = None, [None]*n_curves
	else:
		data = np.loadtxt(answer.format("trajectories"), delimiter=',', skiprows=1) # get the true curves
		XB = data[:, 0]
		if "rhor" in data_header[3].lower():
			YBs = [data[:, 1], data[:, 4], data[:, 3], np.zeros(XB.shape)] # extract the relevant info from them
		elif "rhor" in data_header[1].lower():
			YBs = [data[:, 3], data[:, 4], data[:, 1], np.zeros(XB.shape)]
		elif "dsr" in data_header[2].lower():
			YBs = [data[:, 1], data[:, 3], data[:, 2]*20.4]
		else:
			raise Exception(f"do I have to actually read this format?: {data_header}")

		while np.ptp(XB) < np.ptp(XA)/100:
			XB *= 1e+3
		YBs[0] *= np.sum(YAs[0]*np.gradient(XA))/np.sum(YBs[0]*np.gradient(XB))

else:
	XB, YBs = None, [None]*n

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

fig, host_ax = plt.subplots(figsize=(9.5, 4.5))
fig.subplots_adjust(right=1 - (0.12*(n_plots-1)))
axes = [host_ax]
plots = []
i_ax = 0
for i in range(n_curves):
	if i > 0 and ylabels[i] not in ylabels[:i]:
		i_ax += 1

	if i_ax > 0:
		axes.append(axes[0].twinx())
		axes[i_ax].spines['right'].set_position(('axes', 1 + (i_ax - 1)*.18))
	if i_ax > 1:
		axes[i_ax].set_frame_on(True)
		axes[i_ax].patch.set_visible(False)
		for sp in axes[i_ax].spines.values():
			sp.set_visible(False)
		axes[i_ax].spines['right'].set_visible(True)

	rainge = {
		'Y': (0, None),
		'T': (0, 18.),
		'ρ': (0, 1.0),
		'V': (-100, 100),
		'a': (-1, 1),
		'D': (0, None),
		'S': (0, None),
	}.get(ylabels[i][0], (None, None))
	YAs[i][np.isnan(ΔAs[i])] = np.nan
	plots.append(axes[i_ax].plot(XA, YAs[i], '-o', label=ylabels[i], color=f'C{i}', zorder=20)[0])
	if XB is not None:
		axes[i_ax].plot(XB, YBs[i], '--', color=f'C{i}')
	axes[i_ax].fill_between(XA, YAs[i] - ΔAs[i], YAs[i] + ΔAs[i], color='C' + str(i), alpha=0.3, zorder=10)
	axes[i_ax].set_ylabel(ylabels[i])
	axes[i_ax].set_ylim(*rainge)
	# if "MeV" in xlabel:
	# 	axes[j].set_yscale('symlog', linthresh=max(1, YAs[i].max()/1e3), linscale=1/np.log(10))

	if ylabels[i].startswith('Y'):
		Ymax = YAs[i].max(initial=0, where=np.isfinite(YAs[i]))
		# lims = np.min(XA[YAs[i]/Ymax >= 1e-3]), np.max(XA[YAs[i]/Ymax >= 1e-3])
		lims = (-150, 100)
		if not all(np.isfinite(lims)):
			lims = XA[0], XA[-1]
		axes[0].set_xlim(*lims)
axes[0].xaxis.set_visible(False)
axes[-1].xaxis.set_visible(True)
axes[0].patch.set_visible(False)
axes[-1].patch.set_visible(True)
axes[-1].set_xlabel(xlabel.replace("ns", "ps"))
axes[-1].grid("on")
for i_ax in range(len(axes)):
	axes[i_ax].set_zorder(len(axes) - i_ax)

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

with h5py.File("output/inferred_trajectories.h5", "w") as f:
	for name, value in [("time true", XB), ("time fit", XA),
	                    ("burn true", YBs[0]), ("burn fit", YAs[0]), ("burn error", ΔAs[0]),
	                    ("Tion true", YBs[1]), ("Tion fit", YAs[1]), ("Tion error", ΔAs[1]),
	                    # ("rhoR true", YBs[2]), ("rhoR fit", YAs[2]), ("rhoR error", ΔAs[2]),
	                    ]:
		if value is not None:
			dataset = f.create_dataset(name, value.shape)
			dataset[...] = value
plt.savefig("output/inferred_trajectories.png", dpi=300, transparent=True)
plt.savefig("output/inferred_trajectories.eps")

plt.show()
