import matplotlib.pyplot as plt
import numpy as np
import os


if "python" in os.getcwd():
	os.chdir('../..')

slit_lengths = np.loadtxt("output/focal plane lengths.csv")
slit_widths = np.loadtxt("output/focal plane widths.csv")
slit_positions = np.loadtxt("output/focal plane positions.csv")
energies = np.loadtxt("output/focal plane energies.csv")
particles = np.loadtxt("output/focal plane particles.csv",
                       delimiter=",")
particles = particles.reshape((particles.shape[0], -1, 3))

# go thru the energies and get some moments of the distributions
quantiles = np.linspace(0, 1, 9)
x = []
y = []
for e, Y in zip(energies, particles):
	x.append(Y[:,0].mean()*1e2)
	y.append([])
	for q in quantiles:
		y[-1].append(np.quantile(Y[:,1], q)*1e2)
x = np.array(x)
y = np.array(y)

plt.figure(figsize=(9, 4))

for i in range(quantiles.size):
	plt.plot(x, y[:,i], linewidth=(2-i%4)*1.1, color='C0')

for e, Y in zip(energies[::5], particles[::5]):
	# histogram, x_bins, y_bins = np.histogram2d(Y[:,0]*1e2, Y[:,1]*1e2, bins=20)
	# x_points, y_points = np.meshgrid((x_bins[:-1] + x_bins[1:])/2, (y_bins[:-1] + y_bins[1:])/2, indexing="ij")
	# plt.contour(x_points, y_points, histogram, levels=histogram.max()*np.linspace(0, 1, 21), colors="C2")
	x_min, x_max = Y[:,0].min()*1e2, Y[:,0].max()*1e2
	y_min, y_max = Y[:,1].min()*1e2, Y[:,1].max()*1e2
	plt.plot([x_min, x_min, x_max, x_max, x_min],
	         [y_min, y_max, y_max, y_min, y_min],
	         linewidth=0.8, color='C2')

for x, w, h in zip(slit_positions, slit_lengths, slit_widths):
	plt.plot(np.multiply([x -w/2, x -w/2, x +w/2, x +w/2, x -w/2], 1e2),
	         np.multiply([  -h/2,    h/2,    h/2,   -h/2,   -h/2], 1e2), 'C1')

# plt.plot(
# 100*np.array(
# [-0.25852328457373486, -0.23262669584333442, -0.20837751894256018, -0.18556162751515398, -0.1641414954559396, -0.14396825385559583, -0.12496638330199167, -0.1070535963324801, -0.09014272570712327, -0.07415937627674185, -0.05901801577777054, -0.044679226868282206, -0.031072375412029635, -0.018142456392673315, -0.005817982864627899, 0.005916441531215929, 0.01713761869335928, 0.027890181853226073, 0.03819083039043938, 0.048107810945184674, 0.05767247356247733, 0.06691087775499732, 0.07587672019652382, 0.08456343953911179, 0.09305478387935032, 0.10135956773479679, 0.10949378147473755, 0.1175026786505859, 0.12539752232684992, 0.1332143847630211]
# ), 50*np.array(
# [0.02155079415459317, 0.019882981189509905, 0.01784787546852584, 0.01586133723038026, 0.013814513724120523, 0.011884319375824504, 0.009839432738702796, 0.007865206845313517, 0.005812811557314994, 0.0038660373492532514, 0.0019757989637796437, 3.18422374245041E-4, 0.0020440698791455722, 0.0038074519080877854, 0.005512288629432701, 0.00716757236424969, 0.008829477346068324, 0.010269790901380773, 0.01157425766918338, 0.013003040471176593, 0.014188571776518492, 0.015332041982655795, 0.016369163648262318, 0.017259504960887108, 0.018004987270485998, 0.01863329370205195, 0.019269334906784102, 0.01980349094898133, 0.020095160944690785, 0.020210522988862126]
# ), 'C3-'
# )

# plt.axis('equal')
plt.ylim(-2, 2)
plt.xlabel("x (cm)")
plt.ylabel("y (cm)")
plt.tight_layout()
plt.savefig('temp.png', dpi=300)
plt.show()
