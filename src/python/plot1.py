import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import csv
import sys

xlabel, ylabels, title, n = sys.argv[1:]
ylabels = ylabels.split('\n')
n = int(n)

X = np.genfromtxt('working/{}_x.csv'.format(title), delimiter=',')
Ys = [np.genfromtxt('working/{}_y_{}.csv'.format(title, i), delimiter=',') for i in range(n)]

for i in range(n):
	plt.plot(X, Ys[i], label=ylabels[i])
plt.xlabel(xlabel)
plt.legend()
plt.title(title)

plt.show()
