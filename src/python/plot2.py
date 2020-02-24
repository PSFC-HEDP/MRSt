import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
import csv
import sys

xlabel, ylabel, title = sys.argv[1:]

X = np.genfromtxt('working/{}_x.csv'.format(title), delimiter=',')
Y = np.genfromtxt('working/{}_y.csv'.format(title), delimiter=',')
Z = np.genfromtxt('working/{}_z.csv'.format(title), delimiter=',')

plt.pcolormesh(X, Y, Z, cmap='plasma', norm=colors.SymLogNorm(vmin=0, vmax=Z.max(), linthresh=Z.max()/100, linscale=1))
plt.ylim(12, 16)
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.colorbar()
plt.title(title)

plt.tight_layout()
plt.show()
