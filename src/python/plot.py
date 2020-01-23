import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

xlabel, ylabel, title = sys.argv[1:]

X = np.genfromtxt('working/{}_x.csv'.format(title), delimiter=',')
Y = np.genfromtxt('working/{}_y.csv'.format(title), delimiter=',')
Z = np.genfromtxt('working/{}_z.csv'.format(title), delimiter=',')

plt.pcolormesh(X, Y, Z, cmap='plasma')
plt.xlabel(xlabel)
plt.ylabel(ylabel)
plt.title(title)

plt.show()
