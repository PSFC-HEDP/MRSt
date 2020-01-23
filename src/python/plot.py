import numpy as np
import matplotlib.pyplot as plt
import csv
import sys

filename, width = sys.argv[1], int(sys.argv[2])

X = np.genfromtxt('working/{}_x.csv'.format(filename), delimiter=',')
Y = np.genfromtxt('working/{}_y.csv'.format(filename), delimiter=',')
Z = np.genfromtxt('working/{}_z.csv'.format(filename), delimiter=',')

plt.figure(figsize=(width/150, width/150*0.8), dpi=150)
plt.pcolormesh(X, Y, Z, cmap='plasma')
plt.xlabel("time (ns)")
plt.ylabel("position (cm)")

plt.savefig('working/{}.png'.format(filename))
