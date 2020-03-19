import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

simulations = pd.read_csv('../../working/ensemble.csv')
print(simulations)

plt.figure()
plt.xscale('log')
plt.scatter(simulations['Yield factor'], simulations['Bang time (ns)'])
plt.xlabel("Yield factor")
plt.ylabel("Bang time (ns)")
plt.tight_layout()
plt.show()
