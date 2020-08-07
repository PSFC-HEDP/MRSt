import numpy as np
import matplotlib.pyplot as plt

correct = np.loadtxt("../../data/Yn-rR-Ti_150327_16p26 - Yn-rR-Ti_150327_16p26.csv", skiprows=1, delimiter=',')

wrong_time = np.loadtxt("../../working/data_x.csv") - .1
wrong_yield = np.loadtxt("../../working/Data_y_2.csv")/1e6
wrong_temp = np.loadtxt("../../working/Data_y_0.csv")
wrong_rhoR = np.loadtxt("../../working/Data_y_1.csv")

plt.plot(correct[:,0], correct[:,1], '-C0', label="Yield (10^18)")
plt.plot(correct[:,0], correct[:,3], '-C1', label="Temperature (keV)")
plt.plot(correct[:,0], correct[:,4], '-C2', label="ρR DT (g/cm^2)")
trustworthy = wrong_yield >= wrong_yield.max()/1e3
plt.plot(wrong_time[trustworthy], wrong_yield[trustworthy], '--C0')
plt.plot(wrong_time[trustworthy], wrong_temp[trustworthy], '--C1')
plt.plot(wrong_time[trustworthy], wrong_rhoR[trustworthy], '--C2')
# plt.yscale('log')
plt.legend()
plt.show()
